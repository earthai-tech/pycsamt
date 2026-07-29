# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Electromagnetic utilities.

Aggregates helpers for AMT/CSAMT processing.
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

from ..api.bunch import Bunch
from ..api.typing import EDIO, ZO, ArrayLike, DType, NDArray
from ..compat.aliases import compat_alias
from ..context import nullify_output
from ..exceptions import EdIDataError, EMError
from ..seg.validation import IsEdi
from .arrayops import (
    concat_array_from_list,
    is_iterable,
    reshape,
)
from .cleaner import ismissing
from .conversion import convert_value
from .plot import plot_errorbar
from .stats import get_confidence_ratio, remove_outliers
from .validation import assert_ratio, isinstance_relaxed

__all__ = [
    "check_em_kind",
    "extract_z_list",
    "parse_tensor",
    "compute_qc",
    "full_freq",
    "tensor2d",
    "align_tensor",
    "export_edis",
    "plot_confidence",
    "plot_strike",
    "plot_tensors",
    "plot_station_tensors",
    "wrap_phase",
    "plot_lcurve",
]

_BACKWARD_SINCE = "2.0.0"
_BACKWARD_REMOVE = "2.17.0"


def check_em_kind(objs, /) -> str:  # _assert_z_or_edi_objs
    """
    Validate a collection of EM objects and return their common kind.

    Ensures that *all* elements in ``objs`` are instances of either
    :class:`Edi` *or* :class:`Z`, but not a mix of both. Returns the
    common kind as ``"EDI"`` or ``"Z"``. Raises if the collection is
    empty, contains non-EM objects, or mixes kinds.

    Parameters
    ----------
    objs : iterable
        Iterable of objects expected to be all :class:`Edi` *or* all
        :class:`Z`. Strings and bytes are not valid inputs.

    Returns
    -------
    str
        ``"EDI"`` if all objects are EDI instances, else ``"Z"`` if
        all are Z instances.

    Raises
    ------
    TypeError
        If ``objs`` is not an iterable, or is a string/bytes, or is
        empty.
    EMError
        If an element is neither :class:`Edi` nor :class:`Z`, or if
        the iterable mixes :class:`Edi` and :class:`Z`.

    Notes
    -----
    Uses :func:`is_instance_extended` to be robust to class reloads
    and alternate import paths.

    Examples
    --------
    >>> from pycsamt.utils.em import check_em_kind
    >>> # assuming `eds` is a list of Edi instances
    >>> check_em_kind(eds)
    'EDI'
    >>> # assuming `zs` is a list of Z instances
    >>> check_em_kind(zs)
    'Z'
    """
    # Local imports to avoid import cycles during module init.
    from ..seg.edi import (
        EDIFile as Edi,  # v2: EDIFile replaced Edi
    )
    from ..z.z import Z

    # Basic iterable validation (exclude str/bytes).
    if isinstance(objs, (str, bytes)):
        raise TypeError(
            "Input must be an iterable of Edi or Z objects; got 'str'/'bytes'."
        )
    if not hasattr(objs, "__iter__"):
        raise TypeError("Input must be an iterable of Edi or Z objects.")

    # Materialize once; allows single pass over generators.
    items = list(objs)
    if len(items) == 0:
        raise TypeError("Input iterable is empty.")

    # Determine membership per class with robust isinstance checks.
    is_edi = [bool(isinstance_relaxed(o, Edi)) for o in items]
    is_z = [bool(isinstance_relaxed(o, Z)) for o in items]

    # Any non-EM object?
    if not all(e or z for e, z in zip(is_edi, is_z)):
        raise EMError("All elements must be Edi or Z instances.")

    # Mixed kinds?
    if any(is_edi) and any(is_z):
        raise EMError("Do not mix Edi and Z objects in the same collection.")

    # Uniform kind → return the label.
    return "EDI" if all(is_edi) else "Z"


def extract_z_list(objs, /):  # get_z_from
    """
    Return a list of :class:`Z` objects from EDI or Z inputs.

    If ``objs`` is a collection of :class:`Edi`, each element's
    ``.Z`` attribute is extracted. If ``objs`` is already a
    collection of :class:`Z`, it is returned as a plain list.

    Parameters
    ----------
    objs : iterable
        Iterable of EM objects. Must be all :class:`Edi` or all
        :class:`Z`. See :func:`check_em_kind`.

    Returns
    -------
    list
        List of :class:`Z` instances.

    Raises
    ------
    TypeError
        If ``objs`` is not an iterable or is empty.
    EMError
        If ``objs`` mixes :class:`Edi` and :class:`Z`, contains
        non-EM objects, or an :class:`Edi` is missing the ``.Z``
        attribute.

    Notes
    -----
    Uses :func:`check_em_kind` to validate homogeneity.

    Examples
    --------
    >>> zs = extract_z_list(list_of_edi)  # from EDI inputs
    >>> zs = extract_z_list(list_of_z)  # already Z
    """
    kind = check_em_kind(objs)
    items = list(objs)

    if kind == "Z":
        # Normalize to list; caller may have passed a generator.
        return items

    # kind == "EDI": extract `.Z` from each EDI
    try:
        return [e.Z for e in items]
    except AttributeError as exc:
        raise EMError("An Edi object is missing the '.Z' attribute.") from exc


def parse_tensor(  # validate_tensor
    out: str = "resxy",
    *,
    tensor: str | None = None,
    component: str | None = None,
    kind: str = "complex",
    **kws,
) -> tuple[str, str | None]:
    """
    Parse and validate a tensor request, returning name and component.

    This helper normalizes shorthand like ``'resxy'`` or explicit
    pairs like ``tensor='z', component='xy'`` and validates the
    desired numeric *kind* (e.g., complex, real, imag, modulus).

    Parameters
    ----------
    out : str, default='resxy'
        Compact token specifying the tensor and component, e.g.
        ``'resxy'``, ``'zxy'``, ``'phaseyx'``, or a frequency
        request such as ``'freq'``. When both ``tensor`` and
        ``component`` are given, they override ``out``.
    tensor : str, optional
        Tensor name or alias. Accepted values include:
        ``'z'``, ``'tensor'``, ``'res'``, ``'rho'``, ``'rhoa'``,
        ``'phase'``, ``'phs'``, ``'freq'``, ``'frequency'``.
    component : str, optional
        EM component among ``'xx'``, ``'xy'``, ``'yx'``, ``'yy'``.
        Required for ``'z'``, ``'resistivity'`` and ``'phase'``.
    kind : {'complex','real','imag','modulus'}, default='complex'
        Numeric form of the tensor to be later extracted. Aliases
        are accepted, e.g. ``'re'``→``'real'``, ``'im'``→``'imag'``,
        ``'abs'``/``'mod'``→``'modulus'``, ``'reel'``→``'real'``.
    **kws
        Extra keywords ignored here. Kept for API symmetry with
        callers that also manage frequency expansion.

    Returns
    -------
    (name, comp) : tuple of (str, Optional[str])
        Normalized tensor name and its component. The ``name`` is
        one of ``'z'``, ``'resistivity'``, ``'phase'``, ``'_freq'``,
        or includes the ``'_err'`` suffix for error arrays where
        applicable. ``comp`` is ``'xx'``, ``'xy'``, ``'yx'``,
        ``'yy'`` or ``None`` for frequency.

    Raises
    ------
    EMError
        If only one of ``tensor`` or ``component`` is provided.
    ValueError
        If the parsed tokens are invalid, the component is missing
        for a tensor that requires it, or ``kind`` is unknown.

    Examples
    --------
    >>> parse_tensor("zxy")
    ('z', 'xy')
    >>> parse_tensor(tensor="res", component="yx")
    ('resistivity', 'yx')
    >>> parse_tensor("freq")
    ('_freq', None)
    >>> parse_tensor("resx")
    Traceback (most recent call last):
        ...
    ValueError: 'Resistivity' component is missing...
    """
    # Validate paired arguments; both or none must be given.
    if (tensor and not component) or (component and not tensor):
        raise EMError("Provide both 'tensor' and 'component' or neither.")

    # When both are given, build a compact 'out' token.
    if tensor and component:
        out = f"{str(tensor)}{str(component)}"

    # Normalize tokens
    s_out = str(out).strip().lower()
    s_kind = str(kind).strip().lower()

    # Map kind aliases to canonical names
    kind_map = {
        "imaginary": "imag",
        "im": "imag",
        "re": "real",
        "reel": "real",
        "abs": "modulus",
        "mod": "modulus",
        "amplitude": "modulus",
        "cmp": "complex",
    }
    s_kind = kind_map.get(s_kind, s_kind)

    allowed_kinds = {"modulus", "imag", "real", "complex"}
    if s_kind not in allowed_kinds:
        raise ValueError(
            f"Unacceptable argument {s_kind!r}. Expect 'modulus', "
            f"'imag', 'real', or 'complex'."
        )

    # Regex buckets:
    #  - tensor name or alias
    #  - component token (xx|xy|yx|yy)
    #  - error token (err|error)
    r_name = re.compile(r"(res|rho|rhoa|phase|phs|z|tensor|freq|frequency)")
    r_comp = re.compile(r"(xx|xy|yx|yy)")
    r_err = re.compile(r"(err|error)")

    m_name = r_name.search(s_out)
    m_comp = r_comp.search(s_out)
    m_err = r_err.search(s_out)

    if m_name is None:
        raise ValueError(
            f"{out!r} does not match any of 'resistivity', 'phase', "
            f"'tensor', or 'frequency'."
        )

    name = m_name.group()
    comp = m_comp.group() if m_comp is not None else None
    err = m_err is not None

    # Canonicalize the tensor name
    if name in ("res", "rho", "rhoa"):
        name = "resistivity"
    elif name in ("phase", "phs"):
        name = "phase"
    elif name in ("z", "tensor"):
        name = "z"
    elif name in ("freq", "frequency"):
        name = "_freq"

    # Enforce component presence when required
    if comp is None and name in ("z", "resistivity", "phase"):
        tlabel = "Tensor" if name == "z" else name.title()
        raise ValueError(
            f"{tlabel!r} component is missing. Use e.g. "
            f"'{name}_xy' for 'xy' component."
        )

    # Append error suffix only when meaningful
    if err and name not in ("_freq", "z"):
        name = f"{name}_err"

    # Return normalized (name, component). `s_kind` is validated
    # here for callers that inspect it, but not returned; this
    # function's contract matches the legacy behavior.
    return name, comp


def align_tensor(  # fittensor
    ref_freq: ArrayLike,
    site_freq: ArrayLike,
    z: NDArray[DType[complex]],
    fill_value: float | None = np.nan,
) -> NDArray[DType[complex]]:
    """
    Align a tensor component to a reference frequency grid.

    The reference frequency grid (``ref_freq``) is assumed to be the
    *complete* set of clean frequencies for the survey. Site-level
    measurements (``site_freq``) may be missing some frequencies due
    to interferences. This function maps the provided tensor values
    (``z``) to the reference grid and fills gaps with ``fill_value``.

    Parameters
    ----------
    ref_freq : ArrayLike
        Reference frequency grid collected in the field. It should
        contain *all* survey frequencies.
    site_freq : ArrayLike
        Frequencies measured at a site. All values must be present in
        ``ref_freq`` (i.e., no out-of-grid frequencies).
    z : ndarray of complex
        Tensor component values measured at ``site_freq``. Typically
        the real or imaginary part for one of ``xx``, ``xy``, ``yx``,
        or ``yy``. Length must match ``site_freq``.
    fill_value : float, default=nan
        Value used to fill positions in the reference grid where the
        site tensor is missing.

    Returns
    -------
    ndarray of complex
        Array aligned to ``ref_freq`` with gaps filled by
        ``fill_value``. The dtype matches ``z.dtype``.

    Raises
    ------
    EMError
        If the number of mappable positions inferred from
        ``site_freq`` does not match the length of ``z``.
    ValueError
        If input shapes are inconsistent.

    Notes
    -----
    Internally uses :func:`ismissing` to identify positions in
    ``ref_freq`` that correspond to ``site_freq``. The function
    does *not* interpolate values; it only aligns and fills gaps.

    Examples
    --------
    >>> ref_freq = np.linspace(7e7, 1.0, 20)
    >>> site_freq = np.hstack([ref_freq[:7], ref_freq[12:]])
    >>> z = np.random.randn(len(site_freq)) + 1j * np.random.randn(
    ...     len(site_freq)
    ... )
    >>> z_aligned = align_tensor(ref_freq, site_freq, z)
    >>> np.isnan(z_aligned).sum()  # gaps inserted
    5
    """
    # Validate and normalize reference frequency as a 1-D array.
    ref_freq = np.asarray(ref_freq, dtype=float).ravel()

    # Map site frequencies onto the reference grid.
    # `ismissing` returns (filled_ref, absence_mask) where absence_mask
    # is True where ref_freq[i] is NOT in site_freq (gap positions).
    freq_like_ref, absence_mask = ismissing(
        refarr=ref_freq,
        arr=site_freq,
        return_index="mask",
        fill_value=fill_value,
    )
    # presence_mask is True where site values should be placed
    presence_mask = ~absence_mask

    # Prepare an output array filled with the sentinel value.
    z_aligned = np.full_like(
        freq_like_ref,
        fill_value=fill_value,
        dtype=z.dtype,
    )

    # Ensure shape compatibility between the mask and provided z.
    z_flat = reshape(z)  # robust 1-D view
    if len(z_aligned[presence_mask]) != len(z_flat):
        raise EMError(
            "Inconsistent frequency lengths: tensor values do not "
            "match the number of mapped site frequencies. Frequencies "
            "in `z` must be a subset of the complete `ref_freq`. "
            f"Got {len(z_flat)} values for {presence_mask.sum()} "
            "mapped positions."
        )

    # Place site tensor values at their positions on the reference grid.
    z_aligned[presence_mask] = z_flat

    return z_aligned


def export_edis(  # exportEDIS
    edi_objs: list[EDIO],
    new_z: list[ZO],
    savepath: str | None = None,
    **kws,
) -> None:
    """
    Export new EDI files from a batch of EDI objects and Z tensors.

    Applies updated impedance tensors to each input EDI object and
    writes new EDI files. This is typically used after applying
    corrections or replacements to the impedance tensor.

    Parameters
    ----------
    edi_objs : list of :class:`Edi`
        Collection of EDI objects. All elements must be instances of
        :class:`Edi` (no ``Z`` objects are allowed here).
    new_z : list of ndarray (nfreq, 2, 2), complex
        Collection of impedance tensors matching ``edi_objs`` one-to-
        one. Each tensor is a 3-D complex array with shape
        ``(n_freq, 2, 2)``.
    savepath : str, optional
        Directory to write the new EDI files. If ``None``, the EDI
        writer decides (often the current directory).
    **kws
        Extra arguments forwarded to the EDI writer method
        ``Edi.write_new_edifile`` (e.g., naming options).

    Returns
    -------
    None
        Files are written as a side effect. The underlying writer may
        return paths, but this function does not collect them.

    Raises
    ------
    EdIDataError
        If ``edi_objs`` does not contain exclusively EDI objects.
    ValueError
        If the lengths of ``edi_objs`` and ``new_z`` differ.

    See Also
    --------
    exportedi
        Helper for exporting a single EDI (if available in the API).

    Examples
    --------
    >>> # edi_objs: list[Edi], z_list: list[np.ndarray]
    >>> export_edis(edi_objs, z_list, savepath="out/")
    """
    # Ensure we truly have EDI objects (no Z objects mixed in).
    kind = check_em_kind(edi_objs)
    if kind != "EDI":
        raise EdIDataError("Expected a collection of EDI objects; got non-EDI items.")

    # Normalize containers (accept tuples/generators, exclude strings).
    edi_objs = is_iterable(
        edi_objs,
        exclude_string=True,
        transform=True,
    )
    new_z = is_iterable(
        new_z,
        exclude_string=True,
        transform=True,
    )

    if len(edi_objs) != len(new_z):
        raise ValueError(
            "Length mismatch between 'edi_objs' and 'new_z'. "
            f"Got {len(edi_objs)} EDI(s) and {len(new_z)} tensor(s)."
        )

    # Write each new EDI; forward writer-specific kwargs untouched.
    for e, z in zip(edi_objs, new_z):
        # Each EDI object is assumed to expose this method.
        e.write_new_edifile(
            new_Z=z,
            savepath=savepath,
            **kws,
        )


def plot_confidence(  # plot_confidence_in
    z_or_edis_obj_list: list[EDIO | ZO],
    /,
    tensor: str = "res",
    view: str = "1d",
    drop_outliers: bool = True,
    distance: float | None = None,
    c_line: bool = False,
    view_ci: bool = True,
    figsize: tuple[float, float] = (6.0, 2.0),
    fontsize: float = 4.0,
    dpi: int = 300,
    top_label: str = "Stations",
    rotate_xlabel: float = 90.0,
    fbtw: bool = True,
    savefig: str | None = None,
    **plot_kws: Any,
):
    """
    Plot confidence diagnostics from tensor errors for EM data.

    The default :term:`tensor` for confidence evaluation is the
    resistivity error at TE mode (``'xy'``). This plot helps decide
    which frequencies (and stations) are reliable, recoverable, or
    should be discarded before further processing.

    Three confidence levels are highlighted:

    - **High**: :math:`conf \\ge 0.95`
    - **Soft**: :math:`0.5 \\le conf < 0.95` (often recoverable)
    - **Bad**: :math:`conf < 0.5` (usually discard)

    Parameters
    ----------
    z_or_edis_obj_list : list of EDI or Z
        Collection of :class:`Edi` or :class:`Z` objects.
    tensor : str, default='res'
        Tensor selector. Accepted aliases include resistivity
        (``'res'``, ``'rho'``, ``'rhoa'``), phase (``'phase'``,
        ``'phs'``), or ``'z'``. Error arrays are used for
        resistivity/phase automatically.
    view : {'1d', '2d'}, default='1d'
        Plot as a 1-D profile (by station) or as a 2-D map
        (frequency vs. station).
    drop_outliers : bool, default=True
        If ``True``, suppress outliers in the error tensor before
        plotting (filled with ``nan``).
    distance : float, optional
        Inter-station distance. Used to scale the x-axis in 1-D view.
        If ``None``, a unit spacing of 1 is used.
    c_line : bool, default=False
        If ``True`` and ``view='2d'``, overlay the confidence line.
    view_ci : bool, default=True
        If ``True``, show markers indicating confidence classes.
    figsize : tuple of float, default=(6.0, 2.0)
        Matplotlib figure size.
    fontsize : float, default=4.0
        Base font size used for labels and ticks.
    dpi : int, default=300
        Figure resolution in dots per inch.
    top_label : str, default='Stations'
        Title used for the top x-axis (station labels).
    rotate_xlabel : float, default=90.0
        Rotation angle for station labels on the top x-axis.
    fbtw : bool, default=True
        In 1-D view, fill between the curve and confidence bands.
    savefig : str, optional
        Path to save the figure. If ``None``, the figure is shown.

    Returns
    -------
    matplotlib.axes.Axes
        The Matplotlib Axes with the plotted content.

    Notes
    -----
    Internally, the function computes an error tensor for the chosen
    ``tensor`` (resistivity/phase use error arrays). Confidence is
    aggregated across stations and displayed either as a 1-D line or
    a 2-D image with categorical markers.

    Examples
    --------
    >>> ax = plot_confidence(
    ...     emobj.ediObjs_, distance=20, view="2d", figsize=(6, 2)
    ... )
    >>> ax = plot_confidence(
    ...     emobj.ediObjs_, distance=20, view="1d", figsize=(6, 3), fontsize=5
    ... )
    """
    # Lazy import shared by 1D path only; 2D plot2d imported inside block.
    from .plot import _get_xticks_formatage

    # Normalize options
    view = str(view).strip().lower()
    tkn = str(tensor).strip().lower()

    if view not in {"1d", "2d"}:
        raise ValueError("Invalid 'view'. Expect '1d' or '2d'.")

    # Force error arrays for resistivity/phase requests.
    # Examples of accepted aliases are broadened here.
    if tkn.startswith(("res", "rho", "rhoa")):
        tkn = "res_err"
    elif tkn.startswith(("phase", "phs")):
        tkn = "phase_err"
    # else keep user token (e.g., 'zxy_err', 'z', etc.)

    # Fetch error tensor (2-D) and its frequency vector.
    # rerr.shape = (n_freq, n_stations)
    rerr, freqs = tensor2d(
        z_or_edis_obj_list,
        tensor=tkn,
        return_freqs=True,
    )

    # Confidence ratio per station (across frequency axis).
    ratio = get_confidence_ratio(rerr)

    # Build class masks and display attributes.
    conf_props = dict(
        high_cf=(
            np.where(ratio >= 0.95)[0],
            "#15B01A",
            r"$conf. \geq 0.95$",
        ),
        soft_cf=(
            np.where((ratio < 0.95) & (ratio >= 0.5))[0],
            "#FF81C0",
            r"$0.5 \leq conf. < 0.95$",
        ),
        bad_cf=(
            np.where(ratio < 0.5)[0],
            "#650021",
            r"$conf. < 0.5$",
        ),
    )

    # Station spacing and x-coordinate (1-D distance axis).
    spacing = distance or 1.0
    d = np.arange(rerr.shape[1]) * convert_value(spacing)

    # Label for colorbar / y-axis units.
    clab = (
        r"resistivity ($\Omega\cdot m$)"
        if tkn.startswith("res")
        else (r"phase ($^\circ$)" if tkn.startswith("phase") else tkn)
    )

    if view == "2d":
        from .plot import (
            plot2d,  # only needed for 2D path
        )

        # Optional outlier removal for cleaner 2-D maps.
        ar2d = remove_outliers(rerr, fill_value=np.nan) if drop_outliers else rerr

        ax = plot2d(
            ar2d,
            cmap="binary",
            cb_label=f"Error in {clab}",
            top_label=top_label,
            rotate_xlabel=rotate_xlabel,
            distance=spacing,
            y=np.log10(freqs),
            fig_size=figsize,
            fig_dpi=dpi,
            font_size=fontsize,
        )
    else:
        # 1-D confidence profile along stations.
        fig, ax = plt.subplots(
            figsize=figsize,
            dpi=dpi,
        )
        ax.plot(
            d,
            ratio,
            "ok-",
            markersize=2.0,
            **plot_kws,
        )

        if fbtw:
            # Lower band starts at min between 0.5 and observed min.
            min_ci = 0.5 if ratio.min() <= 0.5 else ratio.min()

            # Define the two confidence reference lines.
            y_high = np.repeat(0.95, len(ratio))
            y_soft = np.repeat(min_ci, len(ratio))

            # Masks for the three regions.
            masks = (
                (ratio >= 0.95),
                (ratio < 0.95) & (ratio >= min_ci),
                (ratio < min_ci),
            )

            # Fill regions with class colors.
            for i, m in enumerate(masks):
                # Skip if mask has no True; avoids warnings.
                if not np.any(m):
                    continue
                ax.fill_between(
                    d,
                    ratio,
                    y_soft if i != 0 else y_high,
                    where=m,
                    facecolor=list(conf_props.values())[i][1],
                    alpha=0.3,
                )
                ax.axhline(
                    y=min_ci if i != 0 else 0.95,
                    color="k",
                    linestyle="--",
                    lw=1.0,
                )

        ax.set_xlabel(
            "Distance (m)",
            fontsize=1.2 * fontsize,
            fontdict={"weight": "bold"},
        )
        ax.set_ylabel(
            "Confidence ratio ×100 (%)",
            fontsize=1.2 * fontsize,
            fontdict={"weight": "bold"},
        )
        ax.tick_params(labelsize=fontsize)
        ax.set_xlim([d.min(), d.max()])

        # Secondary x-axis with station labels on top.
        ax_top = ax.twiny()
        ax_top.set_xticks(range(len(d)), minor=False)
        ax_top.tick_params(labelsize=fontsize)

        _get_xticks_formatage(
            ax_top,
            range(len(d)),
            fmt="E{:02}",
            auto=True,
            rotation=rotate_xlabel,
        )
        ax_top.set_xlabel(
            top_label,
            fontdict={"size": fontsize, "weight": "bold"},
        )

    # Overlay confidence markers / line if requested.
    if view_ci:
        if view == "2d" and c_line:
            # Pull style kwargs if present; keep defaults minimal.
            c = plot_kws.pop("c", "r")
            lw = plot_kws.pop("lw", 0.5)
            ls = plot_kws.pop("ls", "-")

            ax.plot(
                d,
                ratio * np.log10(freqs).max(),
                ls=ls,
                c=c,
                lw=lw,
                label="Confidence line",
            )

        for idx, color, label in conf_props.values():
            if len(idx) == 0:
                continue
            scale = np.log10(freqs).max() if view == "2d" else 1.0
            ax.scatter(
                d[idx],
                ratio[idx] * scale,
                marker="o",
                edgecolors="k",
                color=color,
                label=label,
            )

        ax.legend(
            loc="lower right" if view == "2d" else "best",
            facecolor="white",
            prop=dict(size=fontsize * 2),
        )

    if savefig:
        plt.savefig(savefig, dpi=dpi)
        plt.close()
    else:
        plt.show()

    return ax


@compat_alias(
    "get2dtensor",
    since=_BACKWARD_SINCE,
    remove_in=_BACKWARD_REMOVE,
    export=True,
    extra=(f"Use 'tensor2d'. Removal in v{_BACKWARD_REMOVE}."),
)
def tensor2d(  # get2dtensor
    z_or_edis_obj_list: list[EDIO | ZO],
    /,
    tensor: str = "z",
    component: str = "xy",
    kind: str = "modulus",
    return_freqs: bool = False,
    freqs: ArrayLike | None = None,
    **kws,
):
    """
    Build a 2-D matrix (freq × station) from a tensor collection.

    Converts a collection of :class:`Edi` or :class:`Z` objects
    into a 2-D array where rows are frequencies and columns are
    stations. Missing per-site frequencies are filled with ``NaN``
    (no interpolation).

    Parameters
    ----------
    z_or_edis_obj_list : list of Edi or Z
        Collection of EM objects. All items must be :class:`Edi` or
        all :class:`Z`.
    tensor : str, default='z'
        Tensor name or alias. Examples include ``'z'``,
        ``'res'``/``'rho'``/``'rhoa'``, ``'phase'``/``'phs'``.
        Error arrays are supported (e.g., ``'resistivity_err'``).
    component : {'xx','xy','yx','yy'}, default='xy'
        Component to extract for non-frequency tensors.
    kind : {'modulus','real','imag','complex'}, default='modulus'
        Numeric form for complex ``Z`` tensors. Ignored for real
        arrays (e.g., resistivity/phase).
    return_freqs : bool, default=False
        If ``True``, also return the reference frequency vector.
    freqs : array-like, optional
        Precomputed reference frequency grid. If given, it is used
        directly and ``get_full_frequency`` is not called. The grid
        should include (or supersede) each site's frequencies; any
        missing values will be filled with ``NaN`` during alignment.
    **kws
        Extra keywords forwarded to the frequency extractor when
        ``freqs`` is not provided.

    Returns
    -------
    mat2d : ndarray (n_freq, n_stations)
        2-D matrix of the requested tensor component.
    (mat2d, freqs) : tuple
        If ``return_freqs=True``, also returns the frequency vector.

    Raises
    ------
    EMError
        If inputs are missing, mixed, or illegal (e.g., frequency
        requested as the primary tensor here).
    ValueError
        If an unknown ``kind`` is used for complex ``Z``.

    Notes
    -----
    Each item in the input provides a 3-D tensor of shape
    ``(n_freq, 2, 2)``. Index positions map as::

        xx -> (0, 0)   xy -> (0, 1)
        yx -> (1, 0)   yy -> (1, 1)

    Examples
    --------
    >>> phase_yx = tensor2d(data, tensor="phase", component="yx")
    >>> phase_yx.shape
    (56, 7)
    """
    # Normalize and validate request.
    name, comp = parse_tensor(
        tensor=tensor,
        component=component,
        kind=kind,
        **kws,
    )
    if name == "_freq":
        raise EMError(
            "Frequency as the primary tensor is not supported here. "
            "Use a dedicated frequency-to-2D helper instead."
        )
    if z_or_edis_obj_list is None:
        raise EMError(f"Cannot build {name!r} 2D matrix from a missing collection.")

    # Ensure the collection is homogeneous (all EDI or all Z).
    kind_in = check_em_kind(z_or_edis_obj_list)

    # Reference frequency grid: use provided or compute.
    if freqs is None:
        # Lazy import to avoid circular deps at import time.
        freqs = full_freq(z_or_edis_obj_list, **kws)
    else:
        # Validate provided grid and coerce to a safe 1-D view.
        freqs = np.asarray(freqs, dtype=float).ravel()

    # Component index mapping (rows=freqs).
    idx_map = {
        "xx": (slice(None, len(freqs)), 0, 0),
        "xy": (slice(None, len(freqs)), 0, 1),
        "yx": (slice(None, len(freqs)), 1, 0),
        "yy": (slice(None, len(freqs)), 1, 1),
    }

    # Extract raw per-site vectors at 'comp' from 'name'.
    def _extract_vector(obj):
        container = obj.Z if kind_in == "EDI" else obj
        arr = getattr(container, name)
        return arr[idx_map[comp]]

    cols = [_extract_vector(o) for o in z_or_edis_obj_list]

    # Try a fast stack assuming identical frequency grids.
    try:
        mat2d = np.vstack(cols).T
    except Exception:
        # Fallback: align each site to the reference frequency grid.
        def _site_freq(obj):
            container = obj.Z if kind_in == "EDI" else obj
            # Accept both private and public freq attributes.
            f = getattr(container, "_freq", None)
            if f is None:
                f = getattr(container, "freq", None)
            return f

        aligned = [
            align_tensor(
                freqs,
                _site_freq(o),
                v,
                fill_value=np.nan,
            )
            for o, v in zip(z_or_edis_obj_list, cols)
        ]
        mat2d = concat_array_from_list(aligned, concat_axis=1)

    # If complex Z requested, project to the desired numeric form.
    if "z" in name:
        k = str(kind).strip().lower()
        k = {"abs": "modulus", "mod": "modulus"}.get(k, k)
        if k not in {"modulus", "real", "imag", "complex"}:
            raise ValueError(
                "Unknown 'kind' for complex Z. Expected one of "
                "'modulus', 'real', 'imag', or 'complex'."
            )
        proj = {
            "modulus": np.abs(mat2d),
            "real": mat2d.real,
            "imag": mat2d.imag,
            "complex": mat2d,
        }
        mat2d = proj[k]

    return (mat2d, freqs) if return_freqs else mat2d


@compat_alias(
    "qc",
    since=_BACKWARD_SINCE,
    remove_in=_BACKWARD_REMOVE,
    export=True,
    extra=(f"Use 'compute_qc'. Removal in v{_BACKWARD_REMOVE}."),
)
def compute_qc(
    z_or_edis_obj_list: list[EDIO | ZO],
    /,
    tol: float = 0.5,
    *,
    interpolate_freq: bool = False,
    return_freq: bool = False,
    tensor: str = "res",
    return_data: bool = False,
    to_log10: bool = False,
    return_qco: bool = False,
) -> (
    tuple[float]
    | tuple[float, np.ndarray]
    | tuple[float, np.ndarray, np.ndarray]
    | Bunch
):
    """
    Assess data quality across a collection of EDI/Z objects.

    Computes a global completeness ratio
    :math:`1 - \\#NaN / (N_{freq} N_{sta})` from a 2-D tensor
    (freq × station). Frequencies whose per-row missing-data fraction
    exceeds ``tol`` are dropped. Optionally interpolate the retained
    frequencies and/or return a structured summary object.

    Parameters
    ----------
    z_or_edis_obj_list : list of Edi or Z
        Homogeneous collection of EM objects.
    tol : float, default=0.5
        Tolerance threshold in ``[0, 1]``. A frequency row is
        considered **invalid** and dropped when its fraction of
        missing values exceeds ``tol``.
    interpolate_freq : bool, default=False
        If ``True``, interpolate the retained frequencies on a
        log-spaced grid spanning ``[min, max]`` with the same count.
    return_freq : bool, default=False
        If ``True``, return the retained (possibly interpolated)
        frequency vector.
    tensor : {'z','res','rho','rhoa','phase','phs'}, default='res'
        Tensor family used for QC. The function first attempts the
        **TE** component (``'xy'``); on failure, it falls back to
        **TM** (``'yx'``).
    return_data : bool, default=False
        If ``True``, also return the subset of the tensor data
        corresponding to retained frequencies.
    to_log10 : bool, default=False
        If ``True``, return ``log10`` of the retained frequencies.
        Applied after interpolation if ``interpolate_freq`` is set.
    return_qco : bool, default=False
        If ``True``, return a :class:`Bunch` with the following
        attributes:

        - ``rate_``: global completeness ratio
        - ``component_``: selected component, ``'xy'`` or ``'yx'``
        - ``mode_``: EM mode, ``'TE'`` or ``'TM'``
        - ``freqs_``: retained (optionally interpolated) frequencies
        - ``invalid_freqs_``: frequencies dropped by the QC
        - ``data_``: tensor data at retained frequencies

        Setting this flag forces ``return_freq=True`` and
        ``return_data=True``.

    Returns
    -------
    (rate,) or (rate, freqs) or (rate, freqs, data) or Bunch
        Depending on the flags. ``rate`` is in ``[0, 1]``.

    Notes
    -----
    The 2-D tensor has shape ``(n_freq, n_stations)``. The per-row
    missing-data fraction is ``nan_count / n_stations``.

    Examples
    --------
    >>> (rate,) = compute_qc(data)
    >>> rate, freqs = compute_qc(data, return_freq=True)
    >>> rep = compute_qc(data, return_qco=True)
    >>> rep.rate_, rep.component_, rep.freqs_.shape
    (0.75, 'xy', (56,))
    """
    # Validate tolerance.
    tol = assert_ratio(
        tol,
        bounds=(0.0, 1.0),
        exclude_value="use lower bound",
        name="tolerance",
        in_percent=True,
    )

    # Extract 2-D tensor (TE first, then TM on failure).
    tensor = str(tensor).strip().lower()
    try:
        component, mode = "xy", "TE"
        ar, f = tensor2d(
            z_or_edis_obj_list,
            tensor=tensor,
            component=component,
            return_freqs=True,
        )
    except Exception:
        component, mode = "yx", "TM"
        ar, f = tensor2d(
            z_or_edis_obj_list,
            tensor=tensor,
            component=component,
            return_freqs=True,
        )

    # ar: (n_freq, n_stations); f: (n_freq,)
    # Per-row NaN count and fraction (across stations).
    nan_per_row = np.nansum(np.isnan(ar), axis=1)
    frac_nan = np.around(nan_per_row / ar.shape[1], 2)

    # Global completeness ratio:
    #   1 - (#NaN across all cells) / (n_freq * n_stations)
    rate = 1.0 - float(nan_per_row.sum()) / float(ar.size)
    rate = float(np.around(rate, 2))

    # Identify rows to drop (too many NaNs).
    bad_idx = reshape(np.argwhere(frac_nan > tol))
    keep_mask = np.ones(ar.shape[0], dtype=bool)
    keep_mask[bad_idx] = False

    # Retained frequencies and data.
    f_ret = f[keep_mask]
    ar_ret = ar[keep_mask, :]

    # Ensure descending frequency order (common for EM plots).
    if f_ret[0] < f_ret[-1]:
        f_ret = f_ret[::-1]
        ar_ret = ar_ret[::-1, :]

    # Track invalid (dropped) frequencies.
    invalid_freqs = f[~np.isin(f, f_ret)]

    # Optional interpolation on a log-spaced grid.
    if interpolate_freq and f_ret.size > 1:
        f_log = np.logspace(
            np.log10(f_ret.min()),
            np.log10(f_ret.max()),
            f_ret.size,
        )[::-1]
        f_ret = f_log
        # If user wants natural units, undo log10 later.
        if not to_log10:
            f_ret = np.power(10.0, f_ret)
        # We will handle to_log10 below; ensure correct flag flow.
        to_log10 = False

    # Optional log10 transform of retained frequencies.
    if to_log10:
        if np.any(f_ret <= 0.0):
            raise ValueError("Frequencies must be > 0 for log10 transform.")
        f_ret = np.log10(f_ret)

    # Flatten frequency to a 1-D shape.
    f_ret = reshape(f_ret)

    # If any frequency transform happened, ensure we include freqs.
    if interpolate_freq or to_log10:
        return_freq = True

    # If a report is requested, force returning freqs and data.
    if return_qco:
        return_freq = True
        return_data = True

    # Assemble the return payload.
    out: (
        tuple[float]
        | tuple[float, np.ndarray]
        | tuple[float, np.ndarray, np.ndarray]
        | Bunch
    )

    if return_qco:
        out = Bunch(
            tol=tol,
            tensor=tensor,
            component_=component,
            mode_=mode,
            rate_=rate,
            freqs_=f_ret,
            invalid_freqs_=invalid_freqs,
            data_=ar_ret,
        )
    else:
        parts: list = [rate]
        if return_freq:
            parts.append(f_ret)
        if return_data:
            parts.append(ar_ret)
        out = tuple(parts)  # type: ignore[assignment]

    return out


@compat_alias(
    "get_full_frequency",
    since=_BACKWARD_SINCE,
    remove_in=_BACKWARD_REMOVE,
    export=True,
    extra=(f"Use 'full_freq'. Removal in v{_BACKWARD_REMOVE}."),
)
def full_freq(
    z_or_edis_obj_list: list[EDIO | ZO],
    /,
    to_log10: bool = False,
) -> ArrayLike:
    """
    Return the reference (clean) frequency grid for a collection.

    The *full* frequency grid is taken from the site that contains
    the largest number of frequency samples (i.e., the most complete
    set). This is commonly used as the survey reference grid to which
    per-site tensors are aligned.

    Parameters
    ----------
    z_or_edis_obj_list : list of Edi or Z
        Homogeneous collection of :class:`Edi` or :class:`Z` objects.
    to_log10 : bool, default=False
        If ``True``, return ``log10(freqs)``. Frequencies must be
        strictly positive.

    Returns
    -------
    ndarray of shape (n_freq,)
        The reference frequency vector.

    Raises
    ------
    TypeError
        If the input is empty or not iterable.
    EMError
        If the collection mixes Edi and Z, or any element is missing
        a frequency attribute.
    ValueError
        If ``to_log10=True`` and any frequency is non-positive.

    Notes
    -----
    For each object, the function looks for ``.Z._freq`` / ``.Z.freq``
    (when :class:`Edi`) or ``._freq`` / ``.freq`` (when :class:`Z`),
    using the first available attribute.

    Examples
    --------
    >>> f = full_freq(edi_data)
    >>> f.shape
    (56,)
    >>> flog = full_freq(edi_data, to_log10=True)
    """
    # Ensure homogeneous collection (all EDI or all Z).
    kind = check_em_kind(z_or_edis_obj_list)

    def _get_freq(obj) -> np.ndarray:
        """Extract frequency array from Edi/Z with robust attribute
        lookup; raise meaningful error if missing."""
        container = obj.Z if kind == "EDI" else obj
        f = getattr(container, "_freq", None)
        if f is None:
            f = getattr(container, "freq", None)
        if f is None:
            raise EMError(
                "Frequency attribute not found on object. Expected "
                "'.Z._freq' or '.Z.freq' for Edi, or '._freq'/'freq' "
                "for Z."
            )
        return np.asarray(f, dtype=float)

    # Collect all frequency arrays and pick the longest as reference.
    freq_list = [_get_freq(o) for o in z_or_edis_obj_list]
    lengths = np.array([len(f) for f in freq_list], dtype=int)
    if lengths.size == 0:
        raise TypeError("Input collection is empty.")
    idx = int(np.argmax(lengths))
    f_ref = freq_list[idx]

    if to_log10:
        if np.any(f_ref <= 0.0):
            raise ValueError(
                "Frequencies must be strictly positive to compute log10()."
            )
        f_ref = np.log10(f_ref)

    return f_ref


# XXX TODO:


@compat_alias(
    "plot_tensors2",
    since=_BACKWARD_SINCE,
    remove_in=_BACKWARD_REMOVE,
    export=True,
    extra=(f"Use 'plot_station_tensors'. Removal in v{_BACKWARD_REMOVE}."),
)
def plot_station_tensors(  # plot_tensors2
    z_or_edis_obj_list: list[EDIO | ZO],
    /,
    station: int | str = "S00",
    *,
    plot_z: bool = False,
    show_error_bars: bool = True,
    **kwargs: Any,
) -> ZO:
    """
    Plot tensors for one station: resistivity/phase or Z real/imag.

    By default, the function plots apparent resistivity and phase
    panels for the four components (xx, xy, yx, yy). If ``plot_z`` is
    ``True``, it plots the real and imaginary parts of the impedance
    tensor instead. Error bars can be displayed or hidden.

    Parameters
    ----------
    z_or_edis_obj_list : list of Edi or Z
        Collection containing either :class:`Edi` or :class:`Z`
        objects. When EDI objects are provided, the embedded ``Z``
        is extracted.
    station : int or str, default='S00'
        Target station index or label. Strings like ``'S00'`` are
        parsed to an integer index (0-based).
    plot_z : bool, default=False
        If ``True``, plot real/imag parts of ``Z``. Otherwise, plot
        apparent resistivity and phase.
    show_error_bars : bool, default=True
        Whether to include error bars on each panel.
    **kwargs
        Plot customization such as:
        - ``fig_size`` (tuple, default=(6, 6))
        - ``dpi`` (int, default=300)
        - ``subplot_wspace`` (float, default=0.3)
        - ``phase_limits`` (tuple[min,max] in deg)
        - ``freq_limits`` (tuple[min,max] in Hz)
        - ``period_limits`` (tuple[min,max] in s)
        - ``mod_base`` (int, default=360)
        - style: ``color_mode`` ('color'|'bw'), markers, line widths,
          legend style, font sizes, etc.

    Returns
    -------
    Z
        The :class:`Z` object for the selected station.

    Notes
    -----
    The function expects the station's :class:`Z` object to provide
    (or compute on demand) the following arrays with shape
    ``(n_freq, 2, 2)``:
    - ``resistivity``, ``resistivity_err``
    - ``phase``, ``phase_err``
    - ``z`` (complex), ``z_err`` (complex)
    and a frequency vector ``_freq`` (Hz).

    Examples
    --------
    >>> z = plot_station_tensors(edi_list, station=3)
    >>> z = plot_station_tensors(
    ...     edi_list,
    ...     station="S00",
    ...     plot_z=True,
    ...     show_error_bars=False,
    ...     color_mode="bw",
    ... )
    """
    # Some code paths historically used `zplot=...`; accept alias.
    if "zplot" in kwargs and not kwargs.get("plot_z_set", False):
        plot_z = bool(kwargs.pop("zplot"))
        kwargs["plot_z_set"] = True  # prevent repeated aliasing

    # 1) Resolve station index and extract the Z container.
    idx = _parse_station_index(station)
    obj_type, z_obj = _select_station_obj(z_or_edis_obj_list, idx)

    # 2) Build the 2×4 panel layout (xx, xy, yx, yy × resistivity/phase
    #    or real/imag). Returns (fig, [axes0..axes7]).
    fig, axes = _init_axes(**kwargs)

    # 3) Prepare plot inputs (res/phase vs z.real/z.imag + errors).
    res, res_err, phs, phs_err = _prepare_station_data(z_obj, plot_z)

    # 4) Apply frequency/period and phase filtering.
    freq, res, res_err, phs, phs_err = _apply_filters(
        z_obj._freq,
        res,
        res_err,
        phs,
        phs_err,
        **kwargs,
    )

    # 5) Draw panels and decorate.
    _draw_station_panels(
        fig,
        axes,
        freq,
        res,
        res_err,
        phs,
        phs_err,
        show_error_bars=show_error_bars,
        plot_z=plot_z,
        **kwargs,
    )

    plt.show()
    return z_obj


def _parse_station_index(station: int | str) -> int:
    """
    Extract the 0-based station index from ``station``.

    Accepts either an integer or a string like ``'S00'``. Raises a
    clear error when the string does not contain digits.
    """
    if isinstance(station, int):
        return station
    m = re.search(r"\\d+", str(station), flags=re.IGNORECASE)
    if m is None:
        raise TypeError(
            "Station should be an integer or include a position " "number, e.g., 'S00'."
        )
    return int(m.group())


def _select_station_obj(
    z_or_edis_obj_list: list[EDIO | ZO],
    idx: int,
) -> tuple[str, ZO]:
    """
    Determine object type (EDI or Z) and return the station's Z.
    """
    if idx < 0 or idx >= len(z_or_edis_obj_list):
        raise ValueError(
            "Station index out of range. Only "
            f"{len(z_or_edis_obj_list)} stations available."
        )

    item = z_or_edis_obj_list[idx]
    obj_type = "EDI" if hasattr(item, "Z") else "Z"
    z_obj = item.Z if obj_type == "EDI" else item
    return obj_type, z_obj


def _init_axes(**kwargs: Any) -> tuple[plt.Figure, list[plt.Axes]]:
    """
    Initialize the 2×4 panel layout and return (figure, axes).
    """
    fig_size = kwargs.pop("fig_size", (6.0, 6.0))
    fig_dpi = kwargs.pop("dpi", 300)

    fig = plt.figure(figsize=fig_size, dpi=fig_dpi)
    plt.clf()

    # 2 rows (res/phase or real/imag) × 4 comps (xx, xy, yx, yy).
    gs = gridspec.GridSpec(2, 4, wspace=kwargs.get("subplot_wspace", 0.3))
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(4)]
    return fig, axes


def _prepare_station_data(
    z_obj: ZO,
    plot_z: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Collect arrays for plotting at the station, with errors.
    """
    # Ensure derived quantities exist (resistivity/phase + errors).
    # If `compute_resistivity_phase` is a no-op when present, this
    # still keeps the code safe.
    if hasattr(z_obj, "compute_resistivity_phase"):
        z_obj.compute_resistivity_phase()

    if plot_z:
        # Real/imag parts of Z with corresponding errors.
        res = np.abs(z_obj.z.real)  # use absolute magnitude
        res_err = np.abs(z_obj.z_err.real)
        phs = np.abs(z_obj.z.imag)
        phs_err = np.abs(z_obj.z_err.imag)
    else:
        # Apparent resistivity and phase with errors.
        res = z_obj.resistivity
        res_err = z_obj.resistivity_err
        phs = z_obj.phase
        phs_err = z_obj.phase_err

    return res, res_err, phs, phs_err


def _apply_filters(
    freq: np.ndarray,
    res: np.ndarray,
    res_err: np.ndarray,
    phs: np.ndarray,
    phs_err: np.ndarray,
    **kwargs: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Apply frequency/period limits and (optionally) phase wrapping.
    """
    phase_limits = kwargs.pop("phase_limits", None)
    period_limits = kwargs.pop("period_limits", None)
    freq_limits = kwargs.pop("freq_limits", None)
    mod_base = kwargs.pop("mod_base", 360)

    # If both are provided, prefer freq limits.
    if freq_limits is not None and period_limits is not None:
        print(
            "Both 'freq_limits' and 'period_limits' provided. " "Using 'freq_limits'."
        )
        period_limits = None

    # Convert period limits to frequency if needed.
    if period_limits is not None:
        if (
            not isinstance(period_limits, (list, tuple))
            or len(period_limits) != 2
            or not all(isinstance(x, (int, float)) for x in period_limits)
        ):
            raise ValueError("period_limits must be (min_period, max_period).")
        p_min, p_max = sorted(period_limits)
        # f = 1 / T
        freq_limits = (1.0 / p_max, 1.0 / p_min)

    # Build mask from frequency limits.
    if freq_limits is not None:
        if (
            not isinstance(freq_limits, (list, tuple))
            or len(freq_limits) != 2
            or not all(isinstance(x, (int, float)) for x in freq_limits)
        ):
            raise ValueError("freq_limits must be (min_freq, max_freq).")
        f_min, f_max = sorted(freq_limits)
        m = (freq >= f_min) & (freq <= f_max)
        idx = np.where(m)[0]
        f_new = freq[m]
    else:
        idx = np.arange(len(freq))
        f_new = freq

    # Phase wrapping / clipping if requested.
    if phase_limits is not None:
        phs_new = wrap_phase(  # adjusted phase range
            phs[idx, :, :],
            value_range=phase_limits,
            mod_base=mod_base,
        )
    else:
        phs_new = phs[idx, :, :]

    # Slice data consistently.
    res_new = res[idx, :, :]
    res_err_new = res_err[idx, :, :]
    phs_err_new = phs_err[idx, :, :]

    return f_new, res_new, res_err_new, phs_new, phs_err_new


def _draw_station_panels(
    fig: plt.Figure,
    axes: list[plt.Axes],
    freq: np.ndarray,
    res: np.ndarray,
    res_err: np.ndarray,
    phs: np.ndarray,
    phs_err: np.ndarray,
    *,
    show_error_bars: bool,
    plot_z: bool,
    **kwargs: Any,
) -> None:
    """
    Render the 2×4 panels with styling and legends.
    """

    # Style defaults
    ms = kwargs.pop("ms", 1.5)
    ms_r = kwargs.pop("ms_r", 3.0)  # noqa
    lw = kwargs.pop("lw", 0.5)
    lw_r = kwargs.pop("lw_r", 1.0)  # noqa
    e_capthick = kwargs.pop("e_capthick", 0.5)
    e_capsize = kwargs.pop("e_capsize", 2.0)
    color_mode = kwargs.pop("color_mode", "color")
    leg_style = kwargs.pop("leg_style", "2")
    font_size = kwargs.pop("font_size", 6)
    tick_label_size = kwargs.pop("tick_label_size", 8)

    plt.rcParams["font.size"] = font_size
    fontdict = {"size": font_size + 2, "weight": "bold"}  # noqa

    # Color/marker scheme
    if color_mode == "color":
        cted = kwargs.pop("cted", (0, 0, 1))
        ctmd = kwargs.pop("ctmd", (1, 0, 0))
        mted = kwargs.pop("mted", "s")
        mtmd = kwargs.pop("mtmd", "o")
    elif color_mode == "bw":
        cted = kwargs.pop("cted", (0, 0, 0))
        ctmd = kwargs.pop("ctmd", (0, 0, 0))
        mted = kwargs.pop("mted", "s")
        mtmd = kwargs.pop("mtmd", "o")
    else:
        raise ValueError("color_mode must be 'color' or 'bw'.")

    # Map components to panel indices and styles.
    comp_info = {
        "xx": {"idx": (0, 0), "color": cted, "marker": mted},
        "xy": {"idx": (0, 1), "color": cted, "marker": mted},
        "yx": {"idx": (1, 0), "color": ctmd, "marker": mtmd},
        "yy": {"idx": (1, 1), "color": ctmd, "marker": mtmd},
    }

    #  # --> make key word dictionaries for plotting
    # kw_xx = {'color': cted,
    #          'marker': mted,
    #          'ms': ms,
    #          'ls': ':',
    #          'lw': lw,
    #          'e_capsize': e_capsize,
    #          'e_capthick': e_capthick}

    # kw_yy = {'color': ctmd,
    #          'marker': mtmd,
    #          'ms': ms,
    #          'ls': ':',
    #          'lw': lw,
    #          'e_capsize': e_capsize,
    #          'e_capthick': e_capthick}

    # period = 1 / freq  # Convert frequency to period

    # # Correctly mapping components to their indices in the
    # # 2x2 matrix and respective colors
    # components_info = {
    #     'xx': {'index': (0, 0), 'kw': 'kw_xx'},
    #      # xx and xy share the same color and marker style
    #     'xy': {'index': (0, 1), 'kw': 'kw_xx'},
    #     'yx': {'index': (1, 0), 'kw': 'kw_yy'},
    #      # yx and yy share the same color and marker style
    #     'yy': {'index': (1, 1), 'kw': 'kw_yy'}
    # }
    # res_labels = []
    # phase_labels =[]
    # res_leg_objs =[]
    # phase_leg_objs=[]
    # for i, component in enumerate(components_info):
    #     ax_res = ax_list[i]  # Axes for resistivity plots
    #     ax_phase = ax_list[i + 4]  # Axes for phase plots
    #     index = components_info[component]['index']
    #     kw_arg = components_info[component]['kw']
    #     # Accessing the 3D array for each component
    #     res = plot_res[:, index[0], index[1]]
    #     res_err = plot_res_err[:, index[0], index[1]]
    #     phase = plot_phase[:, index[0], index[1]]
    #     phase_err = plot_phase_err[:, index[0], index[1]]

    #     # Keyword arguments for plotting, differentiating
    #     # between xx/xy and yx/yy components
    #     # Dynamically access the appropriate kwargs based on the component
    #     plot_kwargs = locals()[kw_arg]

    #     # Plot resistivity
    #     er_res=plot_errorbar(
    #         ax_res, period, res, res_err if show_error_bars else None,
    #         **plot_kwargs
    #     )
    #     # Plot phase
    #     er_phase=plot_errorbar(
    #         ax_phase, period, phase, phase_err if show_error_bars else None,
    #         **plot_kwargs
    #     )

    # X axis in period (s) for both rows.
    period = 1.0 / np.asarray(freq, dtype=float)

    res_labels: list[str] = []
    phs_labels: list[str] = []
    res_legs = []
    phs_legs = []

    # Iterate in the fixed order xx, xy, yx, yy.
    for i, comp in enumerate(comp_info):
        ax_r = axes[i]  # top row
        ax_p = axes[i + 4]  # bottom row
        r, c = comp_info[comp]["idx"]
        color = comp_info[comp]["color"]
        marker = comp_info[comp]["marker"]

        # Slice 3-D arrays into 1-D series for the component.
        series_res = res[:, r, c]
        series_res_err = res_err[:, r, c]
        series_phs = phs[:, r, c]
        series_phs_err = phs_err[:, r, c]

        # Common kwargs for errorbar plotting.
        eb_kw = dict(
            color=color,
            marker=marker,
            ms=ms,
            ls=":",
            lw=lw,
            e_capsize=e_capsize,
            e_capthick=e_capthick,
        )

        # Resistivity (or Re[Z]) panel.
        leg_r = plot_errorbar(
            ax_r,
            period,
            series_res,
            series_res_err if show_error_bars else None,
            **eb_kw,
        )
        # Phase (or Im[Z]) panel.
        leg_p = plot_errorbar(
            ax_p,
            period,
            series_phs,
            series_phs_err if show_error_bars else None,
            **eb_kw,
        )

        # Axes scaling and grids.
        ax_r.set_xscale("log")
        ax_r.set_yscale("log")
        ax_p.set_xscale("log")
        ax_r.grid(True, which="both", ls="--", lw=0.5, color="gray", alpha=0.5)
        ax_p.grid(True, which="both", ls="--", lw=0.5, color="gray", alpha=0.5)

        # Data-driven y limits (guard against non-finite).
        def _finite_min_max(a: np.ndarray) -> tuple[float, float]:
            a = a[np.isfinite(a)]
            if a.size == 0:
                return (1e-3, 1.0)
            return (float(a.min()), float(a.max()))

        rmin, rmax = _finite_min_max(series_res)
        pmin, pmax = _finite_min_max(series_phs)

        ax_r.set_ylim([max(rmin / 2.0, 1e-12), rmax * 2.0])
        ax_p.set_ylim([pmin - 5.0, pmax + 5.0])

        # Tick sizes
        ax_r.tick_params(
            axis="both",
            which="major",
            labelsize=kwargs.get("tick_label_size", tick_label_size),
        )
        ax_p.tick_params(
            axis="both",
            which="major",
            labelsize=kwargs.get("tick_label_size", tick_label_size),
        )

        if i == 0:
            ax_r.set_ylabel(
                "Re[Z] (mV/km·nT)" if plot_z else r"App. Res. ($\Omega\cdot m$)",
                fontsize=kwargs.get("font_size", font_size),
            )
            ax_p.set_ylabel(
                "Im[Z] (mV/km·nT)" if plot_z else "Phase (deg.)",
                fontsize=kwargs.get("font_size", font_size),
            )

        # Legend labels with LaTeX formatting.
        res_label = rf"$z_{{{comp}}}$" if plot_z else rf"$\rho_{{{comp}}}$"
        phs_label = rf"$\phi_{{{comp}}}$"
        res_labels.append(res_label)
        phs_labels.append(phs_label)

        res_legs.append(leg_r)
        phs_legs.append(leg_p)

    # Top row (res/Re[Z]) hides x-labels; bottom row shows periods.
    for ax in axes[:4]:
        ax.set_xlabel("")
        ax.tick_params(
            axis="x",
            which="both",
            bottom=False,
            top=False,
            labelbottom=False,
        )
    for ax in axes[4:]:
        ax.set_xlabel(
            "Period (s)",
            fontsize=kwargs.get("font_size", font_size),
        )

    # Legends: style '2' puts a compact legend above each top panel.
    if str(leg_style).lower() == "2":
        for k, ax in enumerate(axes[:4]):
            ax.legend(
                [res_legs[k]],
                [res_labels[k]],
                loc="upper center",
                bbox_to_anchor=(0.5, 1.18),
                markerscale=1.0,
                borderaxespad=0.01,
                labelspacing=0.07,
                handletextpad=0.2,
                borderpad=0.15,
                framealpha=1.0,
                prop={"size": max(font_size, 5)},
            )
    else:
        for i, ax in enumerate(axes[:4]):
            ax.legend(
                [res_labels[i]],
                loc="best",
                fontsize=kwargs.get("font_size", 8),
                frameon=True,
                edgecolor="black",
            )
        for i, ax in enumerate(axes[4:]):
            ax.legend(
                [phs_labels[i]],
                loc="best",
                fontsize=kwargs.get("font_size", 8),
                frameon=True,
                edgecolor="black",
            )

    fig.subplots_adjust(hspace=0.1, wspace=0.3)


def wrap_phase(
    phase: np.ndarray | float | int,
    value_range: tuple[float, float] | list | float | int | None = None,
    mod_base: int = 360,
) -> np.ndarray:
    """
    Wrap phase values to a target range with a given periodic base.

    By default, phases are wrapped into the interval
    ``[0, mod_base)``. When ``value_range`` is provided, the wrapped
    phases are **linearly remapped** from ``[0, mod_base)`` to the
    desired interval.

    Parameters
    ----------
    phase : array-like
        Phase values (any shape), possibly negative or outside the
        target range.
    value_range : {None, scalar, (min, max)}, optional
        Target interval for the output.
        - ``None``: return values in ``[0, mod_base)``.
        - scalar: treated as ``(0, scalar)``.
        - tuple/list ``(min, max)``: custom interval; ``min < max``.
    mod_base : {90, 180, 270, 360}, default=360
        Periodicity used for wrapping (degrees).

    Returns
    -------
    np.ndarray
        Wrapped (and optionally remapped) phase values with the same
        shape as the input and dtype float.

    Notes
    -----
    - For the common symmetric range ``(-180, 180]``, set
      ``value_range=(-180, 180)`` with ``mod_base=360``.
    - The remapping is affine: values in ``[0, mod_base)`` are
      scaled to ``[min, max)``.

    Examples
    --------
    >>> x = np.array([-540, -180, 0, 180, 360, 540])
    >>> # Default: [0, 360)
    >>> wrap_phase(x, mod_base=360)
    array([180., 180.,   0., 180.,   0., 180.])

    >>> # Symmetric range (-180, 180)
    >>> wrap_phase(x, value_range=(-180, 180), mod_base=360)
    array([-180., -180.,    0.,  180.,    0.,  180.])

    >>> # Custom half-range [0, 180)
    >>> wrap_phase(x, value_range=180, mod_base=360)
    array([90., 90.,  0., 90.,  0., 90.])
    """
    # ---- validate inputs -------------------------------------------------
    arr = np.asarray(phase, dtype=float)

    if mod_base not in (90, 180, 270, 360):
        raise ValueError("mod_base must be one of {90, 180, 270, 360}.")

    # ---- wrap to [0, mod_base) ------------------------------------------
    # numpy's remainder yields values in [0, mod_base) for positive base.
    wrapped = np.remainder(arr, float(mod_base))

    # ---- optional linear remap to a custom interval ----------------------
    if value_range is None:
        return wrapped

    # Scalar => (0, scalar)
    if np.isscalar(value_range):
        vmin, vmax = 0.0, float(value_range)
    else:
        if not isinstance(value_range, (tuple, list)) or len(value_range) != 2:
            raise ValueError(
                "value_range must be None, a scalar, or a tuple/list "
                "of two values (min, max)."
            )
        vmin, vmax = float(value_range[0]), float(value_range[1])

    # Ensure vmin < vmax
    if not (vmin < vmax):
        # sort, preserving intent if reversed order was given
        vmin, vmax = min(vmin, vmax), max(vmin, vmax)

    # Affine remap: [0, mod_base) -> [vmin, vmax)
    scale = (vmax - vmin) / float(mod_base)
    out = wrapped * scale + vmin

    return out


def plot_tensors(
    z_or_edis_obj_list: list[EDIO | ZO],
    /,
    station: int | str = "S00",
    zplot: bool = False,
    show_error_bars: bool = False,
    **kwargs: Any,
) -> ZO:
    """
    Plot tensors for one station (compat wrapper).

    This is a compatibility wrapper around
    :func:`plot_station_tensors`. It preserves the legacy
    API—``station``, ``zplot`` (for impedance vs. app. resistivity/
    phase), and ``show_error_bars``—and forwards any additional
    styling options to the underlying plotter.

    Parameters
    ----------
    z_or_edis_obj_list : list of Edi or Z
        Collection of EM objects containing either :class:`Edi`
        (from which the embedded ``Z`` is extracted) or :class:`Z`
        directly.
    station : int or str, default='S00'
        Target station index or label. Strings such as ``'S00'`` are
        parsed to 0-based indices.
    zplot : bool, default=False
        If ``True``, plot real/imag parts of the impedance tensor
        (``Z``). If ``False``, plot apparent resistivity and phase.
    show_error_bars : bool, default=False
        Whether to display error bars.
    **kwargs
        Additional plotting options passed to
        :func:`plot_station_tensors` (e.g., ``color_mode``,
        markers, line widths, legend style, font sizes, etc.).

    Returns
    -------
    Z
        The :class:`Z` object for the selected station.

    Notes
    -----
    The heavy lifting (layout, filtering, styling) is handled by
    :func:`plot_station_tensors`. This wrapper exists to maintain
    source compatibility with v1.x code.

    Examples
    --------
    >>> z = plot_tensors(edi_list, station="S03", zplot=True)
    >>> z = plot_tensors(edi_list, station=0, show_error_bars=True)
    """
    # Reuse the modern implementation; keep legacy param names.
    return plot_station_tensors(
        z_or_edis_obj_list,
        station=station,
        plot_z=zplot,
        show_error_bars=show_error_bars,
        **kwargs,
    )


def plot_strike(
    list_of_edis: str | Iterable[str],
    /,
    kind: int = 2,
    period_tolerance: float = 0.05,
    text_pad: float = 1.65,
    rot_z: float = 0.0,
    **kws,
) -> None:
    """
    Plot strike angles from invariants and phase tensor as rose/polar
    diagrams.

    Accepts a single *.edi* file path, a directory containing *.edi*
    files, or an iterable of *.edi* paths. Files are validated before
    plotting. Output is produced by :class:`mtpy.imaging.plotstrike.
    PlotStrike`, with console output muted.

    Parameters
    ----------
    list_of_edis : str or iterable of str
        Path to an *.edi* file, a directory of *.edi* files, or a
        list/tuple of *.edi* file paths.
    kind : {1, 2}, default=2
        Plot style for :class:`PlotStrike`:
        - ``1``: plot individual decades in one plot.
        - ``2``: plot all period ranges in a polar diagram for each
          strike estimate.
    period_tolerance : float, default=0.05
        Tolerance to match periods across different EDI files.
    text_pad : float, default=1.65
        Padding of the angle label at the bottom of each polar
        diagram.
    rot_z : float, default=0.0
        Clockwise rotation (degrees) applied to the tensor.
    **kws
        Extra keyword arguments forwarded to :class:`PlotStrike`
        (e.g., ``plot_range``, ``plot_tipper``, ``fold``,
        ``plot_orientation``, color settings, etc.).

    Returns
    -------
    None
        Plots are created as a side effect.

    Notes
    -----
    - Files are validated with :class:`IsEdi._assert_edi`.
    - Third-party output is muted via ``nullify_output()``.

    Examples
    --------
    >>> plot_strike("/path/to/edis_dir")
    >>> plot_strike("/path/to/site.edi", kind=1)
    >>> plot_strike(["a.edi", "b.edi"], rot_z=10.0)
    """
    # Lazy imports to avoid hard deps at module import.
    # XXX TODO: Replace with version v-2
    from mtpy.imaging.plotstrike import PlotStrike

    # ---- normalize input into a list of .edi paths ---------------------
    paths: list[str]
    if isinstance(list_of_edis, str):
        p = list_of_edis
        if os.path.isdir(p):
            paths = [
                os.path.join(p, f)
                for f in os.listdir(p)
                if str(f).lower().endswith(".edi")
            ]
        elif os.path.isfile(p):
            paths = [p]
        else:
            raise FileNotFoundError(f"Path does not exist: {p!r}")
    else:
        # Iterable of paths (e.g., list/tuple/generator)
        paths = [str(x) for x in list_of_edis]

    if len(paths) == 0:
        raise ValueError(
            "No .edi files found to plot. Provide a file path, a "
            "directory containing .edi files, or a list of .edi paths."
        )

    # ---- validate each file is a proper EDI ----------------------------
    for f in paths:
        IsEdi._assert_edi(f)

    # ---- ensure PlotStrike receives our key options --------------------
    # (Do not override user values if already provided.)
    kws.setdefault("period_tolerance", period_tolerance)
    kws.setdefault("text_pad", text_pad)
    kws.setdefault("rot_z", rot_z)

    # ---- run the third-party plotter quietly ---------------------------
    with nullify_output():
        PlotStrike(
            fn_list=paths,
            plot_type=kind,
            **kws,
        )


@compat_alias(
    "plot_l_curve",
    since=_BACKWARD_SINCE,
    remove_in=_BACKWARD_REMOVE,
    export=True,
    extra=(f"Use 'plot_lcurve'. Removal in v{_BACKWARD_REMOVE}."),
)
def plot_lcurve(
    rms: Iterable[float | int],
    roughness: Iterable[float | int],
    tau: Iterable[float | int] | None = None,
    hansen_point: str | tuple[float, float] | None = None,
    rms_target: float | int | None = None,
    view_tline: bool = False,
    hpoint_kws: dict | None = None,
    fig_size: tuple[float, float] = (10.0, 4.0),
    ax: plt.Axes | None = None,
    fig: plt.Figure | None = None,
    style: str | None = "classic",
    savefig: str | None = None,
    **plot_kws: Any,
) -> plt.Axes:
    """
    Plot the Hansen L-curve (RMS vs. roughness) with annotations.

    The L-curve criterion helps select a suitable model after running
    several inversions with different :math:`\\tau` values. This
    function plots RMS against model roughness, optionally highlights
    the *Hansen knee point*, displays the RMS target line, and labels
    each point with its :math:`\\tau`.

    Parameters
    ----------
    rms : array-like
        Sequence of RMS values, one per inversion.
    roughness : array-like
        Sequence of roughness values matching ``rms``.
    tau : array-like, optional
        Sequence of :math:`\\tau` values to annotate at each point.
        Length must match ``rms`` and ``roughness``.
    hansen_point : {'auto', (x, y)}, optional
        If ``'auto'``, the knee point is detected automatically.
        Otherwise, pass a 2-tuple ``(roughness, rms)`` to highlight.
    rms_target : float, optional
        Target RMS. If provided and ``view_tline=True``, a horizontal
        line is drawn at this value. If provided and
        ``view_tline=False``, the y-limits are expanded around it.
    view_tline : bool, default=False
        Whether to draw the target RMS horizontal line.
    hpoint_kws : dict, optional
        Matplotlib kwargs to style the Hansen point marker.
    fig_size : tuple of float, default=(10.0, 4.0)
        Figure size used when creating a new figure.
    ax : matplotlib.axes.Axes, optional
        Target axes. If ``None``, a new figure/axes is created.
    fig : matplotlib.figure.Figure, optional
        Figure handle used for saving when ``savefig`` is set. If
        omitted, the figure from ``ax`` is used or a new one created.
    style : str, optional
        Matplotlib style to use within a context (not global).
    savefig : str, optional
        Path to save the figure. If omitted, the plot is shown.
    **plot_kws
        Extra style passed to ``Axes.plot`` for the L-curve line.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the plotted L-curve.

    Notes
    -----
    Knee detection uses a simple triangle-area curvature surrogate.
    See Hansen & O'Leary (1993) for background on the L-curve.

    Examples
    --------
    >>> rough = [0, 50, 100, 150, 200, 250, 300, 350]
    >>> rmse = [3.16, 3.12, 3.10, 3.08, 3.06, 3.04, 3.02, 3.00]
    >>> plot_lcurve(rmse, rough, hansen_point="auto")
    """
    # Normalize inputs to float arrays (1-D).
    rms_arr = np.asarray(
        is_iterable(rms, exclude_string=True, transform=True),
        dtype=float,
    )
    rough_arr = np.asarray(
        is_iterable(roughness, exclude_string=True, transform=True),
        dtype=float,
    )
    if rms_arr.shape != rough_arr.shape:
        raise ValueError(
            "Shapes of 'rms' and 'roughness' must match. "
            f"Got {rms_arr.shape} vs {rough_arr.shape}."
        )

    # Plot style context (does not change global state).
    _ctx = plt.style.context(style) if style is not None else None
    if _ctx is not None:
        _ctx.__enter__()

    try:
        # Create fig/ax if needed.
        created = False
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=fig_size)
            created = True
        if fig is None:
            fig = ax.figure

        # Defaults for the L-curve line, merged with user kwargs.
        plot_kws = _merge_plot_kws(
            plot_kws,
            {"marker": "o", "linestyle": "-", "color": "black"},
        )
        ax.plot(rough_arr, rms_arr, **plot_kws)

        # Resolve Hansen knee point.
        if isinstance(hansen_point, str) and hansen_point.lower() == "auto":
            hansen_point = _hansen_knee(rough_arr, rms_arr)

        if hansen_point is not None:
            if not (isinstance(hansen_point, (tuple, list)) and len(hansen_point) == 2):
                raise ValueError(
                    "Hansen knee point must be a tuple '(roughness, rms)'."
                )
            hx, hy = float(hansen_point[0]), float(hansen_point[1])
            hpoint_kws = _merge_plot_kws(hpoint_kws, {"marker": "o", "color": "red"})
            ax.plot(hx, hy, **hpoint_kws)
            ax.annotate(
                f"{hx:g}",
                (hx, hy),
                textcoords="offset points",
                xytext=(0, 10),
                ha="center",
            )

        # Annotate tau values at each point (skip the highlighted knee).
        if tau is not None:
            tau_arr = np.asarray(is_iterable(tau, exclude_string=True, transform=True))
            if tau_arr.shape[0] != rms_arr.shape[0]:
                raise ValueError(
                    "'tau' must have the same length as 'rms' and "
                    f"'roughness'. Got {tau_arr.shape[0]} vs {rms_arr.shape[0]}."
                )
            # Column-stack roughness and rms into point coordinates.
            pts = np.column_stack([rough_arr, rms_arr])
            for (x, y), tval in zip(pts, tau_arr):
                if hansen_point is not None and np.allclose([x, y], hansen_point):
                    continue
                ax.annotate(
                    str(tval),
                    (float(x), float(y)),
                    textcoords="offset points",
                    xytext=(0, 10),
                    ha="center",
                )

        # Target RMS line and/or y-limits centered on target.
        if rms_target is not None:
            try:
                from .validation import (
                    _assert_all_types as _aat,
                )

                rms_target_val = float(
                    _aat(rms_target, float, int, objname="RMS target")
                )
            except Exception:
                rms_target_val = float(rms_target)

            if view_tline:
                ax.axhline(y=rms_target_val, color="k", linestyle=":")

            # Expand y-limits symmetrically around target.
            extent = abs(rms_target_val - float(np.nanmin(rms_arr)))
            ymin = (rms_target_val - extent) if view_tline else 0.0
            ymax = float(np.nanmax(rms_arr)) + extent
            ax.set_ylim([ymin, ymax])

        # Labels and title.
        ax.set_xlabel("Roughness")
        ax.set_ylabel("RMS")
        ax.set_title("RMS vs. Roughness")

        # Save or show.
        if savefig:
            try:
                from .plotutils import savefigure

                savefigure(fig, savefig, dpi=300)
            except Exception:
                fig.savefig(savefig, dpi=300, bbox_inches="tight")
            plt.close(fig)
        else:
            if created:
                plt.show()

        return ax

    finally:
        if _ctx is not None:
            _ctx.__exit__(None, None, None)


def _hansen_knee(
    roughness: np.ndarray,
    rms: np.ndarray,
) -> tuple[float, float]:
    """
    Estimate the Hansen knee point using a triangle-area surrogate.

    Parameters
    ----------
    roughness : ndarray, shape (N,)
        Roughness vector.
    rms : ndarray, shape (N,)
        RMS vector aligned with ``roughness``.

    Returns
    -------
    (x, y) : tuple of float
        Coordinates of the detected knee point.

    Notes
    -----
    This simple curvature proxy uses the area of successive triangles
    along the polyline. It is fast and works well for monotonic
    L-curves. For noisy curves, consider smoothing first.
    """
    n = roughness.size
    if n < 3:
        # Not enough points to form triangles; fall back to best RMS.
        i_best = int(np.nanargmin(rms))
        return float(roughness[i_best]), float(rms[i_best])

    # Triangle areas via Heron's formula, sliding window of 3 points.
    areas = np.zeros(n - 2, dtype=float)
    for i in range(1, n - 1):
        x1, y1 = roughness[i - 1], rms[i - 1]
        x2, y2 = roughness[i], rms[i]
        x3, y3 = roughness[i + 1], rms[i + 1]

        a = np.hypot(x2 - x1, y2 - y1)
        b = np.hypot(x3 - x2, y3 - y2)
        c = np.hypot(x3 - x1, y3 - y1)
        s = 0.5 * (a + b + c)

        # Guard for degenerate triangles.
        if a == 0 or b == 0 or c == 0:
            area = 0.0
        else:
            # Heron's formula.
            area = max(s * (s - a) * (s - b) * (s - c), 0.0) ** 0.5
        areas[i - 1] = area

    k = int(np.argmax(areas)) + 1
    return float(roughness[k]), float(rms[k])


def _merge_plot_kws(
    user_kws: dict | None,
    defaults: dict,
) -> dict:
    """
    Merge user plotting kwargs with defaults without mutating inputs.
    """
    out = dict(defaults or {})
    if user_kws:
        out.update(user_kws)
    return out


# ============  Backward-compatibility alias(v1.xxx) ==========================

get_full_frequency = full_freq
qc = compute_qc
plot_tensors2 = plot_station_tensors
plot_l_curve = plot_lcurve
get2dtensor = tensor2d
