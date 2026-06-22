"""High-level TEMAVG survey conversion workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..api.property import PyCSAMTObject

from ._base import TEMSounding
from .survey import TEMSurvey, read_temavg_survey

PathLike = str | Path

__all__ = [
    "TEMAVGConversion",
    "read_temavg_soundings",
    "transform_temavg_survey",
]


@dataclass(repr=False)
class TEMAVGConversion(PyCSAMTObject):
    """Result bundle returned by TEMAVG workflow helpers.

    Parameters
    ----------
    survey : TEMSurvey
        Parsed survey folder, including companion files and
        optional coordinate table.
    soundings : list of TEMSounding
        Station-wise soundings built from TEMAVG rows.
    results : list of dict
        Frequency-domain transform result dictionaries.
    collection : object, optional
        EDI collection created when ``return_collection=True``.
    written_paths : list of str
        EDI file paths written when ``savepath`` is supplied.
    """

    survey: TEMSurvey
    soundings: list[TEMSounding]
    results: list[dict[str, Any]] = field(default_factory=list)
    collection: Any | None = None
    written_paths: list[str] = field(default_factory=list)
    verbose: int = 0
    logger: object | None = None

    @property
    def n_soundings(self) -> int:
        """Number of generated station soundings."""
        return len(self.soundings)

    @property
    def n_results(self) -> int:
        """Number of transformed frequency-domain results."""
        return len(self.results)

    @property
    def stations(self) -> list[str]:
        """Station names represented by the generated soundings."""
        return [snd.station_name for snd in self.soundings]


def read_temavg_soundings(
    path: PathLike,
    *,
    pattern: str = "*.AVG",
    coordinate_file: PathLike | None = None,
    stems: list[str] | None = None,
    component: str = "Hz",
    frequency: float | None = None,
    data_column: str = "magnitude",
    magnitude_unit: str = "uV/A",
    data_type: str = "voltage",
    rx_turns: int = 1,
    tx_turns: int = 1,
    min_gates: int = 1,
    verbose: int = 0,
    logger: object | None = None,
) -> list[TEMSounding]:
    """Read a TEMAVG folder and return station soundings.

    This helper performs the first half of the TEMAVG pipeline:
    parse a survey folder, attach coordinates when available, and
    build one :class:`TEMSounding` per station.

    Parameters
    ----------
    path : path-like
        Folder containing TEMAVG ``.AVG`` files.
    pattern : str, default "*.AVG"
        Glob pattern used to select processed average files.
    coordinate_file : path-like, optional
        Coordinate table to attach by ``(profile, station)``.
        If omitted, common coordinate filenames are discovered
        automatically.
    stems : list of str, optional
        File stems to convert, for example ``["TEM100"]``.
        When omitted, all parsed ``.AVG`` files are used.
    component : str, default "Hz"
        TEMAVG component to select.
    frequency : float, optional
        Repetition frequency to select.
    data_column : {"magnitude"}, default "magnitude"
        TEMAVG transient column used for the sounding decay.
    magnitude_unit : str, default "uV/A"
        Unit of the TEMAVG magnitude column.
    data_type : str, default "voltage"
        Data convention stored in each sounding.
    rx_turns, tx_turns : int, default 1
        Receiver and transmitter turn counts.
    min_gates : int, default 1
        Minimum number of selected gates required per station.

    Returns
    -------
    list of TEMSounding
        Station-wise soundings ready for transformation.
    """
    survey = read_temavg_survey(
        path,
        pattern=pattern,
        coordinate_file=coordinate_file,
        verbose=verbose,
        logger=logger,
    )
    return survey.to_soundings(
        stems=stems,
        component=component,
        frequency=frequency,
        data_column=data_column,
        magnitude_unit=magnitude_unit,
        data_type=data_type,
        rx_turns=rx_turns,
        tx_turns=tx_turns,
        min_gates=min_gates,
        verbose=verbose,
        logger=logger,
    )


def transform_temavg_survey(
    path: PathLike,
    *,
    pattern: str = "*.AVG",
    coordinate_file: PathLike | None = None,
    stems: list[str] | None = None,
    component: str = "Hz",
    frequency: float | None = None,
    data_column: str = "magnitude",
    magnitude_unit: str = "uV/A",
    data_type: str = "voltage",
    rx_turns: int = 1,
    tx_turns: int = 1,
    min_gates: int = 1,
    method: str = "late_time",
    freq_convention: str = "skin_depth",
    phase_mode: str = "homogeneous",
    loop_geometry_correction: bool = True,
    return_collection: bool = False,
    savepath: PathLike | None = None,
    verbose: int = 0,
    logger: object | None = None,
) -> TEMAVGConversion:
    """Run a complete TEMAVG-to-frequency workflow.

    The function reads a TEMAVG survey folder, builds station
    soundings, transforms them to pseudo-frequency apparent
    resistivity and impedance, and optionally creates or writes EDI
    output.

    Parameters
    ----------
    path : path-like
        Folder containing TEMAVG ``.AVG`` files.
    pattern : str, default "*.AVG"
        Glob pattern used to select processed average files.
    coordinate_file : path-like, optional
        Coordinate table to attach by profile and station.
    stems : list of str, optional
        File stems to convert. Use this to process one profile
        before scaling to a full survey folder.
    component : str, default "Hz"
        TEMAVG component to select.
    frequency : float, optional
        Repetition frequency to select.
    data_column : {"magnitude"}, default "magnitude"
        TEMAVG transient column used for the sounding decay.
    magnitude_unit : str, default "uV/A"
        Unit of the TEMAVG magnitude column.
    data_type : str, default "voltage"
        Data convention stored in each sounding.
    rx_turns, tx_turns : int, default 1
        Receiver and transmitter turn counts.
    min_gates : int, default 1
        Minimum number of selected gates required per station.
    method : {"late_time", "fourier"}, default "late_time"
        Transformation method used by :class:`TEMtoEDI`.
    freq_convention : str, default "skin_depth"
        Pseudo-frequency convention for late-time transforms.
    phase_mode : str, default "homogeneous"
        Phase model used for synthetic impedance.
    loop_geometry_correction : bool, default True
        Whether to apply in-loop geometry corrections.
    return_collection : bool, default False
        If ``True``, also build an EDI collection in memory.
    savepath : path-like, optional
        Directory where synthetic EDI files are written.

    Returns
    -------
    TEMAVGConversion
        Bundle containing the parsed survey, station soundings,
        transform results, and optional EDI outputs.

    Examples
    --------
    >>> from pycsamt.tdem import transform_temavg_survey
    >>> out = transform_temavg_survey(
    ...     "data/TEMAVG/JIANGSU",
    ...     stems=["TEM100"],
    ... )
    >>> out.n_soundings > 0
    True
    """
    from .transform import FourierTransform, LateTimeTransform, TEMtoEDI

    survey = read_temavg_survey(
        path,
        pattern=pattern,
        coordinate_file=coordinate_file,
        verbose=verbose,
        logger=logger,
    )
    soundings = survey.to_soundings(
        stems=stems,
        component=component,
        frequency=frequency,
        data_column=data_column,
        magnitude_unit=magnitude_unit,
        data_type=data_type,
        rx_turns=rx_turns,
        tx_turns=tx_turns,
        min_gates=min_gates,
        verbose=verbose,
        logger=logger,
    )
    method_name = method.lower()
    if method_name == "late_time":
        transformer = LateTimeTransform(
            freq_convention=freq_convention,
            phase_mode=phase_mode,
            loop_geometry_correction=loop_geometry_correction,
            verbose=verbose,
            logger=logger,
        )
    elif method_name == "fourier":
        transformer = FourierTransform(
            verbose=verbose,
            logger=logger,
        )
    else:
        msg = "method must be 'late_time' or 'fourier'."
        raise ValueError(msg)

    results = transformer.transform_many(soundings)
    converter = None
    if return_collection or savepath is not None:
        converter = TEMtoEDI(
            method=method_name,
            freq_convention=freq_convention,
            phase_mode=phase_mode,
            loop_geometry_correction=loop_geometry_correction,
            out_dir=savepath or "edi_out/tem",
            verbose=verbose,
            logger=logger,
        )
    collection = (
        converter.transform_many(soundings)
        if return_collection and converter is not None
        else None
    )
    written_paths = (
        converter.save(soundings, path=savepath)
        if savepath is not None and converter is not None
        else []
    )
    return TEMAVGConversion(
        survey=survey,
        soundings=soundings,
        results=results,
        collection=collection,
        written_paths=written_paths,
        verbose=verbose,
        logger=logger,
    )
