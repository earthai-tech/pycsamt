# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""High-level inversion workflow entry point.

``pycsamt.inversion.workflow`` is the small orchestration layer that turns an
inversion configuration into a backend run. It keeps application code stable:
callers build an :class:`InversionConfig`, pass it to
:class:`InversionWorkflow`, and receive an
:class:`pycsamt.inversion.results.InversionResult` regardless of whether the
selected backend is built-in, SimPEG, pyGIMLi, Occam2D, or ModEM.

The workflow object does not implement physics itself. It validates and stores
configuration, instantiates the selected backend lazily through the backend
registry, and delegates execution to ``backend.run``.

Examples
--------
>>> from pycsamt.inversion.workflow import run_inversion
>>> result = run_inversion(
...     method="mt",
...     dimension="1d",
...     backend="builtin",
...     data={"freqs": [1.0, 10.0],
...           "rho_a": [100.0, 120.0],
...           "phase": [45.0, 47.0]},
...     max_iter=4,
... )  # doctest: +SKIP

See Also
--------
pycsamt.inversion.config.InversionConfig
    Configuration object consumed by this module.
pycsamt.inversion.results.InversionResult
    Backend-neutral result returned by this module.
"""

from __future__ import annotations

from typing import Any

from ..api.property import PyCSAMTObject
from .config import InversionConfig
from .doc import _inversion_param_docs
from .results import InversionResult

__all__ = ["InversionWorkflow", "run_inversion"]


class InversionWorkflow(PyCSAMTObject):
    def __init__(self, config: InversionConfig | dict[str, Any] | None = None, **kw: Any):
        """Create a workflow and instantiate its backend.

        Parameters
        ----------
        config : InversionConfig or dict, optional
            Complete inversion configuration. Dictionaries are converted to
            :class:`pycsamt.inversion.config.InversionConfig`.
        **kw
            Configuration overrides. When *config* is ``None``, these keyword
            arguments build a new ``InversionConfig``. When *config* is a dict,
            they override matching keys. When *config* is an object, they are
            passed through ``config.clone(**kw)``.

        Attributes
        ----------
        config : InversionConfig
            Normalized configuration used by the backend.
        backend : BaseInversionBackend
            Instantiated backend selected by ``config.backend``.

        Examples
        --------
        >>> from pycsamt.inversion.workflow import InversionWorkflow
        >>> workflow = InversionWorkflow(
        ...     method="mt",
        ...     dimension="1d",
        ...     backend="builtin",
        ...     data={"freqs": [1.0],
        ...           "rho_a": [100.0],
        ...           "phase": [45.0]},
        ... )
        >>> workflow.config.backend
        'builtin'
        """
        if config is None:
            config = InversionConfig(**kw)
        elif isinstance(config, dict):
            values = dict(config)
            values.update(kw)
            config = InversionConfig(**values)
        elif kw:
            config = config.clone(**kw)
        self.config = config
        self.backend = config.to_backend()

    def run(self, data: Any | None = None) -> InversionResult:
        """Run the selected backend.

        Parameters
        ----------
        data : mapping, object, sequence, or path-like, optional
            Optional data override for this call. When omitted, ``config.data``
            is used. Accepted values are backend-dependent but normally pass
            through :class:`pycsamt.inversion.data.EMData.coerce`.

        Returns
        -------
        InversionResult
            Backend-neutral result containing recovered model, RMS/objective
            diagnostics, files, warnings, uncertainty, history, and native
            backend objects when available.

        Examples
        --------
        >>> from pycsamt.inversion.workflow import InversionWorkflow
        >>> workflow = InversionWorkflow(
        ...     method="mt",
        ...     dimension="1d",
        ...     backend="builtin",
        ...     max_iter=4,
        ... )
        >>> result = workflow.run({
        ...     "freqs": [1.0, 10.0],
        ...     "rho_a": [100.0, 120.0],
        ...     "phase": [45.0, 47.0],
        ... })  # doctest: +SKIP
        """
        return self.backend.run(data=data)


def run_inversion(config: InversionConfig | dict[str, Any] | None = None, **kw: Any) -> InversionResult:
    """Run an inversion in one call.

    This convenience function is equivalent to
    ``InversionWorkflow(config, **kw).run()``. It is useful in scripts,
    notebooks, tests, and demos where the workflow object does not need to be
    reused.

    Parameters
    ----------
    config : InversionConfig or dict, optional
        Complete inversion configuration. Dictionaries are converted to
        :class:`pycsamt.inversion.config.InversionConfig`.
    **kw
        Configuration overrides or full configuration fields when *config* is
        omitted.

    Returns
    -------
    InversionResult
        Backend-neutral inversion result.

    Examples
    --------
    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="mt",
    ...     dimension="1d",
    ...     backend="builtin",
    ...     data={"freqs": [1.0, 10.0],
    ...           "rho_a": [100.0, 120.0],
    ...           "phase": [45.0, 47.0]},
    ...     max_iter=4,
    ... )  # doctest: +SKIP
    """
    return InversionWorkflow(config, **kw).run()


InversionWorkflow.__doc__ = rf"""
Execute a configured EM inversion backend.

``InversionWorkflow`` is the object-oriented entry point for running EM
inversions. It owns one normalized
:class:`pycsamt.inversion.config.InversionConfig` and one instantiated backend.
Calling :meth:`run` delegates to that backend and returns a common
:class:`pycsamt.inversion.results.InversionResult`.

Parameters
----------
{_inversion_param_docs.workflow.config}
{_inversion_param_docs.workflow.keyword_overrides}

Returns
-------
{_inversion_param_docs.result.result}

Notes
-----
The workflow does not mutate caller-provided data. A per-call data override can
be supplied to :meth:`run`, which is convenient for applying the same
configuration to multiple soundings or profiles.

{_inversion_param_docs.workflow.examples}

See Also
--------
pycsamt.inversion.config.InversionConfig
    Configuration consumed by the workflow.
pycsamt.inversion.backends.get_backend
    Backend registry used by ``InversionConfig.to_backend``.
pycsamt.inversion.results.InversionResult
    Result returned by :meth:`run`.

{_inversion_param_docs.workflow.references}
"""
