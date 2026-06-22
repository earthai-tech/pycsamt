# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt model zoo — registry and download utilities.

Pre-trained checkpoints are hosted at
https://github.com/earthai-tech/pycsamt-models
and are downloaded on first use to a local cache
(``~/.pycsamt/model_zoo/`` by default).

Model naming convention
-----------------------
``<method>-<arch>-<n_layers>layer-v<version>``

Examples: ``mt1d-resnet-5layer-v1``, ``mt1d-cnn-7layer-v1``,
``joint-drcnn-5layer-v1``.
"""
from __future__ import annotations

import hashlib
import os
import warnings
from pathlib import Path
from typing import Any, Dict, Optional

__all__ = [
    "list_pretrained",
    "get_pretrained_info",
    "download_checkpoint",
]

# ─────────────────────────────────────────────────────────────────────────────
# Registry
# ─────────────────────────────────────────────────────────────────────────────

_GITHUB_RELEASE = (
    "https://github.com/earthai-tech/pycsamt-models/releases/download/v0.1"
)

_MODEL_ZOO: Dict[str, Dict[str, Any]] = {
    "mt1d-resnet-5layer-v1": {
        "url": f"{_GITHUB_RELEASE}/mt1d-resnet-5layer-v1.npz",
        "arch": "resnet",
        "n_layers": 5,
        "solver": "mt1d",
        "n_freqs": 30,
        "include_phase": True,
        "log_thickness": True,
        "description": (
            "MT 1-D ResNet (Liu 2021 style), 5 layers.  "
            "Trained on 50 000 synthetic profiles spanning "
            "ρ ∈ [1, 10 000] Ω·m, depth ≤ 2 km."
        ),
        "md5": None,  # populated when weights are released
    },
    "mt1d-cnn-5layer-v1": {
        "url": f"{_GITHUB_RELEASE}/mt1d-cnn-5layer-v1.npz",
        "arch": "cnn1d",
        "n_layers": 5,
        "solver": "mt1d",
        "n_freqs": 30,
        "include_phase": True,
        "log_thickness": True,
        "description": (
            "MT 1-D CNN (Puzyrev 2021 style), 5 layers.  "
            "Faster inference than ResNet; slightly lower accuracy."
        ),
        "md5": None,
    },
    "mt1d-resnet-7layer-v1": {
        "url": f"{_GITHUB_RELEASE}/mt1d-resnet-7layer-v1.npz",
        "arch": "resnet",
        "n_layers": 7,
        "solver": "mt1d",
        "n_freqs": 30,
        "include_phase": True,
        "log_thickness": True,
        "description": (
            "MT 1-D ResNet, 7 layers.  Higher vertical resolution "
            "for shallow targets."
        ),
        "md5": None,
    },
    "csamt1d-resnet-5layer-v1": {
        "url": f"{_GITHUB_RELEASE}/csamt1d-resnet-5layer-v1.npz",
        "arch": "resnet",
        "n_layers": 5,
        "solver": "csamt1d",
        "n_freqs": 25,
        "include_phase": True,
        "log_thickness": True,
        "description": (
            "CSAMT 1-D ResNet, 5 layers.  Trained with Zonge & Hughes "
            "(1991) near-field correction."
        ),
        "md5": None,
    },
    "tem1d-fcn-5layer-v1": {
        "url": f"{_GITHUB_RELEASE}/tem1d-fcn-5layer-v1.npz",
        "arch": "fcn",
        "n_layers": 5,
        "solver": "tem1d",
        "n_freqs": None,
        "include_phase": False,
        "log_thickness": True,
        "description": (
            "TEM 1-D FCN (Moghadas 2020 style), 5 layers.  "
            "Input: dBz/dt decay curve."
        ),
        "md5": None,
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Public utilities
# ─────────────────────────────────────────────────────────────────────────────

def list_pretrained() -> Dict[str, str]:
    """
    Return a dict of available pre-trained models.

    Returns
    -------
    models : dict
        Mapping ``name → description``.

    Examples
    --------
    >>> from pycsamt.ai._zoo import list_pretrained
    >>> for name, desc in list_pretrained().items():
    ...     print(f"{name:35s}  {desc[:60]}")
    """
    return {k: v["description"] for k, v in _MODEL_ZOO.items()}


def get_pretrained_info(name: str) -> Dict[str, Any]:
    """
    Return metadata for a registered pre-trained model.

    Parameters
    ----------
    name : str
        Model identifier (see :func:`list_pretrained`).

    Returns
    -------
    info : dict

    Raises
    ------
    KeyError
        If *name* is not in the registry.
    """
    if name not in _MODEL_ZOO:
        available = ", ".join(sorted(_MODEL_ZOO))
        raise KeyError(
            f"Unknown model {name!r}.  "
            f"Available models: {available}"
        )
    return dict(_MODEL_ZOO[name])


def _get_cache_dir(cache_dir: Optional[str] = None) -> Path:
    """Return (and create) the local model cache directory."""
    if cache_dir is not None:
        return Path(cache_dir)
    env = os.environ.get("PYCSAMT_MODEL_CACHE")
    if env:
        return Path(env)
    return Path.home() / ".pycsamt" / "model_zoo"


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def download_checkpoint(
    name: str,
    *,
    cache_dir: Optional[str] = None,
    force: bool = False,
    verbose: bool = True,
) -> Path:
    """
    Download a pre-trained checkpoint to the local cache.

    Parameters
    ----------
    name : str
        Model name (see :func:`list_pretrained`).
    cache_dir : str or None
        Override the default cache location
        (``~/.pycsamt/model_zoo/``).
    force : bool
        Re-download even if the file exists locally.
    verbose : bool

    Returns
    -------
    path : Path
        Local path of the downloaded ``.npz`` file.

    Raises
    ------
    KeyError
        If *name* is not in the model zoo.
    RuntimeError
        If the download fails or the MD5 check fails.
    """
    import urllib.request

    info = get_pretrained_info(name)
    url = info["url"]
    cache = _get_cache_dir(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    fname = url.split("/")[-1]
    fpath = cache / fname

    if fpath.exists() and not force:
        if verbose:
            print(f"  Using cached {fname}  ({fpath})")
        return fpath

    if verbose:
        print(f"  Downloading {name} from {url} …")

    try:
        urllib.request.urlretrieve(url, fpath)
    except Exception as exc:
        if fpath.exists():
            fpath.unlink()
        raise RuntimeError(
            f"Failed to download {url!r}: {exc}\n"
            "Pre-trained weights for pycsamt are scheduled for Phase 5 "
            "and may not yet be publicly available.  "
            "Train your own model with EMInverter1D.fit() or check "
            "https://github.com/earthai-tech/pycsamt-models for updates."
        ) from exc

    # Optional MD5 verification
    expected_md5 = info.get("md5")
    if expected_md5:
        actual = _md5_file(fpath)
        if actual != expected_md5:
            fpath.unlink()
            raise RuntimeError(
                f"MD5 mismatch for {fname}: "
                f"expected {expected_md5}, got {actual}."
            )

    if verbose:
        print(f"  Saved to {fpath}")
    return fpath
