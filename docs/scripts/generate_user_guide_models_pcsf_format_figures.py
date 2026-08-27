"""Generate the Occam2D PCSF section figure for the models user guide.

Run from the repository root::

    python docs/scripts/generate_user_guide_models_pcsf_format_figures.py \
        --save-dir docs/source/images/user_guide/models

A single wide Occam2D ``grid2d`` section loaded from a real persisted
PCSF file, cropped to 1.5 km depth and draped over real Tongkeng station
topography, with station name labels along the terrain (smart-thinned to
avoid overlap). A ModEM ``grid3d`` point-cloud view will follow in its
own figure once app/mapview gains direct PCSF/PCSM support.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

#: Occam2D section display depth cap (m) -- matches OccamConfig's own
#: max_depth default (see pycsamt.models.occam2d.config), a sensible
#: near-surface target for CSAMT-scale surveys like Tongkeng.
OCCAM_DISPLAY_MAX_DEPTH_M = 1500.0

#: Lowest elevation (km a.s.l.) shown on the y-axis. A visual crop only
#: -- the draped section is still computed for the full 1.5 km depth
#: below each station; this just keeps the low-lying deepest corner
#: (below the lowest-elevation station) out of the printed figure.
OCCAM_DISPLAY_MIN_ELEV_KM = -0.6


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_tongkeng_topography(root: Path) -> tuple[np.ndarray, np.ndarray]:
    """Real Tongkeng station topography from ``data/avg/k1.stn``.

    ``data/occam2D`` (the mtpy-anonymised S00..S46 profile used by the
    PCSF conversion demo) carries no elevation of its own -- it has no
    source EDI files, only the already-built Occam data/mesh/model. The
    original Tongkeng AVG station list (``k1.stn``, columns
    ``dot, e, n, h``) is the same survey and does carry real elevation
    (``h``, m a.s.l.). Chainage is the cumulative real Easting/Northing
    path distance (not the nominal ``dot`` design stake spacing), which
    tracks the Occam profile's own real station offsets far more closely
    than the nominal 50 m grid does.

    Returns
    -------
    chainage_km, elevation_m : numpy.ndarray
    """
    path = root / "data" / "avg" / "k1.stn"
    e: list[float] = []
    n: list[float] = []
    h: list[float] = []
    with path.open() as fh:
        next(fh)  # header: "dot","e","n","h"
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) != 4:
                continue
            _dot, ei, ni, hi = (float(v) for v in parts)
            e.append(ei)
            n.append(ni)
            h.append(hi)
    e_arr = np.asarray(e, dtype=float)
    n_arr = np.asarray(n, dtype=float)
    elevation_m = np.asarray(h, dtype=float)
    seg = np.sqrt(np.diff(e_arr) ** 2 + np.diff(n_arr) ** 2)
    chainage_km = np.concatenate([[0.0], np.cumsum(seg)]) / 1000.0
    return chainage_km, elevation_m


def make_figure():
    """Return the Occam2D PCSF terrain-draped section figure."""
    from pycsamt.format import read_pcsf
    from pycsamt.models.occam2d.data import OccamData
    from pycsamt.topo import build_topo_section
    from pycsamt.topo.config import TopoConfig
    from pycsamt.topo.overlay import draw_topo_section

    root = _repository_root()
    output = root / "examples" / "pcsf_conversion_demo" / "output"
    occam_model = read_pcsf(output / "occam2d_no_topo.pcsf")

    # occam2d_no_topo.pcsf carries no station elevation of its own (no
    # source EDI, see occam2d_to_pcsf's docstring) -- fall back to the
    # real Tongkeng AVG station topography for the drape.
    occam_geometry = occam_model.geometry
    occam_values_full = np.log10(occam_model.resistivity)
    k1_chainage_km, k1_elevation_m = _load_tongkeng_topography(root)

    # geometry.x carries real station chainage (see occam2d_to_pcsf's
    # docstring) but still includes the mesh's exponentially-widening
    # horizontal padding columns on both ends -- far outside k1.stn's own
    # chainage range, where interp_elev clamps to the boundary elevation
    # and draws a long flat artifact shelf. Crop to the real station
    # extent (with a small margin) before draping.
    occam_data = OccamData.read(root / "data" / "occam2D" / "OccamDataFile.dat")
    station_lo = float(occam_data.offsets.min())
    station_hi = float(occam_data.offsets.max())
    margin = 0.04 * (station_hi - station_lo)
    x_full = np.asarray(occam_geometry.x, dtype=float)
    station_cols = (x_full >= station_lo - margin) & (x_full <= station_hi + margin)
    x_stations = x_full[station_cols]
    values_stations = occam_values_full[:, station_cols]

    occam_topo = build_topo_section(
        (x_stations, occam_geometry.z, values_stations),
        elevation=k1_elevation_m,
        chainage=k1_chainage_km,
        station_x=occam_data.offsets,
        station_names=occam_data.sites,
        model_unit="m",
        depth_max=OCCAM_DISPLAY_MAX_DEPTH_M,
        log_rho=True,
        clip_above_surface=True,
    )

    # Colour range computed after the depth crop, so the scale reflects
    # the displayed 1.5 km window rather than the full 6 km model.
    occam_finite = occam_topo.values[np.isfinite(occam_topo.values)]
    occam_norm = Normalize(*np.nanpercentile(occam_finite, (1.0, 99.0)), clip=True)
    cmap = "viridis"

    plt.rcParams.update(
        {
            "font.size": 9.5,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    # Wide figure: station labels are rotated 90 deg along the terrain
    # (draw_topo_section's built-in smart-thinning already avoids
    # overlap), and a generous width shows more of the 47 real names.
    fig, ax1 = plt.subplots(figsize=(17.0, 6.4), constrained_layout=True)

    mesh = ax1.pcolormesh(
        occam_topo.x_nodes_km,
        occam_topo.z_draped_km,
        occam_topo.values,
        cmap=cmap,
        norm=occam_norm,
        shading="auto",
        rasterized=True,
    )
    ax1.set_title(
        "Occam2D PCSF · grid2d — draped over real topography",
        loc="left",
        fontweight="bold",
    )
    ax1.set_xlabel("Distance along profile (km)")
    ax1.set_ylabel("Elevation (km)")
    if ax1.yaxis_inverted():
        ax1.invert_yaxis()
    _top = float(np.nanmax(occam_topo.surface_km))
    _bottom = max(
        float(np.nanmin(occam_topo.z_draped_km)), OCCAM_DISPLAY_MIN_ELEV_KM
    )
    _margin = 0.04 * max(_top - _bottom, 1e-6)
    ax1.set_ylim(_bottom, _top + _margin)
    ax1.set_xlim(
        float(occam_topo.x_nodes_km.min()), float(occam_topo.x_nodes_km.max())
    )
    ax1.grid(color="0.85", linewidth=0.45, alpha=0.7)
    draw_topo_section(
        ax1,
        occam_topo.chainage_km,
        occam_topo.elev_km * 1000.0,
        occam_topo.station_names,
        station_x_km=occam_topo.station_x_km,
        cfg=TopoConfig(station_pins_at_surface=True, fill_alpha=0.0),
        dark=False,
    )
    ax1.text(
        0.01,
        0.03,
        "Tongkeng CSAMT · 47 stations · depth ≤ 1.5 km\n"
        "topography: data/avg/k1.stn (real station elevation)",
        transform=ax1.transAxes,
        va="bottom",
        ha="left",
        fontsize=8.5,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 3},
    )

    occam_bar = fig.colorbar(
        mesh,
        ax=ax1,
        orientation="vertical",
        fraction=0.028,
        pad=0.012,
        aspect=28,
        extend="both",
    )
    occam_bar.set_label(r"$\log_{10}(\rho\;[\Omega\,\mathrm{m}])$")
    fig.suptitle(
        "Occam2D PCSF section draped over real topography",
        fontsize=14,
        fontweight="bold",
    )
    return fig, (ax1,)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--save-dir",
        type=Path,
        default=None,
        help="Directory to save the figure into (PNG). Shown interactively "
        "when omitted.",
    )
    return parser.parse_args()


def main() -> None:
    args = _arguments()
    fig, _axes = make_figure()
    if args.save_dir is not None:
        args.save_dir.mkdir(parents=True, exist_ok=True)
        out = args.save_dir / "pcsf_occam2d_topo_section.png"
        fig.savefig(out, dpi=220, bbox_inches="tight", facecolor="white")
        print(f"wrote {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
