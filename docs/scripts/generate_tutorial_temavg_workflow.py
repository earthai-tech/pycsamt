"""Reproduce the data, QC, and correction figures for the TEMAVG tutorial.

Everything here is executed against the real, bundled
``data/TEMAVG/JIANGSU`` survey folder -- no synthetic stand-ins. Run this
script directly to regenerate every figure referenced by
``docs/source/tutorials/process_temavg_survey.rst``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pycsamt.tdem import TEMSounding, TEMtoEDI, read_temavg_survey

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "TEMAVG" / "JIANGSU"
OUT = ROOT / "docs" / "source" / "images" / "tutorials" / "process_temavg_survey"
RESULTS = ROOT / "results" / "process_temavg_survey"

PROFILE_STEM = "TEM100"

# Best-effort CRS for the survey's zone-19 "gauss_x"/"gauss_y" coordinates.
# CGCS2000 / Gauss-Kruger zone 19 (EPSG:4497) is the registered CRS whose
# name matches the zone prefix embedded in gauss_y (19,xxx,xxx). Projecting
# through it does NOT land in Jiangsu (lands near 38.8N, 111.1E -- Shanxi),
# so this is a labeled guess, not a verified survey location; see the
# tutorial's CRS provenance note.
GK_EPSG = 4497


def load_profile(stem: str = PROFILE_STEM):
    """Read the JIANGSU survey and build raw soundings for one profile."""
    survey = read_temavg_survey(DATA)
    avg = survey.get(stem)
    soundings = survey.to_soundings(stems=[stem])
    # Sort by along-profile station value; EDICollection/Sites ordering is
    # alphabetical by station name ("TEM100_1000" < "TEM100_120"), which is
    # not the physical along-profile order.
    soundings = sorted(soundings, key=lambda s: float(s.station_name.split("_")[1]))
    return survey, avg, soundings


def attach_geographic_coords(soundings, survey, *, epsg: int = GK_EPSG):
    """Overwrite each sounding's local x/y with projected lon/lat.

    Uses :func:`pycsamt.gis.utils.epsg_project` on the coordinate table's
    ``gauss_x``/``gauss_y`` (already zone-prefixed, so no manual false-
    easting arithmetic is needed). This is a best-effort geographic anchor,
    not a verified one -- see :data:`GK_EPSG`.
    """
    from pycsamt.gis.utils import epsg_project

    for snd in soundings:
        stem, station = snd.station_name.split("_")
        profile = float("".join(ch for ch in stem if ch.isdigit()))
        coord = survey.coordinates.get(profile, float(station))
        if coord is None:
            continue
        lon, lat = epsg_project(coord.gauss_y, coord.gauss_x, epsg, 4326)
        snd.x, snd.y = float(lon), float(lat)
    return soundings


def drop_noise_floor_gates(snd: TEMSounding) -> tuple[TEMSounding, int]:
    """Remove gates whose processed magnitude is negative.

    A negative TEMAVG ``magnitude`` at late time means the stacked
    transient has crossed zero -- the signal has fallen into the noise
    floor and the sign is not physically meaningful for the late-time
    apparent-resistivity formula, which only ever sees ``|dBdt|``. Feeding
    those gates in unfiltered does not raise an error; it silently
    fabricates a resistivity value from noise.
    """
    keep = snd.data > 0.0
    n_dropped = int((~keep).sum())
    cleaned = TEMSounding(
        time_gates=snd.time_gates[keep],
        data=snd.data[keep],
        current=snd.current,
        tx_area=snd.tx_area,
        data_type=snd.data_type,
        tx_turns=snd.tx_turns,
        rx_area=snd.rx_area,
        rx_turns=snd.rx_turns,
        offset=snd.offset,
        loop_shape=snd.loop_shape,
        loop_dims=snd.loop_dims,
        station_name=snd.station_name,
        x=snd.x,
        y=snd.y,
        elevation=snd.elevation,
        error=(snd.error[keep] if snd.error is not None else None),
        waveform=snd.waveform,
    )
    return cleaned, n_dropped


def make_survey_overview(survey, avg, soundings) -> Path:
    """Plot the full station grid and the TEM100 elevation profile."""
    coords = survey.coordinates.to_dataframe()
    line = coords[coords.profile == float(_profile_number(PROFILE_STEM))]
    line = line.sort_values("point")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)

    ax = axes[0]
    ax.scatter(
        coords.gauss_y, coords.gauss_x, s=2, c="#b0b0b0", label="planned grid pegs"
    )
    ax.scatter(
        line.gauss_y, line.gauss_x, s=14, c="#d62728", label=f"{PROFILE_STEM} stations"
    )
    ax.set(
        title=f"JIANGSU grid: {len(survey.avg_files)} surveyed profiles of "
        f"{coords.profile.nunique()} planned",
        # Chinese Gauss-Krüger convention: X = northing, Y = zone-prefixed
        # easting -- the opposite of the usual (easting, northing) reading.
        xlabel="Gauss-Krüger Y (zone-prefixed easting, m)",
        ylabel="Gauss-Krüger X (northing, m)",
    )
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    ax = axes[1]
    ax.plot(line.point, line.elevation, "o-", ms=3.5, color="#1f77b4")
    ax.set(
        title=f"{PROFILE_STEM}: measured topography ({len(line)} stations)",
        xlabel="Station chainage (m)",
        ylabel="Elevation (m)",
    )
    ax.grid(alpha=0.25)

    target = OUT / "jiangsu_survey_overview.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)
    return target


def _profile_number(stem: str) -> float:
    digits = "".join(ch for ch in stem if ch.isdigit())
    return float(digits)


def make_time_domain_qc(avg, soundings) -> tuple[Path, dict]:
    """Plot signed decay curves and the per-window sign-reversal count."""
    n_windows = max(rec.window for rec in avg.records)
    neg_count = np.zeros(n_windows, dtype=int)
    for rec in avg.records:
        if rec.magnitude < 0.0:
            neg_count[rec.window - 1] += 1
    n_stations = len(avg.stations)

    picks = [soundings[0], soundings[len(soundings) // 2], soundings[-1]]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)

    ax = axes[0]
    for snd in picks:
        t_ms = snd.time_gates * 1e3
        pos = snd.data > 0.0
        ax.loglog(t_ms[pos], snd.data[pos], "o-", ms=4, lw=1.2, label=snd.station_name)
        if (~pos).any():
            ax.loglog(
                t_ms[~pos],
                np.abs(snd.data[~pos]),
                "x",
                ms=9,
                mew=2,
                color="black",
                label="_nolegend_",
            )
    ax.set(
        title="Signed decay curves (x = negative magnitude)",
        xlabel="Time (ms)",
        ylabel="|processed magnitude| (V)",
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.bar(np.arange(1, n_windows + 1), neg_count, color="#d62728")
    ax.set(
        title=f"{PROFILE_STEM}: stations with negative magnitude per window",
        xlabel="Time-window number",
        ylabel=f"Station count (of {n_stations})",
    )
    ax.grid(alpha=0.25, axis="y")

    target = OUT / "tem100_time_domain_qc.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)

    stats = {
        "n_windows": n_windows,
        "neg_count_per_window": neg_count,
        "n_stations": n_stations,
        "n_records": len(avg.records),
        "n_negative_total": int(neg_count.sum()),
    }
    return target, stats


def convert_profile_to_edi(soundings, out_dir: Path):
    """Drop noise-floor gates, transform to EDI, write, and reload."""
    from pycsamt.emtools import ensure_sites

    cleaned, n_dropped = [], 0
    for snd in soundings:
        c, n = drop_noise_floor_gates(snd)
        cleaned.append(c)
        n_dropped += n

    conv = TEMtoEDI(method="late_time", phase_mode="homogeneous", out_dir=str(out_dir))
    written = conv.save(cleaned)
    sites = ensure_sites(out_dir, recursive=False).ordered(by="station")

    stats = {
        "n_dropped_total": n_dropped,
        "n_gates_in": sum(s.n_gates for s in soundings),
        "n_written": len(written),
        "n_reloaded": len(sites),
        "n_gates_out_min": min(len(s.freq) for s in sites),
        "n_gates_out_max": max(len(s.freq) for s in sites),
    }
    return sites, stats


def make_edi_processing(raw_sites):
    """Apply Hampel despiking and Hanning static-shift correction."""
    from pycsamt.emtools import correct_static_shift, hampel_filter_freq

    despiked = hampel_filter_freq(
        raw_sites, win=2, nsig=3.0, on="z", domain="magphase", inplace=False
    )
    n_cells_total = 0
    n_cells_changed = 0
    for a, b in zip(raw_sites, despiked):
        za, zb = a.z[..., 0, 1], b.z[..., 0, 1]
        n_cells_total += za.size
        n_cells_changed += int(np.sum(~np.isclose(za, zb, equal_nan=True)))

    corrected = correct_static_shift(
        despiked, window_m=200.0, spacing_m=20.0, comp="xy", inplace=False
    )

    station = np.array([float(s.name.split("_")[1]) for s in raw_sites])
    freq_common = np.unique(np.concatenate([s.freq for s in raw_sites]))
    freq_common.sort()

    def _grid(sites):
        grid = np.full((len(sites), freq_common.size), np.nan)
        for i, s in enumerate(sites):
            rho = s.rho[..., 0, 1]
            idx = np.searchsorted(freq_common, s.freq)
            grid[i, idx] = rho
        return grid

    rho_raw = _grid(raw_sites)
    rho_despiked = _grid(despiked)
    rho_corrected = _grid(corrected)

    with np.errstate(divide="ignore", invalid="ignore"):
        d_hampel = np.log10(rho_despiked) - np.log10(rho_raw)
        d_static = np.log10(rho_corrected) - np.log10(rho_despiked)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    for ax, delta, title in (
        (axes[0], d_hampel, "Hampel despike: Δlog10 ρ_xy"),
        (axes[1], d_static, "Static-shift correction: Δlog10 ρ_xy"),
    ):
        vmax = np.nanmax(np.abs(delta)) or 1e-6
        pc = ax.pcolormesh(
            station,
            freq_common,
            delta.T,
            shading="nearest",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
        )
        ax.set_yscale("log")
        ax.set(title=title, xlabel="Station chainage (m)", ylabel="Pseudo-frequency (Hz)")
        fig.colorbar(pc, ax=ax, label=r"$\Delta\log_{10}\rho_a$")

    target = OUT / "tem100_edi_processing_changes.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)

    log_diff_static = np.abs(d_static[np.isfinite(d_static)])
    stats = {
        "hampel_cells_changed": n_cells_changed,
        "hampel_cells_total": n_cells_total,
        "static_median_abs_log10": round(float(np.median(log_diff_static)), 4),
        "static_max_abs_log10": round(float(np.max(log_diff_static)), 4),
    }
    return target, despiked, corrected, stats


def make_response_panels(raw_sites, corrected_sites) -> Path:
    """Compare raw and corrected rho/phase at three representative stations."""
    picks_idx = [0, len(raw_sites) // 2, len(raw_sites) - 1]
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.2), constrained_layout=True)
    for col, idx in enumerate(picks_idx):
        raw, corr = raw_sites[idx], corrected_sites[idx]
        rho_raw = raw.rho[..., 0, 1]
        rho_corr = corr.rho[..., 0, 1]
        phase_raw = raw.phase[..., 0, 1]
        phase_corr = corr.phase[..., 0, 1]

        ax = axes[0, col]
        ax.loglog(raw.freq, rho_raw, "o-", ms=4, color="#7f7f7f", label="raw")
        ax.loglog(corr.freq, rho_corr, "s-", ms=4, color="#d62728", label="corrected")
        ax.set(title=raw.name, xlabel="", ylabel=r"$\rho_a$ ($\Omega\cdot$m)" if col == 0 else "")
        ax.grid(True, which="both", alpha=0.25)
        if col == 0:
            ax.legend(fontsize=8)

        ax = axes[1, col]
        ax.semilogx(raw.freq, phase_raw, "o-", ms=4, color="#7f7f7f")
        ax.semilogx(corr.freq, phase_corr, "s-", ms=4, color="#d62728")
        ax.set(
            xlabel="Pseudo-frequency (Hz)",
            ylabel=r"$\phi_{xy}$ (deg)" if col == 0 else "",
        )
        ax.grid(True, which="both", alpha=0.25)

    target = OUT / "tem100_three_station_comparison.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)
    return target


def make_mesh(corrected_sites):
    """Build and draw a real graded 2-D triangular mesh from the profile.

    Station chainage and elevation come straight from the corrected EDI
    collection now that :func:`attach_geographic_coords` has given every
    site real, non-NaN ``Site.coords``. ``surface_cell_m`` is set from
    this profile's own skin depth at its highest measured pseudo-
    frequency, not a fixed constant.
    """
    from pycsamt.api.mesh import draw_tri_mesh
    from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh

    station_x = np.array([float(s.name.split("_")[1]) for s in corrected_sites])
    elevation = np.array([s.coords[2] for s in corrected_sites])
    order = np.argsort(station_x)
    station_x, elevation = station_x[order], elevation[order]
    topo_z = elevation.max() - elevation

    freq_max = max(float(s.freq.max()) for s in corrected_sites)
    rho_at_fmax = np.array(
        [s.rho[np.argmax(s.freq), 0, 1] for s in corrected_sites]
    )
    skin_depth_m = 503.0 * np.sqrt(np.nanmedian(rho_at_fmax) / freq_max)
    surface_cell_m = round(float(0.04 * skin_depth_m), -1) or 10.0

    pad = 100.0
    x_range_m = (float(station_x.min() - pad), float(station_x.max() + pad))
    z_range_m = (0.0, 600.0)

    mesh = build_graded_tri_mesh(
        x_range_m,
        z_range_m,
        station_x,
        surface_cell_m=surface_cell_m,
        topo_x_m=station_x,
        topo_z_m=topo_z,
    )

    fig, ax = plt.subplots(figsize=(11.0, 5.0), constrained_layout=True)
    draw_tri_mesh(ax, mesh, preset="diagram")
    ax.plot(station_x, np.interp(station_x, station_x, topo_z), "rv", ms=6, zorder=5)
    ax.set(
        title=(
            f"{PROFILE_STEM}: graded triangular mesh "
            f"({mesh.n_triangles} triangles, surface cell {surface_cell_m:g} m)"
        ),
        xlabel="Chainage (m)",
        ylabel="Depth below highest station (m)",
    )
    ax.invert_yaxis()
    ax.set_xlim(*x_range_m)

    target = OUT / "tem100_triangular_mesh.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)

    stats = {
        "n_nodes": int(mesh.nodes_m.shape[0]),
        "n_triangles": int(mesh.n_triangles),
        "surface_cell_m": surface_cell_m,
        "skin_depth_m_at_fmax": round(float(skin_depth_m), 1),
        "freq_max_hz": round(freq_max, 1),
        "x_range_m": x_range_m,
        "z_range_m": z_range_m,
    }
    return target, mesh, stats


def run_ai_inversion(
    corrected_sites,
    *,
    n_train_profiles: int = 100,
    epochs: int = 100,
    patience: int = 15,
    seed: int = 0,
):
    """Train and gate a real ``Inv2DAgent(physics="mt2d_tri")`` run.

    Not part of the default ``__main__`` regeneration below -- at these
    settings (100 realizations x 15 frequencies x 51 stations, each a real
    :class:`~pycsamt.forward.maxwell.tri_fem2d.TriFEM2DAdapter` solve) this
    takes on the order of 15 minutes. Call it explicitly.

    Restricts to the 15-frequency band (125.0-2729.0 Hz) every one of the
    51 stations shares after the noise-floor gate drop -- some stations
    kept as few as 15 of 25 gates, so this, not an arbitrary round number,
    is the true common band.
    """
    import numpy as np

    from pycsamt.agents import Inv2DAgent
    from pycsamt.forward.maxwell.tri_fem2d import TriFEM2DAdapter
    from pycsamt.site.edit import select_freq_all

    try:
        import torch

        torch.manual_seed(seed)
    except ImportError:
        pass
    np.random.seed(seed)

    usable = select_freq_all(corrected_sites, fmin=125.0)
    freqs_hz = list(list(usable)[0].freq)

    station_x = np.array([float(s.name.split("_")[1]) for s in usable])
    order = np.argsort(station_x)
    station_x = station_x[order]
    elevation = np.array([s.coords[2] for s in usable])[order]
    topo_z = elevation.max() - elevation

    agent = Inv2DAgent(
        physics="mt2d_tri",
        epochs=epochs,
        patience=patience,
        n_freqs=len(freqs_hz),
        depth_max=600.0,
        n_train_profiles=n_train_profiles,
        n_stations_per_profile=len(station_x),
        station_spacing_m=20.0,
        mesh_target_cell_m=20.0,
        field_grid_cell_m=10.0,
        correlation_length_x_m=(100.0, 350.0),
        correlation_length_z_m=(40.0, 150.0),
        topo_x_m=station_x,
        topo_z_m=topo_z,
        mare2dem_adapter=TriFEM2DAdapter(),
    )
    result = agent.execute(
        {
            "sites": usable,
            "freqs": freqs_hz,
            "output_dir": str(RESULTS / "ai2d_tri"),
        }
    )
    if result.status == "success":
        import shutil

        fig_path = result.data["figure_paths"].get("inv2d_tri_section")
        if fig_path:
            shutil.copy(fig_path, OUT / "tem100_ai2d_tri_section.png")
    return result


def write_manifest(corrected_sites):
    """Write the corrected profile EDIs plus a coordinate manifest."""
    from pycsamt.site.export import write_sites

    out = RESULTS / "edi_corrected" / PROFILE_STEM
    paths = write_sites(
        corrected_sites, out, exist_ok=True, manifest_csv=out / "manifest.csv"
    )
    return out, paths


if __name__ == "__main__":
    survey, avg, soundings = load_profile()
    make_survey_overview(survey, avg, soundings)
    make_time_domain_qc(avg, soundings)
    attach_geographic_coords(soundings, survey)
    raw_sites, _ = convert_profile_to_edi(
        soundings, RESULTS / "edi_raw" / PROFILE_STEM
    )
    _, despiked, corrected, _ = make_edi_processing(raw_sites)
    make_response_panels(raw_sites, corrected)
    make_mesh(corrected)
    write_manifest(corrected)
    print("Done.")
