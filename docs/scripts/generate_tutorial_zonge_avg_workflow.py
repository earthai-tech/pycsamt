"""Reproduce the data audit and overview figure for the K1/K2 AVG tutorial."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pycsamt.transformers import AVGtoEDI
from pycsamt.zonge import AVG
from pycsamt.zonge.processing import ASTATIC


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "avg"
OUT = ROOT / "docs" / "source" / "images" / "tutorials" / "zonge_avg_workflow"


def load_line(name: str):
    """Load AVG and STN, project UTM 49N, and materialize EDI objects."""
    avg = AVG.from_file(DATA / f"{name}.AVG")
    avg.add_topography(DATA / f"{name}.stn", epsg=32649)
    measured = np.sort(np.asarray(avg.df["station"].unique(), dtype=float))
    surveyed = np.asarray(avg.topo.frame["station"], dtype=float)
    if not np.all(np.isin(measured, surveyed)):
        # K2 observations are at 25 m receiver midpoints whereas K2.stn
        # contains 50 m pegs.  Interpolate only inside surveyed chainage.
        source = avg.topo.frame.sort_values("station")
        if measured.min() < surveyed.min() or measured.max() > surveyed.max():
            raise ValueError(f"{name}: AVG stations fall outside STN chainage")
        midpoint_topo = pd.DataFrame({"station": measured})
        for column in ("easting", "northing", "elevation"):
            midpoint_topo[column] = np.interp(
                measured, source.station, source[column]
            )
        avg.add_topography(midpoint_topo, epsg=32649)
    avg.topo.convert_coords(to="ll", inplace=True)
    collection = AVGtoEDI().transform(avg)
    return avg, collection


def summarize(name: str, avg, collection) -> str:
    """Return a stable, compact summary suitable for captured documentation."""
    frame = avg.topo.frame
    nfreq = np.array([len(ed.Z.freq) for ed in collection])
    first = collection[0].get_section("head")
    last = collection[-1].get_section("head")
    return (
        f"{name}: {len(collection)} stations; frequencies/station "
        f"{nfreq.min()}-{nfreq.max()}; elevation "
        f"{frame.elevation.min():.1f}-{frame.elevation.max():.1f} m; "
        f"first=({first.lat:.5f}, {first.long:.5f}); "
        f"last=({last.lat:.5f}, {last.long:.5f})"
    )


def make_overview(lines) -> Path:
    """Plot both topographic profiles and raw apparent-resistivity grids."""
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), constrained_layout=True)
    for row, (name, avg, collection) in enumerate(lines):
        topo = avg.topo.frame.sort_values("station")
        axes[row, 0].plot(topo.station, topo.elevation, "o-", ms=3.2, lw=1.2)
        axes[row, 0].set(
            title=f"{name}: measured topography",
            xlabel="Station chainage (m)",
            ylabel="Elevation (m)",
        )
        axes[row, 0].grid(alpha=0.25)

        stations = np.array([float(ed.station.lstrip("S")) for ed in collection])
        # AVG lines have one common frequency grid; sort for pcolormesh.
        freq = np.asarray(collection[0].Z.freq)
        order = np.argsort(freq)
        rho = np.array([np.asarray(ed.Z.res_xy)[order] for ed in collection])
        image = axes[row, 1].pcolormesh(
            stations,
            freq[order],
            np.log10(np.maximum(rho.T, np.finfo(float).tiny)),
            shading="auto",
            cmap="turbo",
        )
        axes[row, 1].set_yscale("log")
        # Keep high frequencies at the top: they are generally more
        # sensitive to shallow structure than the low-frequency response.
        axes[row, 1].set(
            title=f"{name}: raw Ex-Hy response",
            xlabel="Station chainage (m)",
            ylabel="Frequency (Hz)",
        )
        fig.colorbar(image, ax=axes[row, 1], label=r"$\log_{10}\rho_a$ (Zonge field units)")

    target = OUT / "k1_k2_data_overview.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)
    return target


def _pivot(frame, value):
    """Return station, frequency and a sorted station-by-frequency matrix."""
    table = frame.pivot_table(index="station", columns="freq", values=value)
    table = table.sort_index().sort_index(axis=1)
    return (
        table.index.to_numpy(dtype=float),
        table.columns.to_numpy(dtype=float),
        table.to_numpy(dtype=float),
    )


def make_native_qc(lines) -> Path:
    """Plot AVG-native resistivity and phase uncertainties for both lines."""
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), constrained_layout=True)
    for row, (name, avg, _collection) in enumerate(lines):
        for col, (field, title, vmax) in enumerate(
            (("pc_rho", r"apparent-resistivity error (%)", 100.0),
             ("s_phz", "phase scatter (mrad)", 300.0))
        ):
            station, freq, values = _pivot(avg.df, field)
            shown = np.clip(values, 0.0, vmax)
            pc = axes[row, col].pcolormesh(
                station, freq, shown.T, shading="auto", cmap="magma", vmin=0, vmax=vmax
            )
            axes[row, col].set_yscale("log")
            # Keep high frequencies at the top, consistent with make_overview:
            # they are generally more sensitive to shallow structure.
            axes[row, col].set(
                title=f"{name}: {title}", xlabel="Station chainage (m)",
                ylabel="Frequency (Hz)",
            )
            fig.colorbar(pc, ax=axes[row, col], label=f"clipped at {vmax:g}")
    target = OUT / "k1_k2_native_avg_qc.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)
    return target


def make_static_shift_diagnostics(lines):
    """Execute TMA static-shift trials and plot factors plus response changes."""
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), constrained_layout=True)
    summaries = []
    for row, (name, avg_with_topo, _collection) in enumerate(lines):
        # Reload so this diagnostic never mutates the object used for raw EDI export.
        raw = AVG.from_file(DATA / f"{name}.AVG")
        before = raw.df.copy(deep=True)
        reference = float(before.freq.max())
        processor = ASTATIC().read(raw)
        shifts = processor.correct_static_shift(
            reference_freq=reference, filter_method="tma",
            window_size=5, update_components=True,
        )
        after = processor.avg.df

        ax = axes[row, 0]
        ax.semilogy(shifts.station, shifts.rho_original, ".-", label="observed")
        ax.semilogy(shifts.station, shifts.rho_smoothed, "-", lw=2, label="5-point TMA")
        ax.set(title=f"{name}: reference at {reference:g} Hz",
               xlabel="Station chainage (m)", ylabel=r"$\rho_a$ (AVG units)")
        ax.grid(alpha=0.25); ax.legend()

        ax = axes[row, 1]
        ax.semilogy(shifts.station, shifts.shift_factor, "o-", ms=3)
        ax.axhline(1.0, color="black", ls="--", lw=1)
        ax.set(title=f"{name}: trial static-shift factors",
               xlabel="Station chainage (m)", ylabel="multiplicative factor")
        ax.grid(alpha=0.25)

        _, _, rho0 = _pivot(before, "rho")
        _, _, rho1 = _pivot(after, "rho")
        delta = np.log10(np.maximum(rho1, np.finfo(float).tiny)) - np.log10(
            np.maximum(rho0, np.finfo(float).tiny)
        )
        summaries.append(
            (name, reference, float(shifts.shift_factor.min()),
             float(shifts.shift_factor.max()), float(np.nanmedian(shifts.shift_factor)),
             float(np.nanmedian(np.abs(delta))))
        )
    target = OUT / "k1_k2_static_shift_trial.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)
    return target, summaries


def make_ai_two_line(corrected_lines):
    """Train independent K1/K2 models and render topographic audit rows."""
    from pycsamt.agents import Inv2DAgent
    from pycsamt.topo import plot_topo_section

    np.random.seed(20260803)
    try:
        import torch
        torch.manual_seed(20260803)
    except ImportError:
        pass
    freqs = np.array([32.0, 128.0, 512.0, 2048.0, 8192.0])
    results = {}
    for index, name in enumerate(("K1", "K2")):
        np.random.seed(20260803 + index)
        try:
            torch.manual_seed(20260803 + index)
        except (ImportError, NameError):
            pass
        sites = corrected_lines[name]
        agent = Inv2DAgent(
            physics="mt2d", n_depth=8,
            n_stations_per_profile=len(sites), n_train_profiles=32,
            epochs=50, patience=8, depth_max=1600.0,
            station_spacing_m=50.0,
            correlation_length_x_m=(150.0, 500.0),
            correlation_length_z_m=(50.0, 250.0),
            log_resistivity_mean=2.0, log_resistivity_std=0.8,
            lambda_x=0.01, lambda_z=0.005, lambda_tv=0.002,
            mesh_safety_factor=3.0,
        )
        results[name] = agent.execute(
            {"sites": sites, "freqs": freqs, "topography": True}
        )

    fig = plt.figure(figsize=(15.0, 9.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=(1.25, 1.25, 0.9))
    for row, name in enumerate(("K1", "K2")):
        result = results[name]
        pred = np.asarray(result.data["pred_section"], dtype=float)
        chain = np.asarray(result.data["station_chainage_km"], dtype=float)
        elev = np.asarray(result.data["station_elevation_m"], dtype=float)
        depth = np.asarray(result.data["depths_km"], dtype=float)
        station_names = list(result.data["station_names"])
        if pred.shape[0] != chain.size:
            pred = pred.T
        label_step = max(1, len(station_names) // 8)
        sparse_labels = [
            station if i % label_step == 0 or i == len(station_names) - 1 else ""
            for i, station in enumerate(station_names)
        ]
        ax_model = fig.add_subplot(grid[row, :2])
        plot_topo_section(
            {"pred_rho": pred, "depths_km": depth,
             "station_names": station_names},
            ax=ax_model, elevation=elev, chainage=chain,
            station_names=sparse_labels, station_x=chain,
            topo_source="array", model_unit="km", depth_max=float(depth[-1]),
            cmap="turbo", colorbar=True, show_stations=True,
            show_station_names=True,
            title=f"{name}: AI inversion from corrected EDIs",
        )
        ax_audit = fig.add_subplot(grid[row, 2])
        history = result.data["inverter"]._history
        recovery = result.data["mt2d_recovery"]
        epochs = np.arange(1, len(history["train_loss"]) + 1)
        ax_audit.plot(epochs, history["train_loss"], "o-", ms=3,
                      label="training")
        ax_audit.plot(epochs, history["val_loss"], "s-", ms=3,
                      label="validation")
        ax_audit.set(title=f"{name}: learning and inversion audit",
                     xlabel="Epoch", ylabel="Loss")
        ax_audit.grid(alpha=0.25); ax_audit.legend(loc="upper right")
        ax_audit.text(
            0.03, 0.04,
            f"field RMS = {result.data['rms_global']:.3f}\n"
            f"held-out RMSE = {recovery['rmse']:.3f}\n"
            f"held-out R² = {recovery['r2']:.3f}\n"
            f"stopped at epoch {len(history['train_loss'])}/50\n"
            f"32 profiles · 8 depths · 5 frequencies",
            transform=ax_audit.transAxes, va="bottom",
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
        )
    target = OUT / "k1_k2_ai_inversion_audit.png"
    fig.savefig(target, dpi=180)
    plt.close(fig)
    return target, results


def _edi_rho_matrix(sites):
    """Return chainage, common frequency grid, and Ex-Hy rho from EDI Z."""
    from pycsamt.emtools._core import _iter_items

    items = [getattr(item, "edi", item) for item in _iter_items(sites)]
    chain = np.array([float(ed.station.lstrip("S")) for ed in items])
    freq = np.asarray(items[0].Z.freq, dtype=float)
    rho = np.array([
        0.2 * np.abs(np.asarray(ed.Z.z)[:, 0, 1]) ** 2 / np.maximum(ed.Z.freq, 1e-24)
        for ed in items
    ])
    return chain, freq, rho


def process_edi_collection(raw):
    """Apply the selected scalar Ex-Hy EDI conditioning sequence."""
    from pycsamt.emtools import correct_static_shift, hampel_filter_freq

    despiked = hampel_filter_freq(
        raw, win=2, nsig=3.0, on="z", domain="magphase", inplace=False
    )
    return correct_static_shift(
        despiked, window_m=500.0, spacing_m=50.0,
        comp="xy", inplace=False,
    )


def make_edi_processing(lines):
    """Run the scalar-safe EDI processing chain and export corrected EDIs."""
    from pycsamt.emtools import hampel_filter_freq
    from pycsamt.emtools._core import ensure_sites
    from pycsamt.site.export import write_sites

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), constrained_layout=True)
    summaries = []
    corrected_lines = {}
    for row, (name, _avg, raw) in enumerate(lines):
        despiked = hampel_filter_freq(
            raw, win=2, nsig=3.0, on="z", domain="magphase", inplace=False
        )
        corrected = process_edi_collection(raw)
        corrected_lines[name] = corrected
        chain, freq, rho0 = _edi_rho_matrix(raw)
        _, _, rho1 = _edi_rho_matrix(despiked)
        _, _, rho2 = _edi_rho_matrix(corrected)
        tiny = np.finfo(float).tiny
        d_hampel = np.log10(np.maximum(rho1, tiny)) - np.log10(np.maximum(rho0, tiny))
        d_static = np.log10(np.maximum(rho2, tiny)) - np.log10(np.maximum(rho1, tiny))
        for col, (delta, title) in enumerate(
            ((d_hampel, "Hampel frequency despiking"),
             (d_static, "500 m Hanning static shift"))
        ):
            lim = max(0.05, float(np.nanpercentile(np.abs(delta), 98)))
            pc = axes[row, col].pcolormesh(
                chain, freq, delta.T, shading="auto", cmap="RdBu_r",
                vmin=-lim, vmax=lim,
            )
            axes[row, col].set_yscale("log")
            # Keep high frequencies at the top, consistent with make_overview:
            # they are generally more sensitive to shallow structure.
            axes[row, col].set(title=f"{name}: {title}",
                               xlabel="Station chainage (m)", ylabel="Frequency (Hz)")
            fig.colorbar(pc, ax=axes[row, col], label=r"$\Delta\log_{10}\rho_a$")

        outdir = ROOT / "results" / "zonge_avg_tutorial" / "edi_corrected" / name
        paths = write_sites(
            corrected, outdir, exist_ok=True,
            manifest_csv=outdir / "manifest.csv",
        )
        reloaded = ensure_sites(outdir, recursive=False).ordered()
        summaries.append(
            (name, int(np.sum(np.abs(d_hampel) > 1e-12)),
             float(np.nanmedian(np.abs(d_hampel))),
             float(np.nanmedian(np.abs(d_static))), len(paths), len(reloaded))
        )
    target = OUT / "k1_k2_edi_processing_changes.png"
    fig.savefig(target, dpi=180); plt.close(fig)
    return target, corrected_lines, summaries


def make_response_panels(lines, corrected_lines):
    """Use the public API style for three raw and corrected stations per line."""
    from pycsamt.emtools import plot_raw_sites_1d

    outputs = []
    for name, _avg, raw in lines:
        names = [ed.station for ed in raw]
        picks = [names[0], names[len(names) // 2], names[-1]]
        for label, sites, is_raw in (
            ("raw", raw, True), ("corrected", corrected_lines[name], False)
        ):
            fig = plot_raw_sites_1d(
                sites, stations=picks, components=("xy",), raw=is_raw,
                show_error_bars=True, ncols_groups=3,
                figsize_scale=(4.2, 4.2),
                title_group_fmt=f"{name} {{station}} — {label}",
            )
            target = OUT / f"{name.lower()}_three_station_{label}_rho_phase.png"
            fig.subplots_adjust(bottom=0.18, top=0.90)
            fig.savefig(target, dpi=180)
            plt.close(fig); outputs.append(target)
    return outputs


def make_tensor_prerequisite_audit(lines):
    """Plot measured tensor-component availability and geographic line bearing."""
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1), constrained_layout=True)
    summaries = []
    for ax, (name, avg, collection) in zip(axes, lines):
        z = np.asarray(collection[0].Z.z)
        availability = np.any(np.isfinite(z) & (np.abs(z) > 0), axis=0).astype(int)
        ax.imshow(availability, vmin=0, vmax=1, cmap="Greys", origin="upper")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, "measured" if availability[i, j] else "absent",
                        ha="center", va="center",
                        color="white" if availability[i, j] else "black", fontsize=9)
        ax.set_xticks([0, 1], ["Ex", "Ey"]); ax.set_yticks([0, 1], ["Hx", "Hy"])
        frame = avg.topo.frame.sort_values("station")
        lat1, lat2 = np.deg2rad([frame.latitude.iloc[0], frame.latitude.iloc[-1]])
        dlon = np.deg2rad(frame.longitude.iloc[-1] - frame.longitude.iloc[0])
        bearing = np.degrees(np.arctan2(
            np.sin(dlon) * np.cos(lat2),
            np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon),
        )) % 360.0
        ax.set_title(f"{name}: tensor inputs\nprofile bearing {bearing:.1f}°")
        summaries.append((name, bearing, int(availability.sum())))
    target = OUT / "k1_k2_tensor_prerequisite_audit.png"
    fig.savefig(target, dpi=180); plt.close(fig)
    return target, summaries


def main() -> None:
    lines = []
    for name in ("K1", "K2"):
        avg, collection = load_line(name)
        print(summarize(name, avg, collection))
        lines.append((name, avg, collection))
    target = make_overview(lines)
    print(f"figure: {target.relative_to(ROOT).as_posix()}")
    target = make_native_qc(lines)
    print(f"figure: {target.relative_to(ROOT).as_posix()}")
    target, summaries = make_static_shift_diagnostics(lines)
    for name, ref, lo, hi, med, delta in summaries:
        print(
            f"{name} TMA trial: ref={ref:g} Hz; factor={lo:.3f}-{hi:.3f}; "
            f"median={med:.3f}; median |delta log10 rho|={delta:.3f}"
        )
    print(f"figure: {target.relative_to(ROOT).as_posix()}")
    target, _corrected, summaries = make_edi_processing(lines)
    for name, changed, hampel_med, static_med, written, reloaded in summaries:
        print(
            f"{name} EDI processing: Hampel changed={changed} cells; "
            f"median |Hampel delta|={hampel_med:.3f}; "
            f"median |static delta|={static_med:.3f}; "
            f"written/reloaded={written}/{reloaded}"
        )
    print(f"figure: {target.relative_to(ROOT).as_posix()}")
    for target in make_response_panels(lines, _corrected):
        print(f"figure: {target.relative_to(ROOT).as_posix()}")
    target, tensor_summaries = make_tensor_prerequisite_audit(lines)
    for name, bearing, ncomp in tensor_summaries:
        print(
            f"{name} tensor audit: measured components={ncomp}/4; "
            f"profile bearing={bearing:.1f} deg"
        )
    print(f"figure: {target.relative_to(ROOT).as_posix()}")
    target, ai_results = make_ai_two_line(_corrected)
    for name, result in ai_results.items():
        recovery = result.data["mt2d_recovery"]
        print(
            f"{name} AI: status={result.status}; RMS={result.data['rms_global']:.3f}; "
            f"held-out RMSE={recovery['rmse']:.3f}; R2={recovery['r2']:.3f}"
        )
    print(f"figure: {target.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
