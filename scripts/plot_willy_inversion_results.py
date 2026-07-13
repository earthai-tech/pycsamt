"""Plot WILLY ModEM inversion results with pyCSAMT."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.models.modem.plot import (
    PlotMisfit,
    PlotModel3D,
    PlotPseudo,
    PlotResponse,
    PlotSection,
)
from pycsamt.models.modem.results import InversionResult

SOURCE_ROOT = Path(r"C:\Users\Daniel\Downloads\willy-inversion results")
RUNS = {
    "27-frequ-watex-data-02": SOURCE_ROOT
    / "27-frequ-watex-data-02"
    / "27-frequ-watex-data-02",
    "27-frequ-watex-data-03": SOURCE_ROOT
    / "27-frequ-watex-data-03"
    / "27-frequ-watex-data-03",
    "27-freque-watex-data-04": SOURCE_ROOT
    / "27-freque-watex-data-04"
    / "27-freque-watex-data-04",
}
OUT_ROOT = Path("willy_inversion_plots")
DPI = 180


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def _component_names(result: InversionResult) -> list[str]:
    if result.data_obs is None:
        return []
    names = {
        row[5] for block in result.data_obs.blocks for row in block["rows"]
    }
    preferred = ["TE", "TM", "ZXX", "ZXY", "ZYX", "ZYY"]
    return [name for name in preferred if name in names] + sorted(
        names - set(preferred)
    )


def _response_batches(
    result: InversionResult, batch_size: int = 4
) -> list[list[str]]:
    if result.data_obs is None:
        return []
    names = list(result.data_obs.site_names)
    return [
        names[i : i + batch_size] for i in range(0, len(names), batch_size)
    ]


def plot_run(name: str, workdir: Path) -> dict:
    result = InversionResult(workdir)
    out = OUT_ROOT / name
    meta = {
        "source": str(workdir),
        "mode": result.mode,
        "n_models": len(result.models),
        "n_log_records": int(result.n_iter),
        "final_rms": float(result.final_rms),
        "best_rms": float(result.best_rms),
        "best_iteration": int(
            result.iteration_numbers[result.rms_history.argmin()]
        )
        if result.rms_history.size
        else None,
        "n_sites": len(result.data_obs.site_names) if result.data_obs else 0,
        "n_periods": len(result.data_obs.periods) if result.data_obs else 0,
        "components": _component_names(result),
        "plots": [],
        "errors": [],
    }

    def record(kind: str, filename: str, fn) -> None:
        try:
            fig = fn()
            rel = Path(kind) / filename
            _save(fig, out / rel)
            meta["plots"].append(str(rel))
            print(f"saved {name}/{rel}")
        except Exception as exc:  # noqa: BLE001
            msg = f"{kind}/{filename}: {exc}"
            meta["errors"].append(msg)
            print(f"skipped {name}/{msg}")
            plt.close("all")

    record(
        "convergence",
        "rms_misfit.png",
        lambda: PlotMisfit(result=result).plot(),
    )

    for which in ("initial", "final"):
        record(
            "model_slices",
            f"{which}_depth_slices.png",
            lambda which=which: PlotModel3D(
                result=result,
                which=which,
                depths=[0, 100, 250, 500, 1000, 1500, 2500, 4000],
                rho_min=1.0,
                rho_max=5000.0,
                cmap="turbo",
                n_cols=4,
            ).plot(),
        )

    for direction in ("NS", "EW"):
        record(
            "sections",
            f"final_{direction.lower()}_center_section.png",
            lambda direction=direction: PlotSection(
                result=result,
                direction=direction,
                profile_offset=0.0,
                which="final",
                depth_max=5000.0,
                rho_min=1.0,
                rho_max=5000.0,
                cmap="turbo",
                station_tol=800.0,
                title=f"{name} final {direction} center section",
            ).plot(),
        )

    for component in _component_names(result):
        record(
            "pseudosections",
            f"{component.lower()}_observed_pseudosection.png",
            lambda component=component: PlotPseudo(
                result=result,
                component=component,
                rho_min=1.0,
                rho_max=5000.0,
                cmap="turbo",
            ).plot(),
        )

    for i, stations in enumerate(_response_batches(result), start=1):
        label = f"stations_{i:02d}_{stations[0]}_to_{stations[-1]}.png"
        record(
            "responses",
            label,
            lambda stations=stations: PlotResponse(
                result=result,
                stations=stations,
                max_stations=len(stations),
            ).plot(),
        )

    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )
    return meta


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name, workdir in RUNS.items():
        if not workdir.is_dir():
            summary[name] = {
                "source": str(workdir),
                "errors": ["missing source folder"],
            }
            continue
        summary[name] = plot_run(name, workdir)
    (OUT_ROOT / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"wrote {OUT_ROOT.resolve()}")


if __name__ == "__main__":
    main()
