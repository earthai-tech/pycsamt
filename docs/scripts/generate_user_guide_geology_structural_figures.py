"""Generate the structural-geology figures for the geology structural guide."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.geology import (  # noqa: E402
    FaultTrace,
    StructuralMeasurement,
    StructuralModel,
)

IMAGES = ROOT / "docs/source/images/user_guide/geology"

_SENSE_COLOR = {
    "normal": "#1F618D",
    "reverse": "#B9770E",
    "strike_slip": "#6C3483",
    "unknown": "#5D6D7E",
}


def make_structural_measurement_geometry() -> None:
    """Write the strike/dip-direction compass diagram.

    Shows the same ``bedding`` measurement used on the page (strike 45,
    dip direction 135, which ``dip_azimuth_ok``), the two +/-20 degree
    acceptance wedges centred on ``strike +/- 90`` that
    ``dip_direction_deg`` is checked against, and the rejected
    ``dip_direction_deg=200`` example that falls outside both wedges.
    """
    bedding = StructuralMeasurement(
        x=500.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
        dip_direction_deg=135.0,
    )
    tol = bedding.dip_direction_tolerance_deg
    strike = bedding.strike_deg
    dipdir_ok = bedding.dip_direction_deg
    dipdir_bad = 200.0

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_yticklabels([])
    ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
    ax.set_xticklabels(["N", "E", "S", "W"], fontsize=11)

    for center in (strike + 90.0, strike - 90.0):
        lo, hi = np.deg2rad(center - tol), np.deg2rad(center + tol)
        theta = np.linspace(lo, hi, 30)
        ax.fill_between(theta, 0, 1, color="#27AE60", alpha=0.18)

    for ang in (strike, strike + 180.0):
        ax.plot(
            [np.deg2rad(ang), np.deg2rad(ang)], [0, 1],
            color="#2C3E50", linewidth=2.5, solid_capstyle="round", zorder=3,
        )

    ax.annotate(
        "", xy=(np.deg2rad(dipdir_ok), 0.92), xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", color="#1F618D", linewidth=2.5, mutation_scale=22),
    )
    ax.annotate(
        "", xy=(np.deg2rad(dipdir_bad), 0.92), xytext=(0, 0),
        arrowprops=dict(arrowstyle="-|>", color="#C0392B", linewidth=2.5, mutation_scale=22),
    )

    ax.text(np.deg2rad(dipdir_ok), 1.08, f"dip_direction_deg={dipdir_ok:.0f}\n(accepted)",
            ha="center", va="center", fontsize=9, color="#1F618D")
    ax.text(np.deg2rad(dipdir_bad), 1.15, f"dip_direction_deg={dipdir_bad:.0f}\n(rejected)",
            ha="center", va="center", fontsize=9, color="#C0392B")
    ax.text(np.deg2rad(strike) + 0.05, 0.55, f"strike_deg={strike:.0f}",
            fontsize=9, color="#2C3E50")

    ax.set_title(
        "StructuralMeasurement geometry: strike, dip direction, and the\n"
        "+/-20 deg acceptance wedges around strike +/- 90",
        fontsize=11, pad=28,
    )

    fig.tight_layout()
    fig.savefig(IMAGES / "structural_measurement_geometry.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_structural_evidence_section() -> None:
    """Write the profile-section figure for a ``StructuralModel``.

    Draws the two bedding measurements and two fault traces from the
    "Collecting evidence along a profile" example as an actual cross
    section: fault traces as dipping line segments (apparent dip,
    tilted toward ``downthrown_side``) and bedding measurements as
    triangular strike/dip markers at the surface.
    """
    model = StructuralModel()
    model.add_planar(StructuralMeasurement(
        x=200.0, kind="bedding", strike_deg=40.0, dip_deg=25.0, dip_direction_deg=130.0,
    ))
    model.add_planar(StructuralMeasurement(
        x=650.0, kind="bedding", strike_deg=50.0, dip_deg=35.0, dip_direction_deg=140.0,
    ))
    fault1 = FaultTrace(
        x=500.0, dip_deg=70.0, downthrown_side="right", sense="normal",
        throw_m=12.0, evidence="resistivity offset",
    )
    model.add_fault(fault1)
    model.add_fault(FaultTrace(
        x=900.0, dip_deg=60.0, downthrown_side="left", evidence="surface mapping",
    ))

    fig, ax = plt.subplots(figsize=(9, 5))

    z_max = 80.0
    for f in model.faults:
        direction = 1.0 if f.downthrown_side == "right" else -1.0
        z = np.array([0.0, z_max])
        dx = direction * z / np.tan(np.deg2rad(f.dip_deg))
        x_line = f.x + dx
        color = _SENSE_COLOR[f.sense]
        ax.plot(x_line, z, color=color, linewidth=2.5, zorder=3)
        ax.annotate(
            "down", xy=(x_line[1] + 25 * direction, z[1] * 0.35),
            xytext=(x_line[1], z[1] * 0.35),
            color=color, fontsize=8, va="center",
            arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2),
        )
        label = f"{f.sense}\ndip {f.dip_deg:.0f} deg"
        if f.throw_m is not None:
            label += f"\nthrow {f.throw_m:.0f} m"
        ax.text(f.x, -6.0, label, ha="center", va="bottom", fontsize=8, color=color)

    for m in model.planar:
        ax.plot(m.x, 0.0, marker="v", color="#117864", markersize=10, zorder=4)
        ax.text(
            m.x, -6.0, f"{m.kind}\n{m.strike_deg:.0f}/{m.dip_deg:.0f}->{m.dip_direction_deg:.0f}",
            ha="center", va="bottom", fontsize=8, color="#117864",
        )

    ax.set_xlim(0, 1050)
    ax.set_ylim(z_max, -18)
    ax.set_xlabel("Profile position x (m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title("Structural evidence along the profile (StructuralModel)")
    ax.grid(alpha=0.25)
    ax.axhline(0.0, color="0.3", linewidth=0.8)

    handles = [plt.Line2D([0], [0], color=c, lw=2.5, label=s) for s, c in _SENSE_COLOR.items()]
    handles.append(plt.Line2D([0], [0], marker="v", color="#117864", linestyle="none",
                              markersize=9, label="bedding measurement"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              fontsize=8, ncol=4, frameon=False)

    fig.tight_layout()
    fig.savefig(IMAGES / "structural_evidence_section.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_structural_measurement_geometry()
    make_structural_evidence_section()
