r"""
Validate a workflow plan
========================

A :class:`~pycsamt.agents.WorkflowPlan` is the typed hand-off between parsing a
request and executing it.  Validating the plan *before* any files are touched
surfaces problems — an unknown workflow, a missing or non-existent data path —
while they are still cheap to fix.

This example builds one clean plan and one deliberately broken plan, checks
both, shows :func:`~pycsamt.agents.validate_workflow_plan` refusing the broken
one, and renders the two side by side as a validation report card.
"""

# sphinx_gallery_thumbnail_number = 1

# %%
# A plan that passes
# ------------------
# ``data_path`` must point at something that exists on disk, so this example
# uses a real temporary directory.  With a recognised ``workflow_type`` and a
# request string present, :meth:`WorkflowPlan.validation_errors
# <pycsamt.agents.WorkflowPlan.validation_errors>` comes back empty.

import tempfile

from pycsamt.agents import WorkflowPlan, validate_workflow_plan

existing_dir = tempfile.mkdtemp(prefix="pycsamt_survey_")

good = WorkflowPlan(
    request="QC the L22 survey and prepare a quality report",
    workflow_type="qc",
    data_path=existing_dir,
    expected_outputs=["qc_confidence.png", "qc_report.md"],
    provider="offline",
)

print("is_valid          :", good.is_valid())
print("validation_errors :", good.validation_errors())

# %%
# A plan that fails
# -----------------
# This plan names a workflow that does not exist and carries no data path.
# Called with ``raise_on_error=True``,
# :func:`~pycsamt.agents.validate_workflow_plan` raises a ``ValueError`` that
# lists every problem — the guard an application would rely on before spending
# time or money executing the chain.

bad = WorkflowPlan(
    request="Do the thing",
    workflow_type="mega_inversion",  # not a recognised workflow
    data_path="",  # nothing to load
    provider="offline",
)

try:
    validate_workflow_plan(bad, raise_on_error=True)
except ValueError as exc:
    print("Rejected as expected:\n")
    print(exc)

# %%
# A validation report card
# -------------------------
# Rendering both plans together contrasts a ready-to-run plan with one that
# needs attention.  Each plan lists its checks with a pass/fail marker; the
# failing plan spells out exactly what to fix.

import matplotlib.pyplot as plt


def plan_checks(plan):
    """Return (all_ok, [(label, ok, detail), ...]) for a plan."""
    errors = plan.validation_errors()
    joined = " ".join(errors).lower()
    has_path = bool(plan.data_path)
    checks = [
        ("request present", bool(plan.request), ""),
        (
            "known workflow",
            "unknown workflow" not in joined,
            plan.workflow_type,
        ),
        ("data path set", has_path, ""),
        ("data path exists", has_path and "does not exist" not in joined, ""),
    ]
    return not errors, checks, errors


fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
panels = [
    ("Ready to run", good, "#2e7d32"),
    ("Needs attention", bad, "#c62828"),
]

for ax, (title, plan, accent) in zip(axes, panels):
    ok_all, checks, errors = plan_checks(plan)
    ax.axis("off")
    ax.set_title(
        f"{title}\nworkflow_type = {plan.workflow_type!r}",
        fontsize=11,
        color=accent,
        fontweight="bold",
    )

    y = 0.82
    for label, ok, detail in checks:
        mark, mcol = ("PASS", "#2e7d32") if ok else ("FAIL", "#c62828")
        ax.text(
            0.04,
            y,
            mark,
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            color=mcol,
            family="monospace",
        )
        text = label + (f"  ({detail})" if detail else "")
        ax.text(
            0.24, y, text, transform=ax.transAxes, fontsize=9.5, color="0.2"
        )
        y -= 0.13

    if errors:
        ax.text(
            0.04,
            y - 0.02,
            "notes:",
            transform=ax.transAxes,
            fontsize=8.5,
            fontweight="bold",
            color="0.35",
        )
        y -= 0.13
        for err in errors:
            wrapped = err if len(err) <= 52 else err[:49] + "…"
            ax.text(
                0.06,
                y,
                f"• {wrapped}",
                transform=ax.transAxes,
                fontsize=8,
                color="0.35",
            )
            y -= 0.09

    # a coloured frame around each panel
    ax.add_patch(
        plt.Rectangle(
            (0.01, 0.01),
            0.98,
            0.98,
            transform=ax.transAxes,
            fill=False,
            edgecolor=accent,
            linewidth=1.5,
        )
    )

fig.suptitle("validate_workflow_plan — pre-flight check", fontsize=12)
fig.tight_layout()
