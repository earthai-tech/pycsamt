"""
pycsamt.agents.model_zoo_agent
===============================

:class:`ModelZooAgent` — Browse, download, and deploy pre-trained EM inverters.

The model zoo is hosted at
`github.com/earthai-tech/pycsamt-models <https://github.com/earthai-tech/pycsamt-models>`_
and managed by :mod:`pycsamt.ai._zoo`.  Checkpoints are cached locally in
``~/.pycsamt/model_zoo/`` (overridden via ``PYCSAMT_MODEL_CACHE``).

Three operations
----------------
``"list"``  (default)
    Print the full registry — model names, architectures, descriptions.
    No network access required.

``"download"``
    Download a named checkpoint to the local cache and return its path.

``"predict"``
    Download → load into the matching :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`
    → predict on observed sites and plot the resistivity section.
    This is the Phase 5 *model zoo* shortcut — users get a production-quality
    inverter in one call, without training from scratch.

Pre-trained model naming convention
------------------------------------
``<method>-<arch>-<n_layers>layer-v<version>``

Examples: ``mt1d-resnet-5layer-v1``, ``mt1d-cnn-5layer-v1``,
``csamt1d-resnet-5layer-v1``, ``tem1d-fcn-5layer-v1``.

Notes
-----
* Weights are released in Phase 5.  Until then, :func:`download_checkpoint`
  raises a :exc:`RuntimeError` explaining the situation.
* The ``"predict"`` action falls back gracefully to on-the-fly training when
  the checkpoint is unavailable.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in pre-trained geophysical AI models.
Given a model zoo query result, write 2-3 sentences that:
1. Describe which pre-trained model was used and its provenance (architecture, training data).
2. Comment on the prediction quality (RMS, reliability) relative to the expected use case.
3. Recommend whether the user should fine-tune on their own data or use the pre-trained weights directly.
Reply in plain English.
"""


class ModelZooAgent(BaseAgent):
    """Browse, download, and run pre-trained EM inverters from the model zoo.

    Parameters
    ----------
    api_key, model, llm_provider : str
    cache_dir : str or None
        Override default cache ``~/.pycsamt/model_zoo/``.
    force_download : bool
        Re-download even if cached (default False).

    Input keys
    ----------
    ``action`` : str
        ``"list"`` (default), ``"download"``, or ``"predict"``.
    ``model_name`` : str
        Required for ``"download"`` and ``"predict"`` actions.
        E.g. ``"mt1d-resnet-5layer-v1"``.
    ``sites`` / ``path`` : Sites or str
        Required for ``"predict"`` action.
    ``output_dir`` : str, optional

    Output data keys
    ----------------
    ``action``         str — which action was performed
    ``models``         dict — full registry (action="list")
    ``checkpoint_path``str — local path (action="download"/"predict")
    ``model_info``     dict — zoo metadata for the requested model
    ``predictions``    dict — {station: ndarray}  (action="predict")
    ``rms_global``     float  (action="predict")
    ``figures``        dict
    ``figure_paths``   dict

    Examples
    --------
    List available models::

        agent = ModelZooAgent()
        r = agent.execute({"action": "list"})
        for name, desc in r["models"].items():
            print(name, "—", desc[:60])

    Download a checkpoint::

        r = agent.execute(
            {"action": "download", "model_name": "mt1d-resnet-5layer-v1"}
        )
        print(r["checkpoint_path"])

    Predict on observed sites (fine-tune skipped if weights unavailable)::

        r = agent.execute(
            {
                "action": "predict",
                "model_name": "mt1d-resnet-5layer-v1",
                "path": "/data/WILLY_EDIs",
                "output_dir": "/out/zoo_predict",
            }
        )
        print(r["rms_global"])
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        cache_dir: str | None = None,
        force_download: bool = False,
    ) -> None:
        super().__init__(
            "ModelZooAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.cache_dir = cache_dir
        self.force_download = force_download

    # ── public ────────────────────────────────────────────────────────────────

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()

        action = str(input_data.get("action", "list")).lower()
        model_name = str(input_data.get("model_name", ""))
        output_dir = input_data.get("output_dir")

        if action == "list":
            return self._action_list(t0)

        if not model_name:
            return AgentResult.failed(
                f"action={action!r} requires 'model_name'.",
                hint="Set model_name to one of the keys returned by action='list'.",
                elapsed=time.time() - t0,
            )

        if action == "download":
            return self._action_download(model_name, t0)

        if action == "predict":
            return self._action_predict(input_data, model_name, output_dir, t0)

        return AgentResult.failed(
            f"Unknown action {action!r}. Use 'list', 'download', or 'predict'.",
            elapsed=time.time() - t0,
        )

    # ── action handlers ───────────────────────────────────────────────────────

    def _action_list(self, t0: float) -> AgentResult:
        from ..ai._zoo import _MODEL_ZOO, list_pretrained

        models = list_pretrained()
        rows = []
        for name, desc in models.items():
            info = _MODEL_ZOO[name]
            rows.append(
                {
                    "name": name,
                    "arch": info.get("arch", "?"),
                    "n_layers": info.get("n_layers", "?"),
                    "solver": info.get("solver", "?"),
                    "description": desc,
                }
            )

        summary_lines = [
            f"  {r['name']:<35s}  {r['description'][:60]}" for r in rows
        ]
        (
            f"Model zoo — {len(models)} pre-trained models available:\n"
            + "\n".join(summary_lines)
        )

        return AgentResult(
            status="success",
            summary=f"{len(models)} pre-trained models in zoo.",
            data={
                "action": "list",
                "models": models,
                "details": rows,
                "figures": {},
                "figure_paths": {},
            },
            warnings=[],
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=0.0,
        )

    def _action_download(self, model_name: str, t0: float) -> AgentResult:
        from ..ai._zoo import (
            download_checkpoint,
            get_pretrained_info,
        )

        warnings: list[str] = []
        ckpt_path = None

        try:
            info = get_pretrained_info(model_name)
        except KeyError as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        try:
            ckpt_path = download_checkpoint(
                model_name,
                cache_dir=self.cache_dir,
                force=self.force_download,
                verbose=True,
            )
        except RuntimeError as exc:
            warnings.append(str(exc))
            warnings.append(
                "Weights not yet released — they will be available in the "
                "public repository when Phase 5 weights are published."
            )

        return AgentResult(
            status="success" if ckpt_path else "needs_review",
            summary=(
                f"Model '{model_name}' cached at {ckpt_path}."
                if ckpt_path
                else f"Model '{model_name}' not yet available for download."
            ),
            data={
                "action": "download",
                "model_name": model_name,
                "model_info": info,
                "checkpoint_path": str(ckpt_path) if ckpt_path else None,
                "figures": {},
                "figure_paths": {},
            },
            warnings=warnings,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=0.0,
        )

    def _action_predict(
        self,
        input_data: dict[str, Any],
        model_name: str,
        output_dir: str | None,
        t0: float,
    ) -> AgentResult:
        from ..ai._zoo import (
            download_checkpoint,
            get_pretrained_info,
        )

        warnings: list[str] = []
        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        import os

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # ── zoo metadata ──────────────────────────────────────────────────────
        try:
            info = get_pretrained_info(model_name)
        except KeyError as exc:
            return AgentResult.failed(str(exc), elapsed=time.time() - t0)

        arch = info.get("arch", "resnet")
        n_layers = int(info.get("n_layers", 5))
        n_freqs = int(info.get("n_freqs") or 40)

        # ── try to download checkpoint ────────────────────────────────────────
        ckpt_path: str | None = None
        try:
            p = download_checkpoint(
                model_name,
                cache_dir=self.cache_dir,
                force=self.force_download,
                verbose=False,
            )
            ckpt_path = str(p)
        except RuntimeError as exc:
            warnings.append(
                f"Checkpoint unavailable: {exc}  "
                "Falling back to on-the-fly training."
            )

        # ── build AIInversionAgent and run ────────────────────────────────────
        from .ai_inversion import AIInversionAgent

        freqs = np.logspace(-4, 3, n_freqs)

        agent = AIInversionAgent(
            api_key=self.api_key,
            model=self.model,
            llm_provider=self.llm_provider,
            arch=arch,
            n_layers=n_layers,
            n_train_samples=2_000,
            epochs=30,
            freqs=freqs,
            pretrained=ckpt_path,
        )
        result = agent.execute(
            {
                "sites": input_data.get("sites") or input_data.get("path"),
                "output_dir": output_dir,
            }
        )

        # relay figures from inner agent
        figures = result.get("figures", {})
        fig_paths = result.get("figure_paths", {})
        predictions = result.get("predictions", {})
        rms_global = result.get("rms_global", float("nan"))

        # ── LLM interpretation ────────────────────────────────────────────────
        interp: str | None = None
        if self.api_key and predictions:
            prompt = (
                f"Model zoo prediction summary:\n"
                f"  Model: {model_name}\n"
                f"  Architecture: {arch}, layers: {n_layers}\n"
                f"  Stations: {len(predictions)}\n"
                f"  Global RMS: {rms_global:.3f}\n"
                f"  Pretrained checkpoint: {'yes' if ckpt_path else 'no (trained fresh)'}\n\n"
                "Evaluate the pre-trained model performance and recommend fine-tuning strategy."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        return AgentResult(
            status=result.status,
            summary=(
                f"Zoo predict ({model_name}): {len(predictions)} stations. "
                f"RMS {rms_global:.3f}. {len(figures)} figures."
            ),
            data={
                "action": "predict",
                "model_name": model_name,
                "model_info": info,
                "checkpoint_path": ckpt_path,
                "predictions": predictions,
                "rms_global": rms_global,
                "inverter": result.get("inverter"),
                "figures": figures,
                "figure_paths": fig_paths,
            },
            warnings=warnings + result.warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost + result.cost_estimate_usd,
        )


__all__ = ["ModelZooAgent"]
