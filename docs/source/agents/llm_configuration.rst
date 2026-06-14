.. _agents-llm-configuration:

Agent Configuration
===================

Agent configuration controls how pyCSAMT agents use LLM providers, API keys,
models, pricing tables, and session budgets.  The central object is
:data:`pycsamt.agents.AGENT_CONFIG`, a process-local singleton used by every
agent that inherits from :class:`pycsamt.agents.BaseAgent`.

The important design rule is simple: configure once, then instantiate agents.
Agents created without an explicit ``api_key`` inherit the active global
configuration automatically.

What Gets Configured
--------------------

The agent configuration stores:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Setting
     - Meaning
   * - ``provider``
     - Active LLM provider. Supported values are ``"claude"``,
       ``"openai"``, and ``"gemini"``.
   * - ``api_key``
     - Provider API key. Keys can be passed directly, stored per provider, or
       discovered from environment variables.
   * - ``model``
     - Optional model override. If omitted, pyCSAMT uses the provider default.
   * - ``custom rates``
     - Optional USD-per-token pricing overrides used for cost estimates.
   * - ``budget``
     - Optional session spend cap. Once the cap is reached, LLM calls raise
       :class:`pycsamt.agents.BudgetExceededError` before contacting the API.
   * - ``spent_usd``
     - Accumulated estimated LLM spend for the current Python session.

Provider Defaults
-----------------

When no model is supplied, pyCSAMT selects a provider-specific default:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Provider
     - Default model
     - Typical use
   * - ``"claude"``
     - ``"claude-sonnet-4-6"``
     - General workflow planning, interpretation, and report text.
   * - ``"openai"``
     - ``"gpt-4o"``
     - General workflow planning and interpretation.
   * - ``"gemini"``
     - ``"gemini-2.0-flash"``
     - Fast interpretation and lower-cost assistant tasks.

Minimal Global Configuration
----------------------------

Use :func:`pycsamt.agents.configure_agents` for the common case.

.. code-block:: python
   :linenos:

   from pycsamt.agents import configure_agents
   from pycsamt.agents import DataQCAgent, PhaseAnalysisAgent

   configure_agents(
       provider="claude",
       api_key="sk-ant-...",
       model="claude-sonnet-4-6",
   )

   qc_agent = DataQCAgent()
   phase_agent = PhaseAnalysisAgent()

   print(qc_agent.llm_provider)
   print(phase_agent.model)

Both agents inherit the configured Claude key and model because neither agent
was instantiated with an explicit ``api_key``.

Using ``AGENT_CONFIG`` Directly
-------------------------------

The singleton exposes the full configuration API.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG
   from pycsamt.agents import StaticShiftAgent

   AGENT_CONFIG.configure(
       provider="openai",
       api_key="sk-...",
       model="gpt-4o-mini",
   )

   agent = StaticShiftAgent()
   result = agent.execute({
       "path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_static_shift",
   })

   print(result.summary)
   print(result.cost_estimate_usd)

Configuration Resolution Rules
------------------------------

When an agent is instantiated, :class:`pycsamt.agents.BaseAgent` asks
``AGENT_CONFIG`` to resolve the effective provider, key, and model.

Resolution follows these rules:

1. If the agent receives an explicit ``api_key``, the agent uses its own
   ``llm_provider``, ``api_key``, and ``model``.  The global provider is not
   inherited.
2. If the agent does not receive an explicit key and its ``llm_provider`` is
   left at the default ``"claude"``, an active global provider is inherited.
3. The key is resolved as stored key, then environment variable, then
   ``None``.
4. The model is inherited from the global configuration only when the effective
   provider matches the active global provider.
5. If no key is found, the agent runs in no-LLM mode where supported.

This means global configuration works naturally for most workflows, while
individual agents can still opt into a different provider.

Per-Agent Override
------------------

Use an explicit key when one agent in a workflow should use a different
provider or model from the global default.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG
   from pycsamt.agents import DataQCAgent, ReportAgent

   AGENT_CONFIG.configure(
       provider="claude",
       api_key="sk-ant-...",
       model="claude-sonnet-4-6",
   )

   qc_agent = DataQCAgent()

   report_agent = ReportAgent(
       llm_provider="openai",
       api_key="sk-...",
       model="gpt-4o",
   )

   print(qc_agent.llm_provider)      # claude
   print(report_agent.llm_provider)  # openai

Environment Variable Fallback
-----------------------------

API keys do not need to be passed in code.  If no stored key exists for the
active provider, pyCSAMT checks provider-specific environment variables.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Provider
     - Environment variables, checked in order
   * - ``"claude"``
     - ``ANTHROPIC_API_KEY``, ``PYCSAMT_CLAUDE_API_KEY``
   * - ``"openai"``
     - ``OPENAI_API_KEY``, ``PYCSAMT_OPENAI_API_KEY``
   * - ``"gemini"``
     - ``GOOGLE_API_KEY``, ``GOOGLE_GENERATIVEAI_API_KEY``,
       ``PYCSAMT_GEMINI_API_KEY``

Example shell setup:

.. code-block:: bash
   :linenos:

   export ANTHROPIC_API_KEY="sk-ant-..."
   export OPENAI_API_KEY="sk-..."
   export GOOGLE_API_KEY="AIza..."

Then configure only the provider and optional model in Python:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, WorkflowOrchestratorAgent

   AGENT_CONFIG.switch("claude", model="claude-sonnet-4-6")

   agent = WorkflowOrchestratorAgent()
   print(agent.api_key is not None)

Multi-Provider Sessions
-----------------------

Use ``set_key`` to preload keys, then switch providers without passing keys
again.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG
   from pycsamt.agents import DataQCAgent, InterpretationAgent

   AGENT_CONFIG.set_key("claude", "sk-ant-...")
   AGENT_CONFIG.set_key("openai", "sk-...")
   AGENT_CONFIG.set_key("gemini", "AIza...")

   AGENT_CONFIG.switch("claude")
   qc_agent = DataQCAgent()

   AGENT_CONFIG.switch("openai", model="gpt-4o")
   interpretation_agent = InterpretationAgent()

   print(qc_agent.llm_provider)
   print(interpretation_agent.llm_provider)

Temporary Overrides
-------------------

Use ``AGENT_CONFIG.using(...)`` when a provider should apply only inside a
small block.  The original configuration is restored even if an exception is
raised.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, ReportAgent

   AGENT_CONFIG.configure(
       provider="claude",
       api_key="sk-ant-...",
       model="claude-sonnet-4-6",
   )

   with AGENT_CONFIG.using(
       provider="gemini",
       api_key="AIza...",
       model="gemini-2.0-flash",
   ):
       fast_report = ReportAgent()
       print(fast_report.llm_provider)  # gemini

   normal_report = ReportAgent()
   print(normal_report.llm_provider)    # claude

No-LLM Mode
-----------

No-LLM mode is useful for deterministic previews, CI checks, and workflows
where the processing agent can operate without interpretation text.

.. code-block:: python
   :linenos:

   from pycsamt.agents import reset_agents
   from pycsamt.agents import ContextInputAgent, DataQCAgent

   reset_agents()

   context = ContextInputAgent()
   qc_agent = DataQCAgent()

   parsed = context.execute({
       "request": "Load /data/WILLY_EDIs and run QC",
   })

   result = qc_agent.execute({
       "path": parsed["config"]["data_path"],
   })

   print(parsed.llm_interpretation)  # None when fallback parsing is used
   print(result.cost_estimate_usd)   # 0.0 when no LLM call was made

Not every agent has the same no-LLM behavior.  Processing agents generally keep
their deterministic processing path and omit LLM interpretation.  Agents whose
main purpose is language planning or report narration may return simpler
fallback output.

Budget Caps
-----------

Budgets are session-local.  The counter increases after successful LLM calls.
Before each new LLM call, pyCSAMT checks whether the cap has already been
reached.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, BudgetExceededError
   from pycsamt.agents import WorkflowOrchestratorAgent, format_cost

   AGENT_CONFIG.configure(
       provider="claude",
       api_key="sk-ant-...",
       model="claude-sonnet-4-6",
   )
   AGENT_CONFIG.set_budget(usd=1.00)

   agent = WorkflowOrchestratorAgent()

   try:
       result = agent.execute({
           "request": "Load WILLY data, QC, static shift, and report",
           "data_path": "/data/WILLY_EDIs",
           "output_dir": "/out/willy",
       })
   except BudgetExceededError as exc:
       print(exc)
       print(f"Spent: {format_cost(exc.spent_usd)}")
       print(f"Budget: {format_cost(exc.budget_usd)}")

   print(f"Session spent: {format_cost(AGENT_CONFIG.spent_usd)}")
   print(f"Remaining: {AGENT_CONFIG.remaining_usd}")

Reset only the counter:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.reset_budget()

Reset both the counter and the cap:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.reset_budget(cap=True)

Pricing And Cost Estimates
--------------------------

Built-in rates are stored as USD per 1,000,000 input or output tokens.  Cost
estimates use the effective table: custom override, built-in exact match,
built-in prefix match, then provider default.

Estimate one call:

.. code-block:: python
   :linenos:

   from pycsamt.agents import estimate_cost, format_cost

   cost = estimate_cost(
       provider="claude",
       model="claude-sonnet-4-6",
       input_tokens=500,
       output_tokens=150,
   )

   print(format_cost(cost))

Override or add rates:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.set_rate(
       "claude",
       "claude-opus-4-8",
       input=12.0,
       output=60.0,
   )

   AGENT_CONFIG.set_rate(
       "openai",
       "gpt-5-turbo",
       input=5.0,
       output=20.0,
   )

   print(AGENT_CONFIG.get_rate("claude", "claude-opus-4-8"))
   print(AGENT_CONFIG.list_rates("openai"))

Inspecting Active Configuration
-------------------------------

Use ``info()`` for display, logging, or debugging.  API keys are masked.

.. code-block:: python
   :linenos:

   from pprint import pprint
   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.configure(
       provider="claude",
       api_key="sk-ant-example",
       model="claude-sonnet-4-6",
   )
   AGENT_CONFIG.set_budget(usd=2.50)

   pprint(AGENT_CONFIG.info())

Example fields include ``provider``, ``model``, ``has_key``,
``api_key_masked``, ``key_source``, ``stored_providers``,
``custom_rate_models``, ``budget_usd``, ``spent_usd``, and
``remaining_usd``.

Resetting Configuration
-----------------------

Use :func:`pycsamt.agents.reset_agents` to clear the active global
configuration.

.. code-block:: python
   :linenos:

   from pycsamt.agents import reset_agents

   reset_agents()

Keep stored keys but clear the active provider, model, pricing overrides, and
budget:

.. code-block:: python
   :linenos:

   from pycsamt.agents import reset_agents

   reset_agents(keys=False)

Configuration Checklist
-----------------------

Before running an expensive workflow:

1. Choose the provider and model.
2. Confirm the key source with ``AGENT_CONFIG.info()``.
3. Set a budget with ``AGENT_CONFIG.set_budget(usd=...)``.
4. Run a dry-run or preview workflow where available.
5. Inspect ``AgentResult.cost_estimate_usd`` and ``AGENT_CONFIG.spent_usd``.

Related API
-----------

See :data:`pycsamt.agents.AGENT_CONFIG`,
:class:`pycsamt.agents.AgentConfig`,
:func:`pycsamt.agents.configure_agents`,
:func:`pycsamt.agents.reset_agents`,
:func:`pycsamt.agents.estimate_cost`,
:func:`pycsamt.agents.format_cost`, and
:class:`pycsamt.agents.BudgetExceededError`.
