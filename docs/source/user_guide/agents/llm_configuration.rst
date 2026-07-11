.. _agents-llm-configuration:

Agent And LLM Configuration
===========================

The agent layer can run in two modes:

* **deterministic mode**, where agents use pycsamt algorithms, rules, and
  validators without contacting an external LLM;
* **LLM-assisted mode**, where selected agents add natural-language parsing,
  interpretation, reporting, planning, and explanation on top of the same
  pycsamt computations.

The important design point is that LLM configuration is centralized.  You
configure the active provider once through
:data:`pycsamt.agents.AGENT_CONFIG`, and every agent created afterwards
inherits the provider, model, API key, pricing table, and optional session
budget.

This page should be read before using the agent catalogue.  It explains what
is configured, how pycsamt finds credentials, how to keep workflows reproducible,
and how to avoid surprising costs.


Quick decision guide
--------------------

Use the table below to choose the right setup for the situation.

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Situation
     - Recommended setup
     - Why
   * - CI, tests, tutorials, or offline work
     - Do not configure an API key.
     - Agents still run their deterministic pycsamt logic and skip LLM calls.
   * - Personal workstation
     - Store the key in an environment variable and call
       :meth:`AGENT_CONFIG.switch <pycsamt.api.agents.AgentConfig.switch>`.
     - Secrets stay outside notebooks, scripts, and documentation.
   * - Script with an explicit runtime secret
     - Use :func:`pycsamt.agents.configure_agents` or
       :meth:`AGENT_CONFIG.configure <pycsamt.api.agents.AgentConfig.configure>`.
     - The script is explicit about the provider and model for the session.
   * - One agent needs a different model
     - Pass ``llm_provider``, ``model``, and optionally ``api_key`` to that
       agent constructor.
     - Only that agent changes; the global configuration remains intact.
   * - A short experiment needs another provider
     - Use
       :meth:`AGENT_CONFIG.using <pycsamt.api.agents.AgentConfig.using>`
       as a context manager.
     - The original session configuration is restored automatically.
   * - You need a hard cost limit
     - Use
       :meth:`AGENT_CONFIG.set_budget <pycsamt.api.agents.AgentConfig.set_budget>`.
     - pycsamt raises
       :class:`~pycsamt.api.agents.BudgetExceededError` before the next LLM
       call once the cap is reached.


Configuration objects
---------------------

The public configuration API is exposed from :mod:`pycsamt.agents`:

.. code-block:: python
   :linenos:

   from pycsamt.agents import (
       AGENT_CONFIG,
       AgentConfig,
       BudgetExceededError,
       configure_agents,
       reset_agents,
   )

``AGENT_CONFIG`` is a process-local singleton, similar in spirit to the
plotting, section, style, and station-rendering API objects in pycsamt.  It is
not written to disk and is not shared between Python processes.  Configure it
at the beginning of a notebook, script, CLI command, or application session.

``configure_agents(...)`` and ``reset_agents(...)`` are convenience wrappers
around ``AGENT_CONFIG.configure(...)`` and ``AGENT_CONFIG.reset(...)``.


What gets configured
--------------------

.. list-table::
   :header-rows: 1
   :widths: 22 43 35

   * - Field
     - Meaning
     - Used by
   * - ``provider``
     - Active LLM provider: ``"claude"``, ``"openai"``, or ``"gemini"``.
     - Every new agent unless it is explicitly overridden.
   * - ``api_key``
     - Credential resolved from an explicit key or an environment variable.
     - :meth:`BaseAgent.query_llm <pycsamt.agents.BaseAgent.query_llm>`.
   * - ``model``
     - Model identifier for the active provider.
     - All LLM calls made by the agent.
   * - ``custom rates``
     - Optional USD-per-million-token rates for models not in the built-in
       table, or for rates you want to override.
     - Cost estimates and budget accounting.
   * - ``budget_usd``
     - Optional hard session cap.
     - Checked before every LLM call.
   * - ``spent_usd``
     - Accumulated estimated LLM spend in the current Python process.
     - Reports, dashboards, logs, and budget checks.

The API key itself is never returned in full by
:meth:`AGENT_CONFIG.info <pycsamt.api.agents.AgentConfig.info>`; the display
value is masked to the last four characters.


Supported providers and defaults
--------------------------------

pycsamt v2 supports three LLM provider families.  If you select a provider but
do not pass a model, pycsamt uses the default listed here.

.. list-table::
   :header-rows: 1
   :widths: 20 32 48

   * - Provider
     - Default model
     - Environment variables checked
   * - ``"claude"``
     - ``claude-sonnet-4-6``
     - ``ANTHROPIC_API_KEY``, ``PYCSAMT_CLAUDE_API_KEY``
   * - ``"openai"``
     - ``gpt-4o``
     - ``OPENAI_API_KEY``, ``PYCSAMT_OPENAI_API_KEY``
   * - ``"gemini"``
     - ``gemini-2.0-flash``
     - ``GOOGLE_API_KEY``, ``GOOGLE_GENERATIVEAI_API_KEY``,
       ``PYCSAMT_GEMINI_API_KEY``

Install the matching client package only for the provider you intend to use:

.. code-block:: bash
   :linenos:

   pip install anthropic
   pip install openai
   pip install google-generativeai

The imports are lazy.  Importing :mod:`pycsamt.agents` does not require all
three client libraries to be installed.


Recommended setup with environment variables
--------------------------------------------

For day-to-day use, keep secrets outside code and configure the active
provider from the environment.

.. code-block:: bash
   :linenos:

   export ANTHROPIC_API_KEY="sk-ant-..."

Then activate the provider in Python:

.. code-block:: python
   :linenos:

   from pprint import pprint

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.switch("claude", model="claude-sonnet-4-6")

   pprint(AGENT_CONFIG.info())

``switch(...)`` is the right call when the key is already available through
the environment or was previously registered with ``set_key(...)``.  It selects
the provider and optional model without embedding the key in the script.

The resulting information dictionary is safe to log:

.. code-block:: python
   :linenos:

   {
       "provider": "claude",
       "model": "claude-sonnet-4-6",
       "has_key": True,
       "api_key_masked": "...1a2b",
       "key_source": "env",
       "stored_providers": [],
       "custom_rate_models": {},
       "budget_usd": None,
       "spent_usd": 0.0,
       "remaining_usd": None,
   }


Explicit Python setup
---------------------

Use explicit setup when the key is supplied by a secret manager, CLI option,
or application settings layer.

.. code-block:: python
   :linenos:

   from pycsamt.agents import configure_agents

   configure_agents(
       provider="openai",
       api_key=runtime_secret,
       model="gpt-4o-mini",
   )

The same operation can be written with the singleton:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.configure(
       provider="claude",
       api_key=runtime_secret,
       model="claude-sonnet-4-6",
   )

Unlike ``switch(...)``, ``configure(...)`` requires an ``api_key`` argument.
That makes it explicit when code is placing a secret into the in-memory
configuration object.


How agents resolve LLM settings
-------------------------------

Every agent inherits from :class:`~pycsamt.agents.BaseAgent`.  During
initialization, ``BaseAgent`` calls ``AGENT_CONFIG.resolve(...)`` to decide
which provider, key, and model the new agent should use.

The resolution rules are:

1. If the agent constructor receives an explicit ``api_key``, that agent uses
   the constructor's ``llm_provider``, ``api_key``, and ``model``.  The global
   configuration is ignored for that agent.
2. If the constructor does not receive an ``api_key`` and the constructor uses
   the default provider ``"claude"``, the agent inherits the active global
   provider when one has been configured.
3. The key for the effective provider is resolved as: stored key, then
   provider-specific environment variables, then ``None``.
4. The model comes from the constructor if supplied.  Otherwise, when the
   effective provider matches the active global provider, the global model is
   used.  If no model is available, pycsamt uses the provider default.
5. If no key is found, the agent still runs but skips LLM calls.  Returned
   ``AgentResult.llm_interpretation`` values will be ``None`` and LLM cost
   remains zero.

This behavior keeps simple workflows concise:

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, DataQCAgent, PhaseAnalysisAgent

   AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   qc = DataQCAgent()              # inherits OpenAI + gpt-4o-mini
   pt = PhaseAnalysisAgent()       # same provider and model


Per-agent overrides
-------------------

Use a per-agent override when a specific task benefits from a different model
or provider.  For example, you might keep the global configuration on a fast
model, but ask the report agent to use a stronger model for final prose.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, DataQCAgent, ReportAgent

   AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   qc_agent = DataQCAgent()

   report_agent = ReportAgent(
       llm_provider="claude",
       model="claude-sonnet-4-6",
   )

In this example, ``report_agent`` still needs a Claude key.  That key can come
from ``ANTHROPIC_API_KEY``, ``PYCSAMT_CLAUDE_API_KEY``, or from
``AGENT_CONFIG.set_key("claude", "...")``.

If you pass an explicit key to one agent, only that agent uses it:

.. code-block:: python
   :linenos:

   report_agent = ReportAgent(
       llm_provider="gemini",
       api_key=gemini_runtime_secret,
       model="gemini-2.0-flash",
   )


Multi-provider sessions
-----------------------

You can register keys for several providers and switch between them during the
same Python session.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.set_key("claude", claude_key)
   AGENT_CONFIG.set_key("openai", openai_key)
   AGENT_CONFIG.set_key("gemini", gemini_key)

   AGENT_CONFIG.switch("claude")
   # agents created here use Claude

   AGENT_CONFIG.switch("openai", model="gpt-4o-mini")
   # agents created here use OpenAI

Existing agent instances keep the provider, key, and model they resolved at
construction time.  Switch before constructing the agents that should use the
new provider.


Temporary overrides
-------------------

Use ``AGENT_CONFIG.using(...)`` when a block of code should temporarily use a
different provider or model.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, InterpretationAgent

   AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   with AGENT_CONFIG.using(provider="claude", model="claude-sonnet-4-6"):
       interpretation = InterpretationAgent()
       result = interpretation.execute({"model": resistivity_model})

   # The OpenAI configuration is restored here.

The context manager also restores custom rates and budget state if an exception
occurs inside the block.


Running without an LLM
----------------------

LLM support is optional.  If no key is configured, agents skip calls to
external services and return deterministic pycsamt outputs.

.. code-block:: python
   :linenos:

   from pycsamt.agents import DataQCAgent, reset_agents

   reset_agents()

   result = DataQCAgent().execute({"path": "data/edis"})

   assert result.cost_estimate_usd == 0.0
   assert result.llm_interpretation is None

This mode is the recommended default for:

* unit tests;
* documentation examples;
* reproducible batch processing;
* secure environments where network calls are not allowed;
* workflows where only numerical outputs are required.

Agents that are mainly computational, such as QC, tensor rotation, EDI export,
pipeline execution, and inversion preparation, still provide useful results in
deterministic mode.  Agents that are mainly conversational or explanatory, such
as context parsing, interpretation, and report writing, may produce simpler
outputs or mark interpretation fields as unavailable.


Budget caps
-----------

Budget caps protect long workflows from unexpected LLM spend.  The cap is
session-local and is checked before every LLM call.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, BudgetExceededError, ReportAgent

   AGENT_CONFIG.switch("claude")
   AGENT_CONFIG.set_budget(usd=2.00)

   try:
       result = ReportAgent().execute({"survey": survey_summary})
   except BudgetExceededError as exc:
       print(f"Budget ${exc.budget_usd:.2f} reached.")
       print(f"Already spent ${exc.spent_usd:.4f}.")

Inspect the remaining budget at any point:

.. code-block:: python
   :linenos:

   print(AGENT_CONFIG.spent_usd)
   print(AGENT_CONFIG.remaining_usd)

Reset only the spend counter:

.. code-block:: python
   :linenos:

   AGENT_CONFIG.reset_budget()

Reset the spend counter and remove the cap:

.. code-block:: python
   :linenos:

   AGENT_CONFIG.reset_budget(cap=True)


Pricing and rate overrides
--------------------------

pycsamt estimates cost from provider rates expressed in USD per 1,000,000
tokens.  Built-in rates are used when the model is known.  If a model is new,
or if your organization wants to use a different accounting rate, add an
override.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   AGENT_CONFIG.set_rate(
       "openai",
       "gpt-5-turbo",
       input=5.00,
       output=20.00,
   )

   cost = AGENT_CONFIG.estimate_cost(
       "openai",
       "gpt-5-turbo",
       input_tokens=12_000,
       output_tokens=2_000,
   )

   print(f"Estimated cost: ${cost:.4f}")

List the effective rate table:

.. code-block:: python
   :linenos:

   from pprint import pprint

   pprint(AGENT_CONFIG.list_rates("claude"))

Rate resolution order is:

1. exact custom override;
2. exact built-in model rate;
3. prefix match in custom rates or built-in rates;
4. provider-level fallback rate.

Cost estimates are approximate.  They are intended for workflow accounting and
budget control, not invoice reconciliation.


Inspecting active configuration
-------------------------------

Use ``AGENT_CONFIG.info()`` when debugging notebooks, application sessions, or
agent workflows.

.. code-block:: python
   :linenos:

   from pprint import pprint
   from pycsamt.agents import AGENT_CONFIG

   pprint(AGENT_CONFIG.info())

Important fields are:

``provider``
    The active provider, or ``None`` if no provider has been selected.

``model``
    The active model, including provider defaults.

``has_key``
    Whether the active provider has a resolvable key.

``key_source``
    ``"explicit"`` for keys stored in ``AGENT_CONFIG``, ``"env"`` for
    environment variables, or ``"none"``.

``stored_providers``
    Providers with explicit keys stored in the current process.

``custom_rate_models``
    Models whose rates were overridden for this session.

``budget_usd``, ``spent_usd``, ``remaining_usd``
    Current budget state.


Resetting configuration
-----------------------

Reset everything, including stored keys:

.. code-block:: python
   :linenos:

   from pycsamt.agents import reset_agents

   reset_agents()

Reset the active provider/model, custom rates, and budget while keeping stored
keys for later switching:

.. code-block:: python
   :linenos:

   reset_agents(keys=False)

``reset_agents(...)`` does not clear environment variables.  If an environment
variable is still present, a later ``AGENT_CONFIG.switch(provider)`` can resolve
that key again.


Practical recipes
-----------------

Offline QC recipe
~~~~~~~~~~~~~~~~~

Use this pattern for tests, documentation, or secure offline processing.

.. code-block:: python
   :linenos:

   from pycsamt.agents import DataQCAgent, reset_agents

   reset_agents()

   qc = DataQCAgent()
   result = qc.execute({"path": "data/edis"})

   print(result.status)
   print(result.cost_estimate_usd)


Local analyst recipe
~~~~~~~~~~~~~~~~~~~~

Use environment variables, choose a provider, set a small budget, then create
agents normally.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, DataQCAgent, ReportAgent

   AGENT_CONFIG.switch("claude", model="claude-sonnet-4-6")
   AGENT_CONFIG.set_budget(usd=1.50)

   qc = DataQCAgent()
   report = ReportAgent()


Mixed-model recipe
~~~~~~~~~~~~~~~~~~

Use a fast global model for routine steps and a stronger model for the final
interpretation.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG, DataQCAgent, InterpretationAgent

   AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   qc_agent = DataQCAgent()

   with AGENT_CONFIG.using(provider="claude", model="claude-sonnet-4-6"):
       interpretation_agent = InterpretationAgent()


Application session recipe
~~~~~~~~~~~~~~~~~~~~~~~~~~

Desktop or web applications should configure agents once from their settings
layer, then create agents from controllers or service classes.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AGENT_CONFIG

   def configure_agent_session(settings):
       if settings.llm_enabled:
           AGENT_CONFIG.configure(
               provider=settings.provider,
               api_key=settings.api_key,
               model=settings.model,
           )
           if settings.budget_usd:
               AGENT_CONFIG.set_budget(usd=settings.budget_usd)
       else:
           AGENT_CONFIG.reset()


Troubleshooting
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Symptom
     - Likely cause
     - Fix
   * - ``llm_interpretation`` is ``None``
     - No API key was resolved, or the provider call failed.
     - Check ``AGENT_CONFIG.info()`` and confirm ``has_key`` is ``True``.
   * - Agent uses Claude even after configuring OpenAI
     - The agent was constructed before the global switch, or it received an
       explicit constructor override.
     - Switch provider before creating the agent; recreate existing agents.
   * - Environment key is ignored
     - The active provider does not match the environment variable you set.
     - Use ``AGENT_CONFIG.switch("openai")`` for OpenAI keys,
       ``switch("claude")`` for Claude keys, and ``switch("gemini")`` for
       Gemini keys.
   * - ``BudgetExceededError`` appears before a workflow starts
     - The previous workflow already used the budget cap.
     - Call ``AGENT_CONFIG.reset_budget()`` or raise the cap with
       ``AGENT_CONFIG.set_budget(usd=...)``.
   * - Cost estimate looks wrong for a new model
     - The model is not in the built-in rate table and pycsamt used a fallback.
     - Add an explicit rate with ``AGENT_CONFIG.set_rate(...)``.
   * - A provider client import fails
     - The optional provider package is not installed.
     - Install only the needed client library, for example ``pip install openai``.


Configuration checklist
-----------------------

Before running a large agent workflow:

* confirm whether the workflow should use LLM assistance or deterministic mode;
* select the provider and model explicitly;
* keep secrets in environment variables or an application settings layer;
* inspect ``AGENT_CONFIG.info()`` once before creating agents;
* set a budget cap for exploratory notebooks or user-facing applications;
* add custom rates for newly released models;
* construct agents after the desired provider/model has been selected;
* reset configuration at the end of tests to avoid leaking state between cases.


Related API
-----------

* :data:`pycsamt.agents.AGENT_CONFIG`
* :class:`pycsamt.api.agents.AgentConfig`
* :func:`pycsamt.agents.configure_agents`
* :func:`pycsamt.agents.reset_agents`
* :class:`pycsamt.api.agents.BudgetExceededError`
* :class:`pycsamt.agents.BaseAgent`
* :class:`pycsamt.agents.AgentResult`
:html_theme.sidebar_secondary.remove:
