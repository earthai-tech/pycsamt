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
:data:`pycsamt.agents.AGENT_CONFIG`, and every :term:`agent` created
afterwards inherits the provider, model, API key, pricing table, and
optional session budget -- inherits at *construction time* specifically,
which turns out to be the one subtlety worth internalizing before anything
else on this page: an already-built agent never notices a later
``AGENT_CONFIG.switch(...)``, only agents built after it do.

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

.. code-block:: pycon

   >>> from pycsamt.agents import (
   ...     AGENT_CONFIG,
   ...     AgentConfig,
   ...     BudgetExceededError,
   ...     configure_agents,
   ...     reset_agents,
   ... )

``AGENT_CONFIG`` is a process-local singleton, similar in spirit to the
plotting, section, style, and station-rendering API objects in pycsamt.  It is
not written to disk and is not shared between Python processes -- and it really
is the same object everywhere it is imported from, including the lower-level
module it is defined in:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG as via_agents
   >>> from pycsamt.api.agents import AGENT_CONFIG as via_api
   >>> via_agents is via_api
   True

Configure it at the beginning of a notebook, script, CLI command, or
application session.

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
     - Active LLM provider -- see the five supported families below.
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

pycsamt supports five LLM provider families.  If you select a provider but
do not pass a model, pycsamt uses the default listed here.

.. list-table::
   :header-rows: 1
   :widths: 20 28 52

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
     - ``GEMINI_API_KEY``, ``GOOGLE_API_KEY``,
       ``GOOGLE_GENERATIVEAI_API_KEY``, ``PYCSAMT_GEMINI_API_KEY``
   * - ``"deepseek"``
     - ``deepseek-chat``
     - ``DEEPSEEK_API_KEY``, ``PYCSAMT_DEEPSEEK_API_KEY``
   * - ``"minimax"``
     - ``MiniMax-M3``
     - ``MINIMAX_API_KEY``, ``MINIMAX_M3``, ``PYCSAMT_MINIMAX_API_KEY``

DeepSeek and MiniMax are OpenAI-compatible REST APIs reached through the same
``openai`` client with a provider-specific ``base_url``, so no separate
client package is needed for either.  For the other three, install the
matching client package only for the provider you intend to use:

.. code-block:: bash
   :linenos:

   pip install anthropic
   pip install openai
   pip install google-generativeai

The imports are lazy.  Importing :mod:`pycsamt.agents` does not require any
of these client libraries to be installed, and none is imported until an
agent actually calls ``query_llm()`` with a resolved key.


Recommended setup with environment variables
--------------------------------------------

For day-to-day use, keep secrets outside code and configure the active
provider from the environment.

.. code-block:: bash
   :linenos:

   export ANTHROPIC_API_KEY="sk-ant-demo1a2b"

Then activate the provider in Python:

.. code-block:: pycon

   >>> from pprint import pprint
   >>> from pycsamt.agents import AGENT_CONFIG

   >>> AGENT_CONFIG.switch("claude", model="claude-sonnet-4-6")

   >>> pprint(AGENT_CONFIG.info())
   {'api_key_masked': '…1a2b',
    'budget_usd': None,
    'custom_rate_models': {},
    'has_key': True,
    'key_source': 'env',
    'model': 'claude-sonnet-4-6',
    'provider': 'claude',
    'remaining_usd': None,
    'spent_usd': 0.0,
    'stored_providers': []}

``switch(...)`` is the right call when the key is already available through
the environment or was previously registered with ``set_key(...)``.  It selects
the provider and optional model without embedding the key in the script.  The
information dictionary above is safe to log or display: ``api_key_masked``
shows only the last four characters, and ``key_source: 'env'`` confirms the
key came from ``ANTHROPIC_API_KEY`` rather than a value stored in
``AGENT_CONFIG`` itself.


Explicit Python setup
---------------------

Use explicit setup when the key is supplied by a secret manager, CLI option,
or application settings layer.  ``runtime_secret`` below stands in for
whatever your settings layer actually returns:

.. code-block:: pycon

   >>> from pycsamt.agents import configure_agents, AGENT_CONFIG

   >>> runtime_secret = "sk-openai-demo9z9z"
   >>> configure_agents(
   ...     provider="openai",
   ...     api_key=runtime_secret,
   ...     model="gpt-4o-mini",
   ... )

   >>> AGENT_CONFIG.info()["key_source"]
   'explicit'
   >>> AGENT_CONFIG.info()["stored_providers"]
   ['openai']

Unlike the environment-variable path above, ``key_source`` now reads
``'explicit'`` -- the key was handed to pycsamt directly rather than
discovered in the process environment -- and ``stored_providers`` records
that a key is now held for ``"openai"`` specifically, which is what lets a
later :meth:`~pycsamt.api.agents.AgentConfig.switch` back to it skip
re-supplying the secret.  The same operation can be written with the
singleton directly instead of the ``configure_agents`` convenience wrapper:

.. code-block:: pycon

   >>> AGENT_CONFIG.configure(
   ...     provider="claude",
   ...     api_key=runtime_secret,
   ...     model="claude-sonnet-4-6",
   ... )

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
2. If the constructor does not receive an ``api_key`` and its ``llm_provider``
   is ``"claude"`` -- ``BaseAgent``'s default -- the agent inherits the active
   global provider when one has been configured.  This rule cannot tell an
   *explicit* ``llm_provider="claude"`` apart from simply not passing
   ``llm_provider`` at all, since both look identical by the time
   ``resolve()`` sees them; see the callout in
   :ref:`per-agent overrides <assistant-per-agent-overrides>` for the
   consequence.
3. The key for the effective provider is resolved as: stored key, then
   provider-specific environment variables, then ``None``.
4. The model comes from the constructor if supplied.  Otherwise, when the
   effective provider matches the active global provider, the global model is
   used.  If no model is available, pycsamt uses the provider default.
5. If no key is found, the agent still runs but skips LLM calls.  Returned
   ``AgentResult.llm_interpretation`` values will be ``None`` and LLM cost
   remains zero.

This behavior keeps simple workflows concise:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, DataQCAgent, PhaseAnalysisAgent

   >>> AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   >>> qc = DataQCAgent()              # inherits OpenAI + gpt-4o-mini
   >>> pt = PhaseAnalysisAgent()       # same provider and model
   >>> qc.llm_provider, qc.model
   ('openai', 'gpt-4o-mini')
   >>> pt.llm_provider, pt.model
   ('openai', 'gpt-4o-mini')


.. _assistant-per-agent-overrides:

Per-agent overrides
-------------------

Use a per-agent override when a specific task benefits from a different model
or provider.  For example, you might keep the global configuration on a fast
model, but ask the report agent to use a stronger model for final prose.  The
reliable way to do this is to pass an explicit ``api_key`` along with
``llm_provider`` -- rule 1 above then applies and the global configuration is
bypassed completely for that one agent:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, DataQCAgent, ReportAgent

   >>> AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   >>> qc_agent = DataQCAgent()

   >>> gemini_runtime_secret = "AIza-demo"
   >>> report_agent = ReportAgent(
   ...     llm_provider="gemini",
   ...     api_key=gemini_runtime_secret,
   ...     model="gemini-2.0-flash",
   ... )
   >>> report_agent.llm_provider, report_agent.model
   ('gemini', 'gemini-2.0-flash')

Omitting ``api_key`` and passing only ``llm_provider`` and ``model`` looks
like it should do the same thing, and it does -- for every provider *except*
``"claude"``, because ``"claude"`` is also ``BaseAgent``'s constructor
default.  Rule 2 above cannot distinguish "the caller explicitly asked for
Claude" from "the caller did not pass ``llm_provider`` at all", so an
explicit ``llm_provider="claude"`` with no key silently inherits whatever
provider is globally active instead:

.. code-block:: pycon

   >>> broken = ReportAgent(
   ...     llm_provider="claude",
   ...     model="claude-sonnet-4-6",
   ... )
   >>> broken.llm_provider, broken.model
   ('openai', 'claude-sonnet-4-6')

``broken`` now holds an OpenAI provider paired with a Claude model string --
not an agent waiting on a Claude key, as it would be reasonable to assume
from the constructor call alone.  It would run past construction without
error and only fail once ``query_llm()`` actually sent ``"claude-sonnet-4-6"``
to the OpenAI API.  Avoid it the same way the working example above does:
pass an explicit ``api_key`` whenever the target provider is ``"claude"`` and
might differ from whatever is globally active, or wrap the construction in
:meth:`AGENT_CONFIG.using <pycsamt.api.agents.AgentConfig.using>` (below)
instead of relying on constructor arguments alone.


Multi-provider sessions
-----------------------

You can register keys for several providers and switch between them during the
same Python session.

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, DataQCAgent

   >>> claude_key, openai_key, gemini_key = (
   ...     "sk-ant-demo1a2b", "sk-openai-demo9z9z", "AIza-demo",
   ... )
   >>> AGENT_CONFIG.set_key("claude", claude_key)
   >>> AGENT_CONFIG.set_key("openai", openai_key)
   >>> AGENT_CONFIG.set_key("gemini", gemini_key)

   >>> AGENT_CONFIG.switch("claude")
   >>> agent_a = DataQCAgent()         # constructed while Claude is active
   >>> agent_a.llm_provider, agent_a.model
   ('claude', 'claude-sonnet-4-6')

   >>> AGENT_CONFIG.switch("openai", model="gpt-4o-mini")
   >>> agent_b = DataQCAgent()         # constructed while OpenAI is active
   >>> agent_b.llm_provider, agent_b.model
   ('openai', 'gpt-4o-mini')

   >>> agent_a.llm_provider, agent_a.model    # unchanged by the later switch
   ('claude', 'claude-sonnet-4-6')

Existing agent instances keep the provider, key, and model they resolved at
construction time -- ``agent_a`` above never notices the later switch to
OpenAI, exactly the behavior flagged in the introduction.  Switch before
constructing the agents that should use the new provider, not after.


Temporary overrides
-------------------

Use ``AGENT_CONFIG.using(...)`` when a block of code should temporarily use a
different provider or model.  Unlike the per-agent ``llm_provider="claude"``
pitfall above, ``using(...)`` changes what "the active global provider" *is*
for the duration of the block, so the ordinary construction-time inheritance
rules apply correctly inside it -- no special-cased provider name to trip
over:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, InterpretationAgent

   >>> AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   >>> with AGENT_CONFIG.using(provider="claude", model="claude-sonnet-4-6"):
   ...     print("inside:", AGENT_CONFIG.provider, AGENT_CONFIG.model)
   ...     interpretation = InterpretationAgent()
   ...     print("agent: ", interpretation.llm_provider, interpretation.model)
   inside: claude claude-sonnet-4-6
   agent:  claude claude-sonnet-4-6

   >>> print("after: ", AGENT_CONFIG.provider, AGENT_CONFIG.model)
   after:  openai gpt-4o-mini

The OpenAI configuration is restored automatically on exit -- the context
manager also restores custom rates and budget state if an exception occurs
inside the block.  ``interpretation`` itself keeps the Claude provider it
resolved at construction time even after the block ends, following the same
construction-time-inheritance rule as everywhere else on this page; only
agents built *inside* the ``with`` block get the temporary override.


Running without an LLM
----------------------

LLM support is optional.  If no key is configured, agents skip calls to
external services and return deterministic pycsamt outputs.

.. code-block:: pycon

   >>> from pycsamt.agents import DataQCAgent, reset_agents

   >>> reset_agents()

   >>> result = DataQCAgent().execute({"path": "data/AMT/WILLY_DATA/L18PLT"})

   >>> assert result.cost_estimate_usd == 0.0
   >>> assert result.llm_interpretation is None
   >>> result.status
   'success'

Both assertions hold and the QC run still succeeds -- ``reset_agents()``
clears every stored key and the active provider, so this is the same
graceful degradation documented throughout :doc:`agent_catalogue` and
:doc:`foundation_agents`, exercised here explicitly rather than as a side
effect of simply never having configured a key.

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
session-local, and ``remaining_usd`` is a direct read of the two counters
underneath it:

.. math::

   \text{remaining}_{\text{USD}} = \max\bigl(0,\ \text{budget}_{\text{USD}} -
   \text{spent}_{\text{USD}}\bigr),

clamped at zero rather than going negative once the cap is reached.  The
check itself, ``spent_usd >= budget_usd``, runs inside
``BaseAgent.query_llm`` immediately before the network call, so a workflow
fails fast on the call that would exceed the cap rather than after paying
for it.  Simulating "several prior calls already reached the cap" with
``_add_spend`` -- the same accumulator ``query_llm`` updates after every
real response -- triggers the identical, real exception without needing a
network connection or a valid key to demonstrate it:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, BudgetExceededError

   >>> AGENT_CONFIG.switch("claude")
   >>> AGENT_CONFIG.set_budget(usd=2.00)
   >>> AGENT_CONFIG._add_spend(2.00)   # stand-in for prior LLM calls

   >>> try:
   ...     AGENT_CONFIG._check_budget()
   ... except BudgetExceededError as exc:
   ...     print(f"Budget ${exc.budget_usd:.2f} reached.")
   ...     print(f"Already spent ${exc.spent_usd:.4f}.")
   Budget $2.00 reached.
   Already spent $2.0000.

A real workflow never calls ``_check_budget()`` or ``_add_spend()``
directly -- ``query_llm()`` does both automatically; they are shown here
only to make the failure reproducible without live credentials.  Inspect the
remaining budget, and reset it, the same way at any point:

.. code-block:: pycon

   >>> AGENT_CONFIG.spent_usd
   2.0
   >>> AGENT_CONFIG.remaining_usd
   0.0

   >>> AGENT_CONFIG.reset_budget()      # zero the counter, keep the $2.00 cap
   >>> AGENT_CONFIG.spent_usd, AGENT_CONFIG.remaining_usd
   (0.0, 2.0)

   >>> AGENT_CONFIG.reset_budget(cap=True)  # also remove the cap entirely
   >>> print(AGENT_CONFIG.remaining_usd)
   None


Pricing and rate overrides
--------------------------

pycsamt estimates cost from provider rates expressed in USD per 1,000,000
tokens, using exactly the :math:`\text{cost}_{\text{USD}} = (n_{\text{in}}
r_{\text{in}} + n_{\text{out}} r_{\text{out}}) / 10^{6}` formula introduced in
:doc:`foundation_agents`.  Built-in rates are used when the model is known.
If a model is new, or if your organization wants to use a different
accounting rate, add an override:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG

   >>> AGENT_CONFIG.set_rate(
   ...     "openai",
   ...     "gpt-5-turbo",
   ...     input=5.00,
   ...     output=20.00,
   ... )

   >>> cost = AGENT_CONFIG.estimate_cost(
   ...     "openai",
   ...     "gpt-5-turbo",
   ...     input_tokens=12_000,
   ...     output_tokens=2_000,
   ... )

   >>> print(f"Estimated cost: ${cost:.4f}")
   Estimated cost: $0.1000

:math:`(12{,}000 \cdot 5.00 + 2{,}000 \cdot 20.00) / 10^{6} = 0.10` -- the
override took effect immediately, with no need to restart a session or
reconstruct already-built agents, because ``estimate_cost`` and
``query_llm`` both read the rate table live rather than caching it.  List
the effective rate table for a provider to see overrides merged with the
built-in models:

.. code-block:: pycon

   >>> from pprint import pprint

   >>> pprint(AGENT_CONFIG.list_rates("claude"))
   {'claude-3-5-sonnet-20241022': {'input': 3.0, 'output': 15.0},
    'claude-3-haiku-20240307': {'input': 0.25, 'output': 1.25},
    'claude-3-opus-20240229': {'input': 15.0, 'output': 75.0},
    'claude-haiku-4-5': {'input': 0.8, 'output': 4.0},
    'claude-haiku-4-5-20251001': {'input': 0.8, 'output': 4.0},
    'claude-opus-4-6': {'input': 15.0, 'output': 75.0},
    'claude-opus-4-7': {'input': 15.0, 'output': 75.0},
    'claude-opus-4-8': {'input': 15.0, 'output': 75.0},
    'claude-sonnet-4-6': {'input': 3.0, 'output': 15.0}}

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
agent workflows -- it is the same call, and the same dictionary shape, shown
running for real in the "Recommended setup" and "Explicit Python setup"
sections above.  Important fields are:

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

Reset everything, including stored keys, or keep stored keys for later
switching and clear only the active selection, custom rates, and budget:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, reset_agents

   >>> AGENT_CONFIG.set_key("claude", "sk-ant-demo")
   >>> AGENT_CONFIG.switch("claude")
   >>> AGENT_CONFIG.set_budget(usd=1.0)

   >>> reset_agents(keys=False)
   >>> AGENT_CONFIG.info()["stored_providers"]
   ['claude']
   >>> AGENT_CONFIG.info()["provider"], AGENT_CONFIG.info()["budget_usd"]
   (None, None)

   >>> reset_agents()
   >>> AGENT_CONFIG.info()["stored_providers"]
   []

``reset_agents(...)`` does not clear environment variables.  If an environment
variable is still present, a later ``AGENT_CONFIG.switch(provider)`` can resolve
that key again even after a full ``reset_agents()``.


Practical recipes
-----------------

Offline QC recipe
~~~~~~~~~~~~~~~~~

Use this pattern for tests, documentation, or secure offline processing.

.. code-block:: pycon

   >>> from pycsamt.agents import DataQCAgent, reset_agents

   >>> reset_agents()

   >>> qc = DataQCAgent()
   >>> result = qc.execute({"path": "data/AMT/WILLY_DATA/L18PLT"})

   >>> print(result.status)
   success
   >>> print(result.cost_estimate_usd)
   0.0


Local analyst recipe
~~~~~~~~~~~~~~~~~~~~

Use environment variables, choose a provider, set a small budget, then create
agents normally.

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, DataQCAgent, ReportAgent

   >>> AGENT_CONFIG.switch("claude", model="claude-sonnet-4-6")
   >>> AGENT_CONFIG.set_budget(usd=1.50)

   >>> qc = DataQCAgent()
   >>> report = ReportAgent()
   >>> qc.llm_provider, report.llm_provider
   ('claude', 'claude')
   >>> AGENT_CONFIG.remaining_usd
   1.5


Mixed-model recipe
~~~~~~~~~~~~~~~~~~

Use a fast global model for routine steps and a stronger model for the final
interpretation -- ``using(...)`` rather than a bare ``llm_provider="claude"``
constructor argument, per the :ref:`per-agent overrides
<assistant-per-agent-overrides>` caveat above:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG, DataQCAgent, InterpretationAgent

   >>> AGENT_CONFIG.switch("openai", model="gpt-4o-mini")

   >>> qc_agent = DataQCAgent()

   >>> with AGENT_CONFIG.using(provider="claude", model="claude-sonnet-4-6"):
   ...     interpretation_agent = InterpretationAgent()
   ...     print(interpretation_agent.llm_provider, interpretation_agent.model)
   claude claude-sonnet-4-6

   >>> qc_agent.llm_provider, AGENT_CONFIG.provider
   ('openai', 'openai')


Application session recipe
~~~~~~~~~~~~~~~~~~~~~~~~~~

Desktop or web applications should configure agents once from their settings
layer, then create agents from controllers or service classes:

.. code-block:: pycon

   >>> from pycsamt.agents import AGENT_CONFIG

   >>> def configure_agent_session(settings):
   ...     if settings.llm_enabled:
   ...         AGENT_CONFIG.configure(
   ...             provider=settings.provider,
   ...             api_key=settings.api_key,
   ...             model=settings.model,
   ...         )
   ...         if settings.budget_usd:
   ...             AGENT_CONFIG.set_budget(usd=settings.budget_usd)
   ...     else:
   ...         AGENT_CONFIG.reset()

Called once at startup with the application's own settings object --
here stood in for by a minimal stand-in with the same four attributes --
it leaves ``AGENT_CONFIG`` ready for every controller or service class
constructed afterward:

.. code-block:: pycon

   >>> from types import SimpleNamespace

   >>> settings = SimpleNamespace(
   ...     llm_enabled=True,
   ...     provider="openai",
   ...     api_key=runtime_secret,
   ...     model="gpt-4o-mini",
   ...     budget_usd=3.0,
   ... )
   >>> configure_agent_session(settings)
   >>> AGENT_CONFIG.provider, AGENT_CONFIG.model, AGENT_CONFIG.remaining_usd
   ('openai', 'gpt-4o-mini', 3.0)


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
   * - Agent built with ``llm_provider="claude"`` silently uses a different
       provider
     - ``"claude"`` is also ``BaseAgent``'s constructor default, so without an
       explicit ``api_key`` it is treated as "no override" and inherits
       whatever provider is globally active; see :ref:`per-agent overrides
       <assistant-per-agent-overrides>`.
     - Pass an explicit ``api_key`` alongside ``llm_provider="claude"``, or
       construct the agent inside ``AGENT_CONFIG.using(provider="claude")``.
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
* :doc:`foundation_agents` for the ``query_llm`` retry/backoff and cost
  formulas this page builds on
* :doc:`agent_catalogue` for the graceful-degradation behavior every recipe
  above relies on
