.. _api-agent-config:

Agent Configuration
===================

:mod:`pycsamt.api.agents` holds :data:`~pycsamt.api.agents.AGENT_CONFIG`,
the process-local singleton that every :term:`agent` in :mod:`pycsamt.agents`
reads its LLM provider, model, API key, pricing table, and spending cap
from -- configure it once, and every agent constructed afterwards inherits
the setting automatically.

The full guide -- provider setup, credential resolution order, the
``.env.local`` fallback, temporary overrides, budget caps, custom pricing,
and a troubleshooting table -- already lives at
:doc:`../user_guide/agents/llm_configuration`. This page exists only to
place ``AGENT_CONFIG`` in the same family lineup as every other
:mod:`pycsamt.api` singleton (see :doc:`overview`) and give a fast
quick-reference; read the linked page before configuring a real provider.

One thing to notice up front: unlike :data:`~pycsamt.api.style.PYCSAMT_STYLE`
or :data:`~pycsamt.api.interp.PYCSAMT_INTERP`, ``AGENT_CONFIG`` has no
dotted-path ``configure(section__attr=...)`` method. Its settings --
provider, key, model, rates, budget -- are not a nested style tree, so
:meth:`~pycsamt.api.agents.AgentConfig.configure` takes plain keyword
arguments instead:

.. code-block:: pycon

   >>> from pycsamt.api.agents import AGENT_CONFIG, configure_agents, reset_agents

   >>> _ = reset_agents()
   >>> _ = configure_agents(
   ...     provider="claude",
   ...     api_key="sk-ant-demo-0000",
   ...     model="claude-sonnet-4-6",
   ... )
   >>> print(AGENT_CONFIG)
   AgentConfig(provider='claude', model='claude-sonnet-4-6', key='…0000')

``AGENT_CONFIG.info()`` gives the same state as a plain dict, with the key
masked to its last four characters -- safe to log or print without leaking
a credential:

.. code-block:: pycon

   >>> info = AGENT_CONFIG.info()
   >>> info["provider"], info["model"], info["api_key_masked"], info["key_source"]
   ('claude', 'claude-sonnet-4-6', '…0000', 'explicit')

A session budget cap is one call, checked before every LLM request an
agent makes -- see :doc:`../user_guide/agents/llm_configuration` for the
:class:`~pycsamt.api.agents.BudgetExceededError` it raises once the cap is
reached:

.. code-block:: pycon

   >>> AGENT_CONFIG.set_budget(usd=1.5)
   AgentConfig(provider='claude', model='claude-sonnet-4-6', key='…0000', budget=$1.50 (spent=$0.0000))
   >>> AGENT_CONFIG.spent_usd, AGENT_CONFIG.remaining_usd
   (0.0, 1.5)

   >>> _ = reset_agents()

:meth:`~pycsamt.api.agents.AgentConfig.offline` is narrower than it first
sounds: it only suppresses the *environment-variable* fallback for the
current thread, so a stray ``ANTHROPIC_API_KEY`` in the shell can't cause
an unintended LLM call from a script or test run. An explicitly stored key
(from :meth:`~pycsamt.api.agents.AgentConfig.configure` or
:meth:`~pycsamt.api.agents.AgentConfig.set_key`) is unaffected and still
resolves inside the block:

.. code-block:: pycon

   >>> import os
   >>> os.environ["ANTHROPIC_API_KEY"] = "sk-ant-env-fallback-demo"
   >>> AGENT_CONFIG.switch("claude")
   AgentConfig(provider='claude', model='claude-sonnet-4-6', key='…demo')
   >>> AGENT_CONFIG.api_key is not None
   True

   >>> with AGENT_CONFIG.offline():
   ...     print(AGENT_CONFIG.api_key)
   None
   >>> AGENT_CONFIG.api_key is not None
   True

   >>> del os.environ["ANTHROPIC_API_KEY"]
   >>> _ = reset_agents()

Next Steps
----------

* :doc:`../user_guide/agents/llm_configuration` for the full guide: the
  quick decision table, credential resolution order, ``.env.local``,
  :meth:`~pycsamt.api.agents.AgentConfig.using`, custom pricing with
  :meth:`~pycsamt.api.agents.AgentConfig.set_rate`, and troubleshooting.
* :doc:`../user_guide/agents/agent_catalogue` for what each concrete
  agent class does once it has a provider to call.
* :doc:`overview` for how the agents family fits alongside every other
  :mod:`pycsamt.api` configuration family.
