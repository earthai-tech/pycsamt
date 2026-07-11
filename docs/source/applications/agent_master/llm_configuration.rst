.. _applications-agent-master-llm:

LLM Configuration
=================

Agent Master uses a large language model to understand your requests and plan
workflows.  You choose the provider and supply an API key; the key stays on
your machine.  Open the configuration from the **Settings** gear on the command
bar.

The Settings Drawer
-------------------

.. figure:: ../../_static/applications/agent_master/setting-llm-configuration.png
   :alt: The Agent Master Settings drawer with LLM providers, API key, and model
   :class: pycsamt-screenshot

   **Settings**: choose a provider, paste an API key, pick a model, and set the
   active provider and the default figure export.

* **LLM Providers** — **Claude**, **OpenAI**, **Gemini**, or **DeepSeek**.
* **API Key** — paste the key for the selected provider; the eye/key toggle
  reveals it.  **Save API keys** stores it.
* **Model** — the model id to use (for example ``claude-sonnet-4-6``).
* **Active LLM provider** — which configured provider requests actually use.
* **Default figure export** — the format and resolution for figures the agents
  save (for example ``PNG (300 dpi)``).

Choosing A Provider And Model
-----------------------------

1. Pick a **provider** tab.
2. Paste its **API key**.
3. Set the **model** id (a sensible default is pre-filled).
4. Select it under **Active LLM provider**.
5. **Save API keys**.

You can configure more than one provider and switch the active one at any time.
When building AI applications on top of pyCSAMT, prefer the latest and most
capable models each provider offers.

Privacy And Keys
----------------

* API keys are kept **on your machine** and used only to call the provider you
  selected.
* Requests you send are processed by the chosen provider under that provider's
  terms — treat prompts and any pasted data accordingly.
* No key is required to *browse* the interface; a key is required to run the
  LLM-backed planning and chat.

When No Model Is Configured
---------------------------

Without a configured provider and key, the language-driven parts of Agent
Master cannot run — requests that need planning will prompt you to add a key.
Configure a provider in **Settings** first, then send your request again.

.. seealso::

   :doc:`/user_guide/agents/llm_configuration`
       Provider and model configuration for the agents at the library level,
       including environment variables and programmatic setup.

Next Steps
----------

* :doc:`workflows_and_agents` -- what the configured model then lets you run.
* :doc:`troubleshooting` -- provider connection and key problems.
