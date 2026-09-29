.. _agent-master-local-models:

Local language models
=====================

Agent Master supports **Local LLM (Ollama)** alongside its cloud providers.
This runs a downloaded language model on your computer. It is different from
**Offline**, which uses deterministic rules and templates without a model
server. Local inference has no API charge, but uses your CPU/GPU, memory,
storage, and electricity.

Setup
-----

Install Ollama from its official distribution, start its local service, and
download a model that fits your hardware. For example::

    ollama pull qwen2.5-coder:1.5b
    ollama serve

Run ``serve`` only if the Ollama service is not already running. To explicitly
disable Ollama cloud features, set ``OLLAMA_NO_CLOUD=1`` in the server's
environment before starting it. Models must be downloaded before offline use.

In Agent Master settings:

1. Select **Local LLM (Ollama)**.
2. Set the endpoint, normally ``http://127.0.0.1:11434``.
3. Enter the exact installed model name. Custom local model names are accepted.
4. Select **Check connection and model**, then **Save settings**.

No API key is required. Connection checking verifies that the model is
installed; it does not certify answer quality or generation speed. Agent
Master does not download models automatically.

Limits and behavior
-------------------

The default request time limit is 60 seconds, the context window is 8192
tokens, the output limit is 1024 tokens per call, and the request allows up to
four model calls. The total time limit includes retrieval and model loading
in the chat request. Adjust these settings for the hardware and workload.

Local mode uses deterministic routing for clear requests and deterministic
workflow configuration extraction, reserving model generation for answers,
scripts, and scientific interpretations. Ambiguous routing may use the model.
RAG evidence is retained, with smaller context selections for local answers
and code. Arbitrary request fidelity and stronger code verification remain
separate capabilities from the model-provider integration.

The client checks a conservative UTF-8 byte bound before generation; this is
not an exact tokenizer count. If the prompt exceeds that bound, the request
fails explicitly instead of silently dropping evidence. Truncated model
output is also reported as an error. A larger configured context consumes
more memory and can increase latency.

The Stop button interrupts the local HTTP request, including model loading.
Late results cannot overwrite a cancelled chat job. Local token counts and
timings are recorded when Ollama provides them. Dollar-based cloud budgets
do not represent local compute; local requests use time and call limits.

Privacy and failure behavior
----------------------------

This provider accepts loopback HTTP endpoints only, bypasses system proxies,
does not follow redirects, and rejects cloud models and model metadata that
indicates remote hosting. It never falls back to a cloud provider. Its RAG
cache is separated from cloud-capable retrieval, uses lexical retrieval, and
disables cloud embeddings and optional external reranking in local sessions.
Query rewriting remains deterministic.

An unavailable server, missing model, timeout, context overflow, or generation
failure produces an error explaining the next step. Small local models can
still invent APIs or misuse valid functions. The existing syntax/import
validator is not proof that a generated script executes correctly or produces
scientifically correct results. Generating a script does not execute it.

Python use
----------

Use a request-scoped configuration when controlling local settings from Python::

    from pycsamt.agents import LocalSettings, PackageQAAgent, local_session

    settings = LocalSettings(
        model="qwen2.5-coder:1.5b",
        endpoint="http://127.0.0.1:11434",
        timeout=60,
    )
    with local_session(settings):
        result = PackageQAAgent(llm_provider="ollama").execute({
            "question": "What does pycsamt.emtools._core.ensure_sites return?"
        })

For default local settings, ``configure_agents(provider="ollama")`` requires
no key. Use a separate ``local_session`` inside each request thread; context
does not propagate automatically to newly created threads. Deterministic
``AGENT_CONFIG.offline()`` disables local inference as well as cloud inference.

Runtime protocol references: `Ollama chat API <https://docs.ollama.com/api/chat>`_,
`model listing <https://docs.ollama.com/api/tags>`_, and
`error responses <https://docs.ollama.com/api/errors>`_.
