Conversation context and stopping work
======================================

Agent Master keeps active data, survey line, selected parameters, prior scripts
and recorded workflow results with the current conversation. Each background
job receives its own snapshot. Other conversations do not supply fallback data
or appear in the current conversation's **Recent runs** list.

Follow-ups such as ``now run phase tensor analysis`` can reuse the active
dataset. ``Show my previous result`` recalls the recorded workflow summary
without running it again. An explicit code edit can reuse an earlier script
even after several intervening questions. The original script remains intact
when a new version is saved.

The configured project root and line registry identify the project. Inline
registry YAML is identified by a content fingerprint; a registry file uses its path.
Changing that boundary drops earlier project context from subsequent requests.
A previously loaded dataset is not silently carried into the new project;
load the intended dataset or name a line in the selected registry. Changing
the selected line within a project retains the conversation.

**New Chat** cancels the current request and clears conversation memory,
active data, figures, pending actions and workflow parameters. Restoring a
saved conversation also cancels the current request and clears those transient
stores before using the restored conversation's own memory. Code blocks are
restored; figure images are not. Provider settings remain configured.

Long conversations
------------------

Question answering receives at most six recent whole turns within a
6000-character budget, a factual memory summary of at most 2000 characters,
and up to 1000 characters of explicitly abbreviated older user requests.
Code generation retains its existing 8000-character recent-turn budget and
adds a bounded history digest. The current request and a script being edited
are preserved in full; model context limits can still require clarification
or a shorter request. Summaries are deterministic excerpts and recorded facts,
not a claim that every older requirement is still present in model context.

Answers and artifacts
---------------------

Answers show source references when available and distinguish explanations,
generated scripts and recorded workflow runs. A saved-script notice appears
only when the returned script path exists. Generated code includes its
validation report; producing code does not mean that its computation ran.

A narrow guard withholds model answers that explicitly claim an unperformed
execution or saved output. This complements the answer prompt and authoritative
execution labels; it is not a general detector of all model inaccuracies.
Scientific correctness and broader model quality still require evaluation.

Stopping a request
------------------

**Stop** cancels further stages and suppresses late results and memory updates.
The same request signal is checked at retrieval, agent execution, model-call,
retry, repair and script-publication boundaries for all providers. Registered
workflow chains check it before starting another step or saving a checkpoint.
Optional fixture execution also receives the signal.

Cancellation is cooperative: a synchronous provider call or scientific kernel
already running may finish. The UI explains that limitation instead of claiming
that the operating-system work was killed. Files written before cancellation
are retained. New Chat and history restore use the same cancellation behavior.
