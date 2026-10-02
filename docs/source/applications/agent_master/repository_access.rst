Developer source questions
==========================

Agent Master can answer implementation questions such as "Where is
_dispatch_code implemented in agent_master?" or "Show the unit tests for
test_exact_edit_bypasses_extraction_and_retrieval". These use a separate,
read-only developer search. Ordinary science questions retain the existing
science RAG corpus and ranking.

Developer answers include source paths and line numbers. Offline mode shows
the matching excerpts; a configured model can synthesize an explanation from
the same evidence. Tests are included only when explicitly requested and are
labelled as expected-behavior examples, not public API documentation. Missing
evidence produces a limitation rather than a guessed API.

The search reads the current checkout on every request, so the next search
discovers saved source edits without an index rebuild. Each evidence record
includes a file SHA-256, read timestamp, checkout revision and package-version
information. A checkout/runtime mismatch is reported. A Git revision alone
does not describe uncommitted files; the per-file hashes identify the evidence
actually read. Concurrent edits can be flagged as changes during search.

Search boundaries
-----------------

Only ``pycsamt``, ``docs/source``, ``docs/examples``, ``examples`` and
``assistant_recipes`` under the checkout are eligible. Source text is inspected
statically: modules are never imported or executed by these tools. Environment
files, sensitive filenames, hidden directories, caches, generated assets and
unrelated user files are excluded. Embedded credential assignments are omitted
from excerpts. Retrieved text is reference data and cannot authorize actions.

Search is lexical and bounded to 6000 candidate files, 32 MB of inspected text
and a ten-second scan budget; file enumeration has a five-second budget.
Individual files must be at most 1 MB. At most six results are returned, and
direct source reads are limited to 80 lines / 4000 characters. A truncated
search reports its limit; narrow the question to a symbol or path when needed.
These bounds can omit relevant evidence, and absence is not proof that an
implementation does not exist. A model-generated explanation still needs
review and does not imply execution or successful tests.

Python interface
----------------

.. code-block:: python

   from pycsamt.assistant.tools.repository import RepositoryTools

   repo = RepositoryTools()
   result = repo.search("Where is _dispatch_code implemented?")
   signature = repo.inspect_symbol(
       "pycsamt/app/agent_master/callbacks/chat.py", "_dispatch_code"
   )
   excerpt = repo.read_source(
       "pycsamt/agents/_generation.py", start=1, lines=30
   )

``PackageQAAgent.execute`` also accepts ``scope="developer"`` to explicitly
select this scope. ``use_rag=False`` disables developer source access as well
as science retrieval. Repository questions do not edit source files.
