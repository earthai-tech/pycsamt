Generated-code validation and repair
====================================

Agent Master validates scripts before saving them. A script with a detected
failure is returned as an unsaved draft. Earlier script files are preserved;
successful publication uses a new filename and a matching
``.validation.json`` report containing the code hash and repair history.
The chat history retains the report for follow-up work.

Each check has a state: ``passed``, ``failed``, ``not_checked`` or
``unverifiable``. Reports distinguish:

* Syntax: parse and compile, without execution.
* Imports: static existence of supported pyCSAMT modules and symbols in source.
* Arguments: supported function/constructor/method signatures, positional and
  keyword arguments. Dynamic signatures, rebinding, inheritance and expanded
  ``*args`` / ``**kwargs`` can remain unverifiable.
* Scientific literals: finite positive layer resistivities/thicknesses with
  matching dimensions, and literal correction-factor mappings.
* Runtime imports: not checked; source inspection does not prove optional
  dependencies are installed or that module imports will succeed.
* Fixture execution and artifacts: optional Docker checks with an explicitly
  supplied local image ID and output contracts. Unavailable isolation remains
  unverifiable; there is no host execution fallback.

External-library imports are explicitly not checked. Runtime array values,
dataframe contents, types and scientific results are not inferred from syntax.
The compatibility ``ok`` field means no detected errors in performed checks; it is not a
certificate that the script runs or produces valid science.

Bounded repair
--------------

When local generation fails a supported static or request-constraint check,
Agent Master can send the errors, source signatures, request contract and
current script back to the model. It performs at most two repair calls, each
limited to 1024 output tokens, within the original local request deadline and
a maximum 60-second generation/repair window. Model call limits and cancellation
still apply. Empty, unchanged or late responses stop repair; late responses are
discarded. Every changed candidate is validated again before publication.

Automatic repair is currently local-only. Cloud transports lack the shared
deadline required by this path, so cloud drafts still receive validation but
their reports explicitly state that automatic repair was skipped. Offline
templates and exact directory/DPI edits are validated without model repair.
Context overflow or exhausted time leaves an unsaved draft with its report.

Execution boundary
------------------

Generated scripts are never executed on the host by the validator. Passing
``execute_fixture=True`` requests execution checking but returns
``unverifiable`` without an explicit ``fixture_image`` (a trusted, locally
installed Linux image's immutable ``sha256:...`` ID) and working Docker.
The image must contain Python and the script's dependencies. Images are not
downloaded or built automatically. Existing registered workflows remain a
separate execution path.

The disposable container has no network, a read-only root, no capabilities,
an unprivileged user, and CPU, memory, PID, file-size and wall-time limits.
The generated code, the trusted checker and the output contracts travel over
the Docker client's standard input; nothing is mounted from the host, so no
project directories, cloud keys or Docker socket are visible inside the
container. Because no host path is shared, the Docker daemon can run locally,
inside WSL2, or on another host. Generated outputs stay
in a bounded temporary filesystem and are discarded after checking. Cancellation
and timeouts trigger container removal. A failed Docker connection can prevent
cleanup confirmation; this is reported explicitly.

``fixture_outputs`` maps relative output paths under ``/work`` to contracts:
``{"kind": "file"}`` requires a nonempty regular file of at most 1 MiB;
``{"kind": "positive_factors"}`` checks a JSON mapping of finite positive
numbers; ``{"kind": "array", "shape": [2, 2]}`` checks a finite JSON array's
dimensions. Symlink outputs are rejected. Fixtures must use synthetic inputs
and write relative paths; arbitrary survey data is not copied into the container.
The fixture image and contracts are supplied by the caller, not inferred by
the language model. Container reports remain generated-program evidence, not
proof against intentional result fabrication or full scientific correctness.

For example, with a trusted Python image already installed locally::

   report = validate_generated_code(
       "import json\nopen('factors.json', 'w').write(json.dumps({'S1': 1.2}))",
       execute_fixture=True,
       fixture_image="sha256:" + "<64-character-local-image-id>",
       fixture_outputs={"factors.json": {"kind": "positive_factors"}},
   )

Setting up Docker and the fixture image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The validator uses ``docker`` on ``PATH``. Set ``PYCSAMT_DOCKER`` to use a
different client command. On Windows without Docker Desktop, Docker Engine
can run inside a WSL2 distribution, for example Ubuntu::

   wsl -u root -e sh -c "apt-get update && apt-get install -y docker.io && systemctl enable --now docker"

Then point the validator at the WSL client:

.. code-block:: powershell

   $env:PYCSAMT_DOCKER = "wsl.exe -e docker"

``packaging/agent-fixture/build.sh`` builds an image with Python and pyCSAMT's
core dependencies and prints its immutable ID. Only git-tracked package files
at ``HEAD`` enter the build context, so ignored files such as ``.env.local``
cannot reach the image. Run it wherever the Docker client runs::

   wsl -e bash packaging/agent-fixture/build.sh

Live acceptance tests require ``PYCSAMT_FIXTURE_IMAGE`` to contain the printed
``sha256:...`` ID:

.. code-block:: powershell

   $env:PYCSAMT_FIXTURE_IMAGE = "sha256:<64-character-image-id>"
   python -m pytest pycsamt/assistant/tools/tests/test_fixture_execution.py

Without it, live tests are explicitly skipped. Rebuild the image after
changing pyCSAMT so that fixture scripts see the current package.

.. code-block:: python

   from pycsamt.assistant.tools.validation_tools import validate_generated_code

   report = validate_generated_code(
       "from pycsamt.emtools._core import ensure_sites\n"
       "sites = ensure_sites('data/3edis')\n"
   )
   print(report["checks"])

``CodeGenerationAgent.execute`` accepts ``max_repairs`` (clamped to 0–2) and
``execute_fixture``, ``fixture_image`` and ``fixture_outputs``. A saved artifact's full report is returned as
``validation`` and its sidecar path as ``validation_path``. General correctness
and constraint fulfillment still require review.
