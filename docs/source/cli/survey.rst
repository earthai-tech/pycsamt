Survey Commands
===============

``pycsamt survey`` manages the active survey context used by other CLI
commands. It lets you register an EDI file or directory once, cache the loaded
``Sites`` collection, inspect the current context, rebuild stale caches, and
clear or purge cached survey state.

Use this command group when you want commands such as ``pycsamt site``,
``pycsamt pipe``, ``pycsamt info``, or ``pycsamt invert build`` to operate on
the same survey without repeating the input path every time.

What The Survey Context Stores
------------------------------

The active survey context is lightweight metadata stored under the user's
home directory:

.. code-block:: text

   ~/.pycsamt/
   |-- context.json
   `-- cache/
       `-- <12-char-cache-key>/
           |-- sites.pkl
           `-- meta.json

``context.json``
    Stores the active survey path, cache key, timestamp, station count, and
    station names.

``sites.pkl``
    Stores the cached ``Sites`` object for fast reuse by CLI commands.

``meta.json``
    Stores the source path, cache timestamp, source mtime at cache time, and
    station count.

The cache key is a 12-character hash derived from the resolved absolute survey
path. Moving the same EDI folder to a different location creates a different
cache entry.

Command Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Command
     - Purpose
   * - ``pycsamt survey set EDI_DIR``
     - Register an EDI source as the active survey and build or reuse its
       ``Sites`` cache.
   * - ``pycsamt survey show``
     - Print the current active survey summary.
   * - ``pycsamt survey clear``
     - Unset the active survey context while keeping cache files.
   * - ``pycsamt survey rebuild``
     - Force-rebuild the cache for the active survey.
   * - ``pycsamt survey cache list``
     - List cached survey entries.
   * - ``pycsamt survey cache purge``
     - Delete the active survey cache, or all survey caches with ``--all``.

Resolution Rules
----------------

Commands that need EDI data call ``resolve_survey``. The practical resolution
order is:

1. explicit positional path or command-specific ``--survey`` path;
2. active survey context from ``~/.pycsamt/context.json``;
3. usage error with a hint to run ``pycsamt survey set`` or pass ``--survey``.

For an explicit path, pyCSAMT checks whether a valid cache already exists for
that path. If the cache is valid and ``--fresh`` was not requested, the cached
``Sites`` object is reused. If the cache is stale, missing, corrupt, or
``--fresh`` is used, the path is parsed again and the cache is updated.

For the active context, the same cache validity rules apply. If the active
survey cache is stale, pyCSAMT silently rebuilds it before returning the
``Sites`` object.

Cache Validity
--------------

Cache validity is based on source modification time:

* for a single file, the file mtime is used;
* for a directory containing ``*.edi`` files, the newest EDI file mtime is
  used;
* for a directory without EDI files, the directory mtime is used.

If the current source mtime is newer than the cached mtime, the cache is
considered stale. A corrupt or unreadable ``sites.pkl`` is treated as a cache
miss and rebuilt when the survey is resolved.

Use ``--fresh`` on commands that expose it when you want to bypass a valid
cache for one command. Use ``survey rebuild`` when the active survey cache
should be refreshed immediately.

Set A Survey
------------

``survey set`` registers an EDI source as the active survey and pre-builds the
``Sites`` cache.

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt survey set data/AMT/WILLY_DATA/L18PLT --force
   pycsamt survey set data/AMT/WILLY_DATA/L18PLT -v

Options:

``--force``
    Rebuild the cache even when the existing cache is still valid.

``-v, --verbose``
    Print more build progress. Use ``-vv`` for deeper loader verbosity where
    supported.

``--no-color``
    Disable colored terminal output.

After a successful ``set``, commands that accept survey input can run without
an explicit path:

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt site info
   pycsamt pipe run --preset basic_qc --dry-run
   pycsamt invert build --solver occam2d --workdir runs/l18_occam

You can still override the active survey for one command by passing that
command's positional source or ``--survey`` option.

Show The Active Survey
----------------------

``survey show`` prints the active context summary.

.. code-block:: console

   pycsamt survey show
   pycsamt survey show --format json

Text output reports:

* survey path;
* time the context was set;
* station count;
* station names;
* cache status, ``valid`` or ``stale``;
* cache directory path;
* cache build timestamp when available.

JSON output includes keys such as ``survey_path``, ``set_at``,
``n_stations``, ``station_names``, ``cache_key``, ``cache_valid``,
``cache_path``, and, when available, ``cached_at``.

When no active survey is set, the command exits successfully and prints a hint
to run ``pycsamt survey set <edi_dir>``.

Clear The Active Survey
-----------------------

``survey clear`` removes only the active context file. It does not delete
cached ``Sites`` objects.

.. code-block:: console

   pycsamt survey clear
   pycsamt survey clear --yes

Options:

``--yes``
    Skip the confirmation prompt.

``--no-color``
    Disable colored terminal output.

Use ``clear`` when you do not want commands to fall back to the previous
survey. Use ``cache purge`` when you also want to remove cached files.

Rebuild The Active Cache
------------------------

``survey rebuild`` rebuilds the cache for the active survey.

.. code-block:: console

   pycsamt survey rebuild
   pycsamt survey rebuild --force
   pycsamt survey rebuild -v

The current implementation always calls ``set_survey`` with force enabled
inside ``rebuild``. The ``--force`` flag is accepted for clarity and future
compatibility, but the rebuild command already performs an unconditional cache
refresh for the active survey.

Use ``rebuild`` after:

* adding or removing EDI files;
* editing EDI contents outside pyCSAMT;
* seeing stale station names or station counts;
* recovering from a suspected cache issue.

If no active survey exists, ``rebuild`` exits with a usage error and suggests
``pycsamt survey set <edi_dir>``.

List Cached Surveys
-------------------

``survey cache list`` shows the survey cache entries under
``~/.pycsamt/cache``.

.. code-block:: console

   pycsamt survey cache list
   pycsamt survey cache list --format json

Text output marks the active cache entry with ``[active]``. Each row shows the
cache key, source path, station count, pickle size, and cached timestamp when
metadata is available.

JSON output emits one object per cache entry. Objects may include:

``key``
    Cache directory name.

``survey_path``
    Resolved source path stored in cache metadata.

``cached_at``
    Cache creation timestamp.

``cached_at_mtime``
    Source mtime recorded at cache time.

``n_stations``
    Station count stored in cache metadata.

``pkl_size_kb``
    Size of ``sites.pkl`` in KiB.

If no cache directory exists, the command prints ``No cache directory found``.
If the directory exists but has no entries, it prints ``Cache is empty``.

Purge Cache Entries
-------------------

``survey cache purge`` deletes cached ``Sites`` objects.

.. code-block:: console

   pycsamt survey cache purge
   pycsamt survey cache purge --yes
   pycsamt survey cache purge --all
   pycsamt survey cache purge --all --yes

Modes:

``pycsamt survey cache purge``
    Delete the cache for the active survey only.

``pycsamt survey cache purge --all``
    Delete every cached survey under ``~/.pycsamt/cache``.

Options:

``--yes``
    Skip the confirmation prompt.

Important behavior:

* Purging the active cache does not clear ``context.json``.
* After purging the active cache, the next command that resolves the active
  survey rebuilds the cache.
* ``--all`` removes the entire cache root, but does not remove the active
  context file.
* Without an active survey, active-only purge fails and suggests ``--all``.

How Other Commands Use Survey Context
-------------------------------------

``pycsamt site`` commands call the shared resolver through their ``_get_sites``
helper. Their resolution priority is positional ``EDI_DIR``, then
``--survey``, then active context.

``pycsamt pipe run`` also resolves EDI data through survey context. A preset
run can therefore be as short as:

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT
   pycsamt pipe run --preset basic_qc --out results/l18_basic_qc

``pycsamt invert build`` accepts an explicit EDI path or uses active context
when no source is passed. This makes the survey command useful before longer
inversion setup sessions.

Typical Workflow
----------------

For a normal CLI session:

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L18PLT -v
   pycsamt survey show
   pycsamt site info
   pycsamt site select --keep-finite --drop-empty --dry-run
   pycsamt pipe run --preset basic_qc --dry-run
   pycsamt pipe run --preset basic_qc --out results/l18_basic_qc

When the source data changes:

.. code-block:: console

   pycsamt survey show
   pycsamt survey rebuild -v
   pycsamt survey show --format json

When switching surveys:

.. code-block:: console

   pycsamt survey set data/AMT/WILLY_DATA/L22PLT
   pycsamt site info

When cleaning up:

.. code-block:: console

   pycsamt survey clear --yes
   pycsamt survey cache purge --all --yes

Troubleshooting
---------------

No active survey
    Run ``pycsamt survey set EDI_DIR`` or pass an explicit source path to the
    command that needs data.

Command seems to use old station data
    Re-run the command with ``--fresh`` when available, or run
    ``pycsamt survey rebuild`` for the active survey.

``survey show`` says cache is stale
    Run ``pycsamt survey rebuild``. Many data commands also rebuild stale cache
    automatically when they resolve the survey.

Cache list shows old survey paths
    This is normal. Cache entries are path-based and remain until purged. Use
    ``pycsamt survey cache purge --all`` to remove every cached survey.

Active context remains after purge
    ``cache purge`` removes cache files, not ``context.json``. Run
    ``pycsamt survey clear`` to unset the active context.

Rebuild fails after moving data
    The active context stores the old absolute path. Run ``pycsamt survey set``
    with the new path.

Safety Notes
------------

``survey set``
    Writes context and cache files under ``~/.pycsamt``. It does not modify
    input EDI files.

``survey clear``
    Removes only active context metadata. It preserves cache entries.

``survey cache purge``
    Deletes cached ``Sites`` objects. It does not modify the original EDI
    source.

``survey rebuild``
    Re-parses the active survey and replaces its cache entry. It does not
    modify source EDI files.

Python Equivalent
-----------------

The CLI context helpers live in ``pycsamt.cli.survey``:

.. code-block:: python

   from pathlib import Path

   from pycsamt.cli.survey import (
       resolve_survey,
       set_survey,
       survey_summary,
   )

   sites = set_survey(Path("data/AMT/WILLY_DATA/L18PLT"), force=True)
   summary = survey_summary()

   # Later, resolve the active survey.
   active_sites = resolve_survey(None)

For science workflows that do not need persistent CLI context, load sites
directly:

.. code-block:: python

   from pycsamt.emtools._core import ensure_sites

   sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True, verbose=0)

