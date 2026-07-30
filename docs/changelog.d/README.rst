Changelog fragments
====================

This directory holds one small file per user-facing pull request, instead
of everyone editing a single shared ``changelog.rst``. At release time,
``docs/scripts/changelog_release.py`` assembles every fragment here into
the current series page under ``docs/source/changelog/`` and deletes the
fragments it consumed (git history keeps them).

Add a fragment in the same PR as the change it documents.

Filename
--------

``<ref>.<type>.rst``

* ``<ref>`` — the GitHub PR or issue number (``123``). If there is no
  issue/PR number, use ``+<short-slug>`` instead (e.g. ``+ordering-docs``).
* ``<type>`` — one of the types below. An unrecognised type is a hard
  error when the release script runs, so a typo gets caught immediately
  rather than silently dropped.

.. list-table::
   :header-rows: 1
   :widths: 15 20 15

   * - type
     - changelog section
     - badge
   * - ``feature``
     - Added
     - |Feature|
   * - ``new``
     - Added
     - |New|
   * - ``enhancement``
     - Changed
     - |Enhancement|
   * - ``perf``
     - Changed
     - |Perf|
   * - ``fix``
     - Fixed
     - |Fix|
   * - ``api``
     - Changed
     - |API Change|
   * - ``deprecated``
     - Changed
     - |Deprecated|
   * - ``breaking``
     - Changed
     - |Breaking|
   * - ``removed``
     - Removed
     - |Breaking|
   * - ``security``
     - Security
     - |Security|
   * - ``docs``
     - Docs & tooling
     - |Docs|
   * - ``build``
     - Docs & tooling
     - |Build|
   * - ``tests``
     - Docs & tooling
     - |Tests|

Content
-------

Write exactly the bullet body — no leading ``* |Badge|``, the release
script adds that from the filename's ``<type>``. Bold the short name of the
change first, matching the existing changelog style.

``docs/changelog.d/123.fix.rst``::

   **Oblique survey lines** -- order is derived from both geographic
   coordinates and projected chainage instead of longitude, latitude, or
   lexical station names alone.

Cutting a release
------------------

.. code-block:: console

   python docs/scripts/changelog_release.py --version 2.2.0 --date 2026-08-15

This groups every fragment by section, inserts the assembled block at the
top of the matching ``docs/source/changelog/vN.rst`` series file, and
deletes the consumed fragments. It leaves a ``TODO`` placeholder for the
one-paragraph release summary (or pass one with ``--summary``) — write
that by hand, then write the narrative
``docs/source/release_notes/vX.Y.Z.rst`` page as usual.
