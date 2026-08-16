.. _pipeline-observability:

Live Observability
==================

Until this page's features existed, ``Pipeline.run`` was a fully synchronous
black box: progress was ``print``/``tqdm`` text to stdout only, and the only
way to see anything about a run was to wait for the returned
``PipelineResult``. Four additions change that, all opt-in and all able to
be combined freely:

* an :ref:`on_step callback <pipeline-observability-on-step>`, fired once
  per step with its real :class:`~pycsamt.pipeline.StepResult`;
* a :ref:`live-updating terminal table <pipeline-observability-live>`
  (``progress_style="rich"`` / CLI ``--live``);
* :ref:`inline HTML rendering <pipeline-observability-notebook>` of a
  finished ``PipelineResult`` in Jupyter;
* an :ref:`opt-in run history log <pipeline-observability-history>`.

None of these touch how a pipeline actually processes data -- they only
observe it. ``on_step=None`` and ``history=False`` (the defaults) mean zero
behaviour change from every earlier release.

.. _pipeline-observability-on-step:

The on_step Hook
----------------

``Pipeline.run(..., on_step=callback)`` calls *callback* once per step,
immediately after that step's :class:`~pycsamt.pipeline.StepResult` is
built -- for a fresh computation, a :doc:`cache <caching>` hit, or a
warn/skip error alike:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step, configure_pipe
   >>> from pycsamt.emtools import ensure_sites
   >>> configure_pipe(show_progress=False)
   >>> import glob
   >>> paths = sorted(glob.glob("data/AMT/WILLY_DATA/L22PLT/*.edi"))[:2]
   >>> sites = ensure_sites(paths)
   >>> pipe = Pipeline([("notch", Step("NR001")), ("band", Step("FREQ001"))])
   >>> seen = []
   >>> result = pipe.run(
   ...     sites, outdir=None, save_report=False, save_plots=False,
   ...     on_step=lambda sr: seen.append(sr.step_code),
   ... )
   >>> seen
   ['NR001', 'FREQ001']

A step that raises under ``on_step_error="raise"`` never gets a
``StepResult`` at all, so it never fires the callback either -- exactly
consistent with the fact that it never appears in
``result.step_results``. An exception raised *inside* the callback itself
is caught and turned into a warning rather than aborting the run:

.. code-block:: pycon
   :linenos:

   >>> def broken(sr):
   ...     raise RuntimeError("bug in my callback")
   >>> import warnings
   >>> with warnings.catch_warnings(record=True) as caught:
   ...     warnings.simplefilter("always")
   ...     result = pipe.run(
   ...         sites, outdir=None, save_report=False, save_plots=False,
   ...         on_step=broken,
   ...     )
   >>> result.ok
   True
   >>> "bug in my callback" in str(caught[0].message)
   True

This is the hook a notebook cell, a custom script, or a future GUI
integration would use to watch a run as it happens rather than only after
it returns -- see :doc:`caching` for the ``cached`` field this hook already
exposes on every ``StepResult``.

.. _pipeline-observability-live:

Live Terminal Progress
----------------------

``configure_pipe(progress_style="rich")`` -- or the CLI's ``--live`` flag --
replaces the static progress bar with a live-updating
:class:`rich.table.Table` (both ``rich`` and ``tqdm`` are hard pyCSAMT
dependencies, so this needs no optional-import fallback). Every step's row
starts ``pending``, flips to ``running`` right before its transform starts,
then rewrites in place with the final status once it completes:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/AMT/WILLY_DATA/L22PLT \
       --steps NR001,FREQ001,FREQ004 --no-edi --no-report --no-plots --live

Captured against the real 25-station WILLY L22 line:

.. code-block:: text
   :linenos:

                        cli_pipeline  (3 step(s))
   ┌───┬─────────────────┬─────────┬────────┬───────┬────────┬────────┐
   │ # │ Step            │ Code    │ Status │  Time │ Sites  │ Cached │
   ├───┼─────────────────┼─────────┼────────┼───────┼────────┼────────┤
   │ 1 │ notch_powerline │ NR001   │ OK     │ 8.62s │ 25->25 │        │
   │ 2 │ select_band     │ FREQ001 │ OK     │ 0.23s │ 25->25 │        │
   │ 3 │ align_grid      │ FREQ004 │ OK     │ 0.82s │ 25->25 │        │
   └───┴─────────────────┴─────────┴────────┴───────┴────────┴────────┘
   Pipeline 'cli_pipeline' finished  [OK]  9.88s  plots=0

The table is rebuilt from scratch on every repaint rather than mutated in
place (``rich.table.Table`` has no supported in-place cell-update API) --
the standard idiom for driving a changing table through
:class:`rich.live.Live`. This is independent of ``on_step``: supplying both
at once works, and neither interferes with the other.

``progress_style="rich"`` is a new opt-in rendering choice, not a new
default -- ``"bar"`` stays the default, so nobody's existing terminal output
changes without asking. A found-by-testing detail worth knowing if you
extend this table yourself: avoid Unicode arrows (``→``) inside cell text.
``rich``'s legacy-Windows console renderer (a non-UTF-8 codepage, e.g. plain
``cmd.exe``) raises a ``UnicodeEncodeError`` on ``U+2192`` specifically --
the site-count column uses a plain ``->`` for exactly this reason, caught
while building this page.

.. _pipeline-observability-notebook:

Inline Notebook Rendering
-------------------------

A finished ``PipelineResult`` renders as a compact HTML table when it is
the last expression in a Jupyter cell -- Jupyter calls ``_repr_html_()``
automatically:

.. code-block:: pycon
   :linenos:

   >>> "<table>" in result._repr_html_()
   True
   >>> "NR001" in result._repr_html_()
   True

Captured output (real 3-station run, reformatted for readability -- the
actual string has no line breaks):

.. code-block:: html
   :linenos:

   <div class="pycsamt-result">
     <div class="title">PipelineResult 'repr_demo' — <span class="ok">OK</span></div>
     <b>Sites:</b> 3→3 | <b>Steps:</b> 2 (2 ok, 0 err) | <b>Time:</b> 0.05s
     <table>
       <tr><th>#</th><th>Step</th><th>Code</th><th>Status</th><th>Time</th><th>Sites</th><th>Cached</th></tr>
       <tr><td>1</td><td>notch</td><td>NR001</td><td><span class="ok">OK</span></td><td>0.03s</td><td>3→3</td><td>–</td></tr>
       <tr><td>2</td><td>band</td><td>FREQ001</td><td><span class="ok">OK</span></td><td>0.02s</td><td>3→3</td><td>–</td></tr>
     </table>
   </div>

This is a quick-glance view, not a replacement for the full
``report.html`` written by ``save_report=True`` -- no plot thumbnails, no
embedded config YAML. Every CSS rule is scoped under ``.pycsamt-result``
rather than bare selectors like ``table {...}``: the HTML is injected
directly into a notebook output cell's DOM, not rendered as a standalone
document, so an unscoped rule would leak out and restyle every other
table/output on the page -- a real correctness constraint this rendering
path has to respect that a standalone ``report.html`` file does not.

.. _pipeline-observability-history:

Run History
-----------

``Pipeline.run(..., history=True)`` appends one compact JSON line
summarizing the run to a log file (default
``~/.pycsamt/pipeline_history.jsonl``, configurable via
:attr:`~pycsamt.api.pipe.PipelineAPIConfig.history_path`); a path argument
logs there instead. This is deliberately a flat JSON-lines file, not a
database -- it exists to answer "how did my last N runs compare," not to be
a query engine:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import load_history
   >>> import tempfile
   >>> hist_file = tempfile.mkdtemp() + "/history.jsonl"
   >>> _ = pipe.run(sites, outdir=None, save_report=False, save_plots=False, history=hist_file)
   >>> _ = pipe.run(sites, outdir=None, save_report=False, save_plots=False, history=hist_file)
   >>> records = load_history(hist_file)
   >>> len(records)
   2
   >>> [s["code"] for s in records[-1]["steps"]]
   ['NR001', 'FREQ001']

From the CLI, ``--history`` (default location) or ``--history-file PATH``
logs a run; ``pycsamt pipe history`` reads it back:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/AMT/WILLY_DATA/L22PLT \
       --steps NR001,FREQ001 --no-edi --no-report --no-plots \
       --history-file /tmp/pycsamt_history.jsonl

   pycsamt pipe history --file /tmp/pycsamt_history.jsonl

Captured output of the second command:

.. code-block:: text
   :linenos:

   Logged 1 pipeline run(s):
     2026-08-14T10:58:14Z  cli_pipeline         OK                9.75s  sites 25→25

A malformed line -- e.g. from a write interrupted mid-flush -- is skipped
by :func:`~pycsamt.pipeline.load_history` rather than raising, the same
"one bad entry is a non-event" posture :class:`~pycsamt.pipeline.StepCache`
already uses; and a history-write failure at the end of a run is caught and
warned rather than failing the run itself -- a logging feature must not be
able to take down the pipeline it's logging.

Troubleshooting
---------------

``on_step`` never fires
    The step it should have fired for never got a ``StepResult`` -- check
    whether ``on_step_error="raise"`` and this step raised; that step (and
    nothing after it) is skipped for both ``result.step_results`` and
    ``on_step`` alike.

The live table looks garbled or falls back to plain lines
    A non-interactive or very narrow terminal (redirected output, some CI
    log viewers) makes ``rich`` degrade its rendering -- this is ``rich``'s
    own behaviour, not a pyCSAMT bug. Redirect to a file and read it back
    with a real terminal, or use ``progress_style="log"`` for CI logs.

``pycsamt pipe history`` shows nothing
    History logging is opt-in per run -- confirm the run that should have
    logged actually passed ``--history``/``--history-file`` (or
    ``history=`` in Python), and that ``pipe history --file`` points at the
    same location that run used.

Related Pages
-------------

* :doc:`concepts` explains the run lifecycle these hooks observe.
* :doc:`caching` explains ``StepResult.cached``, surfaced by every
  observability feature on this page.
* :doc:`cli_pipe` explains ``pycsamt pipe run`` and ``pycsamt pipe history``
  in full.
