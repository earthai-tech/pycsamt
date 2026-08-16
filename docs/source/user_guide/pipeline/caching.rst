.. _pipeline-caching:

Caching And Resume
==================

Every pipeline step transforms a survey through a well-defined chain:
:math:`S_0 \to S_1 \to \cdots \to S_n`, where each arrow is one registered
step applied with one exact set of parameters. If two runs share the same
:math:`S_0` and the same sequence of `(code, params)` pairs up to step
:math:`i`, step :math:`i` will produce the identical :math:`S_i` both times
-- so there is no reason to recompute it a second time. ``Pipeline.run(...,
cache=True)`` turns that observation into a real feature: a
:term:`step cache` on disk, keyed by a :term:`chain hash` of everything a
step's output actually depends on.

This single mechanism is also what makes an interrupted run resumable.
Rerunning the identical command after a crash replays every
already-completed step from the cache -- instantly -- and only recomputes
whatever didn't finish. There is no separate "resume" API; resuming *is*
caching, applied to a rerun. This is the same idea Make, Docker layer
caching, and Snakemake all build on.

Enabling The Cache
------------------

From Python, ``cache`` accepts ``True`` (use the configured default
location), an explicit path, or ``False`` (the default -- off, identical
behaviour to every earlier pyCSAMT release):

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import Pipeline, Step, configure_pipe
   >>> from pycsamt.emtools import ensure_sites
   >>> configure_pipe(show_progress=False)
   >>> import glob
   >>> paths = sorted(glob.glob("data/AMT/WILLY_DATA/L22PLT/*.edi"))[:2]
   >>> sites = ensure_sites(paths)
   >>> pipe = Pipeline([("notch", Step("NR001")), ("band", Step("FREQ001"))])
   >>> import tempfile
   >>> cache_dir = tempfile.mkdtemp()
   >>> r1 = pipe.run(sites, outdir=None, save_report=False, save_plots=False, cache=cache_dir)
   >>> [sr.cached for sr in r1.step_results]
   [False, False]
   >>> r2 = pipe.run(sites, outdir=None, save_report=False, save_plots=False, cache=cache_dir)
   >>> [sr.cached for sr in r2.step_results]
   [True, True]

From the CLI, ``--cache`` uses the default location
(``~/.pycsamt/pipeline_cache``); ``--cache-dir PATH`` uses an explicit one
(and implies ``--cache``):

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis --preset basic_qc --cache
   pycsamt pipe run data/edis --preset basic_qc --cache-dir /scratch/pycsamt_cache

Configure a project-wide default instead of passing ``--cache-dir`` every
time via :attr:`~pycsamt.api.pipe.PipelineAPIConfig.cache_root`:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import configure_pipe
   >>> configure_pipe(cache_root="/scratch/pycsamt_cache")

What Gets Cached
----------------

A step's cache key is :func:`~pycsamt.pipeline.chain_key` applied to three
things: a content fingerprint of the sites flowing *into* the step, the
step's registry code, and its exact merged parameters. The fingerprint
itself, :func:`~pycsamt.pipeline.fingerprint_sites`, hashes each station's
name plus its frequency and impedance arrays, in order:

.. code-block:: pycon
   :linenos:

   >>> from pycsamt.pipeline import fingerprint_sites
   >>> fingerprint_sites(sites) == fingerprint_sites(ensure_sites(paths))
   True

Reloading the same two EDI files independently produces the identical
fingerprint -- this is the property the whole cache depends on. Because the
key chains forward (step 2's key depends on step 1's *output* fingerprint,
not the original input), changing step 1's parameters automatically
invalidates every step after it, while leaving step 1's own cache entry
untouched for the next run that still wants the old parameters.

Diagnostic steps (``returns_sites=False``, e.g. the QC snapshots or the
opt-in ``AI001`` audit from :doc:`extending`) are cache-participating too --
a hit just skips re-invoking the function, which matters when that function
costs a real torch import and a real computation.

Why Reruns Resume
-----------------

Because caching only ever stores a step's output *after* it succeeds, a
crash mid-step leaves no entry for that step -- the next identical run
recomputes exactly the part that was lost, and nothing more:

.. code-block:: console
   :linenos:

   # First attempt dies partway through, e.g. killed at step 4 of 6
   pycsamt pipe run data/edis --config workflow.yaml --cache
   ^C

   # Rerun the identical command: steps 1-3 replay from cache in a fraction
   # of a second, step 4 onward actually execute
   pycsamt pipe run data/edis --config workflow.yaml --cache

This also means two *different* pipeline configs that happen to share their
first few steps and parameters reuse each other's cached results for that
shared prefix -- there is nothing pipeline-specific baked into the key,
only step code, params, and upstream data.

Real Captured Example
---------------------

Captured against the real 25-station WILLY L22 line, ``--no-plots`` to
isolate transform time from figure-generation time:

.. code-block:: console
   :linenos:

   pycsamt pipe run data/AMT/WILLY_DATA/L22PLT \
       --steps NR001,FREQ001,FREQ004 --cache-dir /tmp/pycsamt_cache \
       --no-plots --format json

First run (cold cache) -- every step computed:

.. code-block:: json
   :linenos:

   {
     "steps": [
       {"code": "NR001",   "elapsed_sec": 11.416, "cached": false},
       {"code": "FREQ001", "elapsed_sec": 0.307,  "cached": false},
       {"code": "FREQ004", "elapsed_sec": 0.664,  "cached": false}
     ]
   }

Identical command run again (warm cache) -- every step replayed:

.. code-block:: json
   :linenos:

   {
     "steps": [
       {"code": "NR001",   "elapsed_sec": 1.106, "cached": true},
       {"code": "FREQ001", "elapsed_sec": 0.308, "cached": true},
       {"code": "FREQ004", "elapsed_sec": 0.188, "cached": true}
     ]
   }

``NR001``'s cache-hit time (1.1s) is still real work -- deserialising the
cached station data with ``joblib`` -- just far less than the 11.4s a real
power-line notch filter costs on 25 stations. ``StepResult.cached`` also
shows up in the text/HTML reports (a ``Cached`` column and badge
respectively) and in ``--format csv``.

What Is Not Safe To Cache
-------------------------

Caching is opt-in, not automatic, for a reason: it can only be correct for
a step whose output is a pure function of its input sites and its params.

A step that is not safe to cache
    Reads wall-clock time, a random seed without a fixed value, an
    environment variable, or a relative path with contents that can change
    between runs. A cache hit would silently replay a stale result -- there
    is no general way to detect this, so it isn't guarded against; disable
    ``--cache`` for a pipeline built from such a step.

Export steps (:doc:`steps`, ``EXPORT00x``)
    Cache-eligible like any other step, but a hit skips
    ``InputBuilder.build()`` -- if the files it would have written were
    deleted or moved outside the pipeline between runs, a cache hit will
    not recreate them. Either don't delete pipeline-managed export output,
    or disable caching for a run where that matters.

No cross-process coordination
    Unlike :class:`pycsamt.forward.maxwell.cache.MaxwellResultCache`, the
    step cache has no file locking and no size-based eviction -- it is a
    disposable, regeneratable local cache aimed at one CLI/script user
    iterating on a pipeline, not a shared multi-worker cache. Delete the
    cache directory to reclaim space; nothing tracks how large it has grown.

Inspecting And Clearing The Cache
---------------------------------

:class:`~pycsamt.pipeline.StepCache` is the same object ``cache=`` builds
internally, usable directly:

.. code-block:: pycon
   :linenos:

   >>> from pathlib import Path
   >>> from pycsamt.pipeline import StepCache
   >>> cache = StepCache(cache_dir)
   >>> cache.root == Path(cache_dir)
   True
   >>> cache.clear()

There is no dedicated ``pycsamt pipe cache`` CLI command -- clearing the
default location is a normal directory removal:

.. code-block:: console
   :linenos:

   rm -rf ~/.pycsamt/pipeline_cache

Troubleshooting
---------------

A rerun isn't hitting the cache
    Confirm ``--cache``/``--cache-dir`` (or ``cache=`` in Python) is passed
    on *both* runs -- caching that was off during the first run has nothing
    to replay. Also confirm the input EDI files and every step's params are
    byte-for-byte the same; even a default that changed between pyCSAMT
    versions changes the key.

A cache entry looks corrupted
    :meth:`StepCache.get` treats an unreadable entry as a miss and warns
    rather than raising -- the run continues, recomputing that step. This
    can happen after an interrupted write from an unrelated process killed
    mid-cache-write; safe to ignore, the next successful write replaces it.

The cache directory keeps growing
    Expected -- there is no automatic eviction (see above). Delete it
    periodically, or point ``--cache-dir`` at a location you already manage
    with your own retention policy.

Related Pages
-------------

* :doc:`concepts` explains the run lifecycle this mechanism hooks into.
* :doc:`extending` explains the opt-in ``AI001`` step, whose real torch +
  computation cost is exactly the kind of thing caching is most worth
  enabling for.
* :doc:`steps` explains the export steps referenced in the caveats above.
* :doc:`observability` explains ``StepResult.cached`` surfaced through
  ``on_step``, the live table, and the notebook HTML repr.
