.. pyCSAMT documentation master file — scikit-learn-style landing page.
   Layout lives in _static/css/pycsamt-home.css, behaviour in
   _static/js/pycsamt-home.js. Sidebars are disabled for this page only
   (html_sidebars in conf.py + the metadata field below).

:html_theme.sidebar_secondary.remove: true

.. title:: pyCSAMT — Electromagnetic geophysics in Python

.. raw:: html

   <div class="pycsamt-home">

   <!-- ================================ HERO ================================ -->
   <section class="pyc-hero">
     <div class="pyc-wrap">
       <div class="pyc-hero-copy">
         <div class="pyc-hero-badges">
           <span class="pyc-chip pyc-chip--accent"><i class="fa-solid fa-tag"></i> v2.0 · LGPL-3.0</span>
           <span class="pyc-chip"><i class="fa-brands fa-python"></i> Python 3.9+</span>
           <span class="pyc-chip"><i class="fa-solid fa-code-branch"></i> Open source</span>
         </div>
         <h1 class="pyc-hero-title">pyCSAMT</h1>
         <p class="pyc-hero-sub">
           Electromagnetic geophysics in Python — built for
           <span class="pyc-rotator"
                 data-words='["CSAMT","AMT","MT","TDEM","CSEM"]'>CSAMT</span>
           surveys.
         </p>
         <p class="pyc-hero-text">
           One coherent, scikit-learn-inspired toolkit for the full EM survey
           lifecycle: multi-format data loading, quality control,
           industry-standard corrections, classical and physics-informed
           neural-network inversion, interpretation, mapping, and
           publication-ready reports.
         </p>
         <div class="pyc-hero-actions">
           <a class="pyc-btn pyc-btn--primary" href="getting_started/index.html">
             <i class="fa-solid fa-rocket"></i> Get started
           </a>
           <a class="pyc-btn pyc-btn--ghost" href="user_guide/index.html">
             <i class="fa-solid fa-book-open"></i> User guide
           </a>
           <a class="pyc-btn pyc-btn--ghost"
              href="https://github.com/earthai-tech/pycsamt"
              target="_blank" rel="noopener">
             <i class="fa-brands fa-github"></i> GitHub
           </a>
         </div>
         <div class="pyc-install">
           <span><span class="pyc-dollar">$</span> <code>pip install pycsamt</code></span>
           <button class="pyc-copy-btn" data-copy="pip install pycsamt"
                   type="button" aria-label="Copy install command">
             <i class="fa-regular fa-copy"></i>
           </button>
         </div>
       </div>

       <div class="pyc-hero-demo">
         <div class="pyc-terminal" role="img"
              aria-label="Terminal session showing a pyCSAMT processing pipeline">
           <div class="pyc-terminal-bar">
             <span></span><span></span><span></span>
             <em>pycsamt — survey L22</em>
           </div>
           <pre><span class="t-prompt">$</span> <span class="t-cmd">pip install pycsamt</span>
   <span class="t-prompt">$</span> <span class="t-cmd">python</span>
   <span class="t-com">>>></span> <span class="t-kw">from</span> pycsamt.pipeline <span class="t-kw">import</span> Pipeline, Step
   <span class="t-com">>>></span> pipe = <span class="t-fn">Pipeline</span>([
   <span class="t-com">...</span>     (<span class="t-str">"notch"</span>,        <span class="t-fn">Step</span>(<span class="t-str">"NR001"</span>, mains_hz=<span class="t-num">50</span>)),
   <span class="t-com">...</span>     (<span class="t-str">"static_shift"</span>, <span class="t-fn">Step</span>(<span class="t-str">"SS001"</span>)),
   <span class="t-com">...</span> ])
   <span class="t-com">>>></span> pipe.<span class="t-fn">run</span>(sites, outdir=<span class="t-str">"outputs/run01/"</span>)
   <span class="t-ok">&#10003;</span> 2 steps · 24 stations · manifest written</pre>
         </div>
       </div>
     </div>

     <div class="pyc-hero-waves" aria-hidden="true">
       <svg class="wave-a" viewBox="0 0 2400 90" preserveAspectRatio="none">
         <path fill="#3e65b0"
               d="M0 60 Q150 5 300 60 T600 60 T900 60 T1200 60 T1500 60 T1800 60 T2100 60 T2400 60 V90 H0 Z"/>
       </svg>
       <svg class="wave-b" viewBox="0 0 2400 90" preserveAspectRatio="none">
         <path fill="#f15a29"
               d="M0 65 Q100 25 200 65 T400 65 T600 65 T800 65 T1000 65 T1200 65 T1400 65 T1600 65 T1800 65 T2000 65 T2200 65 T2400 65 V90 H0 Z"/>
       </svg>
     </div>

     <div class="pyc-stats">
       <div class="pyc-stat">
         <b data-counter="9" data-suffix="+">9+</b>
         <span>data formats</span>
       </div>
       <div class="pyc-stat">
         <b data-counter="25">25</b>
         <span>correction methods</span>
       </div>
       <div class="pyc-stat">
         <b data-counter="3">3</b>
         <span>classical solvers</span>
       </div>
       <div class="pyc-stat">
         <b data-counter="4">4</b>
         <span>interfaces</span>
       </div>
     </div>
   </section>

   <!-- ============================ FEATURE BOXES =========================== -->
   <section class="pyc-section">
     <div class="pyc-wrap">
       <div class="pyc-center pyc-reveal">
         <p class="pyc-kicker">Capabilities</p>
         <h2 class="pyc-h2">Everything an EM survey needs</h2>
         <p class="pyc-lead">
           From raw field files to interpreted resistivity sections —
           each stage is a documented, tested module you can use on its own
           or chain into a reproducible pipeline.
         </p>
       </div>

       <div class="pyc-features">

         <div class="pyc-feature pyc-feature--blue pyc-reveal">
           <div class="pyc-feature-head">
             <span class="pyc-feature-icon"><i class="fa-solid fa-database"></i></span>
             <h3><a href="user_guide/data_loading.html">Data I/O &amp; quality control</a></h3>
           </div>
           <p>
             Load EDI, Zonge AVG, Jones J, TDEM, and MARE2DEM files into one
             site model. Inspect frequency coverage, flag noisy stations, and
             audit usable data before any processing step.
           </p>
           <div class="pyc-feature-tags">
             <a href="getting_started/data_formats.html">Formats</a>
             <a href="site/index.html">Site tools</a>
             <span>EDI</span><span>AVG</span><span>J</span><span>emdata</span>
           </div>
           <a class="pyc-feature-more" href="user_guide/data_loading.html">
             Learn more <i class="fa-solid fa-arrow-right"></i>
           </a>
         </div>

         <div class="pyc-feature pyc-feature--orange pyc-reveal">
           <div class="pyc-feature-head">
             <span class="pyc-feature-icon"><i class="fa-solid fa-sliders"></i></span>
             <h3><a href="user_guide/processing.html">Processing &amp; corrections</a></h3>
           </div>
           <p>
             A catalogue of 25 correction methods in six categories: notch
             filtering, static-shift removal, tensor rotation, phase-tensor
             analysis, and more — each with a stable identifier.
           </p>
           <div class="pyc-feature-tags">
             <a href="user_guide/processing.html">Catalogue</a>
             <span>Static shift</span><span>Notch</span><span>Phase tensor</span>
           </div>
           <a class="pyc-feature-more" href="user_guide/processing.html">
             Learn more <i class="fa-solid fa-arrow-right"></i>
           </a>
         </div>

         <div class="pyc-feature pyc-feature--gold pyc-reveal">
           <div class="pyc-feature-head">
             <span class="pyc-feature-icon"><i class="fa-solid fa-layer-group"></i></span>
             <h3><a href="forward/index.html">Forward modelling</a></h3>
           </div>
           <p>
             Build synthetic layered-earth and 2-D models, compute forward
             responses, add realistic noise, and generate datasets for survey
             design or inverter training.
           </p>
           <div class="pyc-feature-tags">
             <a href="forward/index.html">Synthetics</a>
             <a href="theory/index.html">Theory</a>
             <span>MT</span><span>CSAMT</span><span>TDEM</span>
           </div>
           <a class="pyc-feature-more" href="forward/index.html">
             Learn more <i class="fa-solid fa-arrow-right"></i>
           </a>
         </div>

         <div class="pyc-feature pyc-feature--blue pyc-reveal">
           <div class="pyc-feature-head">
             <span class="pyc-feature-icon"><i class="fa-solid fa-brain"></i></span>
             <h3><a href="user_guide/inversion.html">Inversion — classical &amp; AI</a></h3>
           </div>
           <p>
             Drive Occam2D, ModEM, and MARE2DEM end to end — input builders,
             runners, result loaders — or switch to physics-informed neural
             networks and hybrid deep-learning inverters in 1-D to 3-D.
           </p>
           <div class="pyc-feature-tags">
             <a href="models/index.html">Solvers</a>
             <a href="user_guide/ai_inversion.html">PINN</a>
             <span>Occam2D</span><span>ModEM</span><span>MARE2DEM</span>
           </div>
           <a class="pyc-feature-more" href="user_guide/inversion.html">
             Learn more <i class="fa-solid fa-arrow-right"></i>
           </a>
         </div>

         <div class="pyc-feature pyc-feature--orange pyc-reveal">
           <div class="pyc-feature-head">
             <span class="pyc-feature-icon"><i class="fa-solid fa-map"></i></span>
             <h3><a href="user_guide/interpretation.html">Interpretation &amp; mapping</a></h3>
           </div>
           <p>
             Classify resistivity, derive pseudostratigraphic logs, and render
             station maps, pseudosections, overlays, and 3-D quick-look views
             with the code-first MapView platform.
           </p>
           <div class="pyc-feature-tags">
             <a href="map/index.html">Map tools</a>
             <a href="user_guide/mapping.html">Mapping guide</a>
             <span>Pseudosection</span><span>3-D</span>
           </div>
           <a class="pyc-feature-more" href="user_guide/interpretation.html">
             Learn more <i class="fa-solid fa-arrow-right"></i>
           </a>
         </div>

         <div class="pyc-feature pyc-feature--gold pyc-reveal">
           <div class="pyc-feature-head">
             <span class="pyc-feature-icon"><i class="fa-solid fa-robot"></i></span>
             <h3><a href="pipeline/index.html">Pipelines, agents &amp; apps</a></h3>
           </div>
           <p>
             Define reproducible workflows in YAML, JSON, or Python; let
             LLM-driven agents orchestrate QC, inversion prep, and reporting;
             or work interactively in the web dashboard and desktop GUI.
           </p>
           <div class="pyc-feature-tags">
             <a href="pipeline/index.html">Pipeline</a>
             <a href="agents/index.html">AI agents</a>
             <a href="applications/index.html">Apps</a>
           </div>
           <a class="pyc-feature-more" href="pipeline/index.html">
             Learn more <i class="fa-solid fa-arrow-right"></i>
           </a>
         </div>

       </div>
     </div>
   </section>

   <!-- ============================ CODE IN ACTION ========================== -->
   <section class="pyc-section pyc-section--alt">
     <div class="pyc-wrap pyc-center pyc-reveal">
       <p class="pyc-kicker">Code in action</p>
       <h2 class="pyc-h2">A consistent API, from QC to report</h2>
       <p class="pyc-lead">
         Estimator-style objects, named pipeline steps, and one-call plotting —
         the same patterns across processing, inversion, and interpretation.
       </p>
     </div>
     <div class="pyc-wrap">

.. container:: pyc-code-panel pyc-reveal

   .. tab-set::

      .. tab-item:: Processing pipeline

         .. code-block:: python

            from pycsamt.pipeline import Pipeline, Step

            pipe = Pipeline([
                ("notch",        Step("NR001", mains_hz=50)),
                ("band",         Step("FREQ001")),
                ("static_shift", Step("SS001")),
                ("rotate",       Step("TZ001")),
            ])
            result = pipe.run(sites, outdir="outputs/run01/")
            print(result.summary())

      .. tab-item:: Inversion results

         .. code-block:: python

            from pycsamt.models.mare2dem import (
                InversionResult, PlotResponse, PlotModel,
            )

            result = InversionResult("runs/demo_mt/")
            PlotResponse(result).plot(max_rx=6, savefig="response.pdf")
            PlotModel(result).plot(cmap="turbo_r", savefig="section.pdf")

      .. tab-item:: PINN inversion

         .. code-block:: python

            from pycsamt.ai.inversion import PINN2D

            inv = PINN2D(n_layers=64, epochs=3000, backend="torch")
            model = inv.fit(sites, frequencies=freqs)
            model.plot_section()

      .. tab-item:: AI agent

         .. code-block:: python

            from pycsamt.agents import AgentMaster

            master = AgentMaster(provider="anthropic")
            report = master.run(
                "Load data/edi/, flag stations with RMS > 2, "
                "build an Occam2D input for profile L22, "
                "launch inversion, and produce a PDF report."
            )

      .. tab-item:: CLI

         .. code-block:: bash

            pycsamt survey set data/edi/
            pycsamt edi info data/edi/
            pycsamt invert build data/edi/ --solver occam2d --workdir runs/occam2d/
            pycsamt pipe run --config pipeline.yaml --survey data/edi/ --out outputs/run01/

.. raw:: html

     </div>
   </section>

   <!-- ============================= INTERFACES ============================= -->
   <section class="pyc-section">
     <div class="pyc-wrap">
       <div class="pyc-center pyc-reveal">
         <p class="pyc-kicker">Interfaces</p>
         <h2 class="pyc-h2">Work the way you prefer</h2>
         <p class="pyc-lead">
           The same engine behind four front ends — script it, automate it,
           or point and click.
         </p>
       </div>
       <div class="pyc-interfaces">
         <div class="pyc-interface pyc-reveal">
           <i class="fa-brands fa-python"></i>
           <h3><a href="api/index.html">Python API</a></h3>
           <p>Full programmatic control with a NumPy-documented reference.</p>
           <code>import pycsamt</code>
         </div>
         <div class="pyc-interface pyc-reveal">
           <i class="fa-solid fa-terminal"></i>
           <h3><a href="cli/index.html">Command line</a></h3>
           <p>Survey management, QC, inversion, and pipelines from the shell.</p>
           <code>pycsamt --help</code>
         </div>
         <div class="pyc-interface pyc-reveal">
           <i class="fa-solid fa-chart-line"></i>
           <h3><a href="applications/index.html">Web dashboard</a></h3>
           <p>Interactive Dash app for exploration, agent chat, and monitoring.</p>
           <code>pycsamt-web</code>
         </div>
         <div class="pyc-interface pyc-reveal">
           <i class="fa-solid fa-desktop"></i>
           <h3><a href="applications/index.html">Desktop GUI</a></h3>
           <p>Native point-and-click access to processing and inversion.</p>
           <code>pycsamt-desktop</code>
         </div>
       </div>
     </div>
   </section>

   <!-- =========================== NEWS + COMMUNITY ========================= -->
   <section class="pyc-section pyc-section--alt">
     <div class="pyc-wrap">
       <div class="pyc-duo">

         <div class="pyc-panel pyc-reveal">
           <h3><i class="fa-solid fa-bullhorn"></i> News</h3>
           <ul class="pyc-news-list">
             <li>
               <span class="pyc-news-date">v2.0.0</span>
               <div>
                 <a href="release_notes/v2.0.0.html">pyCSAMT v2 released</a>
                 <p>Ground-up rewrite: pipelines, AI inversion, agents,
                    web and desktop apps.</p>
               </div>
             </li>
             <li>
               <span class="pyc-news-date">Docs</span>
               <div>
                 <a href="release_notes/index.html">Release notes &amp; roadmap</a>
                 <p>v2 documentation is under active construction — APIs may
                    change before the stable release.</p>
               </div>
             </li>
             <li>
               <span class="pyc-news-date">Log</span>
               <div>
                 <a href="changelog.html">Full changelog</a>
                 <p>Every change, version by version.</p>
               </div>
             </li>
           </ul>
         </div>

         <div class="pyc-panel pyc-reveal">
           <h3><i class="fa-solid fa-users"></i> Community &amp; citation</h3>
           <ul class="pyc-community-links">
             <li><a href="https://github.com/earthai-tech/pycsamt"
                    target="_blank" rel="noopener">
               <i class="fa-brands fa-github"></i> Source code</a></li>
             <li><a href="https://github.com/earthai-tech/pycsamt/issues"
                    target="_blank" rel="noopener">
               <i class="fa-solid fa-bug"></i> Issue tracker</a></li>
             <li><a href="contributing.html">
               <i class="fa-solid fa-code-pull-request"></i> Contributing</a></li>
             <li><a href="references.html">
               <i class="fa-solid fa-graduation-cap"></i> References</a></li>
           </ul>
           <div class="pyc-cite-box">
             Using pyCSAMT in published research? Please cite
             <em>Kouadio&nbsp;et&nbsp;al.&nbsp;(2022)</em>, J.&nbsp;Applied
             Geophysics —
             <a href="https://doi.org/10.1016/j.jappgeo.2022.104647"
                target="_blank" rel="noopener">doi:10.1016/j.jappgeo.2022.104647</a>.
             See the <a href="references.html">references page</a> for BibTeX.
           </div>
         </div>

       </div>
     </div>
   </section>

   <!-- ============================ EXPLORE DOCS ============================ -->
   <section class="pyc-section">
     <div class="pyc-wrap">
       <div class="pyc-center pyc-reveal">
         <p class="pyc-kicker">Documentation map</p>
         <h2 class="pyc-h2">Explore the docs</h2>
       </div>
       <div class="pyc-explore pyc-reveal">
         <div>
           <h3>Learn</h3>
           <ul>
             <li><a href="getting_started/index.html">Getting started</a></li>
             <li><a href="getting_started/installation.html">Installation</a></li>
             <li><a href="getting_started/first_survey.html">First survey</a></li>
             <li><a href="user_guide/index.html">User guide</a></li>
             <li><a href="tutorials/index.html">Tutorials</a></li>
             <li><a href="examples/index.html">Examples gallery</a></li>
           </ul>
         </div>
         <div>
           <h3>Reference</h3>
           <ul>
             <li><a href="api/index.html">API reference</a></li>
             <li><a href="cli/index.html">CLI reference</a></li>
             <li><a href="theory/index.html">Scientific background</a></li>
             <li><a href="references.html">Bibliography</a></li>
           </ul>
         </div>
         <div>
           <h3>Systems</h3>
           <ul>
             <li><a href="pipeline/index.html">Pipeline system</a></li>
             <li><a href="agents/index.html">AI agents</a></li>
             <li><a href="forward/index.html">Forward modelling</a></li>
             <li><a href="models/index.html">Model integrations</a></li>
             <li><a href="site/index.html">Site tools</a></li>
             <li><a href="map/index.html">Map tools</a></li>
             <li><a href="applications/index.html">Applications</a></li>
           </ul>
         </div>
         <div>
           <h3>Project</h3>
           <ul>
             <li><a href="release_notes/index.html">Release notes</a></li>
             <li><a href="changelog.html">Changelog</a></li>
             <li><a href="contributing.html">Contributing</a></li>
             <li><a href="development/index.html">Development</a></li>
           </ul>
         </div>
       </div>
     </div>
   </section>

   </div><!-- /.pycsamt-home -->

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting started

   getting_started/index
   installation
   api_guide/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User guide

   user_guide/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Tutorials and examples

   tutorials/index
   examples/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: API reference

   api/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Interfaces

   cli/index
   agents/index
   pipeline/index
   forward/index
   models/index
   site/index
   map/index
   applications/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Scientific background

   theory/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   development/index
   release_notes/index
   changelog
   contributing
   references
