.. _faq:

:html_theme.sidebar_secondary.remove: true

Frequently asked questions
==========================

.. raw:: html

   <div class="pyc-faq" data-faq-root>
     <header class="pyc-faq-hero">
       <span class="pyc-faq-eyebrow"><i class="fa-regular fa-circle-question" aria-hidden="true"></i> Help centre</span>
       <h2>Find the right answer, then keep moving</h2>
       <p>Search by symptom, command, file format, or workflow. Answers lead to the detailed guide when you need more than the short version.</p>
       <label class="pyc-faq-search">
         <span class="visually-hidden">Search frequently asked questions</span>
         <i class="fa-solid fa-magnifying-glass" aria-hidden="true"></i>
         <input type="search" data-faq-search placeholder="Try “EDI”, “static shift”, “Occam2D”, or an error message…" autocomplete="off">
         <kbd>/</kbd>
       </label>
       <div class="pyc-faq-chips" role="group" aria-label="Filter questions by topic">
         <button type="button" class="is-active" data-faq-filter="all" aria-pressed="true">All</button>
         <button type="button" data-faq-filter="start" aria-pressed="false">Getting started</button>
         <button type="button" data-faq-filter="data" aria-pressed="false">Data &amp; QC</button>
         <button type="button" data-faq-filter="processing" aria-pressed="false">Processing</button>
         <button type="button" data-faq-filter="inversion" aria-pressed="false">Inversion</button>
         <button type="button" data-faq-filter="solvers" aria-pressed="false">Solvers &amp; binaries</button>
         <button type="button" data-faq-filter="agents" aria-pressed="false">AI &amp; agents</button>
         <button type="button" data-faq-filter="apps" aria-pressed="false">Apps</button>
         <button type="button" data-faq-filter="errors" aria-pressed="false">Troubleshooting</button>
       </div>
       <p class="pyc-faq-count" data-faq-count aria-live="polite"></p>
     </header>

     <section class="pyc-faq-start" aria-labelledby="faq-start-title">
       <div>
         <span class="pyc-faq-kicker">New to pyCSAMT?</span>
         <h2 id="faq-start-title">Choose your shortest path</h2>
       </div>
       <div class="pyc-faq-paths">
         <a href="getting_started/quickstart.html"><i class="fa-solid fa-rocket" aria-hidden="true"></i><span><b>Run my first survey</b><small>Install, load, inspect, and plot</small></span></a>
         <a href="getting_started/data_formats.html"><i class="fa-solid fa-file-waveform" aria-hidden="true"></i><span><b>Check a data format</b><small>EDI, AVG, J, TDEM, and more</small></span></a>
         <a href="applications/index.html"><i class="fa-solid fa-display" aria-hidden="true"></i><span><b>Choose an interface</b><small>Python, CLI, desktop, or web</small></span></a>
       </div>
     </section>

     <main class="pyc-faq-list" data-faq-list>
       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-flag-checkered" aria-hidden="true"></i> Getting started</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="start" open>
           <summary><span>What is the quickest way to verify that pyCSAMT works?</span><span class="pyc-faq-badge pyc-faq-badge--popular">Popular</span></summary>
           <div class="pyc-faq-answer"><p>Install pyCSAMT in a clean environment, confirm the import, and run the bundled first-survey workflow. This separates installation problems from problems in your own field data.</p><p class="pyc-faq-next"><b>Next:</b> <a href="getting_started/quickstart.html">follow the quickstart</a> or review <a href="getting_started/installation.html">installation options</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="start apps">
           <summary><span>Should I use Python, the CLI, desktop app, or web app?</span><span class="pyc-faq-badge">Choose</span></summary>
           <div class="pyc-faq-answer"><p>Use Python for reproducible research and custom pipelines; the CLI for automation and batch jobs; the desktop app for interactive local work; and the web app for a browser-based workflow. They expose the same project concepts, so you can start visually and move to code later.</p><p class="pyc-faq-next"><b>Compare:</b> <a href="applications/overview.html">interfaces at a glance</a> and the <a href="cli/overview.html">CLI overview</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="start errors">
           <summary><span>Why does an optional feature fail after the base package installs?</span><span class="pyc-faq-badge pyc-faq-badge--setup">Setup</span></summary>
           <div class="pyc-faq-answer"><p>Some plotting, geospatial, AI, and application features use optional dependencies. Read the error for the missing package or extra, install only the feature set you need, then restart the Python process so imports are refreshed.</p><p class="pyc-faq-next"><b>Next:</b> see <a href="getting_started/installation.html">installation and optional dependencies</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="start errors">
           <summary><span>Do I need Numba or joblib installed?</span><span class="pyc-faq-badge">Optional</span></summary>
           <div class="pyc-faq-answer"><p>No. The <code>perf</code> extra (<code>python -m pip install "pycsamt[perf]"</code>) adds optional Numba/joblib acceleration for <code>pycsamt.models.occam1d</code>, batch survey agents, and pipeline execution. Everything that can use them degrades gracefully to a slower pure-Python/NumPy path when they are absent -- <code>perf</code> is a convenience, never a hard requirement, and it is not part of the <code>full</code> extra.</p><p class="pyc-faq-next"><b>Details:</b> <a href="installation.html#optional-feature-groups">the optional feature groups table</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="start data">
           <summary><span>Why does a table-returning function give me an APIFrame instead of a plain pandas DataFrame?</span><span class="pyc-faq-badge pyc-faq-badge--popular">Key concept</span></summary>
           <div class="pyc-faq-answer"><p>Most dataframe-returning functions accept an <code>api</code> keyword that defaults to <code>api=None</code>, which resolves through the package-wide <code>PYCSAMT_API_VIEW</code> switch -- and that switch defaults to <code>"pycsamt"</code>, not <code>"pandas"</code>. A bare call is therefore already wrapped into an <code>APIFrame</code> (or <code>APIResult</code> for multi-table workflows) unless you pass <code>api=False</code> for the raw <code>DataFrame</code>, or disable wrapping globally with <code>configure_api_view(backend=False)</code>.</p><p class="pyc-faq-next"><b>Next:</b> read <a href="api_guide/views.html">API views</a> for the full behaviour, including custom wrappers.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-robot" aria-hidden="true"></i> Direct AI inversion versus inversion agents</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="agents inversion start">
           <summary><span>What is the difference between pycsamt.ai inversion and an inversion agent?</span><span class="pyc-faq-badge pyc-faq-badge--popular">Key concept</span></summary>
           <div class="pyc-faq-answer"><p><code>pycsamt.ai.inversion</code> is the scientific model API: you construct inverter classes, datasets, networks, losses, training loops, checkpoints, and predictions directly. <code>pycsamt.agents</code> is the workflow layer around those classes. An inversion agent can load survey files, build features, train or load the underlying inverter, predict, calculate selected diagnostics, save figures, preserve warnings, and return a standard <code>AgentResult</code>.</p><p class="pyc-faq-next"><b>Compare:</b> <a href="user_guide/ai_inversion/agents.html">AI inversion agents</a> and the <a href="user_guide/agents/overview.html#when-to-use-agents">agent decision guide</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="agents inversion">
           <summary><span>When should I call pycsamt.ai directly?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Research control</span></summary>
           <div class="pyc-faq-answer"><p>Use the direct API when you need control of array contracts, architecture, geological priors, forward operators, loss terms, optimizer, data splits, batching, callbacks, calibration, checkpointing, or an experimental training loop. It is also the better layer for developing and testing a new method because the scientific choices remain explicit instead of being hidden behind a workflow contract.</p><p class="pyc-faq-next"><b>Start:</b> review <a href="user_guide/ai_inversion/index.html">the AI inversion guide</a> and <a href="user_guide/ai_inversion/data_preparation.html">data contracts</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="agents inversion apps">
           <summary><span>When is an inversion agent the better interface?</span><span class="pyc-faq-badge">Workflow</span></summary>
           <div class="pyc-faq-answer"><p>Use an agent when its documented input contract matches your task and you want a repeatable end-to-end step with standard status, warnings, error hints, timing, diagnostics, figures, and output paths. Agents are especially useful in coordinated QC-to-inversion pipelines, application interfaces, batch surveys, and reviewed production runs. They reduce orchestration code; they do not remove the need to record configuration or inspect the science.</p><p class="pyc-faq-next"><b>Learn:</b> <a href="user_guide/agents/overview.html">agent concepts and execution contracts</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="agents inversion processing">
           <summary><span>What does Inv2DAgent do, and when should I use it?</span><span class="pyc-faq-badge">U-Net 2-D</span></summary>
           <div class="pyc-faq-answer"><p><code>Inv2DAgent</code> wraps <code>EMInverter2D</code> in a profile workflow. It consumes the complete station-frequency panel and predicts a depth-by-station section with a U-Net-style model, so it can learn lateral continuity instead of treating every station independently. Use it when the stations form a defensible profile and the synthetic training distribution represents the expected 2-D structures. Its learned smoothness and quick data-space diagnostic are not equivalent to a classical 2-D EM roughness objective or full forward solve.</p><p class="pyc-faq-next"><b>Review:</b> <a href="user_guide/agents/ai_model_zoo_agents.html#inv2dagent">Inv2DAgent assumptions and example</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="agents inversion processing">
           <summary><span>What does Inv3DAgent do, and is it a full 3-D inversion?</span><span class="pyc-faq-badge pyc-faq-badge--warning">Know the scope</span></summary>
           <div class="pyc-faq-answer"><p><code>Inv3DAgent</code> wraps <code>GCNInverter3D</code>. It builds or accepts a graph from station coordinates so neighbouring stations exchange information, predicts layered resistivity at each graph node, and can estimate Monte Carlo dropout spread. It is a spatial graph-based AI inversion, not the same numerical problem as a full 3-D Maxwell inversion such as ModEM. Use it for rapid spatial prediction only after reviewing coordinates, graph edges, training coverage, uncertainty, and a physics-based baseline.</p><p class="pyc-faq-next"><b>Review:</b> <a href="user_guide/agents/ai_model_zoo_agents.html#inv3dagent">Inv3DAgent graph and uncertainty guidance</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="agents inversion errors">
           <summary><span>Do inversion agents require an LLM or make the scientific decision for me?</span><span class="pyc-faq-badge pyc-faq-badge--warning">Human review</span></summary>
           <div class="pyc-faq-answer"><p>No. The numerical workflow can be deterministic without an LLM; a configured language model may add narrative interpretation or help orchestrate a request, but it does not validate the inversion. An agent status of <code>success</code> only means the software workflow completed. Review warnings, input QC, training provenance, out-of-distribution checks, response residuals, uncertainty, dimensionality assumptions, and independent geological evidence before accepting the result.</p><p class="pyc-faq-next"><b>Validate:</b> <a href="user_guide/ai_inversion/validation.html">scientific validation</a> and <a href="user_guide/ai_inversion/agents.html#agentresult">AgentResult semantics</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="agents inversion">
           <summary><span>Will direct AI inversion and an agent always produce the same model?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Reproducibility</span></summary>
           <div class="pyc-faq-answer"><p>Only when both paths use the same survey ordering, frequency grid, features, scaling, architecture, weights or training data, random seeds, hyperparameters, checkpoint, inference mode, and post-processing. Agents may apply documented defaults and preprocessing that differ from a custom notebook. For a fair comparison, export the agent configuration and inverter object metadata, then reproduce the same inputs through the direct API.</p><p class="pyc-faq-next"><b>Trace:</b> see <a href="user_guide/ai_inversion/agents.html#common-execution-contract">the common execution contract</a>.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-microchip" aria-hidden="true"></i> Solvers, compilation, and binaries</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion start">
           <summary><span>Do ModEM, Occam2D, and MARE2DEM come with pyCSAMT?</span><span class="pyc-faq-badge pyc-faq-badge--warning">Important</span></summary>
           <div class="pyc-faq-answer"><p>pyCSAMT provides builders, native-file readers and writers, runners, and result loaders, but it does not ship pre-compiled executables. Occam2D and ModEM are external Fortran programs whose source is vendored for supported builds; MARE2DEM source is downloaded separately under its own licence. You can prepare and validate inputs without running a solver.</p><p class="pyc-faq-next"><b>Start:</b> read <a href="user_guide/models/compilation.html">compiling external solvers</a> and <a href="user_guide/models/overview.html">the model-integration lifecycle</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion">
           <summary><span>How do I compile ModEM 2-D or 3-D?</span><span class="pyc-faq-badge pyc-faq-badge--setup">Fortran</span></summary>
           <div class="pyc-faq-answer"><p>Use <code>pycsamt build modem2d --auto-install</code> for <code>Mod2DMT</code> or <code>pycsamt build modem3d --auto-install</code> for <code>Mod3DMT</code>. The build needs <code>gfortran</code>, <code>make</code>, and LAPACK/BLAS. On Windows, pyCSAMT can create an isolated conda toolchain; Linux and macOS use their normal package managers. MPI and Intel builds are opt-in and require those compilers to be installed already.</p><p class="pyc-faq-next"><b>Commands and platform notes:</b> <a href="user_guide/models/compilation.html#modem">ModEM compilation</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion">
           <summary><span>How do I compile Occam2D?</span><span class="pyc-faq-badge pyc-faq-badge--setup">Fortran</span></summary>
           <div class="pyc-faq-answer"><p>Run <code>pycsamt build occam2d --auto-install -y</code>. The build script selects a modern <code>gfortran</code> instead of the legacy <code>f90</code> name in the original Makefile and places the resulting <code>Occam2D</code> or <code>Occam2D.exe</code> in the vendored source directory. Use <code>--clean</code> when compiler changes or stale objects make a rebuild necessary.</p><p class="pyc-faq-next"><b>Details:</b> <a href="user_guide/models/compilation.html#occam2d">Occam2D compilation</a> and the <a href="user_guide/models/occam2d.html#run-occam2d">runner workflow</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion errors">
           <summary><span>Why is compiling MARE2DEM different?</span><span class="pyc-faq-badge pyc-faq-badge--warning">MPI + MKL</span></summary>
           <div class="pyc-faq-answer"><p>MARE2DEM source is not vendored and its build requires the Intel MPI Fortran/C toolchain plus MKL/ScaLAPACK. Use <code>pycsamt build mare2dem</code> to inspect status, then <code>pycsamt build mare2dem --auto-install -y</code> only after Intel oneAPI is installed and sourced. Native Windows builds are unsupported; use Linux, macOS, WSL2, or an HPC system.</p><p class="pyc-faq-next"><b>Prepare:</b> <a href="user_guide/models/compilation.html#mare2dem">MARE2DEM prerequisites</a> and <a href="user_guide/models/mare2dem.html#source-and-binary-management">source management</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion errors">
           <summary><span>How does pyCSAMT find and run a compiled binary?</span><span class="pyc-faq-badge pyc-faq-badge--popular">Popular</span></summary>
           <div class="pyc-faq-answer"><p>Prefer an explicit executable path in the engine configuration: <code>OccamRunner(binary_path=...)</code>, <code>ModEmConfig(binary_2d=..., binary_3d=...)</code>, or <code>Mare2DEMConfig(binary=...)</code>. The runner launches the program as a subprocess inside the native run directory, captures logs, and leaves native outputs for the result loader. Some runners also search the run directory, system <code>PATH</code>, and pyCSAMT build locations, but an explicit path records better provenance.</p><p class="pyc-faq-next"><b>Configure:</b> <a href="user_guide/models/configuration_and_io.html">configuration and executable I/O</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion errors">
           <summary><span>The binary was compiled, but pyCSAMT cannot run it. What should I verify?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Checklist</span></summary>
           <div class="pyc-faq-answer"><p>Confirm the exact path and executable permission, run the binary directly once, and check architecture plus runtime libraries. On Windows keep the copied MinGW DLLs beside the executable; for MPI verify <code>mpirun</code>, process count, and environment initialization. Also confirm that required native input filenames are present in the working directory. A successful compile does not prove that the runtime environment or model inputs are valid.</p><p class="pyc-faq-next"><b>Diagnose:</b> use the checks in <a href="user_guide/models/compilation.html">the compilation guide</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers processing">
           <summary><span>How should I choose between Occam2D, ModEM, and MARE2DEM?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Decision</span></summary>
           <div class="pyc-faq-answer"><p>Choose Occam2D for a smooth MT/AMT/CSAMT profile where 2-D geology and strike assumptions are defensible. Choose ModEM for native 2-D or 3-D MT workflows, especially area surveys and covariance-controlled 3-D models. Choose MARE2DEM for 2.5-D finite-element MT/CSEM problems requiring adaptive triangular meshes, detailed topography, transmitter geometry, or an established MPI workflow.</p><p class="pyc-faq-next"><b>Compare:</b> <a href="user_guide/models/choosing_backend.html">the backend decision matrix</a>.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-wave-square" aria-hidden="true"></i> Maxwell forward datasets and Python solvers</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers processing data">
           <summary><span>When should I generate a 2-D Maxwell forward dataset?</span><span class="pyc-faq-badge">2-D</span></summary>
           <div class="pyc-faq-answer"><p>Use a 2-D dataset when lateral structure along a profile matters and a layered 1-D response is no longer adequate. Build a <code>Grid2D</code>, run <code>MT2DForward</code> for TE and TM responses, and vary anomaly geometry, resistivity, station positions, frequencies, and noise with controlled seeds. This is appropriate for method tests, survey design, and training data whose assumptions are explicitly 2-D.</p><p class="pyc-faq-next"><b>Build:</b> <a href="user_guide/forward/synthetic_datasets.html#d-and-quasi-3-d-solver-datasets">2-D solver datasets</a> and <a href="user_guide/forward/solvers_and_grids.html#d-mt-solver">the MT2D solver</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers processing data">
           <summary><span>What is the difference between a pseudo-3-D dataset and an MT3DForward dataset?</span><span class="pyc-faq-badge pyc-faq-badge--warning">Know the physics</span></summary>
           <div class="pyc-faq-answer"><p><code>generate_dataset_3d</code> creates multi-station spatial datasets from local 1-D MT responses. <code>MT3DForward</code> uses a 3-D resistivity grid and, in its quasi-3-D mode, solves orthogonal 2-D slices to approximate tensor responses. Neither should be labelled a validated production full-3-D Maxwell result. Use ModEM or another full-3-D engine when full 3-D coupling is part of the scientific claim.</p><p class="pyc-faq-next"><b>Understand:</b> <a href="user_guide/forward/synthetic_datasets.html#pseudo-3-d-survey-datasets">pseudo-3-D datasets</a> and <a href="user_guide/forward/solvers_and_grids.html#quasi-3-d-solver">quasi-3-D limitations</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers processing">
           <summary><span>When should I use MT2DForward versus MT3DForward?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Decision</span></summary>
           <div class="pyc-faq-answer"><p>Use <code>MT2DForward</code> for profile-scale TE/TM physics, lateral anomaly experiments, and comparison with a 2-D inversion. Use <code>MT3DForward</code> for rapid volume sensitivity studies, station-grid survey design, or synthetic AI catalogues where its quasi-3-D approximation is acceptable. Move to a production external solver when off-profile current flow and full 3-D coupling must be resolved quantitatively.</p><p class="pyc-faq-next"><b>Compare:</b> <a href="user_guide/forward/solvers_and_grids.html">solvers and grids</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers data processing">
           <summary><span>What must I save with a 2-D or 3-D synthetic dataset?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Reproducibility</span></summary>
           <div class="pyc-faq-answer"><p>Save the resistivity model or generator parameters, grid including padding and air cells, station geometry, frequency axis, components and array layout, solver and version, boundary settings, random seeds, noise model, units, and train/validation/test split rules. Validate finite values and shapes and plot selected responses before training. A feature matrix without its physics and geometry metadata is not a reusable Maxwell dataset.</p><p class="pyc-faq-next"><b>Audit:</b> <a href="user_guide/forward/synthetic_datasets.html#dataset-quality-checks">dataset quality checks</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion">
           <summary><span>When should I use SimPEG?</span><span class="pyc-faq-badge">Python native</span></summary>
           <div class="pyc-faq-answer"><p>Use SimPEG when you need a Python-native research workflow with explicit meshes, mappings, regularization, optimization, sensitivities, and the ability to customize or differentiate the inverse problem. It is a strong choice for experimentation and integration with scientific Python. Expect to manage optional dependency versions and to validate the chosen electromagnetic simulation against a trusted reference.</p><p class="pyc-faq-next"><b>Decide:</b> <a href="user_guide/models/choosing_backend.html">compare Python and native backends</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion">
           <summary><span>When should I use pyGIMLi?</span><span class="pyc-faq-badge">Python ecosystem</span></summary>
           <div class="pyc-faq-answer"><p>Use pyGIMLi when your workflow benefits from its modelling and inversion framework, especially 1-D or stitched-profile experiments, TDEM support, mesh tools, or integration with other geophysical methods already handled in that ecosystem. Like SimPEG, it is an optional backend: confirm that the installed version supports the method and dimensionality you intend to report.</p><p class="pyc-faq-next"><b>Decide:</b> <a href="user_guide/models/choosing_backend.html">the backend matrix</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="solvers inversion start">
           <summary><span>Should I begin with SimPEG, pyGIMLi, or a compiled Fortran solver?</span><span class="pyc-faq-badge pyc-faq-badge--popular">Rule of thumb</span></summary>
           <div class="pyc-faq-answer"><p>Begin with pyCSAMT's built-in backend for a smoke test. Choose SimPEG or pyGIMLi when Python-level control and extensibility matter most. Choose Occam2D, ModEM, or MARE2DEM when a validated native workflow, established file format, HPC execution, or direct comparison with an existing project matters most. Solver choice does not replace dimensionality analysis, error modelling, or independent validation.</p><p class="pyc-faq-next"><b>Next:</b> <a href="user_guide/models/choosing_backend.html#quick-recommendation">use the quick recommendation table</a>.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-database" aria-hidden="true"></i> Data and quality control</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="data start">
           <summary><span>Which survey file formats can pyCSAMT read?</span><span class="pyc-faq-badge">Reference</span></summary>
           <div class="pyc-faq-answer"><p>pyCSAMT supports the common electromagnetic survey formats documented in the data-format matrix, including EDI, Zonge AVG, Jones J, TDEM, and modelling formats. Check the matrix before converting data: keeping the native format usually preserves more metadata.</p><p class="pyc-faq-next"><b>Check:</b> <a href="getting_started/data_formats.html">supported formats and capabilities</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="data errors">
           <summary><span>My EDI files load, but stations appear in the wrong order. What should I check?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Pro tip</span></summary>
           <div class="pyc-faq-answer"><p>Do not assume filenames define survey order. Verify station identifiers, coordinates, profile geometry, and units; then use the site ordering tools. Plot the station map before processing because a plausible-looking pseudosection can still hide a geometry error.</p><p class="pyc-faq-next"><b>Diagnose:</b> <a href="user_guide/site/index.html">site and station tools</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="data processing">
           <summary><span>What quality checks should I run before correction or inversion?</span><span class="pyc-faq-badge pyc-faq-badge--popular">Popular</span></summary>
           <div class="pyc-faq-answer"><p>Check station geometry, frequency overlap, missing components, apparent-resistivity and phase continuity, error estimates, outliers, and repeatability. Save the unmodified import and record every exclusion: correction should address a diagnosed problem, not make a curve merely look smooth.</p><p class="pyc-faq-next"><b>Workflow:</b> <a href="user_guide/data_loading.html">load and audit data</a>, then use the <a href="user_guide/emtools/index.html">processing catalogue</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="data processing inversion start">
           <summary><span>I have raw Stratagem/Zonge hardware field data, not EDI -- is there a complete worked example?</span><span class="pyc-faq-badge">Tutorial</span></summary>
           <div class="pyc-faq-answer"><p>Yes. A full tutorial covers a real field survey from raw Stratagem hardware output through <code>StratagemRawReader</code>, injecting surveyed coordinates with <code>CoordinateInjector</code>, static-shift correction, frequency filtering, noise removal, QC export, and an Occam2D inversion using the vendored, compilable solver. It also shows cross-checking the correction against an independent tool and against pyCSAMT's own log-frequency smoothing.</p><p class="pyc-faq-next"><b>Follow:</b> <a href="tutorials/process_stratagem_dafang_to_inversion.html">processing Stratagem field data to inversion</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="data start">
           <summary><span>Does pyCSAMT read modern EMTF-XML files, or only classic EDI?</span><span class="pyc-faq-badge">New</span></summary>
           <div class="pyc-faq-answer"><p>Both, through the same boundary. <code>Site</code>/<code>Sites</code> now wrap either format via a lazy dual backend: whichever representation was not natively supplied is materialized from the other on first access and cached, so every existing EDI-only accessor (<code>z</code>, <code>tipper</code>, <code>rho</code>, <code>phase</code>) keeps working unchanged. <code>ensure_sites</code> normalizes a single <code>.xml</code> file, an <code>EMTF</code> object, or a directory mixing <code>*.edi</code> and <code>*.xml</code> files into the same <code>Sites</code> container -- no separate function to learn, and nothing downstream needs to know which format a station arrived in.</p><p class="pyc-faq-next"><b>Read:</b> <a href="user_guide/data_loading.html#data-loading-xml">loading EMTF-XML the same way</a> and <a href="user_guide/site/containers.html">the full Site/Sites XML reference</a>.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-tower-broadcast" aria-hidden="true"></i> Airborne EM (ZTEM, AFMAG, MobileMT)</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="data start">
           <summary><span>Does pyCSAMT support airborne EM surveys like ZTEM, AFMAG, or MobileMT?</span><span class="pyc-faq-badge">New</span></summary>
           <div class="pyc-faq-answer"><p>Yes. <code>pycsamt.airborne</code> is a technology-neutral data model built directly on <code>EMTF</code> documents -- flight lines, datasets, a technology/format registry, and structural QC -- with no EDI bridge and no impedance requirement, since a genuine airborne measurement rarely has one. Three technology subpackages (<code>pycsamt.airborne.afmag</code>, <code>.ztem</code>, <code>.mobilemt</code>) map each system's decoded response onto that shared model.</p><p class="pyc-faq-next"><b>Start:</b> <a href="user_guide/airborne/index.html">the airborne EM guide</a>, beginning with <a href="user_guide/airborne/overview.html">the data-model overview</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="data processing">
           <summary><span>How do I read a directory of airborne EMTF-XML files?</span><span class="pyc-faq-badge">Tool</span></summary>
           <div class="pyc-faq-answer"><p><code>pycsamt.airborne.site.ensure_asites</code> is the airborne counterpart of <code>ensure_sites</code>: point it at a directory of EMTF-XML files and it returns an <code>AirborneSites</code> collection with each station's technology auto-detected from its document, ready for the same kind of selection, mapping, and export operations ground <code>Sites</code> already support.</p><p class="pyc-faq-next"><b>See it read real sample surveys:</b> <a href="user_guide/airborne/site.html">the airborne site view</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="data errors">
           <summary><span>Why is an airborne station's <code>z</code> always <code>None</code>?</span><span class="pyc-faq-badge pyc-faq-badge--warning">Not a bug</span></summary>
           <div class="pyc-faq-answer"><p>ZTEM, AFMAG, and MobileMT have no electric-field channel, so there is no impedance to build. An <code>AirborneSite</code> deliberately never fabricates one -- check <code>has_component("tipper")</code> or <code>"admittance"</code> instead of assuming <code>z</code> will populate. ZTEM/AFMAG carry a tipper or interstation tensor; MobileMT carries a ground-electric-to-airborne-magnetic admittance tensor, a genuinely different physical quantity from either.</p><p class="pyc-faq-next"><b>Understand the four response shapes:</b> <a href="user_guide/airborne/site.html#four-response-families">one container, four response families</a>.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-sliders" aria-hidden="true"></i> Processing and interpretation</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="processing data">
           <summary><span>How do I choose a correction method?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Decision</span></summary>
           <div class="pyc-faq-answer"><p>Start from the observed failure mode and survey physics. Diagnose dimensionality, strike, static shift, cultural noise, and frequency-local outliers separately. Apply the least invasive method that targets the evidence, compare before and after, and retain parameters in a reproducible pipeline.</p><p class="pyc-faq-next"><b>Compare:</b> <a href="user_guide/emtools/index.html">the EM tools catalogue</a> and <a href="theory/index.html">scientific background</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="processing data errors">
           <summary><span>Should I smooth noisy curves before inversion?</span><span class="pyc-faq-badge pyc-faq-badge--warning">Caution</span></summary>
           <div class="pyc-faq-answer"><p>Not automatically. First distinguish isolated outliers, coherent cultural noise, static shift, and genuine geological structure. Over-smoothing can remove useful signal and produce unjustified confidence. Preserve raw data, document exclusions, and propagate realistic error floors into inversion.</p><p class="pyc-faq-next"><b>Learn:</b> <a href="user_guide/emtools/index.html">diagnostics and corrections</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="processing data">
           <summary><span>Which pyCSAMT tool actually performs frequency-domain smoothing, and how do I judge whether it helped?</span><span class="pyc-faq-badge">Tool</span></summary>
           <div class="pyc-faq-answer"><p><code>pycsamt.emtools.remove_noise.smooth_logfreq</code> applies a triangular or box kernel along the log-frequency axis. Judge it by comparing the smoothed curve against both the unmodified import and, where available, an independent processing of the same stations -- agreement between two different corrections is a better signal than a curve that merely looks less noisy.</p><p class="pyc-faq-next"><b>See it applied:</b> <a href="tutorials/process_stratagem_dafang_to_inversion.html#static-shift-frequency-filtering-and-noise-removal">the static shift and smoothing section</a> of the Stratagem tutorial.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="processing inversion">
           <summary><span>How do I know whether an interpretation is defensible?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Expert</span></summary>
           <div class="pyc-faq-answer"><p>A defensible interpretation connects data quality, assumptions, model sensitivity, uncertainty, and independent geological evidence. Compare alternative processing choices and models; report what is resolved and what is inferred. A low misfit alone does not make a model unique.</p><p class="pyc-faq-next"><b>Next:</b> use the <a href="user_guide/interpretation/index.html">interpretation guide</a>.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-chart-line" aria-hidden="true"></i> Inversion and modelling</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="inversion start">
           <summary><span>Which inversion backend should I choose?</span><span class="pyc-faq-badge">Choose</span></summary>
           <div class="pyc-faq-answer"><p>Choose from the survey dimensionality, data type, compute environment, licensing constraints, and the outputs you must defend—not from the solver name alone. Begin with the simplest model compatible with the data, then increase complexity only when diagnostics justify it.</p><p class="pyc-faq-next"><b>Compare:</b> <a href="user_guide/models/index.html">solver integrations</a> and the <a href="user_guide/inversion/index.html">inversion workflow</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="inversion errors">
           <summary><span>Why can pyCSAMT prepare an inversion but not run it?</span><span class="pyc-faq-badge pyc-faq-badge--setup">External tool</span></summary>
           <div class="pyc-faq-answer"><p>Input preparation and result inspection are built into pyCSAMT, while some numerical engines are separate executables. Confirm that the selected engine is installed, licensed where required, and discoverable on your system. The prepared files remain useful even when the executable is unavailable.</p><p class="pyc-faq-next"><b>Next:</b> check <a href="user_guide/models/index.html">backend setup</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="inversion processing">
           <summary><span>When should I use AI inversion?</span><span class="pyc-faq-badge pyc-faq-badge--pro">Advanced</span></summary>
           <div class="pyc-faq-answer"><p>Use AI methods when training coverage, validation design, uncertainty reporting, and deployment constraints are explicit. Always benchmark against a classical method and test domain shift. For a new survey or limited training data, a classical baseline is usually the more interpretable starting point.</p><p class="pyc-faq-next"><b>Evaluate:</b> <a href="user_guide/ai_inversion/index.html">AI inversion guidance</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="inversion solvers start">
           <summary><span>Do I need an external Occam1D binary for 1-D inversion?</span><span class="pyc-faq-badge">New</span></summary>
           <div class="pyc-faq-answer"><p>No, not unless you want one. <code>pycsamt.models.occam1d</code> is a native, pure Python/NumPy 1-D Occam smooth-model engine -- forward model, analytic Jacobian, roughness regularization, and the nonlinear Occam loop are all implemented directly, with no compiled dependency. <code>Occam1DRunner</code> remains available for driving an external <code>Occam1D</code>-compatible executable instead, if you already have one and prefer it.</p><p class="pyc-faq-next"><b>Walk through it:</b> <a href="user_guide/models/occam1d.html">configuration, single-station, and batch inversion</a>.</p></div>
         </details>
       </section>

       <section class="pyc-faq-group" data-faq-group>
         <h2><i class="fa-solid fa-display" aria-hidden="true"></i> Applications and troubleshooting</h2>

         <details class="pyc-faq-item" data-faq-item data-topic="apps data">
           <summary><span>Can I move a project between the desktop, web, and Python interfaces?</span><span class="pyc-faq-badge">Workflow</span></summary>
           <div class="pyc-faq-answer"><p>Yes, when you preserve the project inputs, configuration, and exported outputs rather than relying only on interface state. Treat the project folder as the handoff boundary and record version information when reproducibility matters.</p><p class="pyc-faq-next"><b>See:</b> <a href="applications/index.html">application guides</a>.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="errors apps start">
           <summary><span>What information should I include when reporting a problem?</span><span class="pyc-faq-badge pyc-faq-badge--popular">Best practice</span></summary>
           <div class="pyc-faq-answer"><p>Include the pyCSAMT and Python versions, operating system, interface used, exact command or minimal code, complete traceback, expected and actual behaviour, and a small anonymised input when possible. Remove credentials and confidential coordinates before sharing.</p><p class="pyc-faq-next"><b>Report:</b> search <a href="https://github.com/earthai-tech/pycsamt/issues">existing issues</a> before opening a new one.</p></div>
         </details>

         <details class="pyc-faq-item" data-faq-item data-topic="errors">
           <summary><span>The answer is not here—where should I look next?</span><span class="pyc-faq-badge">Help</span></summary>
           <div class="pyc-faq-answer"><p>Use the global documentation search for the exact exception, class, command, or file suffix. Then check the relevant application troubleshooting page and existing issue reports. If it is reproducible and still unresolved, open an issue with the diagnostic bundle described above.</p><p class="pyc-faq-next"><b>Continue:</b> <a href="applications/index.html">application help</a>, <a href="api/index.html">API reference</a>, or <a href="https://github.com/earthai-tech/pycsamt/issues">issue tracker</a>.</p></div>
         </details>
       </section>

       <div class="pyc-faq-empty" data-faq-empty hidden>
         <i class="fa-solid fa-binoculars" aria-hidden="true"></i>
         <h2>No matching question yet</h2>
         <p>Try fewer words or search the complete documentation.</p>
         <button type="button" data-faq-clear>Clear FAQ search</button>
       </div>
     </main>

     <aside class="pyc-faq-help">
       <div><span class="pyc-faq-kicker">Still blocked?</span><h2>Turn a problem into a useful report</h2><p>Copy the exact error, reduce the workflow to the smallest failing example, and include versions plus a safe sample of the input.</p></div>
       <div class="pyc-faq-help-actions"><a class="pyc-faq-primary" href="https://github.com/earthai-tech/pycsamt/issues"><i class="fa-brands fa-github" aria-hidden="true"></i> Search issues</a><a href="search.html"><i class="fa-solid fa-magnifying-glass" aria-hidden="true"></i> Search all docs</a></div>
     </aside>
   </div>
