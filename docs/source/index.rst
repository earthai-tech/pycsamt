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
     <div class="pyc-hero-bg" aria-hidden="true">
       <div class="pyc-slide pyc-slide--a is-active"></div>
       <div class="pyc-slide pyc-slide--b"></div>
       <div class="pyc-slide pyc-slide--c"></div>
       <div class="pyc-hero-scrim"></div>
     </div>

     <div class="pyc-wrap pyc-hero-grid">
       <div class="pyc-hero-copy">
         <div class="pyc-hero-badges">
           <span class="pyc-chip pyc-chip--accent"><i class="fa-solid fa-tag"></i> v2.0 · LGPL-3.0</span>
           <span class="pyc-chip"><i class="fa-brands fa-python"></i> Python 3.9+</span>
           <span class="pyc-chip"><i class="fa-solid fa-code-branch"></i> Open source</span>
         </div>
         <h1 class="pyc-hero-logo" aria-label="pyCSAMT">
           <svg class="pyc-wm" viewBox="0 0 400.29 66.91"
                xmlns="http://www.w3.org/2000/svg" aria-hidden="true">
             <defs>
               <style>.cls-1{fill:#3e65b0;}.cls-2{fill:#fbb040;stroke:#f15a29;stroke-miterlimit:10;}</style>
               <linearGradient id="wm-sheen-grad" gradientUnits="userSpaceOnUse" x1="0" y1="0" x2="90" y2="55">
                 <stop offset="0" stop-color="#fff" stop-opacity="0"/>
                 <stop offset="0.5" stop-color="#fff" stop-opacity="0.65"/>
                 <stop offset="1" stop-color="#fff" stop-opacity="0"/>
                 <animateTransform attributeName="gradientTransform" type="translate"
                   values="-60 0; 460 0; 460 0" keyTimes="0; 0.6; 1" dur="4.5s" repeatCount="indefinite"/>
               </linearGradient>
             </defs>
             <g class="wm-letters">
               <g class="wm-l wm-l1"><path class="cls-1" d="M362.34,336.93h-11V312.56q0-6.72,3.52-10.93a14.45,14.45,0,0,1,5.14-3.87,15.67,15.67,0,0,1,17.42,3,14.5,14.5,0,0,1,4.48,10.78,14.72,14.72,0,0,1-4.45,10.72,14.34,14.34,0,0,1-10.55,4.52c-.66,0-1.64-.07-2.93-.21V314.46a4.32,4.32,0,0,0,2.7,1.08,3.9,3.9,0,0,0,2.91-1.26,4.15,4.15,0,0,0,1.22-3,3.94,3.94,0,0,0-1.26-2.95,4.2,4.2,0,0,0-3-1.21q-4.21,0-4.22,5.83Z" transform="translate(-351.38 -274.97)"/></g>
               <g class="wm-l wm-l2"><path class="cls-1" d="M390.52,338.92V329.4a10.25,10.25,0,0,0,3.9.94q5.32,0,6.91-6a10.93,10.93,0,0,1-5.62,1.7,8.5,8.5,0,0,1-6.71-3,11.32,11.32,0,0,1-2.64-7.72V296.91h11v15.7c0,1.74.72,2.61,2.17,2.61s2.08-1,2.08-2.87V296.91h11V323q0,7.29-4.11,11.77a15.45,15.45,0,0,1-11.83,5.33A16.69,16.69,0,0,1,390.52,338.92Z" transform="translate(-351.38 -274.97)"/></g>
               <g class="wm-l wm-l3"><path class="cls-2" d="M449,280.78l-.26,31.16.26,26.45a28.45,28.45,0,0,0,3.77.31,27.14,27.14,0,0,0,8.82-1.54,27.92,27.92,0,0,0,8-4.17,30.38,30.38,0,0,0,4.63-4.57,32.19,32.19,0,0,0,3.84-5.59c.56-1.05,1.16-1.58,1.8-1.58,1,0,1.45.41,1.45,1.23s-.73,2.37-2.18,4.46a37.71,37.71,0,0,1-4.76,5.65,30.74,30.74,0,0,1-13.23,7.52,32.7,32.7,0,0,1-8.79,1.27A31.31,31.31,0,0,1,440,338.74a30.55,30.55,0,0,1-15.12-13.22,31.88,31.88,0,0,1,13.71-44.65,32.24,32.24,0,0,1,35.64,5.49c2.93,2.72,4.4,4.69,4.4,5.89a1.19,1.19,0,0,1-.44.88,1.43,1.43,0,0,1-1,.39q-.87,0-2.1-1.71a26.46,26.46,0,0,0-10-8.28,29.16,29.16,0,0,0-13.08-3.06A14.39,14.39,0,0,0,449,280.78Z" transform="translate(-351.38 -274.97)"/></g>
               <g class="wm-l wm-l4"><path class="cls-2" d="M486.23,329.38c.53,0,1.14.4,1.85,1.19a24.55,24.55,0,0,0,8.06,5.78,22.58,22.58,0,0,0,9.43,2.17c3.42,0,5.14-1,5.14-3.12a7.9,7.9,0,0,0-1.8-4.46q-1.8-2.47-8.62-9.82-8.57-9.27-11.42-13.75-3.78-5.94-3.78-11.3,0-8.08,7.2-13.23t17.45-5a30,30,0,0,1,13.32,2.95,17.61,17.61,0,0,1,3.93,2.63c1.1,1,1.65,1.83,1.65,2.51a1.24,1.24,0,0,1-.46.9,1.42,1.42,0,0,1-.95.42c-.5,0-1.23-.46-2.2-1.37a18.61,18.61,0,0,0-6.24-3.71,20.82,20.82,0,0,0-7.42-1.47,6.2,6.2,0,0,0-3.58.94A2.85,2.85,0,0,0,506.4,284a8.53,8.53,0,0,0,1.71,4.4,55.59,55.59,0,0,0,5.8,6.9q12.93,13.76,15.83,18.54a19,19,0,0,1,3,9.76,14.39,14.39,0,0,1-2.15,7.55,17,17,0,0,1-6,5.85q-7.29,4.35-17.8,4.35a29.47,29.47,0,0,1-17-5c-3.42-2.34-5.14-4.17-5.14-5.49a1.48,1.48,0,0,1,.46-1A1.52,1.52,0,0,1,486.23,329.38Z" transform="translate(-351.38 -274.97)"/></g>
               <g class="wm-l wm-l5"><path class="cls-2" d="M570,276l35.14,63.59H574.67l-4.24-8H555.71l-13-.05-2.19,0-2.89,5.37a20.67,20.67,0,0,1-1.78,3,1.61,1.61,0,0,1-1.27.58q-1.5,0-1.5-1.32A5.86,5.86,0,0,1,534,337q.85-1.56,6-10.26l26.6-44.78Zm-13.67,28.87-7.68,12.92-6.14,10.64,13.11.17,13.15-.17Z" transform="translate(-351.38 -274.97)"/></g>
               <g class="wm-l wm-l6"><path class="cls-2" d="M670,277.53l9.57,32.74,4.35,14.5,3.9,13.14.53,1.67H661.66l-8.41-31.07L635.1,335.8l-18.62-25.94q-2.7,9.5-5.88,23.3a39.59,39.59,0,0,1-1.66,6.26,2.05,2.05,0,0,1-2,1c-1.29,0-1.93-.44-1.93-1.32a41.58,41.58,0,0,1,1.23-4.43l2.73-9.36,4.75-16.18L623,277.53,647,312.4Z" transform="translate(-351.38 -274.97)"/></g>
               <g class="wm-l wm-l7"><path class="cls-2" d="M734.13,282.54l-.22,13.75-.14,13,.14,30.24H707.1l.26-35-.13-11.91-.13-10.1-5.8.13-8.31.17-1.23.05c-1.4,0-2.11-.53-2.11-1.58a1.38,1.38,0,0,1,1.23-1.59,17.91,17.91,0,0,1,2.16-.13H747a13.76,13.76,0,0,1,3.34.27,1.44,1.44,0,0,1,.88,1.49c0,1-.75,1.49-2.24,1.49l-3.83-.13L735,282.54Z" transform="translate(-351.38 -274.97)"/></g>
             </g>
             <g class="wm-sheen" fill="url(#wm-sheen-grad)" aria-hidden="true">
               <path d="M362.34,336.93h-11V312.56q0-6.72,3.52-10.93a14.45,14.45,0,0,1,5.14-3.87,15.67,15.67,0,0,1,17.42,3,14.5,14.5,0,0,1,4.48,10.78,14.72,14.72,0,0,1-4.45,10.72,14.34,14.34,0,0,1-10.55,4.52c-.66,0-1.64-.07-2.93-.21V314.46a4.32,4.32,0,0,0,2.7,1.08,3.9,3.9,0,0,0,2.91-1.26,4.15,4.15,0,0,0,1.22-3,3.94,3.94,0,0,0-1.26-2.95,4.2,4.2,0,0,0-3-1.21q-4.21,0-4.22,5.83Z" transform="translate(-351.38 -274.97)"/>
               <path d="M390.52,338.92V329.4a10.25,10.25,0,0,0,3.9.94q5.32,0,6.91-6a10.93,10.93,0,0,1-5.62,1.7,8.5,8.5,0,0,1-6.71-3,11.32,11.32,0,0,1-2.64-7.72V296.91h11v15.7c0,1.74.72,2.61,2.17,2.61s2.08-1,2.08-2.87V296.91h11V323q0,7.29-4.11,11.77a15.45,15.45,0,0,1-11.83,5.33A16.69,16.69,0,0,1,390.52,338.92Z" transform="translate(-351.38 -274.97)"/>
               <path d="M449,280.78l-.26,31.16.26,26.45a28.45,28.45,0,0,0,3.77.31,27.14,27.14,0,0,0,8.82-1.54,27.92,27.92,0,0,0,8-4.17,30.38,30.38,0,0,0,4.63-4.57,32.19,32.19,0,0,0,3.84-5.59c.56-1.05,1.16-1.58,1.8-1.58,1,0,1.45.41,1.45,1.23s-.73,2.37-2.18,4.46a37.71,37.71,0,0,1-4.76,5.65,30.74,30.74,0,0,1-13.23,7.52,32.7,32.7,0,0,1-8.79,1.27A31.31,31.31,0,0,1,440,338.74a30.55,30.55,0,0,1-15.12-13.22,31.88,31.88,0,0,1,13.71-44.65,32.24,32.24,0,0,1,35.64,5.49c2.93,2.72,4.4,4.69,4.4,5.89a1.19,1.19,0,0,1-.44.88,1.43,1.43,0,0,1-1,.39q-.87,0-2.1-1.71a26.46,26.46,0,0,0-10-8.28,29.16,29.16,0,0,0-13.08-3.06A14.39,14.39,0,0,0,449,280.78Z" transform="translate(-351.38 -274.97)"/>
               <path d="M486.23,329.38c.53,0,1.14.4,1.85,1.19a24.55,24.55,0,0,0,8.06,5.78,22.58,22.58,0,0,0,9.43,2.17c3.42,0,5.14-1,5.14-3.12a7.9,7.9,0,0,0-1.8-4.46q-1.8-2.47-8.62-9.82-8.57-9.27-11.42-13.75-3.78-5.94-3.78-11.3,0-8.08,7.2-13.23t17.45-5a30,30,0,0,1,13.32,2.95,17.61,17.61,0,0,1,3.93,2.63c1.1,1,1.65,1.83,1.65,2.51a1.24,1.24,0,0,1-.46.9,1.42,1.42,0,0,1-.95.42c-.5,0-1.23-.46-2.2-1.37a18.61,18.61,0,0,0-6.24-3.71,20.82,20.82,0,0,0-7.42-1.47,6.2,6.2,0,0,0-3.58.94A2.85,2.85,0,0,0,506.4,284a8.53,8.53,0,0,0,1.71,4.4,55.59,55.59,0,0,0,5.8,6.9q12.93,13.76,15.83,18.54a19,19,0,0,1,3,9.76,14.39,14.39,0,0,1-2.15,7.55,17,17,0,0,1-6,5.85q-7.29,4.35-17.8,4.35a29.47,29.47,0,0,1-17-5c-3.42-2.34-5.14-4.17-5.14-5.49a1.48,1.48,0,0,1,.46-1A1.52,1.52,0,0,1,486.23,329.38Z" transform="translate(-351.38 -274.97)"/>
               <path d="M570,276l35.14,63.59H574.67l-4.24-8H555.71l-13-.05-2.19,0-2.89,5.37a20.67,20.67,0,0,1-1.78,3,1.61,1.61,0,0,1-1.27.58q-1.5,0-1.5-1.32A5.86,5.86,0,0,1,534,337q.85-1.56,6-10.26l26.6-44.78Zm-13.67,28.87-7.68,12.92-6.14,10.64,13.11.17,13.15-.17Z" transform="translate(-351.38 -274.97)"/>
               <path d="M670,277.53l9.57,32.74,4.35,14.5,3.9,13.14.53,1.67H661.66l-8.41-31.07L635.1,335.8l-18.62-25.94q-2.7,9.5-5.88,23.3a39.59,39.59,0,0,1-1.66,6.26,2.05,2.05,0,0,1-2,1c-1.29,0-1.93-.44-1.93-1.32a41.58,41.58,0,0,1,1.23-4.43l2.73-9.36,4.75-16.18L623,277.53,647,312.4Z" transform="translate(-351.38 -274.97)"/>
               <path d="M734.13,282.54l-.22,13.75-.14,13,.14,30.24H707.1l.26-35-.13-11.91-.13-10.1-5.8.13-8.31.17-1.23.05c-1.4,0-2.11-.53-2.11-1.58a1.38,1.38,0,0,1,1.23-1.59,17.91,17.91,0,0,1,2.16-.13H747a13.76,13.76,0,0,1,3.34.27,1.44,1.44,0,0,1,.88,1.49c0,1-.75,1.49-2.24,1.49l-3.83-.13L735,282.54Z" transform="translate(-351.38 -274.97)"/>
             </g>
           </svg>
         </h1>
         <p class="pyc-hero-sub">
           Electromagnetic geophysics in Python — built for
           <span class="pyc-rotator"
                 data-words='["CSAMT","AMT","MT","TDEM","CSEM","AFMAG","ZTEM","MobileMT"]'
                 data-airborne-words='["AFMAG","ZTEM","MobileMT"]'>CSAMT</span>
           surveys.
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
       </div>

       <aside class="pyc-hero-points" aria-label="Why pyCSAMT">
         <div class="pyc-point">
           <span class="pyc-point-icon"><i class="fa-solid fa-user-group"></i></span>
           <div class="pyc-point-text">
             <b>Built for researchers and practitioners</b>
             <span>One toolkit from field QC to publication figures.</span>
           </div>
         </div>
         <div class="pyc-point">
           <span class="pyc-point-icon"><i class="fa-solid fa-puzzle-piece"></i></span>
           <div class="pyc-point-text">
             <b>Modular, easy to extend and integrate</b>
             <span>Estimator-style objects you can script and chain.</span>
           </div>
         </div>
         <div class="pyc-point">
           <span class="pyc-point-icon"><i class="fa-solid fa-rocket"></i></span>
           <div class="pyc-point-text">
             <b>From simple scripts to large-scale studies</b>
             <span>The same API drives one site or a whole campaign.</span>
           </div>
         </div>
       </aside>
     </div>

     <div class="pyc-wrap pyc-hero-foot">
       <div class="pyc-hero-foot-row">
         <button type="button" class="pyc-flow-install" data-copy="pip install pycsamt"
                 aria-label="Copy install command: pip install pycsamt">
           <i class="fa-solid fa-terminal" aria-hidden="true"></i>
           <code><span class="pyc-flow-install-prompt">$</span> pip install pycsamt</code>
           <i class="fa-regular fa-copy pyc-flow-install-icon" aria-hidden="true"></i>
           <i class="fa-solid fa-check pyc-flow-install-icon pyc-flow-install-icon--ok" aria-hidden="true"></i>
         </button>

         <nav class="pyc-hero-flow" aria-label="Survey workflow">
           <a class="pyc-flow-step" href="user_guide/data_loading.html">
             <i class="fa-solid fa-wave-square"></i><span>Acquire</span>
           </a>
           <span class="pyc-flow-sep" aria-hidden="true"><i class="fa-solid fa-chevron-right"></i></span>
           <a class="pyc-flow-step" href="user_guide/emtools/index.html">
             <i class="fa-solid fa-sliders"></i><span>Process</span>
           </a>
           <span class="pyc-flow-sep" aria-hidden="true"><i class="fa-solid fa-chevron-right"></i></span>
           <a class="pyc-flow-step" href="user_guide/site/index.html">
             <i class="fa-solid fa-magnifying-glass-chart"></i><span>Analyze</span>
           </a>
           <span class="pyc-flow-sep" aria-hidden="true"><i class="fa-solid fa-chevron-right"></i></span>
           <a class="pyc-flow-step" href="user_guide/inversion/index.html">
             <i class="fa-solid fa-chart-line"></i><span>Invert</span>
           </a>
           <span class="pyc-flow-sep" aria-hidden="true"><i class="fa-solid fa-chevron-right"></i></span>
           <a class="pyc-flow-step" href="user_guide/interpretation/index.html">
             <i class="fa-solid fa-map"></i><span>Interpret</span>
           </a>
         </nav>
       </div>

       <div class="pyc-hero-dots" role="tablist" aria-label="Hero background slides">
         <button type="button" class="is-active" role="tab" aria-selected="true"
                 aria-label="Slide 1"></button>
         <button type="button" role="tab" aria-selected="false" aria-label="Slide 2"></button>
         <button type="button" role="tab" aria-selected="false" aria-label="Slide 3"></button>
       </div>
     </div>
   </section>

   <!-- ===================== ANIMATED STATS (below hero) =================== -->
   <div class="pyc-statsbar">
     <div class="pyc-wrap">
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
     </div>
   </div>

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
           <button type="button" class="pyc-feature-flip" aria-expanded="false">
             <i class="fa-solid fa-arrow-right-arrow-left" aria-hidden="true"></i>
             <span class="sr-only">Show links for Data I/O &amp; quality control</span>
           </button>
           <div class="pyc-feature-inner">
             <div class="pyc-feature-face pyc-feature-front">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-database"></i></span>
                 <h3>Data I/O &amp; quality control</h3>
               </div>
               <p>
                 Load EDI, Zonge AVG, Jones J, TDEM, and MARE2DEM files into one
                 site model. Inspect frequency coverage, flag noisy stations, and
                 audit usable data before any processing step.
               </p>
               <div class="pyc-feature-hero" aria-hidden="true">
                 <img src="_static/images/home/card-dataio.png"
                      alt="" loading="lazy" width="1280" height="600"/>
               </div>
             </div>
             <div class="pyc-feature-face pyc-feature-back">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-database"></i></span>
                 <h3><a href="user_guide/data_loading.html">Data I/O &amp; quality control</a></h3>
               </div>
               <div class="pyc-feature-tags">
                 <a href="getting_started/data_formats.html">Formats</a>
                 <a href="user_guide/site/index.html">Site tools</a>
                 <span>EDI</span><span>AVG</span><span>J</span><span>emdata</span>
               </div>
               <a class="pyc-feature-more" href="user_guide/data_loading.html">
                 Learn more <i class="fa-solid fa-arrow-right"></i>
               </a>
             </div>
           </div>
         </div>

         <div class="pyc-feature pyc-feature--orange pyc-reveal">
           <button type="button" class="pyc-feature-flip" aria-expanded="false">
             <i class="fa-solid fa-arrow-right-arrow-left" aria-hidden="true"></i>
             <span class="sr-only">Show links for Processing &amp; corrections</span>
           </button>
           <div class="pyc-feature-inner">
             <div class="pyc-feature-face pyc-feature-front">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-sliders"></i></span>
                 <h3>Processing &amp; corrections</h3>
               </div>
               <p>
                 A catalogue of 25 correction methods in six categories: notch
                 filtering, static-shift removal, tensor rotation, phase-tensor
                 analysis, and more.
               </p>
               <div class="pyc-feature-hero" aria-hidden="true">
                 <img src="_static/images/home/card-processing.png"
                      alt="" loading="lazy" width="1280" height="600"/>
               </div>
             </div>
             <div class="pyc-feature-face pyc-feature-back">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-sliders"></i></span>
                 <h3><a href="user_guide/emtools/index.html">Processing &amp; corrections</a></h3>
               </div>
               <div class="pyc-feature-tags">
                 <a href="user_guide/emtools/index.html">Catalogue</a>
                 <span>Static shift</span><span>Notch</span><span>Phase tensor</span>
               </div>
               <a class="pyc-feature-more" href="user_guide/emtools/index.html">
                 Learn more <i class="fa-solid fa-arrow-right"></i>
               </a>
             </div>
           </div>
         </div>

         <div class="pyc-feature pyc-feature--gold pyc-reveal">
           <button type="button" class="pyc-feature-flip" aria-expanded="false">
             <i class="fa-solid fa-arrow-right-arrow-left" aria-hidden="true"></i>
             <span class="sr-only">Show links for Forward modelling</span>
           </button>
           <div class="pyc-feature-inner">
             <div class="pyc-feature-face pyc-feature-front">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-layer-group"></i></span>
                 <h3>Forward modelling</h3>
               </div>
               <p>
                 Build synthetic layered-earth and 2-D models, compute forward
                 responses, add realistic noise, and generate datasets for survey
                 design or inverter training.
               </p>
               <div class="pyc-feature-hero" aria-hidden="true">
                 <img src="_static/images/home/card-forward.png"
                      alt="" loading="lazy" width="1280" height="600"/>
               </div>
             </div>
             <div class="pyc-feature-face pyc-feature-back">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-layer-group"></i></span>
                 <h3><a href="user_guide/forward/index.html">Forward modelling</a></h3>
               </div>
               <div class="pyc-feature-tags">
                 <a href="user_guide/forward/index.html">Synthetics</a>
                 <a href="theory/index.html">Theory</a>
                 <span>MT</span><span>CSAMT</span><span>TDEM</span>
               </div>
               <a class="pyc-feature-more" href="user_guide/forward/index.html">
                 Learn more <i class="fa-solid fa-arrow-right"></i>
               </a>
             </div>
           </div>
         </div>

         <div class="pyc-feature pyc-feature--blue pyc-reveal">
           <button type="button" class="pyc-feature-flip" aria-expanded="false">
             <i class="fa-solid fa-arrow-right-arrow-left" aria-hidden="true"></i>
             <span class="sr-only">Show links for Inversion — classical &amp; AI</span>
           </button>
           <div class="pyc-feature-inner">
             <div class="pyc-feature-face pyc-feature-front">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-brain"></i></span>
                 <h3>Inversion — classical &amp; AI</h3>
               </div>
               <p>
                 Drive Occam2D, ModEM, and MARE2DEM end to end — input builders,
                 runners, result loaders — or switch to physics-informed neural
                 networks and hybrid deep-learning inverters in 1-D to 3-D.
               </p>
               <div class="pyc-feature-hero" aria-hidden="true">
                 <img src="_static/images/home/card-inversion.png"
                      alt="" loading="lazy" width="1280" height="600"/>
               </div>
             </div>
             <div class="pyc-feature-face pyc-feature-back">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-brain"></i></span>
                 <h3><a href="user_guide/inversion/index.html">Inversion — classical &amp; AI</a></h3>
               </div>
               <div class="pyc-feature-tags">
                 <a href="user_guide/models/index.html">Solvers</a>
                 <a href="user_guide/ai_inversion/index.html">PINN</a>
                 <span>Occam2D</span><span>ModEM</span><span>MARE2DEM</span>
               </div>
               <a class="pyc-feature-more" href="user_guide/inversion/index.html">
                 Learn more <i class="fa-solid fa-arrow-right"></i>
               </a>
             </div>
           </div>
         </div>

         <div class="pyc-feature pyc-feature--orange pyc-reveal">
           <button type="button" class="pyc-feature-flip" aria-expanded="false">
             <i class="fa-solid fa-arrow-right-arrow-left" aria-hidden="true"></i>
             <span class="sr-only">Show links for Interpretation &amp; mapping</span>
           </button>
           <div class="pyc-feature-inner">
             <div class="pyc-feature-face pyc-feature-front">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-map"></i></span>
                 <h3>Interpretation &amp; mapping</h3>
               </div>
               <p>
                 Classify resistivity, derive pseudostratigraphic logs, and render
                 station maps, pseudosections, overlays, and 3-D quick-look views
                 with the code-first MapView platform, plus uncertainty layers
                 for reporting.
               </p>
               <div class="pyc-feature-hero" aria-hidden="true">
                 <img src="_static/images/home/card-interpretation.png"
                      alt="" loading="lazy" width="1280" height="600"/>
               </div>
             </div>
             <div class="pyc-feature-face pyc-feature-back">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-map"></i></span>
                 <h3><a href="user_guide/interpretation/index.html">Interpretation &amp; mapping</a></h3>
               </div>
               <div class="pyc-feature-tags">
                 <a href="user_guide/map/index.html">Map tools</a>
                 <a href="user_guide/map/index.html">Mapping guide</a>
                 <span>Pseudosection</span><span>3-D</span>
               </div>
               <a class="pyc-feature-more" href="user_guide/interpretation/index.html">
                 Learn more <i class="fa-solid fa-arrow-right"></i>
               </a>
             </div>
           </div>
         </div>

         <div class="pyc-feature pyc-feature--gold pyc-reveal">
           <button type="button" class="pyc-feature-flip" aria-expanded="false">
             <i class="fa-solid fa-arrow-right-arrow-left" aria-hidden="true"></i>
             <span class="sr-only">Show links for Pipelines, agents &amp; apps</span>
           </button>
           <div class="pyc-feature-inner">
             <div class="pyc-feature-face pyc-feature-front">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-robot"></i></span>
                 <h3>Pipelines, agents &amp; apps</h3>
               </div>
               <p>
                 Define reproducible workflows in YAML, JSON, or Python; let
                 LLM-driven agents orchestrate QC, inversion prep, and reporting;
                 or work interactively in the web dashboard and desktop GUI.
               </p>
               <div class="pyc-feature-hero" aria-hidden="true">
                 <img src="_static/images/home/card-pipeline.png"
                      alt="" loading="lazy" width="1280" height="600"/>
               </div>
             </div>
             <div class="pyc-feature-face pyc-feature-back">
               <div class="pyc-feature-head">
                 <span class="pyc-feature-icon"><i class="fa-solid fa-robot"></i></span>
                 <h3><a href="user_guide/pipeline/index.html">Pipelines, agents &amp; apps</a></h3>
               </div>
               <div class="pyc-feature-tags">
                 <a href="user_guide/pipeline/index.html">Pipeline</a>
                 <a href="user_guide/agents/index.html">AI agents</a>
                 <a href="applications/index.html">Apps</a>
               </div>
               <a class="pyc-feature-more" href="user_guide/pipeline/index.html">
                 Learn more <i class="fa-solid fa-arrow-right"></i>
               </a>
             </div>
           </div>
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
            :caption: Load and QC

            from pycsamt.api import read_edis
            from pycsamt.emtools import build_qc_table

            sites = read_edis("data/edi/")
            qc = build_qc_table(sites)     # station coverage, SNR, skew
            print(qc[["station", "n_freq", "snr_med"]])

         .. code-block:: python
            :caption: Named pipeline steps

            from pycsamt.pipeline import Pipeline, Step

            pipe = Pipeline([
                ("notch",        Step("NR001", mains_hz=50)),
                ("band",         Step("FREQ001")),
                ("static_shift", Step("SS001")),
                ("rotate",       Step("TZ001")),
            ])
            result = pipe.run(sites, outdir="outputs/run01/")
            print(result.summary())

         .. code-block:: python
            :caption: Reproducible reruns

            pipe.to_yaml("pipeline.yaml")       # share the exact chain

            from pycsamt.pipeline import Pipeline

            rerun = Pipeline.from_yaml("pipeline.yaml")
            result = rerun.run(sites, outdir="outputs/run02/")
            print(result.ok, f"{result.elapsed_sec:.1f}s")

      .. tab-item:: Inversion results

         .. code-block:: python
            :caption: EDIs to solver input

            from pycsamt.models.mare2dem import make_mt_data_from_edi

            make_mt_data_from_edi(
                "data/edi/L22/", "runs/demo_mt/L22.emdata",
                output_modes="all",
                error_floor_te=0.05, error_floor_tm=0.05,
            )

         .. code-block:: python
            :caption: Response and model plots

            from pycsamt.models.mare2dem import (
                InversionResult, PlotResponse, PlotModel,
            )

            result = InversionResult("runs/demo_mt/")
            PlotResponse(result).plot(max_rx=6, savefig="response.pdf")
            PlotModel(result).plot(cmap="turbo_r", savefig="section.pdf")

      .. tab-item:: AI inversion

         .. code-block:: python
            :caption: Train a 1-D inverter

            from pycsamt.forward.batch import generate_dataset
            from pycsamt.ai.inversion import EMInverter1D

            ds = generate_dataset(n_samples=2_000, n_layers=5, seed=0)
            inv = EMInverter1D(arch="resnet", n_layers=5)
            inv.fit(ds, epochs=30)

         .. code-block:: python
            :caption: PINN 2-D inversion

            from pycsamt.ai.inversion import PINN2D

            inv = PINN2D(n_layers=64, epochs=3000, backend="torch")
            model = inv.fit(sites, frequencies=freqs)
            model.plot_section()

         .. code-block:: python
            :caption: Uncertainty bands

            from pycsamt.ai import EnsembleInverter

            ens = EnsembleInverter(inv, n_estimators=8)
            mean, std = ens.predict_with_uncertainty(ds.X)

      .. tab-item:: AI agents

         .. code-block:: python
            :caption: One line, whole workflow

            from pycsamt.agents import AgentMaster

            master = AgentMaster(provider="anthropic")
            report = master.run(
                "Load data/edi/, flag stations with RMS > 2, "
                "build an Occam2D input for profile L22, launch "
                "inversion, and produce a PDF report."
            )

         .. code-block:: python
            :caption: One agent, one job

            from pycsamt.agents import MTLoaderAgent, DataQCAgent

            loaded = MTLoaderAgent().execute({"path": "data/edi/"})
            qc = DataQCAgent().execute({
                "sites": loaded["sites"],
                "output_dir": "outputs/qc/",
            })
            print(qc.summary)

         .. code-block:: python
            :caption: Plan from plain text

            from pycsamt.agents import WorkflowOrchestratorAgent

            agent = WorkflowOrchestratorAgent()
            plan = agent.execute({
                "request": "QC the EDI files and prepare a short report",
                "data_path": "data/edi/",
                "dry_run": True,        # preview the chain first
            })
            print(plan["workflow_type"], plan["reasoning"])

      .. tab-item:: CLI

         .. code-block:: bash
            :caption: Inspect a survey

            pycsamt survey set data/edi/
            pycsamt edi info data/edi/

         .. code-block:: bash
            :caption: Invert and pipeline

            pycsamt invert build data/edi/ --solver occam2d --workdir runs/occam2d/
            pycsamt pipe run --config pipeline.yaml --survey data/edi/ --out outputs/run01/

      .. tab-item:: IoT telemetry

         .. code-block:: python
            :caption: Edge QC on the recorder

            import numpy as np
            from pycsamt.iot import EdgeProcessor, simulate_amt_station

            station = simulate_amt_station("S01", n_samples=1024, seed=7)
            block = np.column_stack(list(station["data"].values()))
            edge = EdgeProcessor().process(
                block, channel_names=["ex", "ey", "hx", "hy"],
            )
            print(edge.decision, edge.reasons)

         .. code-block:: python
            :caption: Field-session dashboard

            from pycsamt.iot import (
                FieldSession, plot_field_dashboard, simulate_iot_network,
            )

            session = FieldSession("survey-2026", method="amt")
            session.add_packets(simulate_iot_network(n_stations=12, seed=7))

            status = session.assess()             # edge-QC + health roll-up
            plot_field_dashboard(session)         # live acquisition dashboard
            inputs = session.to_pipeline_input()  # hand off to processing

         .. code-block:: python
            :caption: Power budget

            from pycsamt.iot import EnergyConfig, estimate_energy_budget

            budget = estimate_energy_budget(EnergyConfig(
                battery_wh=84.0, active_power_w=1.2, duty_cycle=0.5,
                solar_wh_per_day=20.0, telemetry_power_w=2.5,
                telemetry_seconds_per_day=300.0,
            ))
            print(budget.state, budget.runtime_days)

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
             <li><a href="user_guide/pipeline/index.html">Pipeline system</a></li>
             <li><a href="user_guide/agents/index.html">AI agents</a></li>
             <li><a href="user_guide/forward/index.html">Forward modelling</a></li>
             <li><a href="user_guide/models/index.html">Model integrations</a></li>
             <li><a href="user_guide/site/index.html">Site tools</a></li>
             <li><a href="user_guide/map/index.html">Map tools</a></li>
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
             <li><a href="faq.html">FAQ &amp; help</a></li>
           </ul>
         </div>
       </div>
     </div>
   </section>

   </div><!-- /.pycsamt-home -->

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Install

   Install <install/index>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: API

   API <api_landing>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User guide

   user_guide/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Examples

   Examples <examples/index>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Agents

   user_guide/agents/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Interfaces

   cli/index
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
   resources

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Help

   FAQ <faq>
