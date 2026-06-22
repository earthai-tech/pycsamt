/**
 * analytics_carousel.js
 *
 * Auto-rotating carousel for the 5-slide Station Analytics panel.
 * Fully client-side — no Dash round-trips.
 *
 * Features
 * ─────────
 *  • Slides cycle automatically every HOLD_MS milliseconds
 *  • Each transition: current slide exits left, next slides in from right
 *  • A typewriter effect types a short description as the slide arrives
 *  • Dot indicators reflect the active slide; clicking a dot jumps to it
 *  • window resize event dispatched after slide-in so Plotly repaints correctly
 *  • Carousel starts when the analytics row becomes visible (data loaded)
 *    and resets when it is hidden again (data cleared)
 */
(function () {
  'use strict';

  /* ── Configuration ──────────────────────────────────────────────────── */
  var SLIDE_IDS = [
    'crs-slide-0', 'crs-slide-1', 'crs-slide-2',
    'crs-slide-3', 'crs-slide-4', 'crs-slide-5', 'crs-slide-6'
  ];
  var DOT_IDS = [
    'crs-dot-0', 'crs-dot-1', 'crs-dot-2',
    'crs-dot-3', 'crs-dot-4', 'crs-dot-5', 'crs-dot-6'
  ];
  var LABELS = [
    '1 / 7  ·  Tipper Coverage',
    '2 / 7  ·  Frequencies / Station',
    '3 / 7  ·  Quality Radar',
    '4 / 7  ·  Survey Geometry',
    '5 / 7  ·  Station Ranking',
    '6 / 7  ·  Coverage Heatmap',
    '7 / 7  ·  Per-Line Distribution'
  ];
  var DESCRIPTIONS = [
    'Tipper coverage breakdown. ' +
    'For multi-line surveys each bar shows one profile — ' +
    'stack height splits with-tipper (blue) vs without (grey). ' +
    'Good surveys aim for ≥ 50 % tipper coverage per line.',

    'Frequency count per station. ' +
    'In multi-line mode each line gets its own colour; ' +
    'bars are grouped so you can spot lines with sparse or missing records at a glance. ' +
    'The dashed line marks the overall median.',

    'Quality spider for the selected station ' +
    '(freq range · tipper · coordinates · data fill · elevation rank). ' +
    'With no selection and multiple lines, overlaid spiders compare average line quality. ' +
    'Click any station in the list to drill in.',

    'Spatial layout of the survey. ' +
    'Each profile is drawn in its own colour; dot size scales with frequency count. ' +
    'Lines connect stations within each profile in acquisition order.',

    'Stations ranked by total frequency coverage, best first. ' +
    'Colour bands highlight which survey line each station belongs to. ' +
    'A flat, uniform profile indicates consistent data acquisition.',

    'Coverage heatmap — rows are survey lines, columns are station positions. ' +
    'Colour intensity shows the number of recorded frequencies; ' +
    'dark cells flag missing or under-sampled stations for targeted re-acquisition.',

    'Per-line frequency distribution. ' +
    'Each box shows the median, quartiles, and outliers for one profile. ' +
    'Wide boxes indicate uneven station sampling; outliers flag individual problem stations.'
  ];

  var HOLD_MS   = 5200;   // ms each slide is visible (includes type-out time)
  var IN_MS     = 560;    // must match CSS crs-in animation duration
  var OUT_MS    = 440;    // must match CSS crs-out animation duration
  var TYPE_MS   = 22;     // ms per character typed

  /* ── State ──────────────────────────────────────────────────────────── */
  var current   = 0;
  var slideTimer = null;
  var typeTimer  = null;
  var running   = false;

  /* ── DOM helpers ────────────────────────────────────────────────────── */
  function el(id) { return document.getElementById(id); }

  /* ── Dot indicators ─────────────────────────────────────────────────── */
  function activateDot(idx) {
    DOT_IDS.forEach(function (id, i) {
      var d = el(id);
      if (!d) return;
      d.classList.toggle('crs-dot-active', i === idx);
    });
  }

  /* ── Caption counter ────────────────────────────────────────────────── */
  function fadeLabel(idx) {
    var lbl = el('crs-label');
    if (!lbl) return;
    lbl.classList.add('crs-fade');
    setTimeout(function () {
      lbl.textContent = LABELS[idx];
      lbl.classList.remove('crs-fade');
    }, 190);
  }

  /* ── Typewriter ─────────────────────────────────────────────────────── */
  function stopTyping() {
    if (typeTimer) { clearInterval(typeTimer); typeTimer = null; }
  }

  function startTyping(idx) {
    stopTyping();

    var tw  = el('crs-typewriter');
    var cur = el('crs-cursor');
    if (!tw) return;

    tw.textContent = '';
    if (cur) {
      cur.classList.remove('crs-cursor-dim');
    }

    var text = DESCRIPTIONS[idx] || '';
    var pos  = 0;

    typeTimer = setInterval(function () {
      if (pos < text.length) {
        tw.textContent += text[pos];
        pos++;
      } else {
        stopTyping();
        /* Dim the cursor once typing is complete */
        if (cur) cur.classList.add('crs-cursor-dim');
      }
    }, TYPE_MS);
  }

  /* ── Core slide transition ──────────────────────────────────────────── */
  function goTo(next) {
    var prev   = current;
    current    = next;

    var prevEl = el(SLIDE_IDS[prev]);
    var nextEl = el(SLIDE_IDS[next]);

    /* Exit current slide to the left */
    if (prevEl) {
      prevEl.classList.remove('crs-active');
      prevEl.classList.add('crs-exit');
      setTimeout(function () {
        if (prevEl) prevEl.classList.remove('crs-exit');
      }, OUT_MS);
    }

    /* Enter next slide from the right */
    if (nextEl) {
      nextEl.classList.add('crs-active');
      /*
       * After the animation finishes, fire a resize so Plotly repaints
       * its charts at the correct dimensions inside the newly-visible slide.
       */
      setTimeout(function () {
        window.dispatchEvent(new Event('resize'));
      }, IN_MS + 40);
    }

    activateDot(next);
    fadeLabel(next);

    /* Start typing the description a little after the slide arrives */
    setTimeout(function () { startTyping(next); }, IN_MS + 80);
  }

  /* ── Scheduler ──────────────────────────────────────────────────────── */
  function scheduleNext() {
    slideTimer = setInterval(function () {
      goTo((current + 1) % SLIDE_IDS.length);  /* auto-wraps across all 7 slides */
    }, HOLD_MS);
  }

  /* ── Start ──────────────────────────────────────────────────────────── */
  function start() {
    if (running) return;

    var firstEl = el(SLIDE_IDS[0]);
    if (!firstEl) return;

    current = 0;
    firstEl.classList.add('crs-active');
    activateDot(0);
    fadeLabel(0);

    /* Initial resize + typewriter for slide 0 */
    setTimeout(function () {
      window.dispatchEvent(new Event('resize'));
      startTyping(0);
    }, IN_MS + 80);

    running = true;
    scheduleNext();
  }

  /* ── Stop / Reset ───────────────────────────────────────────────────── */
  function stop() {
    if (slideTimer) { clearInterval(slideTimer); slideTimer = null; }
    stopTyping();
    running = false;
  }

  function reset() {
    stop();
    current = 0;
    SLIDE_IDS.forEach(function (id) {
      var s = el(id);
      if (s) s.classList.remove('crs-active', 'crs-exit');
    });
    var tw = el('crs-typewriter');
    if (tw) tw.textContent = '';
    var cur = el('crs-cursor');
    if (cur) cur.classList.remove('crs-cursor-dim');
    var lbl = el('crs-label');
    if (lbl) lbl.textContent = '';
    activateDot(0);
  }

  /* ── Dot click navigation ───────────────────────────────────────────── */
  function attachDotClicks() {
    DOT_IDS.forEach(function (id, i) {
      var d = el(id);
      if (!d || d._crsReady) return;
      d._crsReady = true;
      d.addEventListener('click', function () {
        if (i === current) return;
        /* Jump to the clicked slide and restart the auto-advance timer */
        if (slideTimer) { clearInterval(slideTimer); slideTimer = null; }
        goTo(i);
        scheduleNext();
      });
    });
  }

  /* ── Observe the analytics row visibility ───────────────────────────── */
  var rowObserver = new MutationObserver(function () {
    var row = el('dash-launch-row');
    if (!row) return;
    var visible = row.style.display && row.style.display !== 'none';
    if (visible && !running) {
      attachDotClicks();
      setTimeout(start, 150);
    } else if (!visible && running) {
      reset();
    }
  });

  /* ── Bootstrap ──────────────────────────────────────────────────────── */
  function init() {
    var row = el('dash-launch-row');
    if (!row) { setTimeout(init, 280); return; }

    /* Watch for Dash setting the style attribute on the row */
    rowObserver.observe(row, {
      attributes: true,
      attributeFilter: ['style']
    });

    /* Also watch the body in case Dash mounts the row after our init */
    var bodyObs = new MutationObserver(function () {
      var r = el('dash-launch-row');
      if (r && r.style.display && r.style.display !== 'none' && !running) {
        attachDotClicks();
        setTimeout(start, 150);
      }
    });
    bodyObs.observe(document.body, { childList: true, subtree: true });

    /* Immediate check — data may already be loaded */
    if (row.style.display && row.style.display !== 'none') {
      attachDotClicks();
      setTimeout(start, 150);
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    setTimeout(init, 220);
  }

}());
