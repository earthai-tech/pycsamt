/**
 * profile_panzoom.js
 *
 * Lightweight pan + zoom for static matplotlib PNG images inside
 * .profile-page-fig and .fig-img-wrap containers.
 *
 * Features
 * ─────────
 *  • Mouse-wheel zoom (centred on cursor position)
 *  • Click-drag to pan
 *  • Double-click to reset
 *  • Keyboard: + zoom in, - zoom out, 0 reset
 *  • Toolbar buttons: zoom-in, zoom-out, reset, fullscreen
 *  • Works on dynamically added images (MutationObserver)
 */
(function () {
  'use strict';

  var MIN_SCALE = 1.0;
  var MAX_SCALE = 8.0;
  var WHEEL_FACTOR = 0.0012;
  var BTN_STEP = 0.25;

  /* ── Per-image state ─────────────────────────────────────────── */
  var _state = new WeakMap();

  function getState(img) {
    if (!_state.has(img)) {
      _state.set(img, { scale: 1, tx: 0, ty: 0, dragging: false, sx: 0, sy: 0, otx: 0, oty: 0 });
    }
    return _state.get(img);
  }

  /* ── Apply CSS transform ──────────────────────────────────────── */
  function applyTransform(img) {
    var s = getState(img);
    img.style.transformOrigin = '0 0';
    img.style.transform = 'translate(' + s.tx + 'px,' + s.ty + 'px) scale(' + s.scale + ')';
    img.style.willChange = 'transform';
    img.style.cursor = s.scale > 1 ? 'grab' : 'default';
  }

  function resetTransform(img) {
    var s = getState(img);
    s.scale = 1; s.tx = 0; s.ty = 0;
    img.style.transform = '';
    img.style.cursor = 'default';
  }

  /* ── Clamp translation so image edge never goes past container ── */
  function clampTranslation(img, s) {
    var rect  = img.parentElement.getBoundingClientRect();
    var iw    = img.naturalWidth  * s.scale;
    var ih    = img.naturalHeight * s.scale;
    var maxTx = Math.max(0, (iw  - rect.width)  / 2);
    var maxTy = Math.max(0, (ih  - rect.height) / 2);
    s.tx = Math.max(-maxTx, Math.min(maxTx, s.tx));
    s.ty = Math.max(-maxTy, Math.min(maxTy, s.ty));
  }

  /* ── Wheel handler ────────────────────────────────────────────── */
  function onWheel(evt) {
    var img  = evt.currentTarget;
    var s    = getState(img);
    var rect = img.getBoundingClientRect();

    evt.preventDefault();

    var delta    = -evt.deltaY * WHEEL_FACTOR;
    var newScale = Math.min(MAX_SCALE, Math.max(MIN_SCALE, s.scale * (1 + delta)));
    if (newScale === s.scale) return;

    /* Zoom centred on cursor */
    var mx = evt.clientX - rect.left;
    var my = evt.clientY - rect.top;
    s.tx   = mx - (mx - s.tx) * (newScale / s.scale);
    s.ty   = my - (my - s.ty) * (newScale / s.scale);
    s.scale = newScale;
    clampTranslation(img, s);
    applyTransform(img);
  }

  /* ── Drag handlers ────────────────────────────────────────────── */
  function onMouseDown(evt) {
    if (evt.button !== 0) return;
    var img = evt.currentTarget;
    var s   = getState(img);
    if (s.scale <= 1) return;
    s.dragging = true;
    s.sx = evt.clientX; s.sy = evt.clientY;
    s.otx = s.tx;       s.oty = s.ty;
    img.style.cursor = 'grabbing';
    evt.preventDefault();
  }

  function onMouseMove(evt) {
    var img = evt.currentTarget;
    var s   = getState(img);
    if (!s.dragging) return;
    s.tx = s.otx + (evt.clientX - s.sx);
    s.ty = s.oty + (evt.clientY - s.sy);
    clampTranslation(img, s);
    applyTransform(img);
  }

  function endDrag(evt) {
    var img = evt.currentTarget;
    var s   = getState(img);
    s.dragging = false;
    img.style.cursor = s.scale > 1 ? 'grab' : 'default';
  }

  /* ── Double-click reset ───────────────────────────────────────── */
  function onDblClick(evt) {
    resetTransform(evt.currentTarget);
  }

  /* ── Toolbar ──────────────────────────────────────────────────── */
  function makeToolbar(img, wrap) {
    var bar = document.createElement('div');
    bar.className = 'pz-toolbar';

    function btn(icon, tip, fn) {
      var b = document.createElement('button');
      b.className = 'pz-btn';
      b.title     = tip;
      b.innerHTML = icon;
      b.addEventListener('click', function (e) { e.stopPropagation(); fn(); });
      bar.appendChild(b);
      return b;
    }

    btn('&#43;', 'Zoom in  (+)', function () {
      var s = getState(img);
      var rect = img.getBoundingClientRect();
      var cx = rect.width / 2, cy = rect.height / 2;
      var ns = Math.min(MAX_SCALE, s.scale + BTN_STEP);
      s.tx = cx - (cx - s.tx) * (ns / s.scale);
      s.ty = cy - (cy - s.ty) * (ns / s.scale);
      s.scale = ns;
      clampTranslation(img, s);
      applyTransform(img);
    });

    btn('&#8722;', 'Zoom out  (−)', function () {
      var s = getState(img);
      var rect = img.getBoundingClientRect();
      var cx = rect.width / 2, cy = rect.height / 2;
      var ns = Math.max(MIN_SCALE, s.scale - BTN_STEP);
      s.tx = cx - (cx - s.tx) * (ns / s.scale);
      s.ty = cy - (cy - s.ty) * (ns / s.scale);
      s.scale = ns;
      if (s.scale <= 1) { s.tx = 0; s.ty = 0; }
      applyTransform(img);
    });

    btn('&#8634;', 'Reset view  (0)', function () {
      resetTransform(img);
    });

    /* Fullscreen */
    btn('&#x26F6;', 'Open fullscreen', function () {
      var dlg = document.getElementById('pz-lightbox');
      if (!dlg) {
        dlg = document.createElement('div');
        dlg.id = 'pz-lightbox';
        dlg.innerHTML = '<div class="pz-lb-inner"><img id="pz-lb-img"/>' +
          '<button id="pz-lb-close" class="pz-lb-close" title="Close (Esc)">&times;</button></div>';
        document.body.appendChild(dlg);
        document.getElementById('pz-lb-close').addEventListener('click', function () {
          dlg.classList.remove('pz-lb-open');
        });
        dlg.addEventListener('click', function (e) {
          if (e.target === dlg) dlg.classList.remove('pz-lb-open');
        });
        document.addEventListener('keydown', function (e) {
          if (e.key === 'Escape') dlg.classList.remove('pz-lb-open');
        });
      }
      document.getElementById('pz-lb-img').src = img.src;
      dlg.classList.add('pz-lb-open');
    });

    /* Export PNG — client-side instant download from img.src */
    (function () {
      var exportBtn = document.createElement('button');
      exportBtn.className = 'pz-btn pz-export-btn';
      exportBtn.title     = 'Download PNG  (Shift+S)';
      exportBtn.innerHTML = '&#11123;';   /* ⬛ down-arrow-to-bar */
      exportBtn.addEventListener('click', function (e) {
        e.stopPropagation();
        if (!img.src || img.src.length < 50) return;
        var id  = img.id || 'figure';
        var a   = document.createElement('a');
        a.href  = img.src;
        a.download = 'pycsamt_' + id.replace(/[^a-z0-9_-]/gi, '_') + '.png';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
      });
      bar.appendChild(exportBtn);
    }());

    wrap.appendChild(bar);
  }

  /* ── Attach to an image ───────────────────────────────────────── */
  function attach(img) {
    if (img._pzReady) return;
    img._pzReady = true;

    var wrap = img.parentElement;
    if (!wrap) return;

    /* Ensure overflow hidden on wrap */
    if (getComputedStyle(wrap).overflow === 'visible') {
      wrap.style.overflow = 'hidden';
    }
    wrap.style.position = 'relative';

    img.addEventListener('wheel',      onWheel,     { passive: false });
    img.addEventListener('mousedown',  onMouseDown);
    img.addEventListener('mousemove',  onMouseMove);
    img.addEventListener('mouseup',    endDrag);
    img.addEventListener('mouseleave', endDrag);
    img.addEventListener('dblclick',   onDblClick);

    makeToolbar(img, wrap);
  }

  /* ── Keyboard shortcuts (global) ─────────────────────────────── */
  document.addEventListener('keydown', function (evt) {
    if (['INPUT', 'TEXTAREA', 'SELECT'].indexOf(document.activeElement.tagName) >= 0) return;
    var key = evt.key;
    /* Find any visible pz-img */
    var imgs = Array.from(document.querySelectorAll('.pz-img'));
    var visible = imgs.filter(function (im) {
      return im.offsetParent !== null && im.src && !im.src.endsWith('=');
    });
    if (!visible.length) return;
    var img = visible[0];
    var s   = getState(img);

    if (key === '+' || key === '=') {
      s.scale = Math.min(MAX_SCALE, s.scale + BTN_STEP);
      applyTransform(img);
    } else if (key === '-') {
      s.scale = Math.max(MIN_SCALE, s.scale - BTN_STEP);
      if (s.scale <= 1) { s.tx = 0; s.ty = 0; }
      applyTransform(img);
    } else if (key === '0') {
      resetTransform(img);
    } else if ((key === 'S' || key === 's') && evt.shiftKey) {
      /* Shift+S — download the visible figure as PNG */
      if (img.src && img.src.length > 50) {
        var a  = document.createElement('a');
        a.href = img.src;
        a.download = 'pycsamt_' + (img.id || 'figure').replace(/[^a-z0-9_-]/gi, '_') + '.png';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
      }
    }
  });

  /* ── Scan and attach ──────────────────────────────────────────── */
  function scanImages() {
    var sel = '.profile-page-fig img, .fig-img-wrap img';
    document.querySelectorAll(sel).forEach(function (img) {
      img.classList.add('pz-img');
      if (img.complete) {
        attach(img);
      } else {
        img.addEventListener('load', function () { attach(img); }, { once: true });
      }
    });
  }

  /* ── Watch for Dash-injected images ─────────────────────────────*/
  var obs = new MutationObserver(function () { scanImages(); });

  function init() {
    scanImages();
    obs.observe(document.body, { childList: true, subtree: true, attributes: true, attributeFilter: ['src'] });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    setTimeout(init, 300);
  }

}());
