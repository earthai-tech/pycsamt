/**
 * popout.js — pop-out / lightbox for pyCSAMT figures.
 *
 * Strategy: MutationObserver watches for image based .fig-area elements and
 * Plotly based .tool-graph elements added by Dash. It injects a small pop-out
 * button (top-right corner, visible on hover). No Dash callbacks involved.
 */
(function () {
    'use strict';

    /* ── Lightbox API ────────────────────────────────────────────── */

    function openLightbox(src, filename) {
        var lb  = document.getElementById('lightbox-overlay');
        var img = document.getElementById('lightbox-img');
        var dl  = document.getElementById('lightbox-download');
        if (!lb || !img) return;

        img.src = src;
        if (dl) {
            dl.href = src;                       // wire Save button
            dl.download = filename || 'pycsamt_figure.png';
        }

        lb.style.display = 'flex';
        lb.classList.remove('lb-hidden');
        document.addEventListener('keydown', _escHandler);
    }

    function _escHandler(e) {
        if (e.key === 'Escape') closeLightbox();
    }

    /* exposed globally so the Redock button's onclick can call it */
    window.pycsamt_closeLightbox = function () {
        var lb = document.getElementById('lightbox-overlay');
        if (!lb) return;
        lb.classList.add('lb-hidden');
        /* wait for exit animation, then hide */
        setTimeout(function () {
            lb.style.display = 'none';
            lb.classList.remove('lb-hidden');
        }, 260);
        document.removeEventListener('keydown', _escHandler);
    };

    /* ── Wire the Redock / close button ─────────────────────────── */

    function wireCloseButton() {
        var btn = document.getElementById('lightbox-close');
        if (btn && !btn._wired) {
            btn._wired = true;
            btn.addEventListener('click', window.pycsamt_closeLightbox);
        }
    }

    function plotTitle(gd) {
        var title = gd && gd._fullLayout && gd._fullLayout.title;
        if (!title) return 'pycsamt_figure';
        if (typeof title === 'string') return title;
        return title.text || 'pycsamt_figure';
    }

    function safeFilename(text) {
        return String(text || 'pycsamt_figure')
            .toLowerCase()
            .replace(/<[^>]+>/g, '')
            .replace(/[^a-z0-9]+/g, '_')
            .replace(/^_+|_+$/g, '')
            .slice(0, 80) || 'pycsamt_figure';
    }

    function openPlotly(area) {
        var gd = area.querySelector('.js-plotly-plot');
        if (!gd || !window.Plotly || !window.Plotly.toImage) return false;
        var filename = safeFilename(plotTitle(gd)) + '.png';
        window.Plotly.toImage(gd, {
            format: 'png',
            scale: 2,
            width: Math.max(900, gd.clientWidth || 900),
            height: Math.max(560, gd.clientHeight || 560)
        }).then(function (src) {
            openLightbox(src, filename);
        }).catch(function (err) {
            console.warn('pyCSAMT plot pop-out failed', err);
        });
        return true;
    }

    /* ── Inject pop-out button into one figure container ─────────── */

    function injectPopout(area) {
        if (area.dataset.popoutDone) return;
        area.dataset.popoutDone = '1';

        var btn = document.createElement('button');
        btn.className = 'popout-btn';
        btn.title     = 'Pop out  ·  Esc to close';
        btn.innerHTML =
            '<img src="/icons/pop-out.svg" class="popout-icon" aria-hidden="true" />';

        btn.addEventListener('click', function (e) {
            e.stopPropagation();
            e.preventDefault();

            if (openPlotly(area)) return;

            /* find the first real <img> inside the fig-area */
            var img = area.querySelector('img:not(.popout-icon)');
            if (img && img.src && img.src.length > 200) {
                openLightbox(img.src);
            }
        });

        area.appendChild(btn);
    }

    /* ── Scan DOM for unprocessed figure containers ──────────────── */

    function scan() {
        document.querySelectorAll(
            '.fig-area:not([data-popout-done]), .tool-graph:not([data-popout-done])'
        ).forEach(injectPopout);
        wireCloseButton();
    }

    /* ── MutationObserver — catches Dash component re-renders ───── */

    var observer = new MutationObserver(function (mutations) {
        var needScan = false;
        mutations.forEach(function (m) {
            m.addedNodes.forEach(function (node) {
                if (node.nodeType !== 1) return;
                if (
                    node.classList &&
                    (node.classList.contains('fig-area') ||
                     node.classList.contains('tool-graph'))
                ) {
                    injectPopout(node);
                } else if (node.querySelectorAll) {
                    node.querySelectorAll(
                        '.fig-area:not([data-popout-done]), .tool-graph:not([data-popout-done])'
                    ).forEach(injectPopout);
                }
                needScan = true;
            });
        });
        if (needScan) wireCloseButton();
    });

    /* ── Bootstrap ───────────────────────────────────────────────── */

    function init() {
        scan();
        observer.observe(document.body, { childList: true, subtree: true });
        /* extra scans for late-rendered Dash pages */
        setTimeout(scan, 800);
        setTimeout(scan, 2500);
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }
})();
