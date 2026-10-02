/*
 * Keep Agent Master's working animations continuous and independent.
 *
 * The thinking row is re-rendered on every 600 ms poll, and a re-created
 * element restarts its CSS animation from zero. Each animated element is
 * therefore given a negative animation-delay equal to the page clock
 * modulo its own period when it appears, so a re-created element resumes
 * exactly where its predecessor was. The logo's signal, its earth layers
 * and the status-text shimmer each follow their own period, so none of
 * them resets or pulses in step with another.
 *
 * MutationObserver callbacks run before the next paint, so there is no
 * visible frame at phase zero. prefers-reduced-motion disables the
 * animations in CSS; the delay is then harmless.
 */
(function () {
    'use strict';

    // Must match the animation durations in master.css.
    var PERIODS_MS = {
        'am-logo-signal': 1600,
        'am-logo-bands': 7000,
        'am-think-lbl': 2600
    };
    var SELECTOR = Object.keys(PERIODS_MS)
        .map(function (cls) { return '.' + cls; })
        .join(',');

    function phase(el) {
        for (var cls in PERIODS_MS) {
            if (el.classList.contains(cls)) {
                var offset = performance.now() % PERIODS_MS[cls];
                el.style.animationDelay = '-' + offset.toFixed(0) + 'ms';
                return;
            }
        }
    }

    function scan(node) {
        if (!node || node.nodeType !== 1) { return; }
        if (node.matches(SELECTOR)) { phase(node); }
        node.querySelectorAll(SELECTOR).forEach(phase);
    }

    new MutationObserver(function (mutations) {
        mutations.forEach(function (m) {
            m.addedNodes.forEach(scan);
        });
    }).observe(document.documentElement, { childList: true, subtree: true });

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', function () {
            scan(document.body);
        });
    } else {
        scan(document.body);
    }
}());
