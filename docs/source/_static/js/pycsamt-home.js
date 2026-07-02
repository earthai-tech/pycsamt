/* pyCSAMT v2 — landing page enhancements.
 * Progressive only: the page is fully usable with JS disabled.
 *  1. Tags the <body> so home-only CSS applies (fallback for :has()).
 *  2. Copy-to-clipboard on the pip-install pill.
 *  3. Scroll-reveal via IntersectionObserver.
 *  4. Count-up animation on the hero stats strip.
 *  5. Rotating survey-method keyword in the hero subtitle.
 *  6. Typewriter animation on the code-in-action tabs.
 * All motion respects prefers-reduced-motion.
 */

(function () {
  "use strict";

  function ready(fn) {
    if (document.readyState !== "loading") fn();
    else document.addEventListener("DOMContentLoaded", fn);
  }

  ready(function () {
    var home = document.querySelector(".pycsamt-home");
    if (!home) return;

    document.body.classList.add("pycsamt-home-page");

    var reduceMotion = window.matchMedia(
      "(prefers-reduced-motion: reduce)"
    ).matches;

    /* ---- 2. copy install command ------------------------------------- */
    document.querySelectorAll(".pyc-copy-btn").forEach(function (btn) {
      btn.addEventListener("click", function () {
        var text = btn.getAttribute("data-copy") || "pip install pycsamt";
        var done = function () {
          if (btn.classList.contains("is-copied")) return;
          btn.classList.add("is-copied");
          // The theme loads Font Awesome in SVG+JS mode, which replaces
          // <i> with <svg>; swap the whole markup and let FA's observer
          // re-convert the restored <i> tag.
          var original = btn.innerHTML;
          btn.innerHTML = '<i class="fa-solid fa-check"></i>';
          window.setTimeout(function () {
            btn.classList.remove("is-copied");
            btn.innerHTML = original;
          }, 1800);
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(text).then(done, done);
        } else {
          var ta = document.createElement("textarea");
          ta.value = text;
          ta.style.position = "fixed";
          ta.style.opacity = "0";
          document.body.appendChild(ta);
          ta.select();
          try {
            document.execCommand("copy");
          } catch (e) {
            /* clipboard unavailable — button feedback still shown */
          }
          document.body.removeChild(ta);
          done();
        }
      });
    });

    /* ---- 3. scroll reveal --------------------------------------------- */
    var revealed = home.querySelectorAll(".pyc-reveal");
    if (reduceMotion || !("IntersectionObserver" in window)) {
      revealed.forEach(function (el) {
        el.classList.add("is-visible");
      });
    } else {
      var io = new IntersectionObserver(
        function (entries) {
          entries.forEach(function (entry) {
            if (entry.isIntersecting) {
              entry.target.classList.add("is-visible");
              io.unobserve(entry.target);
            }
          });
        },
        { threshold: 0.12, rootMargin: "0px 0px -5% 0px" }
      );
      revealed.forEach(function (el) {
        io.observe(el);
      });
    }

    /* ---- 4. stat counters ---------------------------------------------- */
    function animateCounter(el) {
      var target = parseInt(el.getAttribute("data-counter"), 10);
      if (isNaN(target)) return;
      var suffix = el.getAttribute("data-suffix") || "";
      if (reduceMotion) {
        el.textContent = target + suffix;
        return;
      }
      var duration = 1100;
      var start = null;
      function step(ts) {
        if (start === null) start = ts;
        var p = Math.min((ts - start) / duration, 1);
        // ease-out cubic
        var eased = 1 - Math.pow(1 - p, 3);
        el.textContent = Math.round(eased * target) + suffix;
        if (p < 1) window.requestAnimationFrame(step);
      }
      window.requestAnimationFrame(step);
    }

    var counters = home.querySelectorAll("[data-counter]");
    if (!("IntersectionObserver" in window)) {
      counters.forEach(animateCounter);
    } else {
      var cio = new IntersectionObserver(
        function (entries) {
          entries.forEach(function (entry) {
            if (entry.isIntersecting) {
              animateCounter(entry.target);
              cio.unobserve(entry.target);
            }
          });
        },
        { threshold: 0.4 }
      );
      counters.forEach(function (el) {
        cio.observe(el);
      });
    }

    /* ---- 5. rotating hero keyword --------------------------------------- */
    var rotator = home.querySelector(".pyc-rotator");
    if (rotator && !reduceMotion) {
      var words;
      try {
        words = JSON.parse(rotator.getAttribute("data-words") || "[]");
      } catch (e) {
        words = [];
      }
      if (words.length > 1) {
        var idx = 0;
        window.setInterval(function () {
          rotator.classList.add("is-swapping");
          window.setTimeout(function () {
            idx = (idx + 1) % words.length;
            rotator.textContent = words[idx];
            rotator.classList.remove("is-swapping");
          }, 270);
        }, 2600);
      }
    }

    /* ---- 6. typewriter animation on the code-in-action tabs ------------- */
    /* Refills the pygments text nodes character by character, so syntax
     * colouring is preserved while the code "types itself". Runs when a
     * tab is toggled and once when the panel first scrolls into view.
     * Any click inside the panel fast-forwards to the full snippet, which
     * also guarantees sphinx-copybutton always copies complete code. */
    var codePanel = home.querySelector(".pyc-code-panel");
    if (codePanel && !reduceMotion) {
      var finishCurrent = null;

      var typePane = function (pane) {
        if (finishCurrent) finishCurrent();
        var pre = pane.querySelector("div.highlight pre");
        if (!pre) return;

        var walker = document.createTreeWalker(pre, NodeFilter.SHOW_TEXT, null);
        var nodes = [];
        while (walker.nextNode()) nodes.push(walker.currentNode);
        var originals = [];
        var total = 0;
        nodes.forEach(function (n) {
          originals.push(n.nodeValue);
          total += n.nodeValue.length;
        });
        if (!total || total > 1500) return;

        // Reserve the final height so the panel does not grow line by line.
        pre.style.minHeight = pre.offsetHeight + "px";
        nodes.forEach(function (n) {
          n.nodeValue = "";
        });

        var cursor = document.createElement("span");
        cursor.className = "pyc-type-cursor";
        pre.appendChild(cursor);

        var ni = 0;
        var ci = 0;
        var timer = null;
        var done = false;

        var finish = function () {
          if (done) return;
          done = true;
          if (timer) window.clearTimeout(timer);
          nodes.forEach(function (n, k) {
            n.nodeValue = originals[k];
          });
          if (cursor.parentNode) cursor.parentNode.removeChild(cursor);
          pre.style.minHeight = "";
          if (finishCurrent === finish) finishCurrent = null;
        };
        finishCurrent = finish;

        var step = function () {
          if (done) return;
          while (ni < nodes.length && !originals[ni].length) ni++;
          if (ni >= nodes.length) {
            finish();
            return;
          }
          ci++;
          var ch = originals[ni].charAt(ci - 1);
          nodes[ni].nodeValue = originals[ni].slice(0, ci);
          var host = nodes[ni].parentNode;
          if (
            host &&
            (cursor.parentNode !== host ||
              cursor.previousSibling !== nodes[ni])
          ) {
            host.insertBefore(cursor, nodes[ni].nextSibling);
          }
          if (ci >= originals[ni].length) {
            ni++;
            ci = 0;
          }
          var delay = 4 + Math.random() * 10;
          if (ch === "\n") delay = 75;
          else if (ch === " ") delay = 2;
          timer = window.setTimeout(step, delay);
        };
        step();
      };

      var tabInputs = codePanel.querySelectorAll(".sd-tab-set > input");
      var tabPanes = codePanel.querySelectorAll(".sd-tab-content");
      tabInputs.forEach(function (inp, k) {
        inp.addEventListener("change", function () {
          if (inp.checked && tabPanes[k]) typePane(tabPanes[k]);
        });
      });

      // First time the panel scrolls into view, type the active tab.
      var typeActivePane = function () {
        tabInputs.forEach(function (inp, k) {
          if (inp.checked && tabPanes[k]) typePane(tabPanes[k]);
        });
      };
      if ("IntersectionObserver" in window) {
        var pio = new IntersectionObserver(
          function (entries) {
            entries.forEach(function (entry) {
              if (entry.isIntersecting) {
                typeActivePane();
                pio.disconnect();
              }
            });
          },
          { threshold: 0.25 }
        );
        pio.observe(codePanel);
      } else {
        typeActivePane();
      }

      // Capture-phase: fast-forward before copy buttons or labels react.
      codePanel.addEventListener(
        "click",
        function () {
          if (finishCurrent) finishCurrent();
        },
        true
      );
    }
  });
})();
