/* pyCSAMT v2 — landing page enhancements.
 * Progressive only: the page is fully usable with JS disabled.
 *  1. Tags the <body> so home-only CSS applies (fallback for :has()).
 *  2. Scroll-reveal via IntersectionObserver.
 *  3. Count-up animation on the stats bar.
 *  4. Rotating survey-method keyword in the hero subtitle.
 *  5. Hero background carousel (slides + dots, auto-advance).
 *  6. Typewriter animation on the code-in-action tabs.
 *  7. Click-to-copy on the "pip install pycsamt" workflow-strip pill.
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

    /* ---- 2. scroll reveal --------------------------------------------- */
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

    /* ---- 3. stat counters ---------------------------------------------- */
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

    /* ---- 4. rotating hero keyword --------------------------------------- */
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

    /* ---- 5. hero background carousel ------------------------------------ */
    /* Crossfades the .pyc-slide layers and keeps the dots in sync. Dots
     * always work; auto-advance is skipped under prefers-reduced-motion
     * and pauses while the hero is hovered or the tab is hidden. */
    var slides = Array.prototype.slice.call(
      home.querySelectorAll(".pyc-hero .pyc-slide")
    );
    var dots = Array.prototype.slice.call(
      home.querySelectorAll(".pyc-hero-dots button")
    );
    if (slides.length > 1) {
      var current = 0;
      var slideTimer = null;

      var showSlide = function (i) {
        current = ((i % slides.length) + slides.length) % slides.length;
        slides.forEach(function (s, k) {
          s.classList.toggle("is-active", k === current);
        });
        dots.forEach(function (d, k) {
          d.classList.toggle("is-active", k === current);
          d.setAttribute("aria-selected", k === current ? "true" : "false");
        });
      };

      var stopSlides = function () {
        if (slideTimer) {
          window.clearInterval(slideTimer);
          slideTimer = null;
        }
      };

      var playSlides = function () {
        if (reduceMotion || slideTimer || document.hidden) return;
        slideTimer = window.setInterval(function () {
          showSlide(current + 1);
        }, 7000);
      };

      dots.forEach(function (d, k) {
        d.addEventListener("click", function () {
          stopSlides();
          showSlide(k);
          playSlides();
        });
      });

      var hero = home.querySelector(".pyc-hero");
      if (hero) {
        hero.addEventListener("mouseenter", stopSlides);
        hero.addEventListener("mouseleave", playSlides);
      }
      document.addEventListener("visibilitychange", function () {
        if (document.hidden) stopSlides();
        else playSlides();
      });

      playSlides();
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

    /* ---- 7. install-pill click-to-copy ---------------------------------- */
    var installBtn = home.querySelector(".pyc-flow-install");
    if (installBtn) {
      installBtn.addEventListener("click", function () {
        var text = installBtn.getAttribute("data-copy") || "";
        var flash = function () {
          installBtn.classList.add("is-copied");
          window.clearTimeout(installBtn._pycCopyTimer);
          installBtn._pycCopyTimer = window.setTimeout(function () {
            installBtn.classList.remove("is-copied");
          }, 1600);
        };
        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(text).then(flash, flash);
        } else {
          // Fallback for non-secure contexts / older browsers.
          var tmp = document.createElement("textarea");
          tmp.value = text;
          tmp.style.position = "fixed";
          tmp.style.opacity = "0";
          document.body.appendChild(tmp);
          tmp.focus();
          tmp.select();
          try {
            document.execCommand("copy");
          } catch (err) {
            /* clipboard unavailable; the command is still selected/visible */
          }
          document.body.removeChild(tmp);
          flash();
        }
      });
    }
  });
})();
