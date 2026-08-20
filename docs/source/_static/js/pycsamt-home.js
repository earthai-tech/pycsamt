/* pyCSAMT v2 — landing page enhancements.
 * Progressive only: the page is fully usable with JS disabled.
 *  1. Tags the <body> so home-only CSS applies (fallback for :has()).
 *  2. Scroll-reveal via IntersectionObserver.
 *  3. Count-up animation on the stats bar.
 *  4. Rotating survey-method keyword in the hero subtitle.
 *  5. Hero background carousel (slides + dots, auto-advance).
 *  6. Typewriter animation on the code-in-action tabs.
 *  7. Click-to-copy on the "pip install pycsamt" workflow-strip pill.
 *  8. Tap/keyboard toggle for the six flip cards (CSS handles hover
 *     and :focus-within on its own; this is the explicit control).
 *  9. Measures flip-card face heights so the card fits its content
 *     instead of the CSS fallback height.
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
    /* Ground methods (CSAMT/AMT/MT/TDEM/CSEM) keep the default accent
     * colour; airborne EM methods (AFMAG/ZTEM/MobileMT), listed
     * separately via data-airborne-words, get the .is-airborne colour
     * instead, so the rotator visibly flags the category switch. */
    var rotator = home.querySelector(".pyc-rotator");
    if (rotator) {
      var words;
      try {
        words = JSON.parse(rotator.getAttribute("data-words") || "[]");
      } catch (e) {
        words = [];
      }
      var airborneWords;
      try {
        airborneWords = JSON.parse(
          rotator.getAttribute("data-airborne-words") || "[]"
        );
      } catch (e) {
        airborneWords = [];
      }
      var airborneSet = new Set(airborneWords);
      var syncAirborne = function (word) {
        rotator.classList.toggle("is-airborne", airborneSet.has(word));
      };
      syncAirborne(rotator.textContent.trim());

      if (!reduceMotion && words.length > 1) {
        var idx = 0;
        window.setInterval(function () {
          rotator.classList.add("is-swapping");
          window.setTimeout(function () {
            idx = (idx + 1) % words.length;
            rotator.textContent = words[idx];
            syncAirborne(words[idx]);
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

    /* ---- 8. flip-card toggle (touch & keyboard) -------------------------- */
    /* Hover and :focus-within already flip the six capability cards in
     * pure CSS, so this only has to cover the case CSS can't: a tap on a
     * touchscreen. Clicking a card's round button toggles .is-flipped;
     * clicking elsewhere, or Escape, closes whichever card is open. */
    var flipButtons = Array.prototype.slice.call(
      home.querySelectorAll(".pyc-feature-flip")
    );
    if (flipButtons.length) {
      var closeAllFlips = function (except) {
        flipButtons.forEach(function (btn) {
          if (btn === except) return;
          var card = btn.closest(".pyc-feature");
          if (card) card.classList.remove("is-flipped");
          btn.setAttribute("aria-expanded", "false");
        });
      };

      flipButtons.forEach(function (btn) {
        btn.addEventListener("click", function (event) {
          event.stopPropagation();
          var card = btn.closest(".pyc-feature");
          if (!card) return;
          var flipped = card.classList.toggle("is-flipped");
          btn.setAttribute("aria-expanded", flipped ? "true" : "false");
          if (flipped) closeAllFlips(btn);
        });
      });

      document.addEventListener("click", function () {
        closeAllFlips(null);
      });
      document.addEventListener("keydown", function (event) {
        if (event.key === "Escape") closeAllFlips(null);
      });
    }

    /* ---- 9. flip-card height measurement --------------------------------- */
    /* Both faces are position:absolute so the card can't size itself from
     * content; drive --pyc-flip-h from the taller of the two faces
     * instead of relying on the fixed CSS fallback. Re-measured on
     * resize and once each hero image has its real dimensions.
     *
     * scrollHeight can't be read directly off a face: inset:0 pins both
     * its top and bottom to .pyc-feature-inner, so as long as content
     * fits, scrollHeight just reports back whatever height the card
     * currently has (the very value we're trying to compute) instead of
     * how tall the content actually wants to be. Flipping the face to
     * position:static for one synchronous read lets it size to content,
     * then the inline override is removed before the browser paints. */
    var featureCards = Array.prototype.slice.call(
      home.querySelectorAll(".pyc-feature")
    );
    if (featureCards.length) {
      var measureFaceHeight = function (face) {
        var prevPosition = face.style.position;
        face.style.position = "static";
        var h = face.offsetHeight;
        face.style.position = prevPosition;
        return h;
      };

      var sizeFeatureCard = function (card) {
        var front = card.querySelector(".pyc-feature-front");
        var back = card.querySelector(".pyc-feature-back");
        if (!front || !back) return;
        var h = Math.max(measureFaceHeight(front), measureFaceHeight(back));
        if (h) card.style.setProperty("--pyc-flip-h", h + "px");
      };

      var sizeAllFeatureCards = function () {
        featureCards.forEach(sizeFeatureCard);
      };

      sizeAllFeatureCards();

      home.querySelectorAll(".pyc-feature-hero img").forEach(function (img) {
        if (!img.complete) {
          img.addEventListener("load", sizeAllFeatureCards, { once: true });
        }
      });

      var featureResizeTimer = null;
      window.addEventListener("resize", function () {
        window.clearTimeout(featureResizeTimer);
        featureResizeTimer = window.setTimeout(sizeAllFeatureCards, 150);
      });
    }
  });
})();
