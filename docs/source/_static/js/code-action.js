/* pyCSAMT — "Code in action" live-session player (home page).
 *
 * Each tab of the .pyc-code-panel tab-set holds several captioned
 * code-blocks ordered simple -> complex. This script collapses them into
 * ONE fixed stage per tab and plays them like a live coding session:
 * example 1 types itself out (syntax colours intact), rests, is replaced
 * IN PLACE by example 2, and so on in a forward loop. The header shows
 * the running example's subtitle (bright + pulsing dot while typing,
 * dimmed once done); numbered dots allow manual jumps, which show the
 * example instantly for reading.
 *
 * Politeness: hovering finishes the current example immediately (so the
 * code is readable and copyable) and pauses the sequence; hidden browser
 * tabs and off-screen panels don't play; prefers-reduced-motion shows
 * static code with manual dots only. Blocks are hidden with inline
 * styles, so a stale stylesheet can never leave examples stacked.
 */
(function () {
  "use strict";

  var TICK_MS = 16;        /* typing timer interval */
  var CHARS_PER_TICK = 3;  /* characters revealed per tick */
  var REST_MS = 3800;      /* pause on the finished example */
  var FIRST_DELAY = 700;   /* settle time for the tab visible at load */

  var reduceMotion =
    window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  /* All text nodes of a <pre>, in document order, with their full text. */
  function textUnits(pre) {
    var units = [];
    var walker = document.createTreeWalker(pre, NodeFilter.SHOW_TEXT);
    var node;
    while ((node = walker.nextNode())) {
      units.push({ node: node, text: node.textContent });
    }
    return units;
  }

  function buildTab(content, tabIndex) {
    var blocks = Array.prototype.slice.call(
      content.querySelectorAll("div.literal-block-wrapper")
    );
    if (!blocks.length) return;

    /* ---- header: status dot + running subtitle + numbered dots ---- */
    var head = document.createElement("div");
    head.className = "pyc-ca-head";
    var status = document.createElement("span");
    status.className = "pyc-ca-status";
    status.textContent = "▶"; /* small play triangle */
    var title = document.createElement("span");
    title.className = "pyc-ca-title";
    var dotsBox = document.createElement("span");
    dotsBox.className = "pyc-ca-dots";
    head.appendChild(status);
    head.appendChild(title);
    head.appendChild(dotsBox);

    /* ---- one fixed stage holding every example ---- */
    var stage = document.createElement("div");
    stage.className = "pyc-ca-stage";

    var labels = [];
    var dots = [];
    blocks.forEach(function (block, i) {
      var caption = block.querySelector(".code-block-caption");
      var label = "Example " + (i + 1);
      if (caption) {
        var clone = caption.cloneNode(true);
        Array.prototype.forEach.call(
          clone.querySelectorAll("a.headerlink"),
          function (a) { a.parentNode.removeChild(a); }
        );
        label = clone.textContent.trim();
        caption.style.display = "none";
      }
      labels.push(label);

      var pre = block.querySelector("pre");
      block.__units = pre ? textUnits(pre) : [];
      block.classList.add("pyc-ca-block");
      block.style.display = "none";
      stage.appendChild(block);

      var dot = document.createElement("button");
      dot.type = "button";
      dot.className = "pyc-ca-dot";
      dot.textContent = String(i + 1);
      dot.setAttribute("aria-label", "Show example: " + label);
      dot.addEventListener("click", function () {
        play(i, false); /* manual jump: show instantly, keep cycling */
      });
      dotsBox.appendChild(dot);
      dots.push(dot);
    });

    content.insertBefore(stage, content.firstChild);
    content.insertBefore(head, stage);

    /* Fixed height = tallest example, so the stage never jumps. Inactive
     * sphinx-design tabs are display:none at load (heights measure 0),
     * so this also runs lazily on the first play() inside each tab. */
    function measureStage() {
      if (stage.style.minHeight) return;
      /* Never measure while the surrounding tab is hidden — everything
       * reads 0 and would poison the locks. */
      if (!content.offsetHeight) return;
      var maxH = 0;
      blocks.forEach(function (block) {
        var prev = block.style.display;
        block.style.visibility = "hidden";
        block.style.display = "block";
        /* Lock each code box to its full height so it never collapses
         * and regrows while its code is being typed. */
        var hl = block.querySelector(".highlight");
        if (hl && !(parseFloat(hl.style.minHeight) > 0)) {
          var h = hl.offsetHeight;
          if (h) hl.style.minHeight = h + "px";
        }
        maxH = Math.max(maxH, block.offsetHeight);
        block.style.display = prev;
        block.style.visibility = "";
      });
      if (maxH) stage.style.minHeight = maxH + "px";
    }
    measureStage();

    var caret = document.createElement("span");
    caret.className = "pyc-ca-caret";

    var current = 0;
    var typingTimer = null;
    var restTimer = null;
    var hovered = false;
    var onScreen = false;
    var initTime = window.performance ? window.performance.now() : 0;

    function stopTimers() {
      if (typingTimer) { window.clearInterval(typingTimer); typingTimer = null; }
      if (restTimer) { window.clearTimeout(restTimer); restTimer = null; }
    }

    function setHeader(running) {
      title.textContent = labels[current];
      head.classList.toggle("pyc-ca-running", running);
      dots.forEach(function (dot, k) {
        dot.classList.toggle("pyc-ca-dot-on", k === current);
      });
    }

    function fillBlock(block) {
      block.__units.forEach(function (u) { u.node.textContent = u.text; });
      if (caret.parentNode) caret.parentNode.removeChild(caret);
    }

    function showOnly(i) {
      current = ((i % blocks.length) + blocks.length) % blocks.length;
      blocks.forEach(function (block, k) {
        block.style.display = k === current ? "block" : "none";
      });
    }

    function scheduleNext() {
      if (restTimer) window.clearTimeout(restTimer);
      restTimer = window.setTimeout(function () {
        if (hovered || !onScreen || document.hidden || reduceMotion) return;
        play(current + 1, true);
      }, REST_MS);
    }

    /* Play example i. animate=true types it out; false shows it whole. */
    function play(i, animate) {
      stopTimers();
      measureStage();
      showOnly(i);
      var block = blocks[current];

      if (!animate || reduceMotion || document.hidden || !onScreen) {
        fillBlock(block);
        setHeader(false);
        if (!reduceMotion) scheduleNext();
        return;
      }

      /* type it out */
      block.__units.forEach(function (u) { u.node.textContent = ""; });
      setHeader(true);
      var ui = 0;
      var offset = 0;
      typingTimer = window.setInterval(function () {
        var budget = CHARS_PER_TICK;
        while (budget > 0 && ui < block.__units.length) {
          var unit = block.__units[ui];
          var take = Math.min(budget, unit.text.length - offset);
          offset += take;
          budget -= take;
          unit.node.textContent = unit.text.slice(0, offset);
          if (caret.parentNode !== unit.node.parentNode ||
              caret.previousSibling !== unit.node) {
            unit.node.parentNode.insertBefore(caret, unit.node.nextSibling);
          }
          if (offset >= unit.text.length) { ui += 1; offset = 0; }
        }
        if (ui >= block.__units.length) {
          window.clearInterval(typingTimer);
          typingTimer = null;
          setHeader(false);
          window.setTimeout(function () {
            if (caret.parentNode) caret.parentNode.removeChild(caret);
          }, 900);
          scheduleNext();
        }
      }, TICK_MS);
    }

    /* Hover: finish the code instantly so it can be read and copied. */
    function settle() {
      if (typingTimer) {
        window.clearInterval(typingTimer);
        typingTimer = null;
        fillBlock(blocks[current]);
        setHeader(false);
      }
      if (restTimer) { window.clearTimeout(restTimer); restTimer = null; }
    }

    content.addEventListener("pointerenter", function () {
      hovered = true;
      settle();
    });
    content.addEventListener("pointerleave", function () {
      hovered = false;
      scheduleNext();
    });
    content.addEventListener("focusin", function () {
      hovered = true;
      settle();
    });
    content.addEventListener("focusout", function () {
      hovered = false;
      scheduleNext();
    });
    document.addEventListener("visibilitychange", function () {
      if (document.hidden) settle();
      else scheduleNext();
    });

    if ("IntersectionObserver" in window) {
      new IntersectionObserver(function (entries) {
        var was = onScreen;
        onScreen = entries[0].isIntersecting;
        if (onScreen && !was) {
          /* Lock heights BEFORE anything plays, so the first run cannot
           * shift the layout. */
          measureStage();
          /* The tab visible at page load settles briefly; a tab the user
           * just opened starts typing right away, so the code is never
           * seen fully rendered and then wiped. */
          var sinceInit =
            (window.performance ? window.performance.now() : 9e9) - initTime;
          var wait = sinceInit < 1500 ? FIRST_DELAY : 200;
          window.setTimeout(function () {
            if (onScreen && !hovered && !document.hidden) play(current, true);
          }, wait);
        } else if (!onScreen) {
          settle();
        }
      }, { threshold: 0.2 }).observe(content);
    } else {
      onScreen = true;
      window.setTimeout(function () { play(0, true); }, FIRST_DELAY);
    }

    /* initial static state (also the no-play fallback) */
    showOnly(0);
    fillBlock(blocks[0]);
    setHeader(false);

    /* programmatic hook (manual testing) */
    content._pycCa = {
      play: play,
      settle: settle,
      state: function () {
        return {
          current: current,
          typing: !!typingTimer,
          visible: blocks.map(function (b) { return b.style.display; }),
          title: title.textContent,
          running: head.classList.contains("pyc-ca-running"),
        };
      },
    };
  }

  function init() {
    var contents = document.querySelectorAll(
      ".pyc-code-panel .sd-tab-content"
    );
    Array.prototype.forEach.call(contents, buildTab);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
