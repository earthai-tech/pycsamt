/* pyCSAMT — examples hub cards (see docs/source/_ext/gallery_hub.py).
 *
 * Renders one card per gallery section from the inline manifest that
 * gallery_hub.py embeds in examples/index.html
 * (window.__PYCSAMT_GALLERY_HUB__). Each card slides through its
 * section's plot thumbnails right-to-left while the caption types out
 * the current example's title.
 *
 * Politeness: images lazy-load, only on-screen cards animate
 * (IntersectionObserver), hover/focus and hidden tabs pause the timer,
 * and prefers-reduced-motion disables all animation (static first
 * thumbnail, full caption, dots still clickable).
 */
(function () {
  "use strict";

  var ADVANCE_MS = 4200;   /* time each slide rests */
  var STAGGER_MS = 1200;   /* per-card start offset so cards don't sync */
  var TYPE_MS = 24;        /* typewriter speed per character */
  var MAX_DOTS = 10;       /* above this, show a "3 / 22" counter */

  var reduceMotion =
    window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function el(tag, className, parent) {
    var node = document.createElement(tag);
    if (className) node.className = className;
    if (parent) parent.appendChild(node);
    return node;
  }

  function buildCard(section, index) {
    var card = el("article", "pgh-card");

    var head = el("a", "pgh-head", card);
    head.href = section.url;
    var titleBox = el("span", "pgh-titles", head);
    el("span", "pgh-title", titleBox).textContent = section.title;
    if (section.blurb) {
      el("span", "pgh-blurb", titleBox).textContent = section.blurb;
    }
    el("span", "pgh-count", head).textContent = section.count + " examples";

    var view = el("a", "pgh-view", card);
    view.href = section.url;
    view.setAttribute("aria-label", "Open the " + section.title + " gallery");
    var track = el("div", "pgh-track", view);
    section.slides.forEach(function (slide, i) {
      var frame = el("div", "pgh-slide", track);
      var img = el("img", null, frame);
      img.src = slide.thumb;
      img.alt = slide.title;
      img.loading = i === 0 ? "eager" : "lazy";
      img.decoding = "async";
      img.draggable = false;
    });
    /* Clone of the first slide appended after the last: the wrap-around
     * (last -> first) slides forward into the clone and then snaps back,
     * unanimated, to the real first slide — so the loop always moves
     * right-to-left instead of rewinding through every slide. */
    if (section.slides.length > 1) {
      var loopClone = track.children[0].cloneNode(true);
      loopClone.setAttribute("aria-hidden", "true");
      track.appendChild(loopClone);
    }

    var cap = el("div", "pgh-cap", card);
    var capLink = el("a", "pgh-cap-link", cap);
    var capText = el("span", "pgh-cap-text", capLink);
    var caret = el("span", "pgh-caret", cap);
    caret.setAttribute("aria-hidden", "true");

    var foot = el("div", "pgh-foot", card);
    var nav = el("span", "pgh-nav", foot);
    var open = el("a", "pgh-open", foot);
    open.href = section.url;
    open.textContent = "Browse section →";

    var slides = section.slides;
    var current = 0;         /* logical slide index (0..slides.length-1) */
    var pendingSnap = false; /* true while sliding into the loop clone */
    var timer = null;
    var typer = null;
    var hovered = false;
    var onScreen = false;

    var dots = [];
    var counter = null;
    if (slides.length > 1 && slides.length <= MAX_DOTS) {
      slides.forEach(function (slide, i) {
        var dot = el("button", "pgh-dot", nav);
        dot.type = "button";
        dot.setAttribute("aria-label", "Show: " + slide.title);
        dot.addEventListener("click", function (event) {
          event.preventDefault();
          goTo(i);
          schedule();
        });
        dots.push(dot);
      });
    } else if (slides.length > 1) {
      counter = el("span", "pgh-counter", nav);
    }

    function typeCaption(text) {
      if (typer) window.clearInterval(typer);
      if (reduceMotion) {
        capText.textContent = text;
        return;
      }
      capText.textContent = "";
      var shown = 0;
      typer = window.setInterval(function () {
        shown += 1;
        capText.textContent = text.slice(0, shown);
        if (shown >= text.length) window.clearInterval(typer);
      }, TYPE_MS);
    }

    function setTransform(position, animate) {
      if (animate) {
        track.style.transform = "translateX(-" + position * 100 + "%)";
        return;
      }
      track.style.transition = "none";
      track.style.transform = "translateX(-" + position * 100 + "%)";
      void track.offsetWidth; /* flush styles so the next change animates */
      track.style.transition = "";
    }

    /* If the track is parked on (or heading to) the loop clone, rebase it
     * onto the real first slide without animating. */
    function settle() {
      if (!pendingSnap) return;
      pendingSnap = false;
      setTransform(0, false);
    }

    track.addEventListener("transitionend", function (event) {
      if (event.target === track) settle();
    });

    function updateMeta() {
      dots.forEach(function (dot, i) {
        dot.classList.toggle("pgh-on", i === current);
      });
      if (counter) {
        counter.textContent = (current + 1) + " / " + slides.length;
      }
      capLink.href = slides[current].url;
      typeCaption(slides[current].title);
    }

    function goTo(next) {
      if (!slides.length) return;
      settle();
      current = ((next % slides.length) + slides.length) % slides.length;
      setTransform(current, true);
      updateMeta();
    }

    function advance() {
      if (current === slides.length - 1) {
        settle();
        current = 0;
        pendingSnap = true;
        setTransform(slides.length, true); /* slide forward into the clone */
        updateMeta();
      } else {
        goTo(current + 1);
      }
    }

    function stop() {
      if (timer) {
        window.clearInterval(timer);
        timer = null;
      }
    }

    function schedule() {
      stop();
      if (reduceMotion || slides.length < 2) return;
      if (hovered || !onScreen || document.hidden) return;
      timer = window.setInterval(advance, ADVANCE_MS);
    }

    card.addEventListener("pointerenter", function () {
      hovered = true;
      stop();
    });
    card.addEventListener("pointerleave", function () {
      hovered = false;
      schedule();
    });
    card.addEventListener("focusin", function () {
      hovered = true;
      stop();
    });
    card.addEventListener("focusout", function () {
      hovered = false;
      schedule();
    });
    document.addEventListener("visibilitychange", schedule);

    if ("IntersectionObserver" in window) {
      new IntersectionObserver(function (entries) {
        onScreen = entries[0].isIntersecting;
        schedule();
      }, { threshold: 0.15 }).observe(card);
    } else {
      onScreen = true;
    }

    goTo(0);
    if (!reduceMotion && slides.length > 1) {
      window.setTimeout(schedule, 600 + index * STAGGER_MS);
    }
    card._pghAdvance = advance; /* programmatic hook (manual testing) */
    return card;
  }

  function init() {
    var root = document.getElementById("pycsamt-gallery-hub");
    var data = window.__PYCSAMT_GALLERY_HUB__;
    if (!root || !data || !data.length) return;
    root.textContent = "";
    data.forEach(function (section, index) {
      if (section.slides.length) root.appendChild(buildCard(section, index));
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
