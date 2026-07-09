/* pyCSAMT v2 — mega-menu navigation cards.
 * Site-wide (every page shares the same navbar). Converts the top-level
 * nav links named as keys in MENUS from plain links into toggles that
 * reveal a panel of cards instead of navigating directly.
 *
 * To add another menu later: add a "<Nav label>": [ {icon, title, desc,
 * href}, ... ] entry below. href is root-relative (e.g. "api/index.html");
 * the actual per-page prefix is resolved automatically.
 */
(function () {
  "use strict";

  var MENUS = {
    "User guide": [
      {
        icon: "fa-solid fa-compass",
        title: "Guide",
        desc: "Task-based walkthroughs: data loading, processing, inversion, interpretation, mapping.",
        href: "user_guide/index.html"
      },
      {
        icon: "fa-solid fa-graduation-cap",
        title: "Tutorials",
        desc: "Worked, end-to-end notebooks and scripts you can adapt directly.",
        href: "tutorials/index.html"
      },
      {
        icon: "fa-solid fa-bolt",
        title: "EM tools",
        desc: "Module-by-module guide to pycsamt.emtools: QC, dimensionality, strike, corrections.",
        href: "user_guide/emtools/index.html"
      },
      {
        icon: "fa-solid fa-location-dot",
        title: "Site tools",
        desc: "Station-centric tools: EDI sites, survey-line collections, diagnostics.",
        href: "user_guide/site/index.html"
      },
      {
        icon: "fa-solid fa-map",
        title: "Map tools",
        desc: "Station maps, pseudosections, 3-D quick-look views, and overlays.",
        href: "user_guide/map/index.html"
      }
    ]
  };

  function ready(fn) {
    if (document.readyState !== "loading") fn();
    else document.addEventListener("DOMContentLoaded", fn);
  }

  // Sphinx bakes correct relative paths into every page's own asset links
  // at build time; reuse that instead of re-deriving page depth ourselves.
  function rootPrefix() {
    var link = document.querySelector('link[rel="stylesheet"][href*="_static/"]');
    if (!link) return "";
    var href = link.getAttribute("href");
    var idx = href.indexOf("_static/");
    return idx === -1 ? "" : href.slice(0, idx);
  }

  function buildPanel(id, cards, prefix) {
    var panel = document.createElement("div");
    panel.className = "pyc-mega";
    panel.id = id;
    panel.hidden = true;

    var inner = document.createElement("div");
    inner.className = "pyc-mega-inner";

    cards.forEach(function (card) {
      var a = document.createElement("a");
      a.className = "pyc-mega-card";
      a.href = prefix + card.href;
      a.innerHTML =
        '<span class="pyc-mega-card-icon"><i class="' +
        card.icon +
        '" aria-hidden="true"></i></span>' +
        '<span class="pyc-mega-card-body">' +
        '<span class="pyc-mega-card-title">' +
        card.title +
        "</span>" +
        '<span class="pyc-mega-card-desc">' +
        card.desc +
        "</span>" +
        "</span>";
      inner.appendChild(a);
    });

    panel.appendChild(inner);
    return panel;
  }

  ready(function () {
    var header = document.querySelector(".bd-header");
    var navLinks = document.querySelectorAll(".navbar-nav > li > a");
    if (!header || !navLinks.length) return;

    if (window.getComputedStyle(header).position === "static") {
      header.style.position = "relative";
    }

    var prefix = rootPrefix();
    var openPanel = null;
    var openTrigger = null;

    function closeOpen() {
      if (!openPanel) return;
      openPanel.hidden = true;
      openTrigger.setAttribute("aria-expanded", "false");
      openPanel = null;
      openTrigger = null;
    }

    Object.keys(MENUS).forEach(function (name) {
      var link = Array.prototype.filter
        .call(navLinks, function (a) {
          return a.textContent.trim() === name;
        })
        .shift();
      if (!link) return;

      var id = "pyc-mega-" + name.toLowerCase().replace(/\s+/g, "-");
      var panel = buildPanel(id, MENUS[name], prefix);
      header.appendChild(panel);

      link.classList.add("pyc-mega-trigger");
      link.setAttribute("aria-haspopup", "true");
      link.setAttribute("aria-expanded", "false");
      link.setAttribute("aria-controls", id);

      var caret = document.createElement("i");
      caret.className = "fa-solid fa-chevron-down pyc-mega-caret";
      caret.setAttribute("aria-hidden", "true");
      link.appendChild(caret);

      link.addEventListener("click", function (e) {
        e.preventDefault();
        var wasOpen = openPanel === panel;
        closeOpen();
        if (!wasOpen) {
          panel.hidden = false;
          link.setAttribute("aria-expanded", "true");
          openPanel = panel;
          openTrigger = link;
        }
      });
    });

    document.addEventListener("click", function (e) {
      if (
        openTrigger &&
        !openTrigger.contains(e.target) &&
        !openPanel.contains(e.target)
      ) {
        closeOpen();
      }
    });

    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape") closeOpen();
    });
  });
})();
