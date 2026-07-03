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
    Install: [
      {
        icon: "fa-solid fa-flag-checkered",
        title: "Getting Started",
        desc: "A guided first run: install, load a survey, and produce your first result.",
        href: "getting_started/index.html"
      },
      {
        icon: "fa-solid fa-download",
        title: "Installation",
        desc: "pip/conda matrix, optional extras, and backend requirements.",
        href: "installation.html"
      }
    ],
    API: [
      {
        icon: "fa-solid fa-sliders",
        title: "Configure",
        desc: "Set up the 11 configure_* families that drive every pipeline stage.",
        href: "api_guide/index.html"
      },
      {
        icon: "fa-solid fa-book",
        title: "Reference",
        desc: "Full autodoc-generated API reference, organised by subpackage.",
        href: "api/index.html"
      }
    ],
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
      }
    ],
    "AI Agents": [
      {
        icon: "fa-solid fa-robot",
        title: "Agent Master",
        desc: "The orchestrating chat agent — plans and runs multi-step workflows end to end, grounded by our RAG assistant.",
        href: "agents/orchestrator.html"
      },
      {
        icon: "fa-solid fa-users-gear",
        title: "Agents",
        desc: "Browse every specialised agent: foundation, processing, inversion, and more.",
        href: "agents/agent_catalogue.html"
      },
      {
        icon: "fa-solid fa-gear",
        title: "Configure",
        desc: "Point pyCSAMT at your LLM provider and tune agent behaviour.",
        href: "agents/llm_configuration.html"
      }
    ],
    Modeling: [
      {
        icon: "fa-solid fa-layer-group",
        title: "Models",
        desc: "Classical solvers — Occam2D, ModEM, and MARE2DEM — end to end.",
        href: "models/index.html"
      },
      {
        icon: "fa-solid fa-brain",
        title: "AI Inversion",
        desc: "Physics-informed neural networks and hybrid deep-learning inverters.",
        href: "user_guide/ai_inversion.html"
      },
      {
        icon: "fa-solid fa-wave-square",
        title: "Forward Model",
        desc: "Build synthetic models and simulate forward responses for survey design.",
        href: "forward/index.html"
      }
    ],
    // TODO: add an "EM" card (icon fa-solid fa-bolt, href "emtools/index.html")
    // once emtools/index.rst exists for pycsamt/emtools/.
    Tools: [
      {
        icon: "fa-solid fa-location-dot",
        title: "Site",
        desc: "Station-centric tools: EDI sites, survey-line collections, diagnostics.",
        href: "site/index.html"
      },
      {
        icon: "fa-solid fa-map",
        title: "Maps",
        desc: "Station maps, pseudosections, 3-D quick-look views, and overlays.",
        href: "map/index.html"
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
