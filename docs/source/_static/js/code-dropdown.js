/* Progressive enhancements for native pyCSAMT code disclosures. */
(function () {
  "use strict";

  function enhance(panel) {
    if (panel.dataset.pycEnhanced === "true") return;
    panel.dataset.pycEnhanced = "true";

    var summary = panel.querySelector(":scope > .pyc-code-dropdown__summary");
    var hint = summary && summary.querySelector(".pyc-code-dropdown__hint");
    if (!summary) return;

    function updateHint() {
      if (!hint) return;
      hint.textContent = panel.open ? hint.dataset.close : hint.dataset.open;
    }

    panel.addEventListener("toggle", updateHint);
    updateHint();

    var pre = panel.querySelector("pre");
    if (pre && pre.textContent.split("\n").length > 36) {
      var collapse = document.createElement("button");
      collapse.type = "button";
      collapse.className = "pyc-code-dropdown__collapse";
      collapse.innerHTML = '<span aria-hidden="true">&#8593;</span> Collapse source code';
      collapse.addEventListener("click", function () {
        var reduceMotion = window.matchMedia &&
          window.matchMedia("(prefers-reduced-motion: reduce)").matches;
        panel.open = false;
        summary.focus({ preventScroll: true });
        panel.scrollIntoView({
          behavior: reduceMotion ? "auto" : "smooth",
          block: "start"
        });
      });
      panel.appendChild(collapse);
    }
  }

  function init(root) {
    var scope = root || document;
    scope.querySelectorAll("details.pyc-code-dropdown").forEach(enhance);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", function () { init(document); });
  } else {
    init(document);
  }
})();
