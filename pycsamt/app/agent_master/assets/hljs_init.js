/**
 * Auto-highlight code blocks whenever new chat
 * bubbles land in the DOM (MutationObserver).
 * Requires Highlight.js to already be loaded
 * via external_scripts in app.py.
 */
(function () {
  "use strict";

  function highlightNew(node) {
    if (!window.hljs) return;
    node.querySelectorAll(
      "code.language-python:not(.hljs)"
    ).forEach(function (el) {
      hljs.highlightElement(el);
    });
  }

  function init() {
    var chatWin = document.getElementById(
      "am-chat-window"
    );
    if (!chatWin) {
      /* retry if Dash hasn't rendered yet */
      setTimeout(init, 400);
      return;
    }
    /* highlight any blocks already present */
    highlightNew(chatWin);

    var obs = new MutationObserver(
      function (mutations) {
        mutations.forEach(function (m) {
          m.addedNodes.forEach(function (n) {
            if (n.nodeType === 1) {
              highlightNew(n);
            }
          });
        });
      }
    );
    obs.observe(chatWin, {
      childList: true,
      subtree: true,
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener(
      "DOMContentLoaded", init
    );
  } else {
    init();
  }
})();
