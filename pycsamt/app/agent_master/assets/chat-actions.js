/* Chat bubble action handlers (copy / EDI). */
(function () {
  "use strict";

  document.addEventListener(
    "click",
    function (e) {

      /* ── copy button ─────────────────────── */
      /* matches both the message copy (.am-copy-btn) and the
         code-block copy (.am-code-copy-btn) */
      var copyBtn = e.target.closest(
        ".am-copy-btn, .am-code-copy-btn"
      );
      if (copyBtn) {
        /* code block copy: data-code attr */
        var codeAttr = copyBtn.getAttribute(
          "data-code"
        );
        if (codeAttr) {
          if (navigator.clipboard) {
            navigator.clipboard
              .writeText(codeAttr)
              .then(function () {
                var ic =
                  copyBtn.querySelector("i");
                var sp =
                  copyBtn.querySelector("span")
                  || copyBtn.lastChild;
                if (ic) {
                  ic.className =
                    "bi bi-clipboard-check me-1";
                }
                var origText =
                  copyBtn.dataset.origText
                  || "Copy";
                copyBtn.dataset.origText =
                  origText;
                if (
                  copyBtn.childNodes.length > 1
                ) {
                  copyBtn.childNodes[
                    copyBtn.childNodes.length - 1
                  ].textContent = "Copied!";
                }
                setTimeout(function () {
                  if (ic) {
                    ic.className =
                      "bi bi-clipboard me-1";
                  }
                  if (
                    copyBtn.childNodes.length > 1
                  ) {
                    copyBtn.childNodes[
                      copyBtn.childNodes.length - 1
                    ].textContent = origText;
                  }
                }, 2000);
              });
          }
          return;
        }

        /* user bubble copy: text content */
        var row = copyBtn.closest(
          ".am-msg-row"
        );
        var bubble = row
          ? row.querySelector(
              ".am-bubble.user"
            )
          : null;
        var txt = bubble
          ? bubble.textContent.trim()
          : "";
        if (txt && navigator.clipboard) {
          navigator.clipboard
            .writeText(txt)
            .then(function () {
              var ic =
                copyBtn.querySelector("i");
              if (ic) {
                ic.className =
                  "bi bi-clipboard-check";
                setTimeout(function () {
                  ic.className =
                    "bi bi-clipboard";
                }, 1800);
              }
            });
        }
        return;
      }

      /* ── EDI button ──────────────────────── */
      var ediBtn = e.target.closest(
        ".am-edi-msg-btn"
      );
      if (ediBtn) {
        var loadBtn = document.getElementById(
          "am-btn-load-edi"
        );
        if (loadBtn) loadBtn.click();
      }
    },
    false
  );
})();
