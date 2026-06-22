/* Plus-menu toggle & outside-click dismiss. */
(function () {
  "use strict";

  function menu() {
    return document.getElementById(
      "am-plus-popover"
    );
  }
  function btn() {
    return document.getElementById(
      "am-btn-plus"
    );
  }
  function close() {
    var m = menu();
    if (m) m.classList.remove("am-open");
  }

  document.addEventListener(
    "click",
    function (e) {
      var b = btn();
      var m = menu();
      if (!b || !m) return;

      if (b.contains(e.target)) {
        m.classList.toggle("am-open");
        return;
      }
      /* click outside → always close */
      if (!m.contains(e.target)) {
        m.classList.remove("am-open");
        return;
      }
      /* click on action (button or link)
         inside menu → close after tick    */
      if (e.target.closest("button, a")) {
        setTimeout(close, 80);
      }
    },
    false
  );
})();
