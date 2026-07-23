/* pyCSAMT — "hosted instance not open yet" dialog (applications/hosted.rst).
 *
 * The three hosted-app cards on that page are sphinx-design cards with
 * :link: pointing at real subdomains (agent-master/web/mapview.pycsamt.org).
 * Those instances aren't ready for public, per-user-isolated use yet, so
 * this script intercepts the click on each card's stretched link and shows
 * a dialog explaining that instead of following the link. The app name and
 * "install locally" href are read from markup already on the card (the
 * visible title, and a hidden .pyc-hosted-install :doc: link) rather than
 * duplicated here, so they stay correct if the docs ever move.
 */
(function () {
  "use strict";

  var LOGO_SRC = "/_static/logo/pycsamt-v2-symbol.svg";
  var MESSAGE =
    "A secured deployment — authentication and per-user data " +
    "isolation — is in progress for this hosted instance. Until it " +
    "ships, treat it as a shared demo server: don't load sensitive data, " +
    "and check back on this page for the update.";

  var reduceMotion =
    window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  var overlay = null;
  var panel = null;
  var titleEl = null;
  var closeBtn = null;
  var installLink = null;
  var lastFocused = null;

  function focusableItems() {
    return [closeBtn, installLink, panel.querySelector(".phd-btn-secondary")].filter(
      function (node) {
        return node && node.offsetParent !== null;
      }
    );
  }

  function onKeydown(event) {
    if (event.key === "Escape") {
      event.preventDefault();
      close();
      return;
    }
    if (event.key !== "Tab") return;
    var items = focusableItems();
    if (!items.length) return;
    var first = items[0];
    var last = items[items.length - 1];
    if (event.shiftKey && document.activeElement === first) {
      event.preventDefault();
      last.focus();
    } else if (!event.shiftKey && document.activeElement === last) {
      event.preventDefault();
      first.focus();
    }
  }

  function build() {
    overlay = document.createElement("div");
    overlay.className = "phd-overlay";
    overlay.setAttribute("role", "presentation");

    panel = document.createElement("div");
    panel.className = "phd-panel";
    panel.setAttribute("role", "dialog");
    panel.setAttribute("aria-modal", "true");
    panel.setAttribute("aria-labelledby", "phd-title");
    panel.setAttribute("aria-describedby", "phd-body");

    closeBtn = document.createElement("button");
    closeBtn.type = "button";
    closeBtn.className = "phd-close";
    closeBtn.setAttribute("aria-label", "Close");
    closeBtn.innerHTML = "&times;";
    closeBtn.addEventListener("click", close);

    var badge = document.createElement("div");
    badge.className = "phd-badge";
    var ring = document.createElement("div");
    ring.className = "phd-badge-ring";
    ring.setAttribute("aria-hidden", "true");
    var logo = document.createElement("img");
    logo.src = LOGO_SRC;
    logo.alt = "";
    logo.setAttribute("aria-hidden", "true");
    badge.appendChild(ring);
    badge.appendChild(logo);

    titleEl = document.createElement("h2");
    titleEl.className = "phd-title";
    titleEl.id = "phd-title";

    var bodyEl = document.createElement("p");
    bodyEl.className = "phd-body";
    bodyEl.id = "phd-body";
    bodyEl.textContent = MESSAGE;

    var actions = document.createElement("div");
    actions.className = "phd-actions";

    installLink = document.createElement("a");
    installLink.className = "phd-btn phd-btn-primary";

    var gotIt = document.createElement("button");
    gotIt.type = "button";
    gotIt.className = "phd-btn phd-btn-secondary";
    gotIt.textContent = "Got it";
    gotIt.addEventListener("click", close);

    actions.appendChild(installLink);
    actions.appendChild(gotIt);

    panel.appendChild(closeBtn);
    panel.appendChild(badge);
    panel.appendChild(titleEl);
    panel.appendChild(bodyEl);
    panel.appendChild(actions);
    overlay.appendChild(panel);
    document.body.appendChild(overlay);

    overlay.addEventListener("click", function (event) {
      if (event.target === overlay) close();
    });
    panel.addEventListener("click", function (event) {
      event.stopPropagation();
    });
  }

  function open(info) {
    if (!overlay) build();

    titleEl.textContent = info.name + " Isn't Open Yet";
    if (info.installHref) {
      installLink.href = info.installHref;
      installLink.textContent =
        "Install " + info.name + " locally instead";
      installLink.style.display = "";
    } else {
      installLink.style.display = "none";
    }

    lastFocused = document.activeElement;
    overlay.classList.add("phd-open");
    document.addEventListener("keydown", onKeydown, true);
    window.setTimeout(function () {
      closeBtn.focus();
    }, reduceMotion ? 0 : 60);
  }

  function close() {
    if (!overlay || !overlay.classList.contains("phd-open")) return;
    overlay.classList.remove("phd-open");
    document.removeEventListener("keydown", onKeydown, true);
    if (lastFocused && typeof lastFocused.focus === "function") {
      lastFocused.focus();
    }
  }

  function wireCard(card) {
    var link = card.querySelector("a.sd-stretched-link");
    if (!link) return;

    var nameEl = card.querySelector(".sd-card-title");
    var name = nameEl ? nameEl.textContent.trim() : "This app";
    var installEl = card.querySelector(".pyc-hosted-install a");
    var installHref = installEl ? installEl.getAttribute("href") : null;

    link.addEventListener("click", function (event) {
      event.preventDefault();
      open({ name: name, installHref: installHref });
    });
  }

  function init() {
    var grid = document.querySelector(".pyc-hosted-grid");
    if (!grid) return;
    Array.prototype.forEach.call(grid.querySelectorAll(".sd-card"), wireCard);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
