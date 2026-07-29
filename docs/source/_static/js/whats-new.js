/* pyCSAMT v2 — "What's new" badge.
 * Time-boxed (TTL_DAYS from the latest changelog release date) *and*
 * per-visitor dismissal via localStorage: a first-time visitor arriving on
 * day 6 still sees it once, but a returning visitor who already dismissed
 * it -- or who has already read the release note / changelog page -- is
 * not nagged again for the remaining days.
 *
 * Data comes from _static/data/whats_new.json, generated at build time
 * from changelog.rst (see conf.py: _write_whats_new_json). No file (no
 * dated release yet), a fetch error, or blocked localStorage all just mean
 * no badge -- this script never hard-fails the page.
 */
(function () {
  "use strict";

  var TTL_DAYS = 7;
  var STORAGE_KEY = "pycsamt:whats-new:dismissed-version";

  function ready(fn) {
    if (document.readyState !== "loading") fn();
    else document.addEventListener("DOMContentLoaded", fn);
  }

  // Same technique as mega-menu.js: Sphinx bakes the correct relative
  // "_static/" prefix into every page's own asset links at build time.
  function rootPrefix() {
    var link = document.querySelector(
      'link[rel="stylesheet"][href*="_static/"]'
    );
    if (!link) return "";
    var href = link.getAttribute("href");
    var idx = href.indexOf("_static/");
    return idx === -1 ? "" : href.slice(0, idx);
  }

  function getDismissed() {
    try {
      return window.localStorage.getItem(STORAGE_KEY);
    } catch (e) {
      return null;
    }
  }

  function setDismissed(version) {
    try {
      window.localStorage.setItem(STORAGE_KEY, version);
    } catch (e) {
      /* Private mode / storage blocked: badge just keeps showing until
       * the TTL window closes on its own. */
    }
  }

  function currentPageMatchesRelease(info) {
    var path = window.location.pathname;
    var versionPattern = info.version.replace(/\./g, "\\.");
    return (
      /\/changelog\.html$/.test(path) ||
      /\/release_notes\/(index\.html)?$/.test(path) ||
      new RegExp(
        "/release_notes/v" + versionPattern + "\\.html$"
      ).test(path)
    );
  }

  function removeAll() {
    var banner = document.querySelector(".pyc-whats-new-banner");
    if (banner) banner.remove();
    var badges = document.querySelectorAll(".pyc-badge-new");
    for (var i = 0; i < badges.length; i++) badges[i].remove();
  }

  function dismiss(version) {
    setDismissed(version);
    removeAll();
  }

  function addBadges() {
    var links = document.querySelectorAll(
      'a[href$="changelog.html"], a[href$="release_notes/index.html"]'
    );
    links.forEach(function (a) {
      if (a.closest(".pyc-whats-new-banner")) return;
      if (a.querySelector(".pyc-badge-new")) return;
      var badge = document.createElement("span");
      badge.className = "pyc-badge-new";
      badge.textContent = "New";
      a.appendChild(document.createTextNode(" "));
      a.appendChild(badge);
    });
  }

  function addBanner(info) {
    var header = document.querySelector(".bd-header");
    if (!header) return;

    var banner = document.createElement("aside");
    banner.className = "pyc-whats-new-banner";
    banner.setAttribute("aria-label", "Latest release announcement");

    var link = document.createElement("a");
    link.href = rootPrefix() + info.url;
    link.className = "pyc-whats-new-link";
    link.innerHTML =
      '<span class="pyc-badge-new">New</span> pyCSAMT ' +
      info.version +
      " is out — see what changed";
    link.addEventListener("click", function () {
      setDismissed(info.version);
    });

    var dismissBtn = document.createElement("button");
    dismissBtn.type = "button";
    dismissBtn.className = "pyc-whats-new-dismiss";
    dismissBtn.setAttribute("aria-label", "Dismiss");
    dismissBtn.textContent = "×";
    dismissBtn.addEventListener("click", function () {
      dismiss(info.version);
    });

    banner.appendChild(link);
    banner.appendChild(dismissBtn);
    header.insertAdjacentElement("afterend", banner);
  }

  ready(function () {
    var prefix = rootPrefix();
    fetch(prefix + "_static/data/whats_new.json", { cache: "no-store" })
      .then(function (resp) {
        return resp.ok ? resp.json() : null;
      })
      .then(function (info) {
        if (!info || !info.version || !info.date) return;

        if (currentPageMatchesRelease(info)) {
          // Reading the release note / changelog itself counts as "seen".
          setDismissed(info.version);
          return;
        }

        if (getDismissed() === info.version) return;

        var releasedAt = new Date(info.date + "T00:00:00Z").getTime();
        var ageDays = (Date.now() - releasedAt) / 86400000;
        if (ageDays < 0 || ageDays >= TTL_DAYS) return;

        addBadges();
        addBanner(info);
      })
      .catch(function () {
        /* No dated release yet, offline, or blocked request: no badge. */
      });
  });
})();
