/* pcsf_loader.js — native folder picker for the inversion-results tab,
 * mirroring edi_loader.js's pattern: a transparent
 * <input webkitdirectory> is injected inside #mv-btn-inv-browse so the
 * click lands on the OS's native folder dialog instead of an in-app
 * folder navigator.
 *
 * Unlike a raw ModEM run directory (which can be hundreds of MB across
 * many iterations and is never uploaded whole), a PCSF/PCSM result is
 * one small self-contained file per line -- so every matching file
 * found in the picked folder is simply read and staged; there is no
 * file-selection heuristic to mirror here the way modem_loader.js
 * (its predecessor) had to have for InversionResult._scan().
 */
(function () {
    "use strict";

    var MATCH_RE = /\.(pcsf|pcsm|pcsm\.gz)$/i;

    function isMatch(name) {
        return MATCH_RE.test(name);
    }

    function readFileAsDataURL(file) {
        return new Promise(function (resolve, reject) {
            var r = new FileReader();
            r.onload = function (ev) { resolve(ev.target.result); };
            r.onerror = function () { reject(new Error("read error")); };
            r.readAsDataURL(file);
        });
    }

    function showSpinner(msg) {
        var ov = document.getElementById("mv-inv-loader-overlay");
        var lbl = document.getElementById("mv-inv-loader-msg");
        if (lbl) lbl.textContent = msg || "Reading files…";
        if (ov) ov.style.display = "flex";
    }
    function hideSpinner() {
        var ov = document.getElementById("mv-inv-loader-overlay");
        if (ov) ov.style.display = "none";
    }
    function setStatus(msg) {
        var el = document.getElementById("mv-inv-browse-status");
        if (el) el.textContent = msg || "";
    }
    function setCount(n) {
        var el = document.getElementById("mv-inv-file-count");
        if (!el) return;
        el.textContent = n > 0 ? (n + " file" + (n !== 1 ? "s" : "") + " found") : "";
    }

    async function sendToStore(files) {
        if (!files.length) {
            /* Deliberately do NOT write to the store here: a downstream
             * Python callback chain (capture -> classify) re-renders
             * this same status span from the (empty) candidate list,
             * which would immediately blank out the message below if
             * it fired. Setting the DOM text directly and stopping is
             * the same pattern edi_loader.js/modem_loader.js already
             * used for their own "nothing found" case. */
            hideSpinner();
            setStatus("No .pcsf/.pcsm/.pcsm.gz files found in that folder "
                + "— convert your ModEM/Occam2D/MARE2DEM results first.");
            setCount(0);
            return;
        }
        showSpinner("Reading " + files.length + " file"
                    + (files.length !== 1 ? "s" : "") + "…");
        var settled = await Promise.all(files.map(function (f) {
            return readFileAsDataURL(f).then(
                function (b64) { return { ok: true, name: f.name, b64: b64 }; },
                function () { return { ok: false }; }
            );
        }));
        var filenames = [], contents = [];
        settled.forEach(function (r) {
            if (r.ok) { filenames.push(r.name); contents.push(r.b64); }
        });
        setCount(filenames.length);
        setStatus(filenames.length + " candidate(s) found — pick one below.");
        if (window.dash_clientside && window.dash_clientside.set_props) {
            window.dash_clientside.set_props("mv-inv-folder-store", {
                data: { filenames: filenames, contents: contents }
            });
        }
        setTimeout(hideSpinner, 400);
    }

    function makeFolderInput() {
        var inp = document.createElement("input");
        inp.type = "file";
        inp.multiple = true;
        inp.setAttribute("webkitdirectory", "");
        inp.setAttribute("directory", "");
        inp.style.cssText = [
            "position:absolute", "inset:0", "width:100%", "height:100%",
            "opacity:0", "cursor:pointer", "font-size:0",
        ].join(";");
        inp.addEventListener("change", async function () {
            var all = Array.from(inp.files || []);
            inp.value = "";
            /* Top-level files only -- a per-line PCSF/PCSM export is
             * expected directly inside the picked folder, not nested
             * further. */
            var topLevel = all.filter(function (f) {
                var rel = f.webkitRelativePath || f.name;
                return rel.split("/").length === 2;
            });
            showSpinner("Scanning " + topLevel.length + " file"
                        + (topLevel.length !== 1 ? "s" : "") + "…");
            var matches = topLevel.filter(function (f) { return isMatch(f.name); });
            await sendToStore(matches);
        });
        return inp;
    }

    function injectBrowseInput() {
        var btn = document.getElementById("mv-btn-inv-browse");
        if (!btn || btn.dataset.folderInjected) return;
        btn.dataset.folderInjected = "1";
        btn.style.position = "relative";
        btn.style.overflow = "hidden";
        btn.appendChild(makeFolderInput());
    }

    function watchForBrowseButton() {
        injectBrowseInput();
        var obs = new MutationObserver(function () { injectBrowseInput(); });
        obs.observe(document.body, { childList: true, subtree: true });
    }

    function init() { watchForBrowseButton(); }
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();
