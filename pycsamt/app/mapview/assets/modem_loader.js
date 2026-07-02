/* modem_loader.js — native folder picker for the ModEM inversion-results
 * tab, mirroring edi_loader.js's pattern: a transparent
 * <input webkitdirectory> is injected inside #mv-btn-inv-browse so the
 * click lands on the OS's native folder dialog instead of an in-app
 * folder navigator.
 *
 * A full ModEM run directory can be hundreds of MB (many iterations x
 * .rho/.dat/.res/.prm, plus large .log files) — far too much to upload
 * as base64. This mirrors InversionResult._scan()'s own file-selection
 * logic client-side (pycsamt/models/modem/results.py) so only the files
 * the backend would actually read are staged: the final-iteration
 * model/data pair, any small control/covariance/log files, non-numbered
 * .dat files. Everything else (.res, .prm, older iterations, plot
 * output, huge logs) is skipped before it's ever read.
 */
(function () {
    "use strict";

    var KEEP_EXTS = [".rho", ".ws", ".dat", ".log", ".inv", ".cov"];
    var MAX_LOG_BYTES = 5 * 1024 * 1024;   // skip a >5MB log; cosmetic only

    function extOf(name) {
        var i = name.lastIndexOf(".");
        return i === -1 ? "" : name.slice(i).toLowerCase();
    }

    function stemOf(name) {
        var base = name.slice(0, name.length - extOf(name).length);
        return base;
    }

    function firstToken(stem) {
        return stem.split("_")[0];
    }

    function iterNum(stem) {
        var m = /_(\d+)$/.exec(stem);
        return m ? parseInt(m[1], 10) : null;
    }

    var OLD_MODEL_RE = /^m(\d+|i)$/i;

    /* Mirror InversionResult._scan()'s file selection: keep only the
     * final-iteration model/data pair (by extension), old-convention
     * m0/mi/d0/di/pr files, non-numbered .dat files, the alphabetically
     * -first .log (size-capped), and any .inv/.cov control files. */
    function selectRelevantFiles(files) {
        var byExtMaxIter = {};   // ext -> {n, files:[...]}
        var keep = [];
        var logCandidates = [];

        files.forEach(function (f) {
            var ext = extOf(f.name);
            if (KEEP_EXTS.indexOf(ext) === -1) return;
            var stem = stemOf(f.name);
            var tok = firstToken(stem);

            if (OLD_MODEL_RE.test(tok) || tok === "d0" || tok === "di" || tok === "pr") {
                keep.push(f);
                return;
            }
            if (ext === ".log") { logCandidates.push(f); return; }
            if (ext === ".inv" || ext === ".cov") { keep.push(f); return; }

            var n = iterNum(stem);
            if (n === null) {
                if (ext === ".dat") keep.push(f);   // non-numbered data file
                return;
            }
            if (ext === ".rho" || ext === ".ws" || ext === ".dat") {
                var bucket = byExtMaxIter[ext];
                if (!bucket || n > bucket.n) {
                    byExtMaxIter[ext] = { n: n, files: [f] };
                } else if (n === bucket.n) {
                    bucket.files.push(f);
                }
            }
            /* .res / .prm and anything else at a numbered iteration is
             * never read by InversionResult — skip. */
        });

        Object.keys(byExtMaxIter).forEach(function (ext) {
            keep = keep.concat(byExtMaxIter[ext].files);
        });

        if (logCandidates.length) {
            logCandidates.sort(function (a, b) {
                return a.name < b.name ? -1 : a.name > b.name ? 1 : 0;
            });
            /* InversionResult picks the alphabetically-first *.log; if
             * that one is too large to upload, fall back to the next
             * smallest-sorted candidate rather than skipping entirely —
             * losing the exact-parity pick is fine, RMS/iteration count
             * is cosmetic metadata, not required to build the map. */
            var log = logCandidates.find(function (f) { return f.size <= MAX_LOG_BYTES; });
            if (log) keep.push(log);
        }
        return keep;
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
        el.textContent = n > 0 ? (n + " file" + (n !== 1 ? "s" : "") + " ready") : "";
    }

    async function sendToStore(files) {
        if (!files.length) {
            hideSpinner();
            setStatus("No ModEM .rho/.dat files found in that folder.");
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
        setStatus("Sending " + filenames.length + " file(s) to server…");
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
            /* Top-level files only — matches InversionResult._scan()'s
             * non-recursive wd.glob("*.ext"), naturally skipping large
             * subfolders (plotdat/, plotmod/, checkpoints, …). */
            var topLevel = all.filter(function (f) {
                var rel = f.webkitRelativePath || f.name;
                return rel.split("/").length === 2;
            });
            if (!topLevel.length) { setStatus("No files found in that folder."); return; }
            showSpinner("Scanning " + topLevel.length + " file"
                        + (topLevel.length !== 1 ? "s" : "") + "…");
            var relevant = selectRelevantFiles(topLevel);
            await sendToStore(relevant);
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
