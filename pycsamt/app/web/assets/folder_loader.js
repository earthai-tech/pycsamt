/* folder_loader.js
 * Adds recursive folder browse + folder drag-and-drop to the load modal.
 * Results are written to the Dash store "folder-upload-store" via
 * window.dash_clientside.set_props so that the Python callback can load them.
 *
 * Browse Folder approach: a transparent <input type="file" webkitdirectory>
 * is injected *inside* the Browse Folder button, stretched to cover it
 * completely.  The user's click lands directly on the file input so the
 * browser opens the folder picker natively — no programmatic input.click()
 * is needed (which some browsers block for webkitdirectory inputs when called
 * from a delegated event handler).
 *
 * Performance: all FileReader calls run in parallel via Promise.all so that
 * reading 128 files takes ~1× the time of the slowest file instead of 128×.
 */

(function () {
    "use strict";

    var _EXTS = [".edi", ".avg", ".j"];

    /* ── Helpers ─────────────────────────────────────────────────────────── */

    function hasExt(name) {
        var lo = name.toLowerCase();
        return _EXTS.some(function (e) { return lo.endsWith(e); });
    }

    function readFileAsDataURL(file) {
        return new Promise(function (resolve, reject) {
            var r = new FileReader();
            r.onload  = function (ev) { resolve(ev.target.result); };
            r.onerror = function ()   { reject(new Error("FileReader error")); };
            r.readAsDataURL(file);
        });
    }

    function entryToFile(entry) {
        return new Promise(function (resolve, reject) {
            entry.file(resolve, reject);
        });
    }

    /* Read all entries from a DirectoryReader in batches (readEntries
     * returns at most 100 entries per call on Chrome). */
    function readAllEntries(dirReader) {
        return new Promise(function (resolve, reject) {
            var all = [];
            function batch() {
                dirReader.readEntries(function (entries) {
                    if (!entries.length) { resolve(all); return; }
                    all = all.concat(Array.from(entries));
                    batch();
                }, reject);
            }
            batch();
        });
    }

    /* Recursively collect {relativePath, file} pairs for all EDI/AVG/J
     * files under a FileSystemEntry. */
    async function collectEntryFiles(entry, basePath) {
        var relPath = basePath ? basePath + "/" + entry.name : entry.name;
        if (entry.isFile) {
            if (!hasExt(entry.name)) return [];
            var file = await entryToFile(entry);
            return [{ relativePath: relPath, file: file }];
        }
        if (entry.isDirectory) {
            var reader  = entry.createReader();
            var entries = await readAllEntries(reader);
            var results = [];
            for (var i = 0; i < entries.length; i++) {
                var sub = await collectEntryFiles(entries[i], relPath);
                results = results.concat(sub);
            }
            return results;
        }
        return [];
    }

    /* ── Spinner / status helpers ─────────────────────────────────────────── */

    function showSpinner(msg) {
        var ov  = document.getElementById("file-loading-overlay");
        var lbl = document.getElementById("file-loader-msg");
        if (lbl) lbl.textContent = msg || "Reading files…";
        if (ov)  ov.style.display = "flex";

        /* Also hide the drag zone so the spinner takes its place */
        var dz = document.getElementById("upload-drop-wrap");
        if (dz) dz.style.display = "none";
    }

    function hideSpinner() {
        var ov = document.getElementById("file-loading-overlay");
        if (ov) ov.style.display = "none";
        /* Drag zone is hidden by the Python callback once the file manager
         * renders, so we don't need to restore it here. */
    }

    function setStatus(msg, type) {
        var el = document.getElementById("folder-browse-status");
        if (!el) return;
        el.textContent = msg;
        el.className   = "small mt-1 folder-browse-status-" + (type || "info");
    }

    function setCount(n) {
        var badge = document.getElementById("folder-file-count");
        if (!badge) return;
        badge.textContent = n + " file" + (n !== 1 ? "s" : "") + " ready";
        badge.style.display = n > 0 ? "inline-block" : "none";
    }

    /* ── Store update (parallel reads) ───────────────────────────────────── */

    async function sendToStore(pairs) {
        if (!pairs.length) {
            hideSpinner();
            setStatus("No EDI / AVG / J files found.", "warning");
            setCount(0);
            return;
        }

        var n = pairs.length;
        showSpinner("Reading " + n + " file" + (n !== 1 ? "s" : "") + "…");

        /* Read all files in parallel — Promise.all fires every FileReader
         * simultaneously instead of waiting for each one in sequence.
         * On 128 files this is typically 5–10× faster than a for-await loop. */
        var settled = await Promise.all(
            pairs.map(function (pair) {
                return readFileAsDataURL(pair.file).then(
                    function (b64) { return { ok: true,  path: pair.relativePath, b64: b64 }; },
                    function ()    { return { ok: false, path: pair.relativePath }; }
                );
            })
        );

        var filenames = [];
        var contents  = [];
        settled.forEach(function (r) {
            if (r.ok) {
                filenames.push(r.path);
                contents.push(r.b64);
            } else {
                console.warn("folder_loader: could not read", r.path);
            }
        });

        setCount(filenames.length);
        setStatus("Sending " + filenames.length + " file(s) to server…", "info");

        /* Push to Dash store (Dash 2.9+ set_props API).
         * The Python render_upload_manager callback will fire next;
         * the MutationObserver below hides the spinner when it renders. */
        if (window.dash_clientside && window.dash_clientside.set_props) {
            window.dash_clientside.set_props("folder-upload-store", {
                data: { filenames: filenames, contents: contents }
            });
        }
    }

    /* ── Auto-hide spinner when file manager renders ─────────────────────── */

    function watchFileManager() {
        var mgr = document.getElementById("load-upload-file-manager");
        if (!mgr) return;

        var obs = new MutationObserver(function () {
            if (mgr.children.length > 0) {
                hideSpinner();
            }
        });
        obs.observe(mgr, { childList: true });
    }

    /* Re-attach the MutationObserver each time the modal opens (React
     * recreates the DOM subtree on each open). */
    function watchModalForManager() {
        var body = document.body;
        var bodyObs = new MutationObserver(function () {
            var mgr = document.getElementById("load-upload-file-manager");
            if (mgr && !mgr._loaderWatched) {
                mgr._loaderWatched = true;
                watchFileManager();
            }
        });
        bodyObs.observe(body, { childList: true, subtree: true });
    }

    /* ── Browse-folder button: inject a transparent file input ────────────── */

    function makeFolderInput() {
        var inp = document.createElement("input");
        inp.type = "file";
        inp.multiple = true;
        inp.setAttribute("webkitdirectory", "");
        inp.setAttribute("directory", "");
        /* Cover the entire button area; opacity:0 keeps it invisible */
        inp.style.cssText = [
            "position:absolute",
            "inset:0",
            "width:100%",
            "height:100%",
            "opacity:0",
            "cursor:pointer",
            "font-size:0",
        ].join(";");

        inp.addEventListener("change", async function () {
            var files = Array.from(inp.files || []).filter(function (f) {
                return hasExt(f.name);
            });
            inp.value = "";  /* reset so re-selecting same folder fires */

            if (!files.length) {
                setStatus("No EDI / AVG / J files found.", "warning");
                return;
            }

            /* Show spinner immediately — before any async work */
            showSpinner("Scanning " + files.length + " file" + (files.length !== 1 ? "s" : "") + "…");

            var pairs = files.map(function (f) {
                return { relativePath: f.webkitRelativePath || f.name, file: f };
            });
            await sendToStore(pairs);
        });

        return inp;
    }

    function injectBrowseInput() {
        var btn = document.getElementById("btn-browse-folder");
        if (!btn) return false;
        if (btn.dataset.folderInjected) return true;  /* already done */

        btn.dataset.folderInjected = "1";
        btn.style.position = "relative";
        btn.style.overflow = "hidden";
        btn.appendChild(makeFolderInput());
        return true;
    }

    function watchForBrowseButton() {
        /* Try immediately in case the modal is already open */
        injectBrowseInput();

        /* React/Dash re-creates the button each time the modal opens.
         * MutationObserver catches the new element and re-injects. */
        var obs = new MutationObserver(function () {
            var btn = document.getElementById("btn-browse-folder");
            if (!btn || btn.dataset.folderInjected) return;
            injectBrowseInput();
        });
        obs.observe(document.body, { childList: true, subtree: true });
    }

    /* ── Drag-and-drop folder support ─────────────────────────────────────── */

    function setupDragDrop() {
        /* Use capture phase so we intercept before dcc.Upload's own handler */
        document.addEventListener("dragover", function (e) {
            var zone = e.target.closest(".upload-drop");
            if (!zone) return;

            var items = e.dataTransfer && e.dataTransfer.items;
            if (!items) return;

            /* Only intervene if at least one dropped item is a directory */
            var hasDir = false;
            for (var i = 0; i < items.length; i++) {
                if (items[i].kind !== "file") continue;
                var entry = items[i].webkitGetAsEntry && items[i].webkitGetAsEntry();
                if (entry && entry.isDirectory) { hasDir = true; break; }
            }
            if (hasDir) {
                e.preventDefault();
                e.stopPropagation();
                zone.style.borderColor = "#89dceb";
            }
        }, true);

        document.addEventListener("dragleave", function (e) {
            var zone = e.target.closest(".upload-drop");
            if (zone) zone.style.borderColor = "";
        }, true);

        document.addEventListener("drop", function (e) {
            var zone = e.target.closest(".upload-drop");
            if (!zone) return;

            var items = e.dataTransfer && e.dataTransfer.items;
            if (!items) return;

            var entries = [];
            for (var i = 0; i < items.length; i++) {
                if (items[i].kind !== "file") continue;
                var entry = items[i].webkitGetAsEntry && items[i].webkitGetAsEntry();
                if (entry) entries.push(entry);
            }

            var hasDir = entries.some(function (en) { return en.isDirectory; });
            if (!hasDir) {
                /* Plain files — dcc.Upload handles the actual reading.
                 * Show spinner immediately; MutationObserver hides it once
                 * the file manager renders (same path as folder browse). */
                var nDrop = entries.filter(function (en) {
                    return en.isFile && hasExt(en.name);
                }).length;
                if (nDrop > 0) {
                    showSpinner("Reading " + nDrop + " file" + (nDrop !== 1 ? "s" : "") + "…");
                }
                return;
            }

            e.preventDefault();
            e.stopPropagation();
            zone.style.borderColor = "";

            /* Show spinner immediately */
            showSpinner("Scanning dropped folder(s)…");

            (async function () {
                var allPairs = [];
                for (var j = 0; j < entries.length; j++) {
                    var pairs = await collectEntryFiles(entries[j], "");
                    /* Strip a leading slash that collectEntryFiles may produce
                     * when basePath is empty. */
                    pairs = pairs.map(function (p) {
                        return {
                            relativePath: p.relativePath.replace(/^\//, ""),
                            file: p.file
                        };
                    });
                    allPairs = allPairs.concat(pairs);
                }
                await sendToStore(allPairs);
            })();
        }, true);
    }

    /* ── Show spinner when dcc.Upload's file picker is used ─────────────────
     * dcc.Upload renders a hidden <input type="file"> inside #upload-edi.
     * When the user clicks the zone and picks files, that input fires a
     * "change" event.  We intercept it in capture phase (fires before
     * React/dcc.Upload handlers) to show the spinner immediately.         */

    function setupUploadInputWatch() {
        document.addEventListener("change", function (e) {
            if (!e.target || e.target.type !== "file") return;
            /* Only react to file inputs inside the dcc.Upload zone,
             * not the Browse Folder input we inject into #btn-browse-folder */
            if (!e.target.closest("#upload-edi")) return;
            var files = Array.from(e.target.files || []);
            var valid = files.filter(function (f) { return hasExt(f.name); });
            if (valid.length > 0) {
                showSpinner("Reading " + valid.length +
                            " file" + (valid.length !== 1 ? "s" : "") + "…");
            }
        }, true); /* capture phase — fires before React's synthetic handler */
    }

    /* ── Initialise ────────────────────────────────────────────────────────── */

    function init() {
        watchForBrowseButton();
        watchModalForManager();
        setupDragDrop();
        setupUploadInputWatch();
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();
