/* edi_loader.js — folder browse + folder drag-drop for the Map View loader.
 * Adapted from pycsamt/app/web/assets/folder_loader.js.
 *
 * A transparent <input webkitdirectory> is injected inside #mv-btn-browse so
 * the click lands on the native folder picker.  All files are read in parallel
 * (Promise.all) and pushed to the Dash store "mv-folder-store"; the Python
 * callback then decodes + groups them into survey lines.
 */
(function () {
    "use strict";

    var EXTS = [".edi", ".avg", ".j"];

    function hasExt(name) {
        var lo = name.toLowerCase();
        return EXTS.some(function (e) { return lo.endsWith(e); });
    }

    function readFileAsDataURL(file) {
        return new Promise(function (resolve, reject) {
            var r = new FileReader();
            r.onload = function (ev) { resolve(ev.target.result); };
            r.onerror = function () { reject(new Error("read error")); };
            r.readAsDataURL(file);
        });
    }

    function entryToFile(entry) {
        return new Promise(function (resolve, reject) {
            entry.file(resolve, reject);
        });
    }

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

    async function collectEntryFiles(entry, basePath) {
        var rel = basePath ? basePath + "/" + entry.name : entry.name;
        if (entry.isFile) {
            if (!hasExt(entry.name)) return [];
            var file = await entryToFile(entry);
            return [{ relativePath: rel, file: file }];
        }
        if (entry.isDirectory) {
            var entries = await readAllEntries(entry.createReader());
            var out = [];
            for (var i = 0; i < entries.length; i++) {
                out = out.concat(await collectEntryFiles(entries[i], rel));
            }
            return out;
        }
        return [];
    }

    /* ── spinner / status ────────────────────────────────────────────── */

    function showSpinner(msg) {
        var ov = document.getElementById("mv-loader-overlay");
        var lbl = document.getElementById("mv-loader-msg");
        if (lbl) lbl.textContent = msg || "Reading files…";
        if (ov) ov.style.display = "flex";
    }
    function hideSpinner() {
        var ov = document.getElementById("mv-loader-overlay");
        if (ov) ov.style.display = "none";
    }
    function setStatus(msg) {
        var el = document.getElementById("mv-browse-status");
        if (el) el.textContent = msg || "";
    }
    function setCount(n) {
        var el = document.getElementById("mv-file-count");
        if (!el) return;
        el.textContent = n > 0 ? (n + " file" + (n !== 1 ? "s" : "") + " ready") : "";
    }

    /* ── push to Dash store (parallel reads) ─────────────────────────── */

    async function sendToStore(pairs) {
        if (!pairs.length) {
            hideSpinner();
            setStatus("No EDI / AVG / J files found.");
            setCount(0);
            return;
        }
        showSpinner("Reading " + pairs.length + " file"
                    + (pairs.length !== 1 ? "s" : "") + "…");
        var settled = await Promise.all(pairs.map(function (p) {
            return readFileAsDataURL(p.file).then(
                function (b64) { return { ok: true, path: p.relativePath, b64: b64 }; },
                function () { return { ok: false }; }
            );
        }));
        var filenames = [], contents = [];
        settled.forEach(function (r) {
            if (r.ok) { filenames.push(r.path); contents.push(r.b64); }
        });
        setCount(filenames.length);
        setStatus("Sending " + filenames.length + " file(s) to server…");
        if (window.dash_clientside && window.dash_clientside.set_props) {
            window.dash_clientside.set_props("mv-folder-store", {
                data: { filenames: filenames, contents: contents }
            });
        }
        /* hide spinner shortly after; file manager render also covers it */
        setTimeout(hideSpinner, 400);
    }

    /* ── browse-folder input injection ───────────────────────────────── */

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
            var files = Array.from(inp.files || []).filter(function (f) {
                return hasExt(f.name);
            });
            inp.value = "";
            if (!files.length) { setStatus("No EDI / AVG / J files found."); return; }
            showSpinner("Scanning " + files.length + " file"
                        + (files.length !== 1 ? "s" : "") + "…");
            await sendToStore(files.map(function (f) {
                return { relativePath: f.webkitRelativePath || f.name, file: f };
            }));
        });
        return inp;
    }

    function injectBrowseInput() {
        var btn = document.getElementById("mv-btn-browse");
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

    /* ── drag-and-drop folders ───────────────────────────────────────── */

    function setupDragDrop() {
        document.addEventListener("dragover", function (e) {
            var zone = e.target.closest(".mv-upload-drop");
            if (zone) { e.preventDefault(); zone.classList.add("mv-drop-hot"); }
        }, true);
        document.addEventListener("dragleave", function (e) {
            var zone = e.target.closest(".mv-upload-drop");
            if (zone) zone.classList.remove("mv-drop-hot");
        }, true);
        document.addEventListener("drop", function (e) {
            var zone = e.target.closest(".mv-upload-drop");
            if (!zone) return;
            var items = e.dataTransfer && e.dataTransfer.items;
            if (!items) return;
            var entries = [];
            for (var i = 0; i < items.length; i++) {
                if (items[i].kind !== "file") continue;
                var en = items[i].webkitGetAsEntry && items[i].webkitGetAsEntry();
                if (en) entries.push(en);
            }
            var hasDir = entries.some(function (en) { return en.isDirectory; });
            zone.classList.remove("mv-drop-hot");
            if (!hasDir) { return; }  /* plain files → dcc.Upload handles it */
            e.preventDefault();
            e.stopPropagation();
            showSpinner("Scanning dropped folder(s)…");
            (async function () {
                var all = [];
                for (var j = 0; j < entries.length; j++) {
                    var pairs = await collectEntryFiles(entries[j], "");
                    all = all.concat(pairs.map(function (p) {
                        return {
                            relativePath: p.relativePath.replace(/^\//, ""),
                            file: p.file
                        };
                    }));
                }
                await sendToStore(all);
            })();
        }, true);
    }

    function init() {
        watchForBrowseButton();
        setupDragDrop();
    }
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();
