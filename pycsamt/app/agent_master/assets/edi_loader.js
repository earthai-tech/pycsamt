/* edi_loader.js
 * Folder browse + drag-and-drop for Agent Master.
 * Adapted from pycsamt/app/web/assets/folder_loader.js.
 *
 * Browse Folder: transparent <input webkitdirectory>
 * injected inside #am-btn-browse so a click lands
 * directly on the native file input — no programmatic
 * .click() needed.
 *
 * Results are pushed to Dash store "am-folder-store"
 * via window.dash_clientside.set_props.
 */
(function () {
    "use strict";

    var EXTS = [".edi"];

    function hasExt(name) {
        var lo = name.toLowerCase();
        return EXTS.some(function (e) {
            return lo.endsWith(e);
        });
    }

    function readFileAsDataURL(file) {
        return new Promise(function (resolve, reject) {
            var r = new FileReader();
            r.onload  = function (ev) {
                resolve(ev.target.result);
            };
            r.onerror = function () {
                reject(new Error("FileReader error"));
            };
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
                    if (!entries.length) {
                        resolve(all);
                        return;
                    }
                    all = all.concat(Array.from(entries));
                    batch();
                }, reject);
            }
            batch();
        });
    }

    async function collectEntryFiles(entry, basePath) {
        var relPath = basePath
            ? basePath + "/" + entry.name
            : entry.name;
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
                var sub = await collectEntryFiles(
                    entries[i], relPath
                );
                results = results.concat(sub);
            }
            return results;
        }
        return [];
    }

    /* ── Status helpers ───────────────────────────── */

    function setStatus(msg) {
        var el = document.getElementById(
            "am-edi-load-status"
        );
        if (el) el.textContent = msg;
    }

    /* ── Push to Dash store ───────────────────────── */

    async function sendToStore(pairs) {
        if (!pairs.length) {
            setStatus("No EDI files found.");
            return;
        }
        var n = pairs.length;
        setStatus("Reading " + n + " file(s)...");

        var settled = await Promise.all(
            pairs.map(function (pair) {
                return readFileAsDataURL(pair.file)
                    .then(
                        function (b64) {
                            return {
                                ok: true,
                                path: pair.relativePath,
                                b64: b64
                            };
                        },
                        function () {
                            return {
                                ok: false,
                                path: pair.relativePath
                            };
                        }
                    );
            })
        );

        var filenames = [];
        var contents  = [];
        settled.forEach(function (r) {
            if (r.ok) {
                filenames.push(r.path);
                contents.push(r.b64);
            }
        });

        setStatus(
            filenames.length + " file(s) read — loading..."
        );

        if (
            window.dash_clientside
            && window.dash_clientside.set_props
        ) {
            window.dash_clientside.set_props(
                "am-folder-store",
                { data: {
                    filenames: filenames,
                    contents: contents
                }}
            );
        }
    }

    /* ── Browse-folder button injection ───────────── */

    function makeFolderInput() {
        var inp = document.createElement("input");
        inp.type = "file";
        inp.multiple = true;
        inp.setAttribute("webkitdirectory", "");
        inp.setAttribute("directory", "");
        inp.style.cssText = [
            "position:absolute",
            "inset:0",
            "width:100%",
            "height:100%",
            "opacity:0",
            "cursor:pointer",
            "font-size:0"
        ].join(";");

        inp.addEventListener(
            "change",
            async function () {
                var files = Array.from(
                    inp.files || []
                ).filter(function (f) {
                    return hasExt(f.name);
                });
                inp.value = "";

                if (!files.length) {
                    setStatus("No EDI files found.");
                    return;
                }

                var pairs = files.map(function (f) {
                    return {
                        relativePath:
                            f.webkitRelativePath
                            || f.name,
                        file: f
                    };
                });
                await sendToStore(pairs);
            }
        );
        return inp;
    }

    function injectBrowseInput() {
        var btn = document.getElementById(
            "am-btn-browse"
        );
        if (!btn) return false;
        if (btn.dataset.folderInjected) return true;

        btn.dataset.folderInjected = "1";
        btn.style.position = "relative";
        btn.style.overflow = "hidden";
        btn.appendChild(makeFolderInput());
        return true;
    }

    function watchForBrowseButton() {
        injectBrowseInput();
        var obs = new MutationObserver(function () {
            var btn = document.getElementById(
                "am-btn-browse"
            );
            if (!btn || btn.dataset.folderInjected) {
                return;
            }
            injectBrowseInput();
        });
        obs.observe(document.body, {
            childList: true,
            subtree: true
        });
    }

    /* ── Drag-and-drop folder support ─────────────── */

    function setupDragDrop() {
        document.addEventListener(
            "dragover",
            function (e) {
                var zone = e.target.closest(
                    ".am-upload-zone"
                );
                if (!zone) return;
                var items =
                    e.dataTransfer && e.dataTransfer.items;
                if (!items) return;
                var hasDir = false;
                for (var i = 0; i < items.length; i++) {
                    if (items[i].kind !== "file") continue;
                    var entry =
                        items[i].webkitGetAsEntry
                        && items[i].webkitGetAsEntry();
                    if (entry && entry.isDirectory) {
                        hasDir = true;
                        break;
                    }
                }
                if (hasDir) {
                    e.preventDefault();
                    e.stopPropagation();
                    zone.style.borderColor = "#89b4fa";
                }
            },
            true
        );

        document.addEventListener(
            "dragleave",
            function (e) {
                var zone = e.target.closest(
                    ".am-upload-zone"
                );
                if (zone) zone.style.borderColor = "";
            },
            true
        );

        document.addEventListener(
            "drop",
            function (e) {
                var zone = e.target.closest(
                    ".am-upload-zone"
                );
                if (!zone) return;
                var items =
                    e.dataTransfer && e.dataTransfer.items;
                if (!items) return;

                var entries = [];
                for (var i = 0; i < items.length; i++) {
                    if (items[i].kind !== "file") continue;
                    var entry =
                        items[i].webkitGetAsEntry
                        && items[i].webkitGetAsEntry();
                    if (entry) entries.push(entry);
                }

                var hasDir = entries.some(function (en) {
                    return en.isDirectory;
                });
                if (!hasDir) return; /* dcc.Upload handles */

                e.preventDefault();
                e.stopPropagation();
                zone.style.borderColor = "";

                (async function () {
                    var allPairs = [];
                    for (var j = 0; j < entries.length; j++) {
                        var pairs = await collectEntryFiles(
                            entries[j], ""
                        );
                        pairs = pairs.map(function (p) {
                            return {
                                relativePath:
                                    p.relativePath
                                    .replace(/^\//, ""),
                                file: p.file
                            };
                        });
                        allPairs = allPairs.concat(pairs);
                    }
                    await sendToStore(allPairs);
                })();
            },
            true
        );
    }

    /* ── Init ─────────────────────────────────────── */

    function init() {
        watchForBrowseButton();
        setupDragDrop();
    }

    if (document.readyState === "loading") {
        document.addEventListener(
            "DOMContentLoaded", init
        );
    } else {
        init();
    }
})();
