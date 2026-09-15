# Packaging pyCSAMT Format Studio

`pycsamt_converter.spec` builds a standalone binary of the format
converter app (`pycsamt.app.converter`) with
[PyInstaller](https://pyinstaller.org/). This is the only frozen app in
the pyCSAMT repo — the full `pycsamt-desktop` suite, the web app, Map
View, and the agent app are not packaged this way.

The converter app never needs a binary to run: `pycsamt-converter` /
`python -m pycsamt.app.converter` work directly once
`pip install 'pycsamt[app]'` is done. Packaging exists for handing a
double-click-able build to someone who doesn't want a Python
environment at all.

## Build

1. Use (or create) a Python environment that has `pycsamt` installed
   together with its `app` extra (PySide6) — e.g. a conda env with
   `pip install -e '.[app]'` run from the repo root — plus
   `pyinstaller` itself (`pip install pyinstaller`).
2. From that environment:

   - Windows: `powershell -ExecutionPolicy Bypass -File packaging\pyinstaller\build_converter.ps1`
   - Linux/macOS: `bash packaging/pyinstaller/build_converter.sh`

   Both scripts just wrap `pyinstaller --noconfirm --clean
   packaging/pyinstaller/pycsamt_converter.spec` — run that directly if
   you want more control over PyInstaller's own flags.
3. The build is a **onedir** build (not onefile): the app lands in
   `dist/pycsamt-converter/` as `pycsamt-converter.exe` (or
   `pycsamt-converter` on Linux/macOS) plus an `_internal/` folder of
   its dependencies. Onedir starts faster and is easier to debug than
   onefile's self-extracting-to-temp approach; the folder is what you
   zip up and hand to someone, not just the exe.

## Size and what's bundled

The frozen build is roughly 600 MB — that's the scientific Python
stack (NumPy, SciPy, pandas, Matplotlib, h5py, pyproj, openpyxl) plus
Qt/PySide6, not the converter app's own code. `pycsamt_converter.spec`
explicitly `excludes=` PyTorch/TensorFlow and the other pyCSAMT apps
(`pycsamt.ai`, `pycsamt.app.desktop`/`web`/`mapview`/`agent_master`) —
none of those are reachable from `pycsamt.app.converter.jobs` (see that
module's docstring), so leaving them out doesn't remove functionality.
If a future page needs one of them, drop it from the `excludes` list in
the spec rather than fighting a `ModuleNotFoundError` at runtime.

`hiddenimports` collects whole submodule trees (`pycsamt.format`,
`pycsamt.emtf`, `pycsamt.geology`, the three solver-result packages
under `pycsamt.models`, `pycsamt.metadata`, `pycsamt.gis`,
`pycsamt.seg`) rather than hand-listing leaf modules, because several
of those do function-local (`from pycsamt.X import Y` inside a
function body) imports that PyInstaller's static analysis does catch,
but a hand-maintained list would be easy to let drift out of sync as
the format package grows. Their own `.tests` subpackages are filtered
out of that collection (never needed at runtime, and a couple import
`pytest`, which is excluded).

## Smoke-testing a build

Before handing out a build, at minimum:

1. Launch `dist/pycsamt-converter/pycsamt-converter.exe` on a clean
   machine (or the same machine, but note it proves nothing about
   missing DLLs if built and run on the same box).
2. Convert one real file on each page — an inversion result folder or
   `.npz` on the Inversion page, a `.pcsf`↔`.pcsm` transcode, a real
   `.edi` on the EDI↔XML page, a CSV on the PCBH/PCGL/PCGS/PCPT
   builder pages, and a mixed batch queue.
3. Watch for a `ModuleNotFoundError` in a console — rebuild with
   `console=True` in the `.spec`'s `EXE(...)` call temporarily if a
   page silently does nothing, since the windowed
   (`console=False`) build swallows uncaught tracebacks.

`pycsamt/app/converter/jobs.py`'s own test suite
(`pytest pycsamt/app/converter/tests/test_jobs.py`) already covers
every conversion function's correctness in-process; the frozen-build
smoke test above is only checking that PyInstaller's dependency
bundling didn't drop something the tests don't exercise (import
resolution, Qt plugin loading, data files).

## Known caveats

- **Qt plugin path.** PySide6's platform plugins (`qwindows.dll`,
  etc.) ship inside `_internal/PySide6/plugins/`; PyInstaller's
  PySide6 hook wires this up automatically, so this should never need
  manual `QT_PLUGIN_PATH` configuration — if a frozen build opens a
  blank/no window, that hook not finding the plugins is the first
  thing to check.
- **Unsigned-executable warnings.** The build is not code-signed, so
  Windows SmartScreen / macOS Gatekeeper will warn on first launch.
  Signing is out of scope here; mention this to whoever receives the
  build.
- **First launch is slow.** A onedir PySide6 + NumPy/SciPy build takes
  a few seconds to start even after this warm-up, longer on the very
  first launch while the OS caches the DLLs — this is normal, not a
  hang.
- **Rebuild after any dependency bump.** If `pycsamt`'s own
  dependencies change (a new format/emtf submodule, a new third-party
  package), rebuild rather than assume the old bundle still covers it
  -- there is no CI job producing this binary automatically yet.
