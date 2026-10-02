# Packaging pyCSAMT's frozen apps

[PyInstaller](https://pyinstaller.org/) specs live here for two of pyCSAMT's
apps:

- `pycsamt_converter.spec` — the standalone format converter app
  (`pycsamt.app.converter`, "pyCSAMT Format Studio"). Small, narrow import
  graph, ~600 MB.
- `pycsamt_desktop.spec` — the full `pycsamt-desktop` suite
  (`pycsamt.app.desktop`). Reaches nearly everything in the library,
  including torch/tensorflow-backed AI features — substantially larger,
  see that spec's own docstring for the full accounting of what it reaches
  and why it can't be trimmed the way the converter's is.

Neither app *needs* a binary to run: `pycsamt-converter` / `pycsamt-desktop`
(or `python -m pycsamt.app.converter` / `python -m pycsamt.app.desktop`)
work directly once `pip install 'pycsamt[app]'` (or `[desktop,agents]` for
every desktop feature) is done. Packaging exists for handing a
double-click-able build to someone who doesn't want a Python environment at
all — and, for the desktop app, for wrapping that build into a real
installer (Windows Setup.exe / Linux self-extracting script — see
`packaging/inno/` and `packaging/linux/`).

## Build: the converter app

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

### Size and what's bundled (converter)

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

## Build: the full desktop suite

1. Use (or create) a Python environment with `pycsamt` installed
   together with **every extra this build reaches** —
   `pip install -e '.[desktop,agents,agent-master,torch,geo,patterns]'`
   at minimum, `+tensorflow` too if you want the CAE Anomaly
   Detection / Denoising agents to actually work in the frozen build
   (see "Size and what's bundled" below) — plus `pyinstaller`.
2. From that environment:

   - Windows: `powershell -ExecutionPolicy Bypass -File packaging\pyinstaller\build_desktop.ps1`
   - Linux/macOS: `bash packaging/pyinstaller/build_desktop.sh`

   Output lands in `dist/pycsamt-desktop/` (onedir, same rationale as
   the converter's build).
3. To wrap that onedir build into a real installer instead of handing
   someone a bare folder:

   - Windows: `powershell -ExecutionPolicy Bypass -File packaging\inno\build_installer.ps1`
     (requires [Inno Setup 6](https://jrsoftware.org/isinfo.php)) →
     `dist/installer/pycsamt-desktop-setup-<version>.exe` with Start
     Menu / Desktop shortcuts and a real uninstaller.
   - Linux: `bash packaging/linux/build_installer.sh` →
     `dist/installer/pycsamt-desktop-<version>-linux-x86_64.sh`, a
     self-extracting installer (installs to
     `~/.local/share/pycsamt-desktop`, adds a launcher on `PATH` and a
     `.desktop` entry; `--uninstall` removes them). AppImage remains a
     stretch goal, not built here.
   - A manual (`workflow_dispatch`-only) GitHub Action,
     `.github/workflows/build-desktop.yml`, runs the Windows onedir +
     installer build in CI on demand — no CI currently builds either
     frozen app automatically on every push (true even for the
     converter), and this build is large/slow enough that it should
     stay an intentional, on-demand action.

### Size and what's bundled (desktop)

Unlike the converter, this build **cannot exclude torch/tensorflow** —
both are genuinely reachable via deep, function-local imports throughout
`pycsamt.ai`/`pycsamt.agents` (confirmed by grep; see
`pycsamt_desktop.spec`'s own docstring for the exact call sites). Excluding
either would silently break real, already-shipped features the moment a
user clicks them:

- **torch** backs AI Inversion (2-D/3-D) and most of the AI-processing
  agent suite (Imputer, Uncertainty Calibration, Distortion
  Classification, Time-Series Denoising).
- **tensorflow** specifically backs the CAE-based Anomaly Detection and
  Denoising agents (`pycsamt.ai.processing.anomaly`/`denoise`'s
  `tensorflow.keras` path) — a genuinely separate dependency the
  project's own `full` extra doesn't even include (it "prefers torch
  backend"). If the build environment doesn't have `tensorflow`
  installed, PyInstaller warns and skips it rather than failing the
  whole build — the frozen app still runs, just with those two specific
  agents raising a clear `ImportError` instead of working, rather than
  silently doing nothing.

`pycsamt.app.converter` (embedded in-process — Phase 9's "Open in Format
Studio…") and `pycsamt.app.agent_master` (launched as a **subprocess**, not
in-process — see `packaging/pyinstaller/entry_desktop.py` and
`pycsamt/app/desktop/agent_master_bridge.py` for the frozen-build-specific
re-exec this needed) are both bundled into the same build for exactly that
reason. `pycsamt.app.web` and `pycsamt.app.mapview` are confirmed
unreachable (no import, no launcher, no subprocess call anywhere under
`pycsamt/app/desktop/`) and excluded.

Expect several hundred MB to a few GB, well above the converter's ~600 MB
— that's torch (and tensorflow, if included) on top of the same
NumPy/SciPy/Matplotlib/Qt stack the converter already carries.

## Smoke-testing a build

Before handing out a build, at minimum:

1. Launch the frozen exe on a clean machine (or the same machine, but
   note it proves nothing about missing DLLs if built and run on the
   same box).
2. **Converter:** convert one real file on each page — an inversion
   result folder or `.npz` on the Inversion page, a `.pcsf`↔`.pcsm`
   transcode, a real `.edi` on the EDI↔XML page, a CSV on the
   PCBH/PCGL/PCGS/PCPT builder pages, and a mixed batch queue.
   **Desktop:** open Profile/Map/QC once each against real bundled
   data, run one AI-processing agent that needs torch (e.g. Uncertainty
   Calibration) and, if tensorflow was included, one that needs it
   (Anomaly Detection), launch "Open in Format Studio…" (embedded
   converter) and "Agents" (Agent Master subprocess — confirms the
   frozen-build re-exec path in `entry_desktop.py` actually works, not
   just the desktop UI itself), and confirm Preferences ▸ License and
   the trial-expiry gate (Phase 11/12) render.
3. Watch for a `ModuleNotFoundError` in a console — rebuild with
   `console=True` in the `.spec`'s `EXE(...)` call temporarily if a
   page silently does nothing, since the windowed
   (`console=False`) build swallows uncaught tracebacks.

`pycsamt/app/converter/jobs.py`'s own test suite
(`pytest pycsamt/app/converter/tests/test_jobs.py`) and the desktop app's
own extensive suite (`pytest pycsamt/app/desktop/tests`) already cover
correctness in-process; the frozen-build smoke test above is only checking
that PyInstaller's dependency bundling didn't drop something the tests
don't exercise (import resolution, Qt plugin loading, data files, the
Agent Master subprocess re-exec).

## Known caveats

- **Qt plugin path.** PySide6's platform plugins (`qwindows.dll`,
  etc.) ship inside `_internal/PySide6/plugins/`; PyInstaller's
  PySide6 hook wires this up automatically, so this should never need
  manual `QT_PLUGIN_PATH` configuration — if a frozen build opens a
  blank/no window, that hook not finding the plugins is the first
  thing to check. The desktop build additionally uses
  `QtWebEngineWidgets` (`pycsamt/app/desktop/widgets/plotly_view.py`'s
  `PlotlyView`, Phase 2) — PyInstaller's PySide6 hook bundles
  `QtWebEngineProcess.exe` and its resources/locales automatically too,
  but this is a heavier, less commonly exercised hook path than plain
  `QtWidgets`; if the PCSF 3-D viewer opens blank, check for a missing
  `QtWebEngineProcess.exe` or `resources/` folder next to it first.
- **Unsigned-executable warnings.** The build is not code-signed, so
  Windows SmartScreen / macOS Gatekeeper will warn on first launch.
  Signing is out of scope here (an explicit open decision in
  `PYCSAMT-DESKTOP-V2.6-MODERNIZATION-PLAN.md` Section 4); mention this
  to whoever receives the build.
- **First launch is slow.** A onedir PySide6 + NumPy/SciPy (+ torch,
  for the desktop build) build takes a few seconds to start even after
  warm-up, longer on the very first launch while the OS caches the
  DLLs — this is normal, not a hang.
- **Rebuild after any dependency bump.** If `pycsamt`'s own
  dependencies change (a new format/emtf submodule, a new third-party
  package), rebuild rather than assume the old bundle still covers it
  -- there is no CI job producing either binary automatically except
  the desktop app's manual `workflow_dispatch` action.
- **Trial/license state survives reinstall on purpose.** Neither the
  Inno Setup uninstaller nor the Linux installer's `--uninstall` touch
  `~/.pycsamt` or the `TrialTracker`/`OfflineLicenseManager` QSettings
  state (registry on Windows, `~/.config/earthai-tech/` on Linux) — see
  each installer's own comment. Deleting that on uninstall would let a
  reinstall silently reset a user's trial clock or drop their activated
  license key.
