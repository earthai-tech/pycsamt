# pycsamt.inversion Roadmap

This file is the working checklist for making `pycsamt.inversion` a robust
physics-based EM inversion package. Move to the next task only after the
current task has code, tests, and a clear user-facing API.

## Milestone 1 - Real Optional Backends

- [x] True SimPEG backend
  - Natural-source MT/AMT/CSAMT 1-D inversion API wired.
  - Stitched profile mode for 2-D sections.
  - Natural-source MT/AMT/CSAMT 3-D inversion API wired through
    `Simulation3DPrimarySecondary`.
  - Missing-dependency behavior tested.
  - Still needs live validation in an environment with SimPEG installed.

- [x] True pyGIMLi backend
  - Implemented 1-D MT/AMT/CSAMT through `pygimli.physics.em.MT1d*Modelling`.
  - Implemented 1-D TDEM through `pygimli.physics.em.TDEM*Modelling`.
  - Added stitched 2-D profile mode using station-by-station 1-D inversions.
  - Added missing-dependency and observation-packing tests.
  - Live-validate in an environment with pyGIMLi installed.

## Milestone 2 - Built-In Physics

- [ ] True 2-D physics inversion inside `builtin`
  - Replace "stitched 1-D" with a real 2-D EM forward/inverse path where feasible.
  - Keep stitched mode as a fast screening option.

- [ ] Mesh-building abstraction used by real physics engines
  - Centralize 1-D/2-D/3-D mesh configuration.
  - Support depth padding, lateral padding, topography hooks, and active cells.

- [ ] Regularization actually wired into all solvers
  - Smooth, damped, blocky options.
  - Reference model support.
  - Consistent parameters across built-in, SimPEG, pyGIMLi, Occam2D, ModEM.

## Milestone 3 - External Engines

- [ ] True 3-D inversion execution
  - Robust ModEM runner lifecycle.
  - SimPEG 3-D natural-source path where practical.
  - Clear compute/runtime warnings.

- [ ] Real Occam2D/ModEM runner lifecycle validation
  - Prepare, run, monitor, load, validate outputs.
  - Detect missing executables and failed iterations cleanly.

- [ ] Real backend compatibility tests
  - Test missing dependencies.
  - Test installed-backend smoke runs behind optional markers.
  - Test result conversion for each backend.

## Milestone 4 - Data and Errors

- [ ] Stronger data readers
  - `Sites`
  - EDI collections
  - TDEM survey objects
  - Existing `pycsamt.tdem` workflows

- [ ] Better error models per component
  - Rho/phase floors.
  - Impedance real/imag errors.
  - TDEM voltage or dB/dt relative/absolute floors.
  - Station/component masks.

- [ ] Robust uncertainty outputs
  - Posterior/proxy uncertainty where available.
  - Sensitivity/Jacobian diagnostics.
  - Station and depth confidence maps.

## Milestone 5 - Products

- [ ] Convergence histories
  - Iteration, beta/lambda, phi_d, phi_m, RMS/chi2.
  - Common output structure across backends.

- [ ] Richer exports
  - VTK
  - GeoJSON
  - GeoTIFF where raster grids are available
  - Backend-native archive snapshots

- [ ] Better documentation/examples
  - Minimal 1-D MT example.
  - Minimal 1-D TDEM example.
  - Stitched 2-D profile example.
  - Occam2D and ModEM adapter examples.
  - Hydrogeophysical interpretation after inversion.
