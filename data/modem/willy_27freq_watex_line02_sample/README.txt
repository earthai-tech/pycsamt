Willy 27-frequency WATEX ModEM sample
=====================================

This folder contains a compact subset of one ModEM/NLCG inversion-result
directory from the Willy project.  It is intended for documentation-gallery
examples and lightweight tests that need real ModEM artefacts without bundling
the complete inversion run.

Source directory
----------------

The files were selected from a local result folder:

    C:\Users\Daniel\Downloads\willy-inversion results\
    27-frequ-watex-data-02\27-frequ-watex-data-02

Only a representative line/run subset was copied.  The full source directory
contains many iteration snapshots and large solver logs; those are deliberately
not vendored here.

Included files
--------------

Run controls and metadata
~~~~~~~~~~~~~~~~~~~~~~~~~

``inv.ctrl``
    ModEM inversion-control file.

``fwd.ctrl``
    Forward-control file used by the run.

``run.slurm``
    Original batch-launch note.  Useful for documenting how the production run
    was submitted, but gallery examples should not execute it.

``fort.2000``
    Small auxiliary solver output.

``CSUr2.err``
    Error stream captured during the run.

``Modular_NLCG.log``
    Main nonlinear conjugate-gradient convergence log.  This is the primary
    file for convergence-monitoring examples.

Initial/input artefacts
~~~~~~~~~~~~~~~~~~~~~~~

``27-freq-run-watex01.cov``
    Covariance/smoothing file.

``27-freq-run-watex01.dat``
    Input data file used by the inversion.

``27-freq-run-watex01.rho``
    Starting/reference resistivity model.

Representative iteration snapshots
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Modular_NLCG_000.dat`` and ``Modular_NLCG_000.rho``
    Initial predicted data/model snapshot.

``Modular_NLCG_030.dat``, ``Modular_NLCG_030.res`` and
``Modular_NLCG_030.rho``
    Mid-run response/residual/model snapshot.

``Modular_NLCG_073.dat``, ``Modular_NLCG_073.res`` and
``Modular_NLCG_073.rho``
    Late/final response/residual/model snapshot copied for result-loading and
    model-comparison examples.

Notes for examples
------------------

Use this folder for read-only documentation workflows:

* parse ``Modular_NLCG.log`` for RMS/objective/model-norm convergence;
* load representative ModEM result files;
* compare early, middle, and late iteration models/responses;
* demonstrate checks for missing, zero-byte, or incomplete solver artefacts.

Do not assume this folder is a complete restartable inversion directory.  It is
a compact teaching dataset, not the full production run.
