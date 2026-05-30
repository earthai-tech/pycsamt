#!/usr/bin/env bash
# install_sources.sh — refresh the bundled ModEM Fortran sources
#
# The ModEM v6.2.6 sources are already committed under _source/2D/ and
# _source/3D/.  You only need this script when you want to update them
# from a newer ModEM release.
#
# Usage:
#   bash install_sources.sh                   # auto-detect ModEMv626/ at repo root
#   bash install_sources.sh /path/to/ModEM    # explicit path to ModEM root
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Walk up from this script to the repository root (5 levels: _source → modem → models → pycsamt → pycsamt → repo)
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../../../.." && pwd)"

if [[ $# -ge 1 ]]; then
    MODEM_ROOT="$1"
else
    MODEM_ROOT="${REPO_ROOT}/ModEMv626/ModEM"
fi

F90="${MODEM_ROOT}/f90"

if [[ ! -d "${F90}" ]]; then
    echo "ERROR: ModEM f90 source directory not found: ${F90}"
    echo ""
    echo "  The ModEM v6.2.6 sources are already bundled — you do NOT need"
    echo "  to run this script unless upgrading to a newer ModEM release."
    echo ""
    echo "  To upgrade, provide the path to your new ModEM distribution:"
    echo "    bash install_sources.sh /path/to/ModEM"
    exit 1
fi

echo "Using ModEM source: ${F90}"

# ---- 2D ----
DIR2D="${SCRIPT_DIR}/2D"
echo "Refreshing ${DIR2D} ..."
for sub in UTILS INV SENS LAPACK MPI; do
    rm -rf "${DIR2D:?}/${sub}"
    cp -r  "${F90}/${sub}" "${DIR2D}/"
done
rm -rf "${DIR2D:?}/2D_MT"
cp -r "${F90}/2D_MT"        "${DIR2D}/"
cp    "${F90}/UserCtrl.f90" "${DIR2D}/"
cp    "${F90}/Mod2DMT.f90"  "${DIR2D}/"
echo "  2D sources refreshed."

# ---- 3D ----
DIR3D="${SCRIPT_DIR}/3D"
echo "Refreshing ${DIR3D} ..."
for sub in UTILS INV SENS LAPACK MPI FIELDS; do
    rm -rf "${DIR3D:?}/${sub}"
    cp -r  "${F90}/${sub}" "${DIR3D}/"
done
rm -rf "${DIR3D:?}/3D_MT"
cp -r "${F90}/3D_MT"        "${DIR3D}/"
cp    "${F90}/UserCtrl.f90" "${DIR3D}/"
cp    "${F90}/Mod3DMT.f90"  "${DIR3D}/"
echo "  3D sources refreshed."

echo ""
echo "Done.  Now compile:"
echo "  cd ${DIR2D} && make"
echo "  cd ${DIR3D} && make        # serial"
echo "  cd ${DIR3D} && make mpi    # MPI-parallel (recommended for 3D)"
