# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build the pycsamt-desktop standalone binary (Windows).
#
# Usage (from anywhere):
#   powershell -ExecutionPolicy Bypass -File packaging\pyinstaller\build_desktop.ps1
#
# Requires: the environment that has pycsamt installed together with its
# `desktop` AND `agents` extras (PySide6, pyqtgraph, contextily, torch,
# scikit-learn, the LLM provider clients) -- this build cannot exclude
# torch/tensorflow the way the converter build does, see
# pycsamt_desktop.spec's own docstring -- plus `pyinstaller` itself
# (pip install pyinstaller). See packaging/pyinstaller/README.md for
# details and troubleshooting.

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$SpecFile = Join-Path $PSScriptRoot "pycsamt_desktop.spec"

Push-Location $RepoRoot
try
{
    python -c "import PyInstaller" 2>$null
    if ($LASTEXITCODE -ne 0)
    {
        Write-Error "PyInstaller is not installed in this Python environment. Run: pip install pyinstaller"
        exit 1
    }

    Write-Host "Cleaning previous build artifacts..."
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue (Join-Path $RepoRoot "build\pycsamt-desktop")
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue (Join-Path $RepoRoot "dist\pycsamt-desktop")

    Write-Host "Running PyInstaller (this build is substantially larger than the converter's -- expect several minutes)..."
    pyinstaller --noconfirm --clean $SpecFile

    $ExePath = Join-Path $RepoRoot "dist\pycsamt-desktop\pycsamt-desktop.exe"
    if (Test-Path $ExePath)
    {
        Write-Host ""
        Write-Host "Build succeeded: $ExePath"
        Write-Host "Next: wrap it with Inno Setup (packaging/inno/pycsamt_desktop.iss) for a real installer."
    }
    else
    {
        Write-Error "Build finished but $ExePath was not found."
        exit 1
    }
}
finally
{
    Pop-Location
}
