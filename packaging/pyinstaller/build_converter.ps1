# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build the pyCSAMT Format Studio standalone binary (Windows).
#
# Usage (from anywhere):
#   powershell -ExecutionPolicy Bypass -File packaging\pyinstaller\build_converter.ps1
#
# Requires: the environment that has pycsamt + PySide6 installed also
# has `pyinstaller` (pip install pyinstaller). See
# packaging/pyinstaller/README.md for details and troubleshooting.

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$SpecFile = Join-Path $PSScriptRoot "pycsamt_converter.spec"

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
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue (Join-Path $RepoRoot "build\pycsamt-converter")
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue (Join-Path $RepoRoot "dist\pycsamt-converter")

    Write-Host "Running PyInstaller..."
    pyinstaller --noconfirm --clean $SpecFile

    $ExePath = Join-Path $RepoRoot "dist\pycsamt-converter\pycsamt-converter.exe"
    if (Test-Path $ExePath)
    {
        Write-Host ""
        Write-Host "Build succeeded: $ExePath"
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
