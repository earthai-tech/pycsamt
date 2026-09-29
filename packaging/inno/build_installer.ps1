# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
#
# Build the pycsamt-desktop Windows installer (Setup.exe) from the
# PyInstaller onedir build, via Inno Setup.
#
# Usage (from anywhere):
#   powershell -ExecutionPolicy Bypass -File packaging\inno\build_installer.ps1
#
# Requires:
#   - Inno Setup 6 installed (https://jrsoftware.org/isinfo.php), with
#     ISCC.exe on PATH or at the default install location this script
#     checks below.
#   - dist\pycsamt-desktop\ already built (this script builds it first via
#     build_desktop.ps1 if it's missing -- but does NOT rebuild it if it's
#     already there, since that build takes several minutes; delete
#     dist\pycsamt-desktop\ yourself first to force a rebuild).

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$DistDir = Join-Path $RepoRoot "dist\pycsamt-desktop"
$IssFile = Join-Path $PSScriptRoot "pycsamt_desktop.iss"
$BuildDesktopScript = Join-Path $RepoRoot "packaging\pyinstaller\build_desktop.ps1"

Push-Location $RepoRoot
try
{
    if (-not (Test-Path (Join-Path $DistDir "pycsamt-desktop.exe")))
    {
        Write-Host "dist\pycsamt-desktop\pycsamt-desktop.exe not found -- building it first..."
        # -NoProfile: a user profile running `conda init`'s hook re-activates
        # conda's *base* env in the child shell, silently swapping out the
        # build env this script was launched from (base has no PyInstaller).
        & powershell -NoProfile -ExecutionPolicy Bypass -File $BuildDesktopScript
        if ($LASTEXITCODE -ne 0 -or -not (Test-Path (Join-Path $DistDir "pycsamt-desktop.exe")))
        {
            # Stop here rather than letting ISCC fail later with a misleading
            # "No files found matching dist\pycsamt-desktop\*".
            Write-Error "build_desktop.ps1 failed (exit $LASTEXITCODE) -- not running Inno Setup."
            exit 1
        }
    }

    $Iscc = Get-Command "ISCC.exe" -ErrorAction SilentlyContinue
    if ($Iscc)
    {
        $IsccPath = $Iscc.Source
    }
    else
    {
        # Two install locations in practice: the machine-wide default
        # (Program Files (x86) -- what `choco install innosetup` and the CI
        # workflow use) and winget's default non-elevated, per-user install
        # (%LocalAppData%\Programs) -- check both before giving up.
        $Candidates = @(
            "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
            "${env:LocalAppData}\Programs\Inno Setup 6\ISCC.exe"
        )
        $IsccPath = $Candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
        if (-not $IsccPath)
        {
            Write-Error "ISCC.exe (Inno Setup 6) not found on PATH or at either default install location. Install it from https://jrsoftware.org/isinfo.php"
            exit 1
        }
    }

    $Version = python -c "import pycsamt; print(pycsamt.__version__)"
    if ($LASTEXITCODE -ne 0 -or -not $Version)
    {
        Write-Error "Could not read pycsamt.__version__ from the active Python environment."
        exit 1
    }
    Write-Host "Building installer for version $Version..."

    & $IsccPath "/DMyAppVersion=$Version" $IssFile

    $SetupExe = Join-Path $RepoRoot "dist\installer\pycsamt-desktop-setup-$Version.exe"
    if (Test-Path $SetupExe)
    {
        Write-Host ""
        Write-Host "Installer built: $SetupExe"
    }
    else
    {
        Write-Error "ISCC finished but $SetupExe was not found."
        exit 1
    }
}
finally
{
    Pop-Location
}
