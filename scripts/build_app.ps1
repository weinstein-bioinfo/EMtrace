# Build the EMtrace standalone Windows bundle using PyInstaller.
#
# Usage (PowerShell):
#   .\scripts\build_app.ps1
#
# Output:
#   dist\EMtrace-x86_64\EMtrace.exe
#   dist\EMtrace-win-x86_64-{version}.zip
#
# Prerequisites:
#   - Python 3.10-3.12 installed and on PATH (3.12 recommended: widest wheel
#     coverage for PyQt6 / pysam / kaleido)
#   - Git (to clone the repo on a build machine)
#
# The venv is created and populated automatically.

$ErrorActionPreference = "Stop"
Push-Location (Split-Path -Parent $PSScriptRoot)

Write-Host "=== EMtrace Windows App Builder ==="
Write-Host ""

# ------------------------------------------------------------------
# 1. Check Python
# ------------------------------------------------------------------
try {
    $pyVer = python --version 2>&1
    Write-Host "Python: $pyVer"
} catch {
    Write-Host "ERROR: Python not found on PATH. Install Python 3.10-3.12 and retry."
    exit 1
}

# ------------------------------------------------------------------
# 2. Create / activate venv
# ------------------------------------------------------------------
$venvDir = ".venv"
if (-not (Test-Path "$venvDir\Scripts\python.exe")) {
    Write-Host "Creating virtual environment in $venvDir ..."
    python -m venv $venvDir
}
Write-Host "Activating venv ..."
& "$venvDir\Scripts\Activate.ps1"

# ------------------------------------------------------------------
# 3. Install dependencies
# ------------------------------------------------------------------
Write-Host ""
Write-Host "Installing / updating dependencies ..."
pip install --upgrade pip --quiet

# pysam has no Windows wheels and does not build on MSVC. The pipeline only
# imports it inside the coverage-plot helper (steps/plots.py), which degrades
# gracefully, so the Windows bundle ships without it.
$constraints = "scripts\win_constraints.txt"
pip install -e ".[gui,excel]" --quiet --constraint $constraints
pip install matplotlib plotly kaleido seaborn PyComplexHeatmap --quiet --constraint $constraints
pip install pyinstaller pillow --quiet
Write-Host "Dependencies ready."

# ------------------------------------------------------------------
# 4. Ensure the app icon exists
# ------------------------------------------------------------------
$iconIco = "src\emtrace\gui\resources\icons\app.ico"
if (-not (Test-Path $iconIco)) {
    Write-Host ""
    Write-Host "App icon not found - generating ..."
    python scripts\generate_icon.py
}

# ------------------------------------------------------------------
# 5. Build
# ------------------------------------------------------------------
Write-Host ""
Write-Host "Building EMtrace with PyInstaller ..."
Write-Host ""
pyinstaller EMtrace.spec --clean --noconfirm

# ------------------------------------------------------------------
# 6. Verify output
# ------------------------------------------------------------------
Write-Host ""
Write-Host "=== Build complete ==="

# Read the version from pyproject.toml
$version = (Select-String -Path "pyproject.toml" -Pattern '^version\s*=\s*"(.+)"' |
    ForEach-Object { $_.Matches[0].Groups[1].Value })

$arch = "x86_64"
$exePath = "dist\EMtrace-$arch\EMtrace.exe"
$zipPath = "dist\EMtrace-win-$arch-$version.zip"

if (Test-Path $exePath) {
    $exeSize = (Get-Item $exePath).Length / 1MB
    Write-Host ("  EXE: {0} ({1:N1} MB)" -f $exePath, $exeSize)
} else {
    Write-Host "ERROR: $exePath was not created."
    exit 1
}

if (Test-Path $zipPath) {
    $zipSize = (Get-Item $zipPath).Length / 1MB
    Write-Host ("  ZIP: {0} ({1:N1} MB)" -f $zipPath, $zipSize)
}

Write-Host ""
Write-Host "To run:         .\$exePath"
Write-Host "CLI:            .\$exePath --help"
Write-Host "To distribute:  share $zipPath"

Pop-Location
