param(
  # Path to the `the-omega/docs` folder (a git working tree).
  [string]$TheOmegaDocsPath = (Join-Path $PSScriptRoot "..\\..\\the-omega\\docs"),

  # Destination legacy docs folder inside this repo.
  [string]$DestLegacyDocsPath = (Join-Path $PSScriptRoot "..\legacy-docs"),

  # If set, mirror directories (delete files that no longer exist in source).
  # Default is safer: copy/update only, keep extra files in destination.
  [switch]$Mirror,

  # Include PDFs (disabled by default).
  [switch]$IncludePdf,

  # Only sync source documents + scripts (recommended / default).
  # This keeps `legacy-docs/` lightweight: tex/md/bib + small scripts, no images/data.
  [bool]$SourcesOnly = $true
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Assert-PathExists {
  param([string]$Path, [string]$Label)
  if (-not (Test-Path -LiteralPath $Path)) {
    throw "$Label does not exist: $Path"
  }
}

function Prune-NonSources {
  param([Parameter(Mandatory=$true)][string]$Root)

  if (-not $SourcesOnly) { return }

  function Remove-TreeForce {
    param([Parameter(Mandatory=$true)][string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) { return }
    # PowerShell Remove-Item can be slow / flaky on very deep trees; use cmd rmdir.
    cmd /c rmdir /s /q "$Path" | Out-Null
  }

  $allowedExt = @(
    ".tex", ".md", ".bib", ".txt", ".rst",
    ".py", ".sh", ".ps1", ".bat", ".cmd",
    ".toml", ".ini", ".yaml", ".yml", ".json"
  )
  if ($IncludePdf) { $allowedExt += ".pdf" }

  # Fast path: remove top-level folders we never want to import.
  Remove-TreeForce -Path (Join-Path $Root "books")
  Remove-TreeForce -Path (Join-Path $Root "assets")

  # Remove known data/output directories regardless of file type.
  $dropDirNames = @(
    "data", "datasets", "dist", "artifacts", "outputs", "figures", "images", "img",
    "books", "assets"
  )
  foreach ($name in $dropDirNames) {
    $matches = Get-ChildItem -LiteralPath $Root -Directory -Recurse -Force -ErrorAction SilentlyContinue |
      Where-Object { $_.Name -ieq $name }
    foreach ($d in $matches) {
      Remove-Item -LiteralPath $d.FullName -Recurse -Force -ErrorAction SilentlyContinue
    }
  }

  # Remove non-source files by extension.
  Get-ChildItem -LiteralPath $Root -Recurse -File -Force -ErrorAction SilentlyContinue |
    Where-Object {
      $ext = $_.Extension
      -not ($allowedExt -contains $ext)
    } |
    ForEach-Object {
      Remove-Item -LiteralPath $_.FullName -Force -ErrorAction SilentlyContinue
    }
}

function Run-Robocopy {
  param(
    [Parameter(Mandatory=$true)][string]$Source,
    [Parameter(Mandatory=$true)][string]$Dest
  )

  # /S copies subdirs but skips empty dirs (keeps output tidy when SourcesOnly is on).
  $modeArgs = @("/S")
  if ($Mirror) { $modeArgs = @("/MIR") }

  $fileSpecs = @("*")
  if ($SourcesOnly) {
    # Source docs + scripts only (no figures/assets/data).
    $fileSpecs = @(
      "*.tex",
      "*.md",
      "*.bib",
      "*.txt",
      "*.rst",
      "*.py",
      "*.sh",
      "*.ps1",
      "*.bat",
      "*.cmd",
      "requirements.txt",
      "pyproject.toml"
    )
  }

  # Exclude common build/cache artifacts that should not be migrated.
  $excludeFiles = @(
    "*.aux", "*.bbl", "*.blg", "*.log", "*.out", "*.toc", "*.lof", "*.lot",
    "*.fls", "*.fdb_latexmk", "*.synctex.gz",
    "*.nav", "*.snm", "*.vrb",
    "*.pyc",
    "*.pkl", "*.pkl.gz", "*.pickle",
    "*.npy", "*.npz",
    "*.pt", "*.pth", "*.ckpt"
  )
  if (-not $IncludePdf) {
    $excludeFiles += "*.pdf"
  }
  if ($SourcesOnly) {
    # Non-source assets / binaries that we don't want in this repo.
    $excludeFiles += @(
      "*.png", "*.jpg", "*.jpeg", "*.gif", "*.webp", "*.svg",
      "*.epub",
      "*.zip", "*.7z", "*.rar",
      "*.docx", "*.pptx",
      "*.mp4", "*.mov"
    )
  }
  $excludeDirs = @(
    "__pycache__", ".ipynb_checkpoints", ".venv", "node_modules", ".git",
    "_cache", ".cache", "cache"
  )
  if ($SourcesOnly) {
    # Data / generated outputs inside the-omega docs tree.
    $excludeDirs += @(
      "data", "datasets", "dist", "artifacts", "outputs", "figures", "images", "img",
      # Do not import non-repo books/assets by default.
      "books", "assets"
    )
  }

  # IMPORTANT: Build a flat argument list; otherwise PowerShell will pass
  # nested arrays as "System.Object[]" to robocopy.
  $args = @(
    $Source,
    $Dest
  ) + $fileSpecs + $modeArgs + @(
    "/R:2",
    "/W:1",
    "/COPY:DAT",
    "/DCOPY:DAT",
    "/NP"
  )

  if ($excludeFiles.Count -gt 0) {
    $args += @("/XF") + $excludeFiles
  }
  if ($excludeDirs.Count -gt 0) {
    $args += @("/XD") + $excludeDirs
  }

  Write-Host "robocopy $Source -> $Dest" -ForegroundColor Cyan
  & robocopy @args | Out-Host
  $rc = $LASTEXITCODE

  # Robocopy exit codes: 0-7 are "success" variants.
  if ($rc -ge 8) {
    throw "robocopy failed with exit code $rc"
  }
}

function Copy-File {
  param(
    [Parameter(Mandatory=$true)][string]$Source,
    [Parameter(Mandatory=$true)][string]$Dest
  )
  $destDir = Split-Path -Parent $Dest
  if (-not (Test-Path -LiteralPath $destDir)) {
    New-Item -ItemType Directory -Force -Path $destDir | Out-Null
  }
  Copy-Item -LiteralPath $Source -Destination $Dest -Force
}

Assert-PathExists -Path $TheOmegaDocsPath -Label "TheOmegaDocsPath"
Assert-PathExists -Path $DestLegacyDocsPath -Label "DestLegacyDocsPath"

Write-Host "==> Pulling latest in: $TheOmegaDocsPath" -ForegroundColor Green
& git -C $TheOmegaDocsPath pull | Out-Host
if ($LASTEXITCODE -ne 0) { throw "git pull failed (exit $LASTEXITCODE)" }

Write-Host "==> Syncing to: $DestLegacyDocsPath" -ForegroundColor Green

#
# IMPORTANT: This repo lives under OneDrive for many users. Copying directly from
# the working tree can fail with ERROR 426 (cloud file hydration timeout).
# To be robust, export the tracked `docs/` content from git objects first, then
# sync from a local temp folder.
#
$theOmegaRepoPath = Split-Path -Parent $TheOmegaDocsPath
Assert-PathExists -Path $theOmegaRepoPath -Label "TheOmega repo root"

$tempRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("the-omega-docs-export-" + [System.Guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null

$zipPath = Join-Path $tempRoot "the-omega-docs.zip"
Write-Host "==> Exporting git-tracked docs to temp: $tempRoot" -ForegroundColor Green
& git -C $theOmegaRepoPath archive --format=zip --output $zipPath "HEAD:docs" | Out-Host
if ($LASTEXITCODE -ne 0) { throw "git archive failed (exit $LASTEXITCODE)" }

$extractDir = Join-Path $tempRoot "extracted"
New-Item -ItemType Directory -Force -Path $extractDir | Out-Null

# Use tar (bsdtar) to avoid Expand-Archive edge cases with directory entries.
& tar -xf $zipPath -C $extractDir
if ($LASTEXITCODE -ne 0) { throw "tar extract failed (exit $LASTEXITCODE)" }
Remove-Item -LiteralPath $zipPath -Force -ErrorAction SilentlyContinue

Run-Robocopy -Source $extractDir -Dest $DestLegacyDocsPath
Prune-NonSources -Root $DestLegacyDocsPath
Remove-Item -LiteralPath $tempRoot -Recurse -Force -ErrorAction SilentlyContinue

Write-Host "==> Done." -ForegroundColor Green
