param(
    [Parameter(Mandatory = $true)][string]$WorkDir,
    [Parameter(Mandatory = $true)][string]$GctxPath,
    [string]$Python = 'python'
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$phase171A = Join-Path $repo 'data\frozen_phase171a'
$scripts = Join-Path $repo 'src\phase171b'

New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
$env:CLEANROOM_PHASE171A = (Resolve-Path -LiteralPath $phase171A).Path
$env:CLEANROOM_PHASE171B = (Resolve-Path -LiteralPath $WorkDir).Path
$env:CLEANROOM_GCTX = (Resolve-Path -LiteralPath $GctxPath).Path
$env:PYTHONPATH = $scripts

& $Python (Join-Path $repo 'scripts\preflight.py')
if ($LASTEXITCODE -ne 0) { throw "preflight failed: $LASTEXITCODE" }
& $Python (Join-Path $scripts 'prepare_expression.py')
if ($LASTEXITCODE -ne 0) { throw "prepare_expression failed: $LASTEXITCODE" }
& $Python (Join-Path $scripts 'observed_scores.py')
if ($LASTEXITCODE -ne 0) { throw "observed_scores failed: $LASTEXITCODE" }
