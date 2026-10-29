param([Parameter(Mandatory=$true)][string]$Phase)
$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
$taskRepo = (Get-Location).Path
$taskPreceding = Join-Path $taskRepo 'tools/benchmarks/tmb_marginal_likelihood/functional_precision_20261028'
. (Join-Path $taskStudy 'resource-controller.ps1')
if (Test-Path -LiteralPath (Join-Path $taskStudy "$Phase-process.json")) { throw 'Retain existing unit attempt.' }
Assert-TaskAdmission
& 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskPreceding 'verify-candidate.py')
if ($LASTEXITCODE -ne 0) { throw 'Preceding package gate/source is incomplete or changed.' }
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
$taskLibraries = @(Get-Content -LiteralPath (Join-Path $taskPreceding 'compatibility01.json') -Raw | ConvertFrom-Json).libraries
$env:R_LIBS = $taskLibraries -join ';'
$env:NOT_CRAN = 'true'
$env:MAKEFLAGS = '-j1'
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
Invoke-TaskPhase $Phase 'C:/Program Files/R/R-4.6.1/bin/x64/Rscript.exe' @('--vanilla',('"' + (Join-Path $taskStudy 'validate-adapter.R') + '"'),$Phase) $taskRepo
