# Resume the unfinished validation phases against the unchanged built archive.
$ErrorActionPreference = 'Stop'
$taskRepo = (Get-Location).Path
$taskStudy = $PSScriptRoot
$taskRAMReserve = 4GB
. (Join-Path $taskStudy 'resource-controller.ps1')
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
$taskSnapshot = Get-Content -LiteralPath (Join-Path $taskStudy 'snapshot01.json') -Raw | ConvertFrom-Json
$taskRoot = $taskSnapshot.workspace
$taskLibraries = @((Join-Path $taskRepo 'tools/benchmarks/joint_scale_assessment/phase31/library'), 'C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library', 'C:/Users/peterman.73/AppData/Local/R/win-library/4.6', 'C:/Program Files/R/R-4.6.1/library')
$env:R_LIBS = $taskLibraries -join ';'
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
$env:PATH = $env:RSTUDIO_PANDOC + ';' + $env:PATH
$env:MAKEFLAGS = '-j1'
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
$env:NOT_CRAN = 'true'
$taskRetryRoot = Join-Path $taskRoot 'check02'
if (Test-Path -LiteralPath $taskRetryRoot) { throw 'Retain existing retry workspace.' }
New-Item -ItemType Directory -Path $taskRetryRoot | Out-Null
$taskArchive = Join-Path $taskRoot 'multiScaleR_0.7.5.tar.gz'
Copy-Item -LiteralPath $taskArchive -Destination $taskRetryRoot
Invoke-TaskPhase 'check02' 'C:/Program Files/R/R-4.6.1/bin/x64/Rcmd.exe' @('check','multiScaleR_0.7.5.tar.gz','--as-cran','--no-manual','--no-tests') $taskRetryRoot
$taskCheckLog = Join-Path $taskRetryRoot 'multiScaleR.Rcheck/00check.log'
Copy-Item -LiteralPath $taskCheckLog -Destination (Join-Path $taskStudy 'R-CMD-check02.log')
if ((Get-Content -LiteralPath $taskCheckLog -Raw) -notmatch '(?m)^Status: OK\s*$') { throw 'Package check is not clean.' }
$taskTests = Join-Path $taskRoot 'exact-archive-tests'
& 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'extract-tests.py') $taskArchive $taskTests
if ($LASTEXITCODE -ne 0) { throw 'Exact archive test extraction failed.' }
$env:R_LIBS = (@((Join-Path $taskRetryRoot 'multiScaleR.Rcheck')) + $taskLibraries) -join ';'
[Environment]::SetEnvironmentVariable('_R_CHECK_PACKAGE_NAME_','multiScaleR','Process')
Invoke-TaskPhase 'installed01' 'C:/Program Files/R/R-4.6.1/bin/x64/Rterm.exe' @('--vanilla','--no-echo','--no-save','--no-restore','-f','testthat.R') $taskTests
Write-TaskJson @{passed=$true;check_status_OK=$true;check_phase='check02';tests_complete=$true;
  NOT_CRAN=$true;tests_from_exact_archive=$true;workspace=$taskRoot;completed_at=(Get-Date -Format o);
  concurrency_exception=$true;RAM_reserve_bytes=$taskRAMReserve;timing_comparison_eligible=$false} (Join-Path $taskStudy 'validation01.json')
