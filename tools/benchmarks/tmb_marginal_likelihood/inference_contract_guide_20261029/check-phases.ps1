$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
$taskRepo = (Get-Location).Path
. (Join-Path $taskStudy 'resource-controller.ps1')
if (Test-Path -LiteralPath (Join-Path $taskStudy 'check01-input.json')) { throw 'Retain existing check attempt.' }
Assert-TaskAdmission
$taskSnapshot = Get-Content -LiteralPath (Join-Path $taskStudy 'snapshot01-input.json') -Raw | ConvertFrom-Json
foreach ($taskFile in $taskSnapshot.source_sha256.PSObject.Properties) {
  if ((Get-FileHash -LiteralPath (Join-Path $taskRepo $taskFile.Name) -Algorithm SHA256).Hash.ToLowerInvariant() -ne $taskFile.Value) { throw "Source drift: $($taskFile.Name)" }
}
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
$taskLibraries = @(Get-Content -LiteralPath (Join-Path $taskStudy 'render02.json') -Raw | ConvertFrom-Json).libraries
$env:R_LIBS = $taskLibraries -join ';'
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
if (-not (Test-Path -LiteralPath (Join-Path $env:RSTUDIO_PANDOC 'pandoc.exe'))) { throw 'Pandoc is unavailable.' }
$env:PATH = $env:RSTUDIO_PANDOC + ';' + $env:PATH
$env:MAKEFLAGS = '-j1'
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
$taskRoot = $taskSnapshot.workspace
Write-TaskJson @{workspace=$taskRoot;source_files=$taskSnapshot.source_files;installed_tests_rerun=$false;
  Rscript_forwarder_SHA256=$taskRscriptShimHash;R_BATCH_forwarder_SHA256=$taskRexeHash;
  NOT_CRAN=$true;started_at=(Get-Date -Format o)} (Join-Path $taskStudy 'check01-input.json')
Invoke-TaskPhase 'build01' 'C:/Program Files/R/R-4.6.1/bin/x64/Rcmd.exe' @('build','multiScaleR') $taskRoot
$taskArchive = Join-Path $taskRoot 'multiScaleR_0.7.4.tar.gz'
$taskSHA = (Get-FileHash -LiteralPath $taskArchive -Algorithm SHA256).Hash
Write-TaskJson @{archive=$taskArchive;SHA256=$taskSHA;build_complete=$true} (Join-Path $taskStudy 'build01.json')
Invoke-TaskPhase 'check01' 'C:/Program Files/R/R-4.6.1/bin/x64/Rcmd.exe' @('check','multiScaleR_0.7.4.tar.gz','--as-cran','--no-manual','--no-tests') $taskRoot
$taskCheckLog = Join-Path $taskRoot 'multiScaleR.Rcheck/00check.log'
Copy-Item -LiteralPath $taskCheckLog -Destination (Join-Path $taskStudy 'R-CMD-check.log')
if ((Get-Content -LiteralPath $taskCheckLog -Raw) -notmatch '(?m)^Status: OK\s*$') { throw 'Package check is not clean.' }
$taskStable = $true
foreach ($taskFile in $taskSnapshot.source_sha256.PSObject.Properties) {
  if ((Get-FileHash -LiteralPath (Join-Path $taskRepo $taskFile.Name) -Algorithm SHA256).Hash.ToLowerInvariant() -ne $taskFile.Value) { $taskStable = $false }
}
Write-TaskJson @{passed=$taskStable;source_stable=$taskStable;check_status_OK=$true;installed_tests_rerun=$false;
  retained_complete_test_evidence='spatial_source_contract_20261029/final-source-evidence.json';
  archive_SHA256=$taskSHA;workspace=$taskRoot;completed_at=(Get-Date -Format o)} (Join-Path $taskStudy 'check01.json')
if (-not $taskStable) { throw 'Live source drifted during validation.' }
