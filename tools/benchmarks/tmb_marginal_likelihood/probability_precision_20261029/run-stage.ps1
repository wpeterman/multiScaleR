param([Parameter(Mandatory=$true)][string]$Phase)
$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
$taskRepo = (Get-Location).Path
. (Join-Path $taskStudy 'resource-controller.ps1')
if (Test-Path -LiteralPath (Join-Path $taskStudy "$Phase-process.json")) {throw 'Retain existing attempt.'}
Assert-TaskAdmission
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
$env:NOT_CRAN = 'true'
$env:MAKEFLAGS = '-j1'
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
$env:PATH = $env:RSTUDIO_PANDOC + ';' + $env:PATH
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
$taskPhaseLimit = 600
$taskDriver = if ($Phase -eq 'reproduction01') {'reproduce.R'} elseif ($Phase -eq 'documentation01') {'document.R'} elseif ($Phase -match '^parity[0-9]+$') {'parity.R'} else {'validate.R'}
Invoke-TaskPhase $Phase 'C:/Program Files/R/R-4.6.1/bin/x64/Rscript.exe' @('--vanilla',('"'+(Join-Path $taskStudy $taskDriver)+'"'),$Phase) $taskRepo
