$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
$taskRepo = (Get-Location).Path
if (Test-Path -LiteralPath (Join-Path $taskStudy 'supplement05-process.json')) { throw 'Inspect existing supplemental attempt.' }
# Reuse only our existing controller function definitions, without executing its
# top-level check launcher or modifying package/test code.
$taskTokens = $null
$taskParseErrors = $null
$taskAst = [System.Management.Automation.Language.Parser]::ParseFile((Join-Path $taskStudy 'check-phases-v2.ps1'), [ref]$taskTokens, [ref]$taskParseErrors)
if ($taskParseErrors.Count) { throw 'Controller parsing failed.' }
foreach ($taskDefinition in @($taskAst.EndBlock.Statements | Where-Object { $_ -is [System.Management.Automation.Language.FunctionDefinitionAst] })) {
  . ([ScriptBlock]::Create($taskDefinition.Extent.Text))
}
Assert-TaskAdmission
$taskCompleted = Get-Content -LiteralPath (Join-Path $taskStudy 'check04.json') -Raw | ConvertFrom-Json
if (-not $taskCompleted.passed) { throw 'Original check/test phases have not passed.' }
$taskTestDir = $taskCompleted.installed_test_directory
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
$taskLibraries = @(Get-Content -LiteralPath (Join-Path $taskStudy 'compatibility02.json') -Raw | ConvertFrom-Json).libraries
$env:R_LIBS = (@((Join-Path $taskCompleted.check_root 'multiScaleR.Rcheck')) + @($taskLibraries)) -join ';'
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
$env:PATH = $env:RSTUDIO_PANDOC + ';' + $env:PATH
$env:NOT_CRAN = 'true'
$env:MSR_CONTRACT_STUDY = $taskStudy
$env:MSR_CONTRACT_INSTALL = Join-Path $taskCompleted.check_root 'multiScaleR.Rcheck/multiScaleR'
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
Copy-Item -LiteralPath (Join-Path $taskStudy 'supplement.R') -Destination (Join-Path $taskTestDir 'supplement-driver.R')
Invoke-TaskPhase 'supplement05' 'C:/Program Files/R/R-4.6.1/bin/x64/Rterm.exe' @('--vanilla','--no-echo','--no-save','--no-restore','-f','supplement-driver.R') $taskTestDir
