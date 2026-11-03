param([Parameter(Mandatory=$true)][string]$Plan, [Parameter(Mandatory=$true)][string]$Phase)
$ErrorActionPreference = 'Stop'
$taskPlan = Get-Content -LiteralPath $Plan -Raw | ConvertFrom-Json
$taskPlanHash = (Get-FileHash -LiteralPath $Plan -Algorithm SHA256).Hash
$taskControl = Get-Content -LiteralPath $taskPlan.control_file -Raw | ConvertFrom-Json
$taskRepo = $taskPlan.repo
$taskStudy = $taskPlan.output
$taskRAMReserve = 8GB
$taskPhase = @($taskPlan.phases | Where-Object { $_.name -eq $Phase })
if ($taskPhase.Count -ne 1) { throw 'Declare one exact phase.' }
$taskPhase = $taskPhase[0]
$taskPhaseIsLong = [bool]$taskPhase.long
. (Join-Path $PSScriptRoot 'resource_control.ps1')
$taskRscriptShimHashes = @('C:/Program Files/R/R-4.6.0/bin/Rscript.exe','C:/Program Files/R/R-4.6.1/bin/Rscript.exe') | Where-Object { Test-Path -LiteralPath $_ } | ForEach-Object { (Get-FileHash -LiteralPath $_ -Algorithm SHA256).Hash }
$taskRexeHashes = @('C:/Program Files/R/R-4.6.0/bin/x64/R.exe','C:/Program Files/R/R-4.6.1/bin/x64/R.exe') | Where-Object { Test-Path -LiteralPath $_ } | ForEach-Object { (Get-FileHash -LiteralPath $_ -Algorithm SHA256).Hash }
$env:R_LIBS = $taskPlan.library_paths -join ';'
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
$env:PATH = $env:RSTUDIO_PANDOC + ';' + $env:PATH
$env:MAKEFLAGS = '-j1'
foreach ($taskName in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskName, '1', 'Process') }
foreach ($taskName in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskName, $null, 'Process') }
$env:NOT_CRAN = 'true'
if ($taskPhase.installed_suite) { [Environment]::SetEnvironmentVariable('_R_CHECK_PACKAGE_NAME_', 'multiScaleR', 'Process') }
try {
  Invoke-TaskPhase $Phase $taskPhase.executable @($taskPhase.argument_line) $taskPhase.working_directory
} catch {
  # A controller exception must not orphan its own R workers. A fresh snapshot
  # must agree with the recorded root birth/path/command before descendants count.
  $taskRecordPath = Join-Path $taskStudy "$Phase-process.json"
  if (Test-Path -LiteralPath $taskRecordPath) {
    $taskRootRecord = Get-Content -LiteralPath $taskRecordPath -Raw | ConvertFrom-Json
    $taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
    $taskRootLive = @($taskAll | Where-Object { $_.ProcessId -eq $taskRootRecord.ProcessId })
    if ($taskRootLive.Count -eq 1 -and $taskRootLive[0].CreationDate -eq [datetime]$taskRootRecord.CreationDate -and $taskRootLive[0].ExecutablePath -eq $taskRootRecord.ExecutablePath -and $taskRootLive[0].CommandLine -eq $taskRootRecord.CommandLine) {
      $taskIds = @($taskRootLive[0].ProcessId)
      do {
        $taskN = $taskIds.Count
        $taskIds = @($taskIds + @($taskAll | Where-Object { $_.ParentProcessId -in $taskIds }).ProcessId | Sort-Object -Unique)
      } while ($taskIds.Count -gt $taskN)
      foreach ($taskOld in @($taskAll | Where-Object { $_.ProcessId -in $taskIds } | Sort-Object ProcessId -Descending)) {
        $taskLive = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskOld.ProcessId)" -ErrorAction SilentlyContinue
        if ($taskLive -and $taskLive.CreationDate -eq $taskOld.CreationDate -and $taskLive.ExecutablePath -eq $taskOld.ExecutablePath -and $taskLive.CommandLine -eq $taskOld.CommandLine) {
          Stop-Process -Id $taskLive.ProcessId -Force -ErrorAction SilentlyContinue
        }
      }
    }
  }
  throw
}
if ((Get-FileHash -LiteralPath $Plan -Algorithm SHA256).Hash -ne $taskPlanHash) { throw 'Frozen phase plan changed.' }
