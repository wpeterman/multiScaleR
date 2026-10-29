param([Parameter(Mandatory=$true)][string]$Mode)
$ErrorActionPreference = 'Stop'
$taskNow = [DateTimeOffset]::UtcNow
if ($taskNow -ge [DateTimeOffset]'2026-10-29T04:45:00-04:00') { throw 'Authorized work window is closed.' }
if ($Mode -eq 'check01' -and $taskNow -ge [DateTimeOffset]'2026-10-29T04:00:00-04:00') {
  throw 'No long check launch after 04:00 Eastern.'
}
if ($taskNow -ge [DateTimeOffset]'2026-10-29T04:15:00-04:00') { throw 'Checkpoint-only window; no numerical launches.' }
$taskStudy = Join-Path $PSScriptRoot '.'
if (Test-Path -LiteralPath (Join-Path $taskStudy "$Mode.log")) { throw 'Retain earlier attempts; choose a new mode.' }
$taskProcesses = @(Get-CimInstance Win32_Process -ErrorAction Stop)
$taskR = @($taskProcesses | Where-Object { $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$' })
if (@($taskR | Where-Object { -not $_.CreationDate -or -not $_.CommandLine -or -not $_.ExecutablePath }).Count) {
  throw 'Incomplete R process inventory; numerical admission closed.'
}
if ($taskR.Count -ge 4) { throw 'Global R interpreter budget is full.' }
if ($Mode -eq 'check01' -and $taskR.Count) { throw 'Full check needs a foreign-R-free window.' }
if (@($taskR | Where-Object { $_.CommandLine -like '*prediction_contracts_20261028*' }).Count) {
  throw 'A validation worker is already active.'
}
$taskMemory = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
$taskDisk = Get-PSDrive C -ErrorAction Stop
if (-not $taskMemory.FreePhysicalMemory -or -not $taskDisk.Free) { throw 'Unknown resource inventory.' }
if ($taskMemory.FreePhysicalMemory * 1KB -le 12GB -or $taskDisk.Free -le 5GB) { throw 'Resource reserve unavailable.' }
@{mode=$Mode; foreign_R=@($taskR | Select-Object ProcessId,ParentProcessId,CommandLine,CreationDate);
  available_RAM_bytes=$taskMemory.FreePhysicalMemory*1KB; free_disk_bytes=$taskDisk.Free} |
  ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $taskStudy "$Mode-admission.json")
foreach ($taskName in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) {
  [Environment]::SetEnvironmentVariable($taskName, '1', 'Process')
}
foreach ($taskLocale in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) {
  [Environment]::SetEnvironmentVariable($taskLocale, $null, 'Process')
}
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
$taskScript = if ($Mode -eq 'parity01') { 'parity.R' } elseif ($Mode -eq 'documentation01') { 'document.R' } elseif ($Mode -eq 'check01') { 'check.R' } else { 'validate.R' }
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' --vanilla (Join-Path $taskStudy $taskScript) $Mode *> (Join-Path $taskStudy "$Mode.log")
$taskExit = $LASTEXITCODE
@{mode=$Mode;exit_code=$taskExit;completed_at=(Get-Date -Format o)} | ConvertTo-Json |
  Set-Content -LiteralPath (Join-Path $taskStudy "$Mode-exit.json")
exit $taskExit
