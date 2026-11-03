param([Parameter(Mandatory=$true)][string]$Output, [Parameter(Mandatory=$true)][string]$Phase,
      [Parameter(Mandatory=$true)][string]$TargetName)
$ErrorActionPreference = 'Stop'
if ($Phase -notmatch '^[a-z][a-z0-9_-]*$' -or $TargetName -notmatch '^[a-z][a-z0-9_-]*[.]json$') { throw 'Use literal phase and record names.' }
$taskTarget = Join-Path $Output $TargetName
if (Test-Path -LiteralPath $taskTarget) { throw 'Retain the existing ownership seal.' }
$taskExit = Get-Content -LiteralPath (Join-Path $Output "$Phase-exit.json") -Raw | ConvertFrom-Json
$taskRecords = @(foreach ($taskLine in Get-Content -LiteralPath (Join-Path $Output "$Phase-resources.jsonl")) {
  $taskValue = $taskLine | ConvertFrom-Json
  foreach ($taskState in $taskValue.states) { if ($taskState.owned) { $taskState } }
})
$taskProcess = Get-Content -LiteralPath (Join-Path $Output "$Phase-process.json") -Raw | ConvertFrom-Json
$taskHasher = [Security.Cryptography.SHA256]::Create()
try { $taskHash = [BitConverter]::ToString($taskHasher.ComputeHash([Text.Encoding]::UTF8.GetBytes([string]$taskProcess.CommandLine))).Replace('-','').ToLowerInvariant() } finally { $taskHasher.Dispose() }
$taskRecords += @{PID=$taskProcess.ProcessId;creation=$taskProcess.CreationDate;path=$taskProcess.ExecutablePath;command_sha256=$taskHash}
$taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
$taskLive = @(foreach ($taskOld in $taskRecords) {
  foreach ($taskCurrent in $taskAll) {
    if ($taskCurrent.ProcessId -ne $taskOld.PID) { continue }
    if ($null -eq $taskCurrent.CreationDate) { throw 'Unknown matching-PID creation identity.' }
    if ($taskCurrent.CreationDate.ToUniversalTime().Ticks -ne ([DateTimeOffset]$taskOld.creation).UtcTicks) { continue }
    if (-not $taskCurrent.ExecutablePath) { throw 'Unknown matching-birth path identity.' }
    if ($taskCurrent.ExecutablePath -ne $taskOld.path) { continue }
    $taskHasher = [Security.Cryptography.SHA256]::Create()
    try { $taskHash = [BitConverter]::ToString($taskHasher.ComputeHash([Text.Encoding]::UTF8.GetBytes([string]$taskCurrent.CommandLine))).Replace('-','').ToLowerInvariant() } finally { $taskHasher.Dispose() }
    if ($taskHash -ne $taskOld.command_sha256) { throw 'Matching-birth process command changed.' }
    $taskCurrent.ProcessId
  }
})
$taskLive = @($taskLive | Sort-Object -Unique)
if ($taskLive.Count) { throw ('This phase still owns live processes: ' + ($taskLive -join ',')) }
$taskValue = @{phase=$Phase;checked_at=(Get-Date -Format o);verified_phase_owned_processes_remaining=0;IDs=$taskLive;exit_code=$taskExit.exit_code;other_phases_and_foreign_jobs_untouched=$true}
[IO.File]::WriteAllText(($taskTarget+'.tmp'),($taskValue | ConvertTo-Json -Depth 5),[Text.UTF8Encoding]::new($false))
Move-Item -LiteralPath ($taskTarget+'.tmp') -Destination $taskTarget
$taskValue | ConvertTo-Json -Depth 5
