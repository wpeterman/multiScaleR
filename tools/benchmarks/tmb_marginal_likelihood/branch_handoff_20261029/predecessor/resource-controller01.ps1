# Morning user authorization permits concurrent foreign R work. Preserve identity
# checks and scientific gates; record this separately from the overnight limits.
function Write-TaskJson($Value, $Path) {
  if (Test-Path -LiteralPath $Path) { throw "Retain existing record: $Path" }
  $taskTemp = $Path + '.tmp'
  $Value | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskTemp
  Move-Item -LiteralPath $taskTemp -Destination $Path
}
function Get-TaskInventory($OwnedIds) {
  $taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
  $taskStates = @()
  foreach ($taskR in @($taskAll | Where-Object { $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$' })) {
    $taskState = 'unknown'
    try {
      $taskModules = @(Get-Process -Id $taskR.ProcessId -Module -ErrorAction Stop)
      if (@($taskModules | Where-Object { $_.ModuleName -ieq 'R.dll' }).Count) {
        $taskState = 'R_interpreter'
      } elseif ($taskR.ExecutablePath) {
        $taskHash = (Get-FileHash -LiteralPath $taskR.ExecutablePath -Algorithm SHA256).Hash
        if ($taskR.Name -ieq 'Rscript.exe' -and $taskHash -eq $taskRscriptShimHash) {
          $taskState = 'verified_installed_forwarder'
        } elseif ($taskR.Name -ieq 'R.exe' -and $taskHash -eq $taskRexeHash -and $taskR.CommandLine -match '(?i)\bCMD\s+BATCH\b') {
          $taskChildren = @($taskAll | Where-Object { $_.ParentProcessId -eq $taskR.ProcessId -and $_.Name -ieq 'Rterm.exe' })
          if ($taskChildren.Count -eq 1) { $taskState = 'verified_BATCH_forwarder' }
        }
      }
    } catch {
      if (-not (Get-Process -Id $taskR.ProcessId -ErrorAction SilentlyContinue)) { $taskState = 'exited_during_inventory' }
    }
    # Commands may contain foreign RStudio launcher credentials. Store a digest.
    $taskHasher = [Security.Cryptography.SHA256]::Create()
    $taskCommandHash = [BitConverter]::ToString($taskHasher.ComputeHash([Text.Encoding]::UTF8.GetBytes([string]$taskR.CommandLine))).Replace('-','').ToLowerInvariant()
    $taskHasher.Dispose()
    $taskStates += @{PID=$taskR.ProcessId;parent=$taskR.ParentProcessId;creation=$taskR.CreationDate;
      path=$taskR.ExecutablePath;command_sha256=$taskCommandHash;state=$taskState;owned=($taskR.ProcessId -in $OwnedIds)}
  }
  $taskM = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
  $taskD = Get-PSDrive C -ErrorAction Stop
  return @{time=(Get-Date -Format o);states=$taskStates;
    actual_interpreters=@($taskStates | Where-Object { $_.state -eq 'R_interpreter' }).Count;
    owned_interpreters=@($taskStates | Where-Object { $_.owned -and $_.state -eq 'R_interpreter' }).Count;
    unknown=@($taskStates | Where-Object { $_.state -eq 'unknown' }).Count;
    foreign=@($taskStates | Where-Object { -not $_.owned -and $_.state -eq 'R_interpreter' }).Count;
    RAM_bytes=$taskM.FreePhysicalMemory*1KB;disk_bytes=$taskD.Free;
    user_concurrency_exception=$true;one_owned_validation_pipeline=$true;threads=1;
    RAM_reserve_bytes=8GB;disk_reserve_bytes=5GB}
}
function Assert-TaskAdmission {
  $taskInventory = Get-TaskInventory @()
  if ($taskInventory.unknown -or $taskInventory.RAM_bytes -le 8GB -or $taskInventory.disk_bytes -le 5GB) {
    $taskInventory | ConvertTo-Json -Depth 6 | Write-Output
    throw 'Unknown inventory or morning resource reserve unavailable.'
  }
  return $taskInventory
}
function Invoke-TaskPhase($Phase, $Executable, $Arguments, $Directory) {
  if (Test-Path -LiteralPath (Join-Path $taskStudy "$Phase-process.json")) { throw 'Retain existing attempt.' }
  $taskAdmission = Assert-TaskAdmission
  Write-TaskJson $taskAdmission (Join-Path $taskStudy "$Phase-admission.json")
  $taskClock = [Diagnostics.Stopwatch]::StartNew()
  $taskRoot = Start-Process -FilePath $Executable -ArgumentList $Arguments -WorkingDirectory $Directory -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskStudy "$Phase.log") -RedirectStandardError (Join-Path $taskStudy "$Phase.err")
  $taskIdentity = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskRoot.Id)"
  Write-TaskJson ($taskIdentity | Select-Object ProcessId,ParentProcessId,CreationDate,ExecutablePath,CommandLine) (Join-Path $taskStudy "$Phase-process.json")
  $taskGood = $true
  while (-not $taskRoot.HasExited) {
    $taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
    $taskIds = @($taskRoot.Id)
    do {
      $taskN = $taskIds.Count
      $taskIds = @($taskIds + @($taskAll | Where-Object { $_.ParentProcessId -in $taskIds }).ProcessId | Sort-Object -Unique)
    } while ($taskIds.Count -gt $taskN)
    $taskInventory = Get-TaskInventory $taskIds
    $taskInventory.phase = $Phase
    $taskInventory.elapsed_seconds = $taskClock.Elapsed.TotalSeconds
    $taskInventory | ConvertTo-Json -Depth 8 -Compress | Add-Content -LiteralPath (Join-Path $taskStudy "$Phase-resources.jsonl")
    if ($taskInventory.unknown -or $taskInventory.owned_interpreters -gt 4 -or $taskInventory.RAM_bytes -le 8GB -or $taskInventory.disk_bytes -le 5GB) {
      $taskGood = $false
      foreach ($taskOld in @($taskAll | Where-Object { $_.ProcessId -in $taskIds })) {
        $taskFresh = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskOld.ProcessId)" -ErrorAction SilentlyContinue
        if ($taskFresh -and $taskFresh.CreationDate -eq $taskOld.CreationDate -and $taskFresh.ExecutablePath -eq $taskOld.ExecutablePath -and $taskFresh.CommandLine -eq $taskOld.CommandLine) {
          Stop-Process -Id $taskFresh.ProcessId -Force -ErrorAction SilentlyContinue
        }
      }
      break
    }
    Start-Sleep -Seconds 15
    $taskRoot.Refresh()
  }
  $taskRoot.WaitForExit()
  Write-TaskJson @{operational_passed=$taskGood;exit_code=$taskRoot.ExitCode;
    elapsed_seconds=$taskClock.Elapsed.TotalSeconds;completed_at=(Get-Date -Format o)} (Join-Path $taskStudy "$Phase-exit.json")
  if (-not $taskGood -or $taskRoot.ExitCode -ne 0) { throw "Validation phase failed: $Phase" }
}
