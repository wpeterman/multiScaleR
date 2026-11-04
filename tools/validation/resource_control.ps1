# Central resource engine. Deadlines come from a frozen phase plan; reviewed caps stay fixed.
function Write-TaskJson($Value, $Path) {
  if (Test-Path -LiteralPath $Path) { throw "Retain existing record: $Path" }
  $taskTemp = $Path + '.tmp'
  $Value | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskTemp
  Move-Item -LiteralPath $taskTemp -Destination $Path
}
function Get-TaskInventory($OwnedIds) {
  $taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
  # Installed architecture/CMD dispatchers are not R interpreters. Verify the
  # executable bytes and CMD mode; actual R.dll-bearing children count below.
  if (-not $script:taskVerifiedCMDHashes) {
    $script:taskVerifiedCMDHashes = @(
      'C:/Program Files/R/R-4.6.0/bin/R.exe',
      'C:/Program Files/R/R-4.6.0/bin/x64/R.exe',
      'C:/Program Files/R/R-4.6.1/bin/R.exe',
      'C:/Program Files/R/R-4.6.1/bin/x64/R.exe'
    ) | Where-Object { Test-Path -LiteralPath $_ } | ForEach-Object {
      (Get-FileHash -LiteralPath $_ -Algorithm SHA256).Hash
    }
  }
  # Enumerate descendants from this same process snapshot. A child can appear
  # between the caller's inventory and this one.
  if ($OwnedIds.Count) {
    do {
      $taskN = $OwnedIds.Count
      $OwnedIds = @($OwnedIds + @($taskAll | Where-Object { $_.ParentProcessId -in $OwnedIds }).ProcessId | Sort-Object -Unique)
    } while ($OwnedIds.Count -gt $taskN)
  }
  $taskStates = @()
  foreach ($taskR in @($taskAll | Where-Object { $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$' })) {
    $taskState = 'unknown'
    try {
      $taskModules = @(Get-Process -Id $taskR.ProcessId -Module -ErrorAction Stop)
      if (@($taskModules | Where-Object { $_.ModuleName -ieq 'R.dll' }).Count) {
        $taskState = 'R_interpreter'
      } elseif ($taskR.ExecutablePath) {
        $taskHash = (Get-FileHash -LiteralPath $taskR.ExecutablePath -Algorithm SHA256).Hash
        if ($taskR.Name -ieq 'Rscript.exe' -and $taskHash -in $taskRscriptShimHashes) {
          $taskState = 'verified_installed_forwarder'
        } elseif ($taskR.Name -ieq 'R.exe' -and $taskHash -in $script:taskVerifiedCMDHashes -and $taskR.CommandLine -match '(?i)\bCMD\s+\S+') {
          $taskState = 'verified_installed_CMD_forwarder'
        } elseif ($taskR.Name -ieq 'R.exe' -and $taskHash -in $script:taskVerifiedCMDHashes -and $taskR.CommandLine -match '--(no-|vanilla|slave|silent)') {
          # Verified installed R frontends also delegate startup options through
          # Windows cmd.exe. R.dll-bearing processes always take precedence.
          $taskCmdHashes = @('C:/Windows/System32/cmd.exe','C:/Windows/SysWOW64/cmd.exe') | Where-Object { Test-Path -LiteralPath $_ } | ForEach-Object { (Get-FileHash -LiteralPath $_).Hash }
          $taskCmdChildren = @($taskAll | Where-Object { $_.ParentProcessId -eq $taskR.ProcessId -and $_.Name -ieq 'cmd.exe' })
          if ($taskCmdChildren.Count -eq 1 -and $taskCmdChildren[0].ExecutablePath -and (Get-FileHash -LiteralPath $taskCmdChildren[0].ExecutablePath).Hash -in $taskCmdHashes) {
            $taskState = 'verified_installed_startup_forwarder'
          }
        } elseif ($taskR.Name -ieq 'R.exe' -and $taskHash -in $taskRexeHashes -and $taskR.CommandLine -match '(?i)\bCMD\s+BATCH\b') {
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
    RAM_reserve_bytes=8GB;RAM_override_authorized=$true;disk_reserve_bytes=5GB}
}
function Assert-TaskAdmission {
  if ([DateTimeOffset]::UtcNow -ge [DateTimeOffset]::Parse([string]$taskControl.launch_until)) {
    throw 'Resumed development launch window closed; checkpoint only.'
  }
  $taskInventory = Get-TaskInventory @()
  if ($taskInventory.unknown -or $taskInventory.disk_bytes -le 20GB -or $taskInventory.RAM_bytes -lt 8GB -or ($taskInventory.foreign -gt 2 -and -not $taskControl.RAM_and_known_foreign_R_exceptions_authorized)) {
    Write-TaskJson $taskInventory (Join-Path $taskStudy "$Phase-admission-blocked.json")
    throw 'Unknown inventory or morning resource reserve unavailable.'
  }
  $taskCPU=Get-CimInstance Win32_PerfFormattedData_PerfOS_Processor -Filter "Name='_Total'" -ErrorAction Stop
  $taskInventory.global_CPU_percent=[int]$taskCPU.PercentProcessorTime
  if ($taskInventory.global_CPU_percent -gt 75) { Write-TaskJson $taskInventory (Join-Path $taskStudy "$Phase-admission-blocked.json"); throw 'Other work is using the CPU; defer this new producer.' }
  return $taskInventory
}
function Invoke-TaskPhase($Phase, $Executable, $Arguments, $Directory) {
  if (($taskPhaseIsLong -or $Phase -match 'build|check|install|full') -and [DateTimeOffset]::UtcNow -ge [DateTimeOffset]::Parse([string]$taskControl.long_launch_until)) {
    throw 'Long validation launch window closed; checkpoint only.'
  }
  if (Test-Path -LiteralPath (Join-Path $taskStudy "$Phase-process.json")) { throw 'Retain existing attempt.' }
  $taskAdmission = Assert-TaskAdmission
  Write-TaskJson $taskAdmission (Join-Path $taskStudy "$Phase-admission.json")
  $taskClock = [Diagnostics.Stopwatch]::StartNew()
  $taskRoot = Start-Process -FilePath $Executable -ArgumentList $Arguments -WorkingDirectory $Directory -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskStudy "$Phase.log") -RedirectStandardError (Join-Path $taskStudy "$Phase.err")
  $taskIdentity = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskRoot.Id)"
  Write-TaskJson ($taskIdentity | Select-Object ProcessId,ParentProcessId,CreationDate,ExecutablePath,CommandLine) (Join-Path $taskStudy "$Phase-process.json")
  $taskGood = $true
  # Fast phases can finish before this first poll. Retain at least one actual
  # inventory so their completion has the same resource/identity record.
  do {
    $taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
    $taskIds = @($taskRoot.Id)
    do {
      $taskN = $taskIds.Count
      $taskIds = @($taskIds + @($taskAll | Where-Object { $_.ParentProcessId -in $taskIds }).ProcessId | Sort-Object -Unique)
    } while ($taskIds.Count -gt $taskN)
    $taskInventory = Get-TaskInventory $taskIds
    $taskPriorityFailed=((Get-FileHash -LiteralPath $Plan -Algorithm SHA256).Hash -ne $taskPlanHash)
    # Set normal priority only on identity-verified owned interpreters; leave every foreign job alone.
    foreach ($taskOwned in @($taskInventory.states | Where-Object { $_.owned -and $_.state -eq 'R_interpreter' })) {
      $taskFresh=Get-CimInstance Win32_Process -Filter ("ProcessId="+$taskOwned.PID) -ErrorAction SilentlyContinue
      if ($taskFresh -and $taskFresh.CreationDate -eq $taskOwned.creation -and $taskFresh.ExecutablePath -eq $taskOwned.path) {
        $taskHasher=[Security.Cryptography.SHA256]::Create()
        try { $taskDigest=[BitConverter]::ToString($taskHasher.ComputeHash([Text.Encoding]::UTF8.GetBytes([string]$taskFresh.CommandLine))).Replace('-','').ToLowerInvariant() } finally { $taskHasher.Dispose() }
        if ($taskDigest -eq $taskOwned.command_sha256) {
          $taskHandle=Get-Process -Id $taskOwned.PID -ErrorAction SilentlyContinue
          if ($taskHandle -and -not $taskHandle.HasExited) {
            try { $taskHandle.PriorityClass='Normal' }
            catch {
              $taskStillLive=Get-CimInstance Win32_Process -Filter ("ProcessId="+$taskOwned.PID) -ErrorAction SilentlyContinue
              if ($taskStillLive -and $taskStillLive.CreationDate -eq $taskFresh.CreationDate) { $taskPriorityFailed=$true }
            }
          }
        }
      }
    }
    $taskInventory.owned_priority='Normal'
    $taskInventory.phase = $Phase
    $taskInventory.elapsed_seconds = $taskClock.Elapsed.TotalSeconds
    $taskInventory | ConvertTo-Json -Depth 8 -Compress | Add-Content -LiteralPath (Join-Path $taskStudy "$Phase-resources.jsonl")
    # Allow for the 15-second polling interval before the user's hard cutoff.
    if ($taskPriorityFailed -or $taskInventory.unknown -or $taskInventory.owned_interpreters -gt 3 -or $taskInventory.disk_bytes -le 5GB -or $taskInventory.RAM_bytes -lt 8GB -or [DateTimeOffset]::UtcNow -ge [DateTimeOffset]::Parse([string]$taskControl.worker_stop_at)) {
      $taskGood = $false
      $taskOwnedSnapshot = @($taskInventory.states | Where-Object { $_.owned })
      foreach ($taskOld in $taskOwnedSnapshot) {
        $taskFresh = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskOld.PID)" -ErrorAction SilentlyContinue
        if (-not $taskFresh) { continue }
        $taskHasher = [Security.Cryptography.SHA256]::Create()
        $taskCommandHash = [BitConverter]::ToString($taskHasher.ComputeHash([Text.Encoding]::UTF8.GetBytes([string]$taskFresh.CommandLine))).Replace('-','').ToLowerInvariant()
        $taskHasher.Dispose()
        if ($taskFresh.CreationDate -eq $taskOld.creation -and $taskFresh.ExecutablePath -eq $taskOld.path -and $taskCommandHash -eq $taskOld.command_sha256) {
          Stop-Process -Id $taskFresh.ProcessId -Force -ErrorAction SilentlyContinue
        }
      }
      $taskFreshRoot = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskRoot.Id)" -ErrorAction SilentlyContinue
      if ($taskFreshRoot -and $taskIdentity -and $taskFreshRoot.CreationDate -eq $taskIdentity.CreationDate -and $taskFreshRoot.ExecutablePath -eq $taskIdentity.ExecutablePath -and $taskFreshRoot.CommandLine -eq $taskIdentity.CommandLine) {
        Stop-Process -Id $taskRoot.Id -Force -ErrorAction SilentlyContinue
      }
      break
    }
    if ($taskRoot.HasExited) { break }
    Start-Sleep -Seconds 15
    $taskRoot.Refresh()
  } while (-not $taskRoot.HasExited)
  $taskRoot.WaitForExit()
  Write-TaskJson @{operational_passed=$taskGood;exit_code=$taskRoot.ExitCode;
    elapsed_seconds=$taskClock.Elapsed.TotalSeconds;completed_at=(Get-Date -Format o)} (Join-Path $taskStudy "$Phase-exit.json")
  if (-not $taskGood -or $taskRoot.ExitCode -ne 0) { throw "Validation phase failed: $Phase" }
}
