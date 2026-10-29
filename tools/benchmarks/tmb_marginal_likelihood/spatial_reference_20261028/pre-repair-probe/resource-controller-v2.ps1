function Write-TaskJson($Value, $Path) {
  $taskTemp = $Path + '.tmp'
  $Value | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $taskTemp
  if (Test-Path -LiteralPath $Path) {
    [IO.File]::Replace($taskTemp, $Path, $null)
  } else {
    [IO.File]::Move($taskTemp, $Path)
  }
}
function Assert-TaskAdmission {
  if ([DateTimeOffset]::UtcNow -ge [DateTimeOffset]'2026-10-29T04:00:00-04:00') { throw 'Long-launch window closed.' }
  $taskExisting = @(Get-CimInstance Win32_Process -ErrorAction Stop | Where-Object {
    $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$'
  })
  if ($taskExisting.Count) { throw 'Validation phase needs an R-free window.' }
  $taskM = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
  $taskD = Get-PSDrive C -ErrorAction Stop
  if (-not $taskM.FreePhysicalMemory -or -not $taskD.Free -or
      $taskM.FreePhysicalMemory*1KB -le 12GB -or $taskD.Free -le 5GB) { throw 'Resource reserve unavailable.' }
}
function Invoke-TaskPhase($Phase, $Executable, $Arguments, $Directory) {
  Assert-TaskAdmission
  $taskPhaseClock = [Diagnostics.Stopwatch]::StartNew()
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
    $taskStates = @()
    foreach ($taskR in @($taskAll | Where-Object { $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$' })) {
      $taskState = 'unknown'
      try {
        $taskModules = @(Get-Process -Id $taskR.ProcessId -Module -ErrorAction Stop)
        if (@($taskModules | Where-Object { $_.ModuleName -ieq 'R.dll' }).Count) {
          $taskState = 'R_interpreter'
        } elseif ($taskR.ProcessId -in $taskIds -and $taskR.ExecutablePath) {
          $taskBinary = (Get-FileHash -LiteralPath $taskR.ExecutablePath -Algorithm SHA256).Hash
          if ($taskR.Name -ieq 'Rscript.exe' -and $taskBinary -eq $taskRscriptShimHash) {
            $taskState = 'verified_installed_forwarder'
          } elseif ($taskR.Name -ieq 'R.exe' -and $taskBinary -eq $taskRexeHash -and $taskR.CommandLine -match '(?i)\bCMD\s+BATCH\b') {
            $taskChild = @($taskAll | Where-Object { $_.ParentProcessId -eq $taskR.ProcessId -and $_.Name -ieq 'Rterm.exe' })
            if ($taskChild.Count -eq 1 -and @((Get-Process -Id $taskChild[0].ProcessId -Module -ErrorAction Stop) | Where-Object { $_.ModuleName -ieq 'R.dll' }).Count) {
              $taskState = 'verified_BATCH_forwarder'
            }
          }
        }
      } catch {
        $taskFreshProcess = Get-Process -Id $taskR.ProcessId -ErrorAction SilentlyContinue
        if (-not $taskFreshProcess -or $taskFreshProcess.HasExited) { $taskState = 'exited_during_inventory' }
      }
      $taskStates += @{PID=$taskR.ProcessId;parent=$taskR.ParentProcessId;creation=$taskR.CreationDate;
        path=$taskR.ExecutablePath;command=$taskR.CommandLine;state=$taskState;owned=($taskR.ProcessId -in $taskIds)}
    }
    $taskCount = @($taskStates | Where-Object { $_.state -eq 'R_interpreter' }).Count
    $taskUnknown = @($taskStates | Where-Object { $_.state -eq 'unknown' }).Count
    $taskForeign = @($taskStates | Where-Object { -not $_.owned -and $_.state -ne 'exited_during_inventory' }).Count
    $taskM = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
    $taskD = Get-PSDrive C -ErrorAction Stop
    @{time=(Get-Date -Format o);phase=$Phase;actual_interpreters=$taskCount;unknown=$taskUnknown;foreign=$taskForeign;
      states=$taskStates;RAM_bytes=$taskM.FreePhysicalMemory*1KB;disk_bytes=$taskD.Free;
      elapsed_seconds=$taskPhaseClock.Elapsed.TotalSeconds;phase_limit_seconds=$taskPhaseLimit} |
      ConvertTo-Json -Depth 6 -Compress | Add-Content -LiteralPath (Join-Path $taskStudy "$Phase-resources.jsonl")
    if (($taskPhaseLimit -gt 0 -and $taskPhaseClock.Elapsed.TotalSeconds -ge $taskPhaseLimit) -or
        $taskUnknown -or $taskForeign -or $taskCount -gt 4 -or -not $taskM.FreePhysicalMemory -or -not $taskD.Free -or
        $taskM.FreePhysicalMemory*1KB -le 12GB -or $taskD.Free -le 5GB -or [DateTimeOffset]::UtcNow -ge [DateTimeOffset]'2026-10-29T04:35:00-04:00') {
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
  $taskResult = @{operational_passed=$taskGood;exit_code=$taskRoot.ExitCode;completed_at=(Get-Date -Format o);
    elapsed_seconds=$taskPhaseClock.Elapsed.TotalSeconds;phase_limit_seconds=$taskPhaseLimit}
  Write-TaskJson $taskResult (Join-Path $taskStudy "$Phase-exit.json")
  if (-not $taskGood -or $taskRoot.ExitCode -ne 0) { throw "Validation phase failed: $Phase" }
}
