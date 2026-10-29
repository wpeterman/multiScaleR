$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
if ([DateTimeOffset]::UtcNow -ge [DateTimeOffset]'2026-10-29T04:00:00-04:00') { throw 'Long-launch window closed.' }
if (Test-Path -LiteralPath (Join-Path $taskStudy 'check02.json')) { throw 'Retain completed attempt.' }
if (Test-Path -LiteralPath (Join-Path $taskStudy 'check02-process.json')) { throw 'Inspect existing attempt before reuse.' }
$taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
$taskR = @($taskAll | Where-Object { $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$' })
if ($taskR.Count) { throw 'Direct full check requires an R-free admission window.' }
$taskMemory = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
$taskDisk = Get-PSDrive C -ErrorAction Stop
if (-not $taskMemory.FreePhysicalMemory -or -not $taskDisk.Free -or
    $taskMemory.FreePhysicalMemory*1KB -le 12GB -or $taskDisk.Free -le 5GB) { throw 'Resource reserve unavailable.' }
$taskRepo = (Get-Location).Path
$taskManifest = Import-Csv -LiteralPath (Join-Path $taskStudy 'check01/checked-source-md5.csv')
foreach ($taskFile in $taskManifest) {
  if ((Get-FileHash -LiteralPath (Join-Path $taskRepo $taskFile.file) -Algorithm MD5).Hash.ToLowerInvariant() -ne $taskFile.md5) {
    throw "Source drift before archive reuse: $($taskFile.file)"
  }
}
$taskPrevious = (Get-Content -LiteralPath (Join-Path $taskStudy 'check01/check-location.txt')).Trim()
$taskArchive = Join-Path $taskPrevious 'check/multiScaleR_0.7.4.tar.gz'
if (-not (Test-Path -LiteralPath $taskArchive)) { throw 'Completed build archive missing.' }
$taskArchiveHash = (Get-FileHash -LiteralPath $taskArchive -Algorithm SHA256).Hash
$taskCheckRoot = Join-Path $env:TEMP ('msr-contract-direct-' + [Guid]::NewGuid().ToString('N').Substring(0,12))
New-Item -ItemType Directory -Path $taskCheckRoot | Out-Null
$taskCopiedArchive = Join-Path $taskCheckRoot 'multiScaleR_0.7.4.tar.gz'
Copy-Item -LiteralPath $taskArchive -Destination $taskCopiedArchive
if ((Get-FileHash -LiteralPath $taskCopiedArchive -Algorithm SHA256).Hash -ne $taskArchiveHash) { throw 'Archive copy mismatch.' }
@{source_manifest=(Join-Path $taskStudy 'check01/checked-source-md5.csv'); archive=$taskArchive;
  SHA256=$taskArchiveHash; check_root=$taskCheckRoot; reused_completed_build=$true;
  no_supervising_R_interpreter=$true; initial_foreign_R=0; initial_RAM_bytes=$taskMemory.FreePhysicalMemory*1KB;
  initial_disk_bytes=$taskDisk.Free} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $taskStudy 'check02-input.json')
foreach ($taskName in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) {
  [Environment]::SetEnvironmentVariable($taskName,'1','Process')
}
foreach ($taskName in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) {
  [Environment]::SetEnvironmentVariable($taskName,$null,'Process')
}
$taskLibraries = @(Get-Content -LiteralPath (Join-Path $taskStudy 'compatibility02.json') -Raw | ConvertFrom-Json).libraries
$env:R_LIBS = $taskLibraries -join ';'
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
$env:MAKEFLAGS = '-j1'
$taskRoot = Start-Process -FilePath 'C:/Program Files/R/R-4.6.1/bin/x64/Rcmd.exe' -ArgumentList @('check','multiScaleR_0.7.4.tar.gz','--as-cran','--no-manual') -WorkingDirectory $taskCheckRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskStudy 'check02.log') -RedirectStandardError (Join-Path $taskStudy 'check02.err')
$taskRootIdentity = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskRoot.Id)"
$taskRootIdentity | Select-Object ProcessId,ParentProcessId,CreationDate,ExecutablePath,CommandLine | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $taskStudy 'check02-process.json')
$taskOperational = $true
$taskReason = ''
while (-not $taskRoot.HasExited) {
  $taskAll = @(Get-CimInstance Win32_Process -ErrorAction Stop)
  $taskOwned = @($taskRoot.Id)
  do {
    $taskN = $taskOwned.Count
    $taskOwned = @($taskOwned + @($taskAll | Where-Object { $_.ParentProcessId -in $taskOwned }).ProcessId | Sort-Object -Unique)
  } while ($taskOwned.Count -gt $taskN)
  $taskR = @($taskAll | Where-Object { $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$' })
  $taskStates = @()
  foreach ($taskProcess in $taskR) {
    $taskState = 'unknown'
    try {
      $taskLoaded = @(Get-Process -Id $taskProcess.ProcessId -Module -ErrorAction Stop)
      if (@($taskLoaded | Where-Object { $_.ModuleName -ieq 'R.dll' }).Count) {
        $taskState = 'R_interpreter'
      } elseif ($taskProcess.ExecutablePath -ieq 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' -and
                $taskProcess.ProcessId -in $taskOwned) {
        # The installed front executable forwards to bin/x64/Rscript.exe.
        # Absence of R.dll verifies that this known owned shim is not an R host.
        $taskState = 'verified_Rscript_forwarder'
      }
    } catch {
      if (-not (Get-Process -Id $taskProcess.ProcessId -ErrorAction SilentlyContinue)) {
        $taskState = 'exited_during_inventory'
      }
    }
    $taskStates += @{PID=$taskProcess.ProcessId;parent=$taskProcess.ParentProcessId;creation=$taskProcess.CreationDate;
      path=$taskProcess.ExecutablePath;command=$taskProcess.CommandLine;state=$taskState;owned=($taskProcess.ProcessId -in $taskOwned)}
  }
  $taskActual = @($taskStates | Where-Object { $_.state -eq 'R_interpreter' }).Count
  $taskUnknown = @($taskStates | Where-Object { $_.state -eq 'unknown' }).Count
  $taskForeign = @($taskStates | Where-Object { -not $_.owned -and $_.state -ne 'exited_during_inventory' }).Count
  $taskMemory = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
  $taskDisk = Get-PSDrive C -ErrorAction Stop
  @{time=(Get-Date -Format o); actual_interpreters=$taskActual; unknown=$taskUnknown;
    foreign=$taskForeign; raw_R_processes=$taskR.Count; states=$taskStates;
    RAM_bytes=$taskMemory.FreePhysicalMemory*1KB; disk_bytes=$taskDisk.Free} |
    ConvertTo-Json -Depth 5 -Compress | Add-Content -LiteralPath (Join-Path $taskStudy 'check02-resources.jsonl')
  if ($taskUnknown -or $taskForeign -or $taskActual -gt 4 -or
      -not $taskMemory.FreePhysicalMemory -or -not $taskDisk.Free -or
      $taskMemory.FreePhysicalMemory*1KB -le 12GB -or $taskDisk.Free -le 5GB -or
      [DateTimeOffset]::UtcNow -ge [DateTimeOffset]'2026-10-29T04:35:00-04:00') {
    $taskOperational = $false
    $taskReason = 'Stopped on unknown/foreign inventory, interpreter/resource budget, or deadline guard.'
    # The root and descendants were observed from this exact registered tree.
    foreach ($taskOld in @($taskAll | Where-Object { $_.ProcessId -in $taskOwned })) {
      $taskFresh = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskOld.ProcessId)" -ErrorAction SilentlyContinue
      if ($taskFresh -and $taskFresh.CreationDate -eq $taskOld.CreationDate -and
          $taskFresh.ExecutablePath -eq $taskOld.ExecutablePath -and $taskFresh.CommandLine -eq $taskOld.CommandLine) {
        Stop-Process -Id $taskFresh.ProcessId -Force -ErrorAction SilentlyContinue
      }
    }
    break
  }
  Start-Sleep -Seconds 15
  $taskRoot.Refresh()
}
$taskRoot.WaitForExit()
$taskExit = $taskRoot.ExitCode
$taskOut = Join-Path $taskStudy 'check02'
New-Item -ItemType Directory -Path $taskOut | Out-Null
Copy-Item -LiteralPath (Join-Path $taskStudy 'check01/checked-source-md5.csv') -Destination (Join-Path $taskOut 'checked-source-md5.csv')
$taskLog = Join-Path $taskCheckRoot 'multiScaleR.Rcheck/00check.log'
if (Test-Path -LiteralPath $taskLog) { Copy-Item -LiteralPath $taskLog -Destination (Join-Path $taskOut 'R-CMD-check.log') }
foreach ($taskName in @('testthat.Rout','testthat.Rout.fail')) {
  $taskTest = Join-Path $taskCheckRoot "multiScaleR.Rcheck/tests/$taskName"
  if (Test-Path -LiteralPath $taskTest) { Copy-Item -LiteralPath $taskTest -Destination (Join-Path $taskOut "installed-$taskName") }
}
$taskStable = $true
foreach ($taskFile in $taskManifest) {
  if ((Get-FileHash -LiteralPath (Join-Path $taskRepo $taskFile.file) -Algorithm MD5).Hash.ToLowerInvariant() -ne $taskFile.md5) { $taskStable = $false }
}
$taskOK = (Test-Path -LiteralPath $taskLog) -and ((Get-Content -LiteralPath $taskLog -Raw) -match '(?m)^Status: OK\s*$')
$taskPassed = $taskExit -eq 0 -and $taskOK -and $taskStable -and $taskOperational
@{passed=$taskPassed;exit_code=$taskExit;check_status_OK=$taskOK;source_stable=$taskStable;
  operational_passed=$taskOperational;reason=$taskReason;check_root=$taskCheckRoot;
  archive_SHA256=$taskArchiveHash;completed_at=(Get-Date -Format o)} |
  ConvertTo-Json | Set-Content -LiteralPath (Join-Path $taskStudy 'check02.json')
if (-not $taskPassed) { exit 1 }
