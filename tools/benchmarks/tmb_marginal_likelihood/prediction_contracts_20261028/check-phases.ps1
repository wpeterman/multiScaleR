$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
$taskRepo = (Get-Location).Path
$taskName = 'check03'
function Write-TaskJson($Value, $Path) {
  $taskTemp = $Path + '.tmp'
  $Value | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $taskTemp
  Move-Item -LiteralPath $taskTemp -Destination $Path
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
      states=$taskStates;RAM_bytes=$taskM.FreePhysicalMemory*1KB;disk_bytes=$taskD.Free} |
      ConvertTo-Json -Depth 6 -Compress | Add-Content -LiteralPath (Join-Path $taskStudy "$Phase-resources.jsonl")
    if ($taskUnknown -or $taskForeign -or $taskCount -gt 4 -or -not $taskM.FreePhysicalMemory -or -not $taskD.Free -or
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
  $taskResult = @{operational_passed=$taskGood;exit_code=$taskRoot.ExitCode;completed_at=(Get-Date -Format o)}
  Write-TaskJson $taskResult (Join-Path $taskStudy "$Phase-exit.json")
  if (-not $taskGood -or $taskRoot.ExitCode -ne 0) { throw "Validation phase failed: $Phase" }
}
if (Test-Path -LiteralPath (Join-Path $taskStudy "$taskName-input.json")) { throw 'Retain existing attempt.' }
Assert-TaskAdmission
$taskManifest = Import-Csv -LiteralPath (Join-Path $taskStudy 'check01/checked-source-md5.csv')
foreach ($taskFile in $taskManifest) {
  if ((Get-FileHash -LiteralPath (Join-Path $taskRepo $taskFile.file) -Algorithm MD5).Hash.ToLowerInvariant() -ne $taskFile.md5) { throw "Source drift: $($taskFile.file)" }
}
$taskPrevious = (Get-Content -LiteralPath (Join-Path $taskStudy 'check01/check-location.txt')).Trim()
$taskArchive = Join-Path $taskPrevious 'check/multiScaleR_0.7.4.tar.gz'
$taskSHA = (Get-FileHash -LiteralPath $taskArchive -Algorithm SHA256).Hash
$taskCheckRoot = Join-Path $env:TEMP ('msr-contract-phases-' + [Guid]::NewGuid().ToString('N').Substring(0,12))
New-Item -ItemType Directory -Path $taskCheckRoot | Out-Null
Copy-Item -LiteralPath $taskArchive -Destination (Join-Path $taskCheckRoot 'multiScaleR_0.7.4.tar.gz')
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
Write-TaskJson @{archive=$taskArchive;SHA256=$taskSHA;check_root=$taskCheckRoot;reused_completed_build=$true;
  tests_run_separately=$true;Rscript_forwarder_SHA256=$taskRscriptShimHash;R_BATCH_forwarder_SHA256=$taskRexeHash;
  coverage='Standard vignette-enabled R CMD check without nested tests, followed by all unchanged archive tests against its installed package'} (Join-Path $taskStudy "$taskName-input.json")
$taskLibraries = @(Get-Content -LiteralPath (Join-Path $taskStudy 'compatibility02.json') -Raw | ConvertFrom-Json).libraries
$env:R_LIBS = $taskLibraries -join ';'
$env:RSTUDIO_PANDOC = 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
$env:MAKEFLAGS = '-j1'
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
Invoke-TaskPhase $taskName 'C:/Program Files/R/R-4.6.1/bin/x64/Rcmd.exe' @('check','multiScaleR_0.7.4.tar.gz','--as-cran','--no-manual','--no-tests') $taskCheckRoot
$taskCheckLog = Join-Path $taskCheckRoot 'multiScaleR.Rcheck/00check.log'
if ((Get-Content -LiteralPath $taskCheckLog -Raw) -notmatch '(?m)^Status: OK\s*$') { throw 'Package check is not clean.' }
$taskTestDir = Join-Path $taskCheckRoot 'exact-archive-tests'
& 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'extract-tests.py') $taskArchive $taskTestDir
if ($LASTEXITCODE -ne 0) { throw 'Exact test extraction failed.' }
$env:R_LIBS = (@((Join-Path $taskCheckRoot 'multiScaleR.Rcheck')) + @($taskLibraries)) -join ';'
$env:NOT_CRAN = 'false'
[Environment]::SetEnvironmentVariable('_R_CHECK_PACKAGE_NAME_','multiScaleR','Process')
Invoke-TaskPhase 'installed03' 'C:/Program Files/R/R-4.6.1/bin/x64/Rterm.exe' @('--vanilla','--no-echo','--no-save','--no-restore','-f','testthat.R') $taskTestDir
$taskStable = $true
foreach ($taskFile in $taskManifest) {
  if ((Get-FileHash -LiteralPath (Join-Path $taskRepo $taskFile.file) -Algorithm MD5).Hash.ToLowerInvariant() -ne $taskFile.md5) { $taskStable = $false }
}
$taskOut = Join-Path $taskStudy $taskName
New-Item -ItemType Directory -Path $taskOut | Out-Null
Copy-Item -LiteralPath $taskCheckLog -Destination (Join-Path $taskOut 'R-CMD-check.log')
Copy-Item -LiteralPath (Join-Path $taskStudy 'installed03.log') -Destination (Join-Path $taskOut 'installed-tests.log')
Copy-Item -LiteralPath (Join-Path $taskStudy 'check01/checked-source-md5.csv') -Destination (Join-Path $taskOut 'checked-source-md5.csv')
Write-TaskJson @{passed=$taskStable;source_stable=$taskStable;check_status_OK=$true;tests_complete=$true;tests_run_separately=$true;
  archive_SHA256=$taskSHA;check_root=$taskCheckRoot;installed_test_directory=$taskTestDir;completed_at=(Get-Date -Format o)} (Join-Path $taskStudy "$taskName.json")
if (-not $taskStable) { throw 'Live source drifted during validation.' }
