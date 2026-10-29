$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
$taskRepo = (Get-Location).Path
. (Join-Path $taskStudy 'resource-controller-v2.ps1')
if (Test-Path -LiteralPath (Join-Path $taskStudy 'panel-v3-started.json')) { throw 'Continuation already admitted.' }
Assert-TaskAdmission
$taskPlan = Get-Content -LiteralPath (Join-Path $taskStudy 'RUNLIST.json') -Raw | ConvertFrom-Json
$taskContinuation = Get-Content -LiteralPath (Join-Path $taskStudy 'CONTINUATION_BASELINE.json') -Raw | ConvertFrom-Json
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
$env:R_LIBS = (@(Get-Content -LiteralPath (Join-Path $taskStudy '../functional_precision_20261028/compatibility01.json') -Raw | ConvertFrom-Json).libraries) -join ';'
$env:MAKEFLAGS = '-j1'
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
& 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'continuation.py')
if ($LASTEXITCODE -ne 0) { throw 'Continuation identities changed.' }
$taskController = Get-CimInstance Win32_Process -Filter "ProcessId=$PID" -ErrorAction Stop
Write-TaskJson @{started_at=(Get-Date -Format o);planned_total=4;only_pending=$taskContinuation.remaining;worker_limit=1;
  controller=($taskController | Select-Object ProcessId,ParentProcessId,CreationDate,ExecutablePath,CommandLine)} (Join-Path $taskStudy 'panel-v3-started.json')
$taskUsed = $taskContinuation.prior_wall_seconds
$taskResults = @(@{id='null-a';completed=$true;qualified=$false;elapsed_seconds=796.8090446})
$taskSequence = 0
try {
  foreach ($taskRun in @($taskPlan.runs | Where-Object {$_.id -in $taskContinuation.remaining})) {
    $taskAttempt = $taskRun.id + '-v2'
    if (@(Get-ChildItem -LiteralPath $taskStudy -Filter ($taskAttempt+'*') -File).Count) { throw 'Pending phase already has records; inspect first.' }
    $taskPhaseLimit = [Math]::Min($taskPlan.fit_wall_limit_seconds, $taskPlan.cumulative_wall_limit_seconds - $taskUsed)
    if ($taskPhaseLimit -le 0) { throw 'Cumulative wall budget exhausted.' }
    & 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'continuation.py')
    if ($LASTEXITCODE -ne 0) { throw 'Source/runtime/producer/evidence identity drift.' }
    $taskSequence++
    Write-TaskJson @{time=(Get-Date -Format o);state='running';active_run=$taskRun.id;completed=$taskResults;
      cumulative_seconds=$taskUsed;planned_runs=4} (Join-Path $taskStudy ('checkpoint-v3-{0:D2}.json' -f $taskSequence))
    Invoke-TaskPhase $taskAttempt 'C:/Program Files/R/R-4.6.1/bin/x64/Rscript.exe' @('--vanilla',('"' + (Join-Path $taskStudy 'run-fit-v2.R') + '"'),$taskRun.id) $taskRepo
    & 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'continuation.py')
    if ($LASTEXITCODE -ne 0) { throw 'Post-run identities drifted.' }
    $taskExit = Get-Content -LiteralPath (Join-Path $taskStudy ($taskAttempt+'-exit.json')) -Raw | ConvertFrom-Json
    $taskUsed += $taskExit.elapsed_seconds
    $taskOutput = Get-Content -LiteralPath (Join-Path $taskStudy ($taskAttempt+'.json')) -Raw | ConvertFrom-Json
    $taskResults += @{id=$taskRun.id;completed=$taskOutput.completed;qualified=$taskOutput.reference_qualified;elapsed_seconds=$taskExit.elapsed_seconds}
    $taskSequence++
    Write-TaskJson @{time=(Get-Date -Format o);state='phase_completed';completed=$taskResults;
      cumulative_seconds=$taskUsed;planned_runs=4} (Join-Path $taskStudy ('checkpoint-v3-{0:D2}.json' -f $taskSequence))
  }
  Write-TaskJson @{time=(Get-Date -Format o);state='all_four_completed';completed=$taskResults;
    cumulative_seconds=$taskUsed;planned_runs=4;paired_agreement='pending analysis'} (Join-Path $taskStudy 'panel-v3-complete.json')
} catch {
  $taskSequence++
  Write-TaskJson @{time=(Get-Date -Format o);state='attention_required';error=$_.Exception.Message;
    active_run=$taskRun.id;completed=$taskResults;cumulative_seconds=$taskUsed;planned_runs=4} (Join-Path $taskStudy ('checkpoint-v3-{0:D2}.json' -f $taskSequence))
  throw
}
