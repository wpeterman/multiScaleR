$ErrorActionPreference = 'Stop'
$taskStudy = $PSScriptRoot
$taskRepo = (Get-Location).Path
. (Join-Path $taskStudy 'resource-controller.ps1')
if (Test-Path -LiteralPath (Join-Path $taskStudy 'panel-started.json')) { throw 'Panel already admitted; inspect existing identities/outcomes.' }
Assert-TaskAdmission
$taskPlan = Get-Content -LiteralPath (Join-Path $taskStudy 'RUNLIST.json') -Raw | ConvertFrom-Json
$taskRscriptShimHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' -Algorithm SHA256).Hash
$taskRexeHash = (Get-FileHash -LiteralPath 'C:/Program Files/R/R-4.6.1/bin/x64/R.exe' -Algorithm SHA256).Hash
$env:R_LIBS = (@(Get-Content -LiteralPath (Join-Path $taskStudy '../functional_precision_20261028/compatibility01.json') -Raw | ConvertFrom-Json).libraries) -join ';'
$env:MAKEFLAGS = '-j1'
foreach ($taskEnv in @('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($taskEnv,'1','Process') }
foreach ($taskEnv in @('LANG','LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME')) { [Environment]::SetEnvironmentVariable($taskEnv,$null,'Process') }
& 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'verify-study.py')
if ($LASTEXITCODE -ne 0) { throw 'Frozen identities did not pass.' }
$taskController = Get-CimInstance Win32_Process -Filter "ProcessId=$PID" -ErrorAction Stop
Write-TaskJson @{started_at=(Get-Date -Format o);planned_runs=4;worker_limit=1;
  controller=($taskController | Select-Object ProcessId,ParentProcessId,CreationDate,ExecutablePath,CommandLine);
  protocol_sha256=$taskPlan.protocol_sha256} (Join-Path $taskStudy 'panel-started.json')
$taskUsed = 0.0
$taskResults = @()
try {
  foreach ($taskRun in $taskPlan.runs) {
    $taskPhaseLimit = [Math]::Min($taskPlan.fit_wall_limit_seconds, $taskPlan.cumulative_wall_limit_seconds - $taskUsed)
    if ($taskPhaseLimit -le 0) { throw 'Cumulative wall budget exhausted.' }
    & 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'verify-study.py')
    if ($LASTEXITCODE -ne 0) { throw 'Source/runtime/producer/evidence identity drift.' }
    Write-TaskJson @{time=(Get-Date -Format o);state='running';active_run=$taskRun.id;completed=$taskResults;
      cumulative_seconds=$taskUsed;planned_runs=4;remaining=(@($taskPlan.runs.id | Where-Object {$_ -notin $taskResults.id}))} (Join-Path $taskStudy 'checkpoint-running.json')
    Invoke-TaskPhase $taskRun.id 'C:/Program Files/R/R-4.6.1/bin/x64/Rscript.exe' @('--vanilla',('"' + (Join-Path $taskStudy 'run-fit.R') + '"'),$taskRun.id) $taskRepo
    & 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' (Join-Path $taskStudy 'verify-study.py')
    if ($LASTEXITCODE -ne 0) { throw 'Post-run source/runtime/producer/evidence identity drift.' }
    $taskExit = Get-Content -LiteralPath (Join-Path $taskStudy ($taskRun.id + '-exit.json')) -Raw | ConvertFrom-Json
    $taskUsed += $taskExit.elapsed_seconds
    $taskOutput = Get-Content -LiteralPath (Join-Path $taskStudy ($taskRun.id + '.json')) -Raw | ConvertFrom-Json
    $taskResults += @{id=$taskRun.id;completed=$taskOutput.completed;qualified=$taskOutput.reference_qualified;elapsed_seconds=$taskExit.elapsed_seconds}
  }
  Write-TaskJson @{time=(Get-Date -Format o);state='all_four_completed';completed=$taskResults;planned_runs=4;
    cumulative_seconds=$taskUsed;paired_agreement='pending analysis'} (Join-Path $taskStudy 'panel-complete.json')
  Write-TaskJson @{time=(Get-Date -Format o);state='all_four_completed';completed=$taskResults;planned_runs=4;
    cumulative_seconds=$taskUsed;remaining=@()} (Join-Path $taskStudy 'checkpoint-running.json')
} catch {
  Write-TaskJson @{time=(Get-Date -Format o);state='attention_required';error=$_.Exception.Message;completed=$taskResults;
    cumulative_seconds=$taskUsed;planned_runs=4;active_run=$taskRun.id;
    remaining=(@($taskPlan.runs.id | Where-Object {$_ -notin $taskResults.id}))} (Join-Path $taskStudy 'checkpoint-running.json')
  throw
}
