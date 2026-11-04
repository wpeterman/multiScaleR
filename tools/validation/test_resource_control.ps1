param([string]$Engine, [string]$Output, [string]$Case)
$ErrorActionPreference='Stop'
. $Engine
$taskStudy=$Output
$Plan=Join-Path $Output 'plan.txt'
[IO.File]::WriteAllText($Plan,'non-numerical process fixture')
$taskPlanHash=(Get-FileHash -LiteralPath $Plan -Algorithm SHA256).Hash
$taskControl=@{launch_until=[DateTimeOffset]::UtcNow.AddHours(1).ToString('o');long_launch_until=[DateTimeOffset]::UtcNow.AddHours(1).ToString('o');worker_stop_at=[DateTimeOffset]::UtcNow.AddHours(1).ToString('o')}
$taskPhaseIsLong=$false
$script:fixtureRoot=[pscustomobject]@{Id=123456789;HasExited=($Case -ne 'running');ExitCode=0}
$script:fixtureRoot | Add-Member ScriptMethod Refresh {$this.HasExited=$true}
$script:fixtureRoot | Add-Member ScriptMethod WaitForExit {}
$script:fixtureIdentity=[pscustomobject]@{ProcessId=123456789;ParentProcessId=1;CreationDate=(Get-Date);ExecutablePath='fixture_Rscript.exe';CommandLine='never executed'}
$script:fixtureStops=@()
function Start-Process {param($FilePath,$ArgumentList,$WorkingDirectory,$WindowStyle,[switch]$PassThru,$RedirectStandardOutput,$RedirectStandardError) $script:fixtureRoot}
function Get-CimInstance {param($ClassName,$Filter,$ErrorAction) if ($Filter) {$script:fixtureIdentity}}
function Get-TaskInventory {param($OwnedIds) @{states=@();unknown=[int]($Case -eq 'unknown');owned_interpreters=0;actual_interpreters=0;foreign=0;disk_bytes=40GB;RAM_bytes=16GB;threads=1}}
function Assert-TaskAdmission {Get-TaskInventory @()}
function Start-Sleep {param($Seconds) if ($Case -ne 'running') {throw 'Fast phase must not sleep'}}
function Stop-Process {param($Id,[switch]$Force,$ErrorAction) if ($Id -ne 123456789) {throw 'Foreign process stop attempted'}; $script:fixtureStops += $Id}
$taskFailed=$false
try {Invoke-TaskPhase 'fixture' 'fixture_Rscript.exe' @('never-executed') $Output} catch {$taskFailed=$true;if ($Case -ne 'unknown') {throw}}
$taskResources=@(Get-Content -LiteralPath (Join-Path $Output 'fixture-resources.jsonl'))
$taskExit=Get-Content -LiteralPath (Join-Path $Output 'fixture-exit.json') -Raw | ConvertFrom-Json
if ($taskResources.Count -ne 1) {throw 'Expected exactly one real inventory call in the process fixture'}
if ($Case -eq 'unknown') {
 if (-not $taskFailed -or $taskExit.operational_passed -or $script:fixtureStops.Count -ne 1) {throw 'Unknown inventory must retain its existing owned-root stop and failure'}
} elseif ($taskFailed -or -not $taskExit.operational_passed -or $script:fixtureStops.Count) {throw 'Completed known phase must pass without stopping any process'}
Write-Output "$Case process-recording fixture passed; no R process launched."
