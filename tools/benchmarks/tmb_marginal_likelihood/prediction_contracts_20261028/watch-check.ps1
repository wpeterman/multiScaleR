$ErrorActionPreference = 'Stop'
$taskRootRecords = @(Get-Content -LiteralPath (Join-Path $PSScriptRoot 'check01-process-v2.json') -Raw | ConvertFrom-Json)
$taskLog = Join-Path $PSScriptRoot 'check01-resources.jsonl'
if (Test-Path -LiteralPath $taskLog) { throw 'Retain existing resource monitor output.' }
while ($true) {
  try {
    $taskProcesses = @(Get-CimInstance Win32_Process -ErrorAction Stop)
    $taskRoots = @($taskProcesses | Where-Object {
      $taskProcess = $_
      @($taskRootRecords | Where-Object {
        $_.ProcessId -eq $taskProcess.ProcessId -and $_.CommandLine -eq $taskProcess.CommandLine -and
        ([DateTime]$_.CreationDate).ToUniversalTime() -eq $taskProcess.CreationDate.ToUniversalTime()
      }).Count -gt 0
    })
    if (-not $taskRoots.Count) { break }
    $taskOwnedIds = @($taskRoots.ProcessId)
    do {
      $taskOldCount = $taskOwnedIds.Count
      $taskOwnedIds = @($taskOwnedIds + @($taskProcesses | Where-Object {
        $_.ParentProcessId -in $taskOwnedIds
      }).ProcessId | Sort-Object -Unique)
    } while ($taskOwnedIds.Count -gt $taskOldCount)
    $taskR = @($taskProcesses | Where-Object { $_.Name -match '^(Rscript|Rterm|R|Rgui|Rserve|rsession.*)\.exe$' })
    $taskForeign = @($taskR | Where-Object { $_.ProcessId -notin $taskOwnedIds })
    $taskUnknown = @($taskR | Where-Object { -not $_.CreationDate -or -not $_.CommandLine -or -not $_.ExecutablePath })
    $taskMemory = Get-CimInstance Win32_OperatingSystem -ErrorAction Stop
    $taskDisk = Get-PSDrive C -ErrorAction Stop
    $taskRecord = @{time=(Get-Date -Format o); actual_R_count=$taskR.Count;
      foreign_R=@($taskForeign | Select-Object ProcessId,ParentProcessId,CreationDate,CommandLine);
      inventory_complete=($taskUnknown.Count -eq 0 -and $taskMemory.FreePhysicalMemory -gt 0 -and $taskDisk.Free -gt 0);
      RAM_bytes=$taskMemory.FreePhysicalMemory*1KB; disk_bytes=$taskDisk.Free}
    $taskRecord | ConvertTo-Json -Depth 4 -Compress | Add-Content -LiteralPath $taskLog
  } catch {
    @{time=(Get-Date -Format o); inventory_complete=$false; error=$_.Exception.Message} |
      ConvertTo-Json -Compress | Add-Content -LiteralPath $taskLog
    break
  }
  Start-Sleep -Seconds 15
}
