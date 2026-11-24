# 会话记录备份：D:\PI -> F:\PI_sessions
# 会话 = pi 的 *.jsonl 对话记录（含当前长期会话）；配置 = D:\PI\agent\models*.json
# 用法： powershell -ExecutionPolicy Bypass -File tools\sync_sessions_to_F.ps1 [-DryRun]
param([switch]$DryRun)

$ErrorActionPreference = 'Continue'
$Sessions = 'D:\PI'
$Dst      = 'F:\PI_sessions'
$CfgSrc   = 'D:\PI\agent'
$CfgDst   = 'F:\PI_sessions\_agent_config'

if (-not (Test-Path $Sessions)) { throw "missing: $Sessions" }
New-Item -ItemType Directory -Force -Path $Dst | Out-Null
New-Item -ItemType Directory -Force -Path $CfgDst | Out-Null

Write-Host "== sessions -> $Dst"
$argsS = @($Sessions, $Dst, '*.jsonl', '/R:1', '/W:1', '/NFL', '/NDL', '/NP')
if ($DryRun) { $argsS += '/L' }
& robocopy @argsS | Select-Object -Last 6

Write-Host "== agent config (excluding auth/secrets) -> $CfgDst"
Copy-Item -Force (Join-Path $CfgSrc 'models.json')        $CfgDst -ErrorAction SilentlyContinue
Copy-Item -Force (Join-Path $CfgSrc 'models-store.json')  $CfgDst -ErrorAction SilentlyContinue
Copy-Item -Force (Join-Path $CfgSrc 'office-adapter.json') $CfgDst -ErrorAction SilentlyContinue
Get-ChildItem (Join-Path $CfgSrc 'models.json.bak*') -ErrorAction SilentlyContinue |
    ForEach-Object { if (-not $DryRun) { Copy-Item -Force $_.FullName $CfgDst } }

Write-Host "== manifest"
if (-not $DryRun) {
    $files = Get-ChildItem $Dst -Filter *.jsonl | Sort-Object Length -Descending
    $lines = @("session backup $(Get-Date -Format o)", "source $Sessions", "count $($files.Count)", "")
    foreach ($f in $files) { $lines += ("{0,12}  {1}" -f $f.Length, $f.Name) }
    $lines | Set-Content -Encoding UTF8 (Join-Path $Dst 'SESSION_MANIFEST.txt')
    Write-Host "   $($files.Count) sessions listed"
}
Write-Host "done. (auth.json with secrets is NOT backed up)"
