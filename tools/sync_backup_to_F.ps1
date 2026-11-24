# RFAM -> F: 备份同步（Windows PowerShell）
#
# 策略（三级）：
#   1) 代码 + 文档 + 小产物：镜像到 F:\realFluidAllMach（排除 build*、.git、大场文件）
#   2) 关键大场（交付网格 + 最新冷流场）：白名单复制到 F:\realFluidAllMach_fields
#   3) 每次运行写一份清单到 F:\realFluidAllMach_fields\BACKUP_MANIFEST.txt
#
# 用法： powershell -ExecutionPolicy Bypass -File tools\sync_backup_to_F.ps1
#        powershell -ExecutionPolicy Bypass -File tools\sync_backup_to_F.ps1 -DryRun

param([switch]$DryRun)

$ErrorActionPreference = 'Stop'
$Src = 'D:\realFluidAllMach'
$Dst = 'F:\realFluidAllMach'
$FieldDst = 'F:\realFluidAllMach_fields'

if (-not (Test-Path $Src)) { throw "source missing: $Src" }
# 老版 PowerShell（无 Get-FileHash）也用得上；robocopy 非零退出码是“有文件被复制”，不是错误。

Write-Host "== 1/3 code + docs + small artifacts -> $Dst"
# /MIR 镜像；排除构建目录、git 和大场文件
robocopy $Src $Dst /MIR /R:1 /W:1 /NFL /NDL /NP /XD build build-phase-foundation build-a406-tetra build-diag .git work\fields /XF *.rfamtetf *.vtu *.o *.a *.so @robo | Out-Null

Write-Host "== 2/3 key large fields -> $FieldDst"
New-Item -ItemType Directory -Force -Path $FieldDst | Out-Null
$whitelist = @(
    'work\mesh_rebuild\delivery\mesh\fluid_volume.rfamtet',
    'work\coldflow_B_verified\benchmark_mesh\fluid_volume.rfamtet',
    'work\coldflow_B_verified\mesh\fluid_volume.rfamtet',
    'work\a406_rootcause\property_table_op2.bin',
    'work\a406_rootcause\property_table_b.bin',
    'work\a406_rootcause\aroute_a44_final.rfamtetf',
    'work\a406_rootcause\aroute_a45_s90.rfamtetf',
    'work\a406_rootcause\aroute_a44_s40.rfamtetf',
    'work\a406_rootcause\aroute_a45_s90.vtu',
    'work\a406_rootcause\aroute_a45_s90.pvsm',
    'work\a406_rootcause\aroute_a45_s90_mach.png'
)
foreach ($rel in $whitelist) {
    $s = Join-Path $Src $rel
    if (-not (Test-Path $s)) { Write-Host "   skip (missing): $rel"; continue }
    $d = Join-Path $FieldDst $rel
    New-Item -ItemType Directory -Force -Path (Split-Path $d) | Out-Null
    if (-not $DryRun) { Copy-Item -Force $s $d }
    Write-Host "   copied: $rel"
}

Write-Host "== 3/3 manifest"
if (-not $DryRun) {
    $lines = @("backup $(Get-Date -Format o)", "source $Src", "mirror $Dst", "fields $FieldDst", "")
    foreach ($rel in $whitelist) {
        $d = Join-Path $FieldDst $rel
        if (Test-Path $d) {
            # certutil exists on every Windows; fall back to a placeholder.
            $h = (& certutil -hashfile $d SHA256 2>$null | Select-Object -Index 1)
            if (-not $h) { $h = 'HASH_UNAVAILABLE' }
            $h = $h -replace '\s', ''
            $n = (Get-Item $d).Length
            $lines += ("{0}  {1,12}  {2}" -f $h, $n, $rel)
        }
    }
    $lines | Set-Content -Encoding UTF8 (Join-Path $FieldDst 'BACKUP_MANIFEST.txt')
    Write-Host "   manifest written"
}
Write-Host "done."
