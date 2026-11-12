<#
.SYNOPSIS
    FluxCell 一键构建安装脚本（等效 make install）。

.DESCRIPTION
    流程：CMake 配置（preset: windows-default，幂等）→ 构建 → 安装到 dist/。
    产物：dist/bin（FluxCell.exe、spm_cli.exe、运行时 DLL）+ dist/params（默认参数）。

.PARAMETER Config
    构建配置：Release（默认）或 Debug。

.PARAMETER Prefix
    安装目录（相对仓库根），默认 dist/。

.PARAMETER Test
    安装前运行完整单元测试（ctest），失败则中止。

.PARAMETER Zip
    额外用 CPack 生成便携 ZIP 包（build/FluxCell-*-win64.zip）。

.EXAMPLE
    pwsh scripts/install.ps1              # 构建 Release 并安装到 dist/
    pwsh scripts/install.ps1 -Test        # 先跑测试再安装
    pwsh scripts/install.ps1 -Zip         # 安装 + CPack ZIP
#>
param(
    [ValidateSet("Release", "Debug")]
    [string]$Config = "Release",
    [string]$Prefix = "dist",
    [switch]$Test,
    [switch]$Zip
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$buildDir = Join-Path $root "build"

function Invoke-Step {
    param([string]$Desc, [scriptblock]$Action)
    Write-Host "==> $Desc" -ForegroundColor Cyan
    & $Action
    if ($LASTEXITCODE -ne 0) {
        Write-Host "==> $Desc 失败（退出码 $LASTEXITCODE）" -ForegroundColor Red
        exit $LASTEXITCODE
    }
}

# ── 1. 配置（已配置则幂等刷新缓存）──
Invoke-Step "配置（preset: windows-default）" {
    Push-Location $root
    try { cmake --preset windows-default }
    finally { Pop-Location }
}

# ── 2. 构建 ──
Invoke-Step "构建（$Config）" {
    cmake --build $buildDir --config $Config
}

# ── 3.（可选）单元测试 ──
if ($Test) {
    Invoke-Step "单元测试（ctest $Config）" {
        ctest --test-dir $buildDir -C $Config --output-on-failure
    }
}

# ── 4. 安装 ──
$prefixPath = Join-Path $root $Prefix
Invoke-Step "安装到 $Prefix/" {
    cmake --install $buildDir --config $Config --prefix $prefixPath
}

# ── 5.（可选）CPack 便携 ZIP ──
if ($Zip) {
    Invoke-Step "打包（CPack ZIP $Config）" {
        Push-Location $buildDir
        try { cpack -C $Config }
        finally { Pop-Location }
    }
}

# ── 摘要 ──
Write-Host ""
Write-Host "安装完成：$prefixPath" -ForegroundColor Green
if (Test-Path (Join-Path $prefixPath "bin")) {
    Get-ChildItem (Join-Path $prefixPath "bin") |
        Format-Table @{ Label = "bin/"; Expression = { $_.Name } },
                     @{ Label = "大小"; Expression = { "{0:N0} KB" -f ($_.Length / 1KB) } }
}
if (Test-Path (Join-Path $prefixPath "params")) {
    $n = (Get-ChildItem (Join-Path $prefixPath "params")).Count
    Write-Host "params/：$n 个参数文件"
}
