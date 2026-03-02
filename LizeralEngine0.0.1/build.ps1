# LizeralEngine Build Script
# Usage: .\build.ps1

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "LizeralEngine Build Script" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# Check if we're in the right directory
$currentDir = Get-Location
$expectedDir = "C:\Lizeral Engine\LizeralEngine0.0.1"
if ($currentDir.Path -ne $expectedDir) {
    Write-Host "Warning: Recommended to run from $expectedDir" -ForegroundColor Yellow
    Write-Host "Current directory: $($currentDir.Path)" -ForegroundColor Yellow
    $confirm = Read-Host "Continue? (y/n)"
    if ($confirm -ne 'y') {
        exit
    }
}

# Step 1: Build the project
Write-Host "`n[1/3] Building project..." -ForegroundColor Green
cd "engine\build"
cmake --build . --config Release

# Check if build succeeded
if ($LASTEXITCODE -eq 0) {
    Write-Host "Build successful!" -ForegroundColor Green
} else {
    Write-Host "Build failed!" -ForegroundColor Red
    exit 1
}

# Step 2: Check executables (同时检查两个引擎程序)
Write-Host "`n[2/3] Checking executables..." -ForegroundColor Green
$allPassed = $true

# 检查 Sandbox
$physicsExePath = "..\bin\LizeralPhysicsSandbox.exe"
if (Test-Path $physicsExePath) {
    $physicsFileInfo = Get-Item $physicsExePath
    Write-Host "✓ Physics Sandbox created: $physicsExePath" -ForegroundColor Green
} else {
    Write-Host "✗ Physics Sandbox not found: $physicsExePath" -ForegroundColor Red
    $allPassed = $false
}

# 检查 Editor
$editorExePath = "..\bin\LizeralEditor.exe"
if (Test-Path $editorExePath) {
    $editorFileInfo = Get-Item $editorExePath
    Write-Host "✓ Lizeral Editor created:  $editorExePath" -ForegroundColor Green
    Write-Host "  File size: $([math]::Round($editorFileInfo.Length/1KB, 2)) KB" -ForegroundColor Gray
} else {
    Write-Host "✗ Lizeral Editor not found: $editorExePath" -ForegroundColor Red
    $allPassed = $false
}

if (-not $allPassed) {
    exit 1
}

# Step 3: Run the Editor (默认运行编辑器)
Write-Host "`n[3/3] Running Lizeral Editor..." -ForegroundColor Green
Write-Host "----------------------------------------" -ForegroundColor DarkGray
cd "..\bin"
Write-Host "Starting Editor..." -ForegroundColor Cyan

# 【关键修改】：这里改为运行 Editor
.\LizeralEditor.exe

$editorExitCode = $LASTEXITCODE
Write-Host "----------------------------------------" -ForegroundColor DarkGray

# Summary
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Build Summary:" -ForegroundColor Cyan
Write-Host "  Editor exit code: $editorExitCode" -ForegroundColor Gray

if ($editorExitCode -eq 0) {
    Write-Host "Lizeral Editor closed cleanly!" -ForegroundColor Green
} else {
    Write-Host "Lizeral Editor crashed or had issues." -ForegroundColor Red
}
Write-Host "========================================" -ForegroundColor Cyan

Write-Host "`nNext time, just run: .\build.ps1"