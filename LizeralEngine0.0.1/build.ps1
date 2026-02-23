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
Write-Host "`n[1/4] Building project..." -ForegroundColor Green
cd "engine\build"
cmake --build . --config Release

# Check if build succeeded
if ($LASTEXITCODE -eq 0) {
    Write-Host "Build successful!" -ForegroundColor Green
} else {
    Write-Host "Build failed!" -ForegroundColor Red
    exit 1
}

# Step 2: Check executables
Write-Host "`n[2/4] Checking executables..." -ForegroundColor Green

# Check Physics Sandbox
$physicsExePath = "..\bin\LizeralPhysicsSandbox.exe"
if (Test-Path $physicsExePath) {
    $physicsFileInfo = Get-Item $physicsExePath
    Write-Host "✓ Physics Sandbox created: $physicsExePath" -ForegroundColor Green
    Write-Host "  File size: $([math]::Round($physicsFileInfo.Length/1KB, 2)) KB" -ForegroundColor Gray
} else {
    Write-Host "✗ Physics Sandbox not found: $physicsExePath" -ForegroundColor Red
    exit 1
}

# Check Rendering Test
$renderingExePath = "..\bin\LizeralRenderingTest.exe"
if (Test-Path $renderingExePath) {
    $renderingFileInfo = Get-Item $renderingExePath
    Write-Host "✓ Rendering Test created: $renderingExePath" -ForegroundColor Green
    Write-Host "  File size: $([math]::Round($renderingFileInfo.Length/1KB, 2)) KB" -ForegroundColor Gray
} else {
    Write-Host "✗ Rendering Test not found: $renderingExePath" -ForegroundColor Red
    exit 1
}

# Step 3: Run Physics Sandbox
Write-Host "`n[3/4] Running Physics Sandbox..." -ForegroundColor Green
Write-Host "----------------------------------------" -ForegroundColor DarkGray
cd "..\bin"
Write-Host "Starting Physics Sandbox..." -ForegroundColor Cyan
.\LizeralPhysicsSandbox.exe
$physicsExitCode = $LASTEXITCODE
Write-Host "----------------------------------------" -ForegroundColor DarkGray

# Step 4: Run Rendering Test
Write-Host "`n[4/4] Running Rendering Test..." -ForegroundColor Green
Write-Host "----------------------------------------" -ForegroundColor DarkGray
Write-Host "Starting Rendering Test..." -ForegroundColor Cyan
#.\LizeralRenderingTest.exe
$renderingExitCode = $LASTEXITCODE
Write-Host "----------------------------------------" -ForegroundColor DarkGray

# Summary
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Build Summary:" -ForegroundColor Cyan
Write-Host "  Physics Sandbox exit code: $physicsExitCode" -ForegroundColor Gray
Write-Host "  Rendering Test exit code: $renderingExitCode" -ForegroundColor Gray

if ($physicsExitCode -eq 0 -and $renderingExitCode -eq 0) {
    Write-Host "Both programs completed successfully!" -ForegroundColor Green
} elseif ($physicsExitCode -eq 0) {
    Write-Host "Physics Sandbox completed, Rendering Test had issues" -ForegroundColor Yellow
} elseif ($renderingExitCode -eq 0) {
    Write-Host "Rendering Test completed, Physics Sandbox had issues" -ForegroundColor Yellow
} else {
    Write-Host "Both programs had issues" -ForegroundColor Red
}
Write-Host "========================================" -ForegroundColor Cyan

Write-Host "`nNext time, just run: .\build.ps1"
