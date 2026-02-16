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

# Step 2: Check executable
Write-Host "`n[2/3] Checking executable..." -ForegroundColor Green
$exePath = "..\bin\LizeralPhysicsSandbox.exe"
if (Test-Path $exePath) {
    $fileInfo = Get-Item $exePath
    Write-Host "Executable created: $exePath" -ForegroundColor Green
    Write-Host "File size: $([math]::Round($fileInfo.Length/1KB, 2)) KB" -ForegroundColor Gray
} else {
    Write-Host "Executable not found: $exePath" -ForegroundColor Red
    exit 1
}

# Step 3: Run the program
Write-Host "`n[3/3] Running program..." -ForegroundColor Green
Write-Host "----------------------------------------" -ForegroundColor DarkGray
cd "..\bin"
.\LizeralPhysicsSandbox.exe
$runExitCode = $LASTEXITCODE
Write-Host "----------------------------------------" -ForegroundColor DarkGray

# Summary
Write-Host "`n========================================" -ForegroundColor Cyan
if ($runExitCode -eq 0) {
    Write-Host "Build and run completed!" -ForegroundColor Green
} else {
    Write-Host "Program exited with code: $runExitCode" -ForegroundColor Yellow
}
Write-Host "========================================" -ForegroundColor Cyan

Write-Host "`nNext time, just run: .\build.ps1"