param(
    [ValidateSet("Editor", "Sandbox")]
    [string]$RunTarget = "Editor",

    [switch]$SkipRun
)

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "LizeralEngine Build Script" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

$currentDir = Get-Location
$expectedDir = "C:\Lizeral Engine\LizeralEngine0.0.1"
if ($currentDir.Path -ne $expectedDir) {
    Write-Host "Warning: Recommended to run from $expectedDir" -ForegroundColor Yellow
    Write-Host "Current directory: $($currentDir.Path)" -ForegroundColor Yellow
    $confirm = Read-Host "Continue? (y/n)"
    if ($confirm -ne "y") {
        exit
    }
}

$engineRoot = Join-Path $expectedDir "engine"
$buildDir = Join-Path $engineRoot "build"
$binDir = Join-Path $engineRoot "bin"

if (-not (Test-Path $buildDir)) {
    New-Item -ItemType Directory -Path $buildDir | Out-Null
}

if (-not (Test-Path (Join-Path $buildDir "CMakeCache.txt"))) {
    Write-Host "`n[0/3] Configuring CMake..." -ForegroundColor Green
    cmake -S $engineRoot -B $buildDir
    if ($LASTEXITCODE -ne 0) {
        Write-Host "CMake configure failed!" -ForegroundColor Red
        exit 1
    }
}

Write-Host "`n[1/3] Building project..." -ForegroundColor Green
    cmake --build $buildDir

if ($LASTEXITCODE -eq 0) {
    Write-Host "Build successful!" -ForegroundColor Green
} else {
    Write-Host "Build failed!" -ForegroundColor Red
    exit 1
}

Write-Host "`n[2/3] Checking executables..." -ForegroundColor Green
$allPassed = $true

$physicsExePath = Join-Path $binDir "LizeralPhysicsSandbox.exe"
if (Test-Path $physicsExePath) {
    $physicsFileInfo = Get-Item $physicsExePath
    Write-Host "Physics Sandbox created: $physicsExePath" -ForegroundColor Green
    Write-Host "  File size: $([math]::Round($physicsFileInfo.Length / 1KB, 2)) KB" -ForegroundColor Gray
} else {
    Write-Host "Physics Sandbox not found: $physicsExePath" -ForegroundColor Red
    $allPassed = $false
}

$editorExePath = Join-Path $binDir "LizeralEditor.exe"
if (Test-Path $editorExePath) {
    $editorFileInfo = Get-Item $editorExePath
    Write-Host "Lizeral Editor created: $editorExePath" -ForegroundColor Green
    Write-Host "  File size: $([math]::Round($editorFileInfo.Length / 1KB, 2)) KB" -ForegroundColor Gray
} else {
    Write-Host "Lizeral Editor not found: $editorExePath" -ForegroundColor Red
    $allPassed = $false
}

if (-not $allPassed) {
    exit 1
}

if ($SkipRun) {
    Write-Host "`n[3/3] Run step skipped." -ForegroundColor Yellow
    exit 0
}

$exeToRun = if ($RunTarget -eq "Sandbox") { $physicsExePath } else { $editorExePath }

Write-Host "`n[3/3] Running $RunTarget..." -ForegroundColor Green
Write-Host "----------------------------------------" -ForegroundColor DarkGray
Set-Location $binDir
& $exeToRun
$exitCode = $LASTEXITCODE
Write-Host "----------------------------------------" -ForegroundColor DarkGray

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "Build Summary:" -ForegroundColor Cyan
Write-Host "  Target: $RunTarget" -ForegroundColor Gray
Write-Host "  Exit code: $exitCode" -ForegroundColor Gray

if ($exitCode -eq 0) {
    Write-Host "$RunTarget closed cleanly!" -ForegroundColor Green
} else {
    Write-Host "$RunTarget exited with errors." -ForegroundColor Red
}
Write-Host "========================================" -ForegroundColor Cyan

Write-Host "`nExamples:" -ForegroundColor Cyan
Write-Host "  .\build.ps1 -RunTarget Editor" -ForegroundColor Gray
Write-Host "  .\build.ps1 -RunTarget Sandbox" -ForegroundColor Gray
Write-Host "  .\build.ps1 -SkipRun" -ForegroundColor Gray
