param(
    [ValidateSet("Editor", "Sandbox")]
    [string]$RunTarget = "Editor",

    [switch]$SkipRun,

    [switch]$Clean
)

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "LizeralEngine Build Script" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

$projectRoot = $PSScriptRoot
$engineRoot = Join-Path $projectRoot "engine"
$buildDir = Join-Path $engineRoot "build-fast"
$binDir = Join-Path $engineRoot "bin"
$toolsRoot = Join-Path $projectRoot ".tools"

$cmakeExe = Join-Path $toolsRoot "cmake-3.31.6-windows-x86_64/bin/cmake.exe"
$ninjaExe = Join-Path $toolsRoot "ninja/ninja.exe"
$qtRoot = Join-Path $toolsRoot "Qt/6.10.0/mingw_64"
$mingwRoot = Join-Path $toolsRoot "Qt/Tools/mingw1310_64"
$gccExe = Join-Path $mingwRoot "bin/gcc.exe"
$gxxExe = Join-Path $mingwRoot "bin/g++.exe"
$vulkanRoot = Join-Path $toolsRoot "VulkanSDK"
$vulkanInclude = Join-Path $vulkanRoot "Include"
$vulkanLibrary = "C:/Windows/System32/vulkan-1.dll"
$windeployqtExe = Join-Path $qtRoot "bin/windeployqt.exe"

$requiredTools = @($cmakeExe, $ninjaExe, $gccExe, $gxxExe, $windeployqtExe)
foreach ($tool in $requiredTools) {
    if (-not (Test-Path $tool)) {
        Write-Host "Missing local build tool: $tool" -ForegroundColor Red
        Write-Host "Please let CodeBuddy run the dependency setup again." -ForegroundColor Yellow
        exit 1
    }
}

$env:Path = "$(Join-Path $qtRoot 'bin');$(Join-Path $mingwRoot 'bin');$(Split-Path $cmakeExe);$(Split-Path $ninjaExe);$env:Path"
$env:VULKAN_SDK = $vulkanRoot

if ($Clean -and (Test-Path $buildDir)) {
    Write-Host "`n[0/4] Cleaning build directory..." -ForegroundColor Green
    Remove-Item -Recurse -Force $buildDir
}

Write-Host "`n[1/4] Configuring CMake..." -ForegroundColor Green
& $cmakeExe `
    -S $engineRoot `
    -B $buildDir `
    -G Ninja `
    -DCMAKE_BUILD_TYPE=Debug `
    "-DCMAKE_PREFIX_PATH=$qtRoot" `
    "-DCMAKE_C_COMPILER=$gccExe" `
    "-DCMAKE_CXX_COMPILER=$gxxExe" `
    "-DVulkan_INCLUDE_DIR=$vulkanInclude" `
    "-DVulkan_LIBRARY=$vulkanLibrary"

if ($LASTEXITCODE -ne 0) {
    Write-Host "CMake configure failed!" -ForegroundColor Red
    exit 1
}

Write-Host "`n[2/4] Building project..." -ForegroundColor Green
& $cmakeExe --build $buildDir --config Debug -j 8
if ($LASTEXITCODE -ne 0) {
    Write-Host "Build failed!" -ForegroundColor Red
    exit 1
}
Write-Host "Build successful!" -ForegroundColor Green

Write-Host "`n[3/4] Checking executables and deploying Qt runtime..." -ForegroundColor Green
$physicsExePath = Join-Path $binDir "LizeralPhysicsSandbox.exe"
$editorExePath = Join-Path $binDir "LizeralEditor.exe"

$allPassed = $true
foreach ($exePath in @($physicsExePath, $editorExePath)) {
    if (Test-Path $exePath) {
        $fileInfo = Get-Item $exePath
        Write-Host "Created: $exePath" -ForegroundColor Green
        Write-Host "  File size: $([math]::Round($fileInfo.Length / 1KB, 2)) KB" -ForegroundColor Gray
    } else {
        Write-Host "Missing executable: $exePath" -ForegroundColor Red
        $allPassed = $false
    }
}

if (-not $allPassed) {
    exit 1
}

& $windeployqtExe $editorExePath --dir $binDir --compiler-runtime | Out-Host
if ($LASTEXITCODE -ne 0) {
    Write-Host "windeployqt failed!" -ForegroundColor Red
    exit 1
}

if ($SkipRun) {
    Write-Host "`n[4/4] Run step skipped." -ForegroundColor Yellow
    exit 0
}

$exeToRun = if ($RunTarget -eq "Sandbox") { $physicsExePath } else { $editorExePath }

Write-Host "`n[4/4] Running $RunTarget..." -ForegroundColor Green
Write-Host "----------------------------------------" -ForegroundColor DarkGray
Push-Location $binDir
try {
    & $exeToRun
    $exitCode = $LASTEXITCODE
} finally {
    Pop-Location
}
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
Write-Host "  .\build.ps1 -Clean -SkipRun" -ForegroundColor Gray

exit $exitCode
