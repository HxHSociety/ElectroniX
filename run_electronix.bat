@echo off
setlocal enabledelayedexpansion
REM ============================================================
REM  ElectroniX - build everything, then launch the application
REM
REM  Usage:
REM    run_electronix.bat            Dev mode (debug build + hot-reload UI)
REM    run_electronix.bat release    Optimised build, launches the .exe
REM    run_electronix.bat clean      Remove staged sidecars, then dev mode
REM
REM  Steps:
REM    1. Check prerequisites (cargo, node, npm)
REM    2. Build helper CLIs: gltf_convertor, rpim_pc, rpim_solver
REM    3. Stage them in src-tauri\binaries with the target-triple suffix
REM       (required by tauri.conf.json -> bundle.externalBin)
REM    4. npm install for the frontend (only when node_modules is missing)
REM    5. Launch the app (tauri dev, or build + run the release exe)
REM ============================================================

cd /d "%~dp0"
set "ROOT=%CD%"
set "MODE=dev"
if /i "%~1"=="release" set "MODE=release"
if /i "%~1"=="clean"   set "MODE=clean"

echo.
echo   ElectroniX build + launch
echo   Mode : %MODE%
echo   Root : %ROOT%
echo.

REM ---- 1. Prerequisites --------------------------------------
where cargo >nul 2>nul || ( echo [ERROR] cargo not found. Install Rust from https://rustup.rs & goto :fail )
where node  >nul 2>nul || ( echo [ERROR] node not found. Install Node.js 20+ from https://nodejs.org & goto :fail )
where npm   >nul 2>nul || ( echo [ERROR] npm not found. & goto :fail )

REM Host target triple, e.g. x86_64-pc-windows-msvc
for /f "tokens=2" %%T in ('rustc -vV ^| findstr /b "host:"') do set "TRIPLE=%%T"
if "%TRIPLE%"=="" ( echo [ERROR] Could not detect Rust host triple. & goto :fail )
echo   Triple: %TRIPLE%

if "%MODE%"=="clean" (
    if exist "src-tauri\binaries" rmdir /s /q "src-tauri\binaries"
    set "MODE=dev"
)

set "PROFILE=debug"
set "CARGO_FLAG="
if "%MODE%"=="release" (
    set "PROFILE=release"
    set "CARGO_FLAG=--release"
)

REM ---- 2. Build helper CLIs ----------------------------------
echo.
echo == [1/4] Building helper binaries (%PROFILE%) ==
for %%P in (gltf_convertor rpim_pc rpim_solver) do (
    echo   cargo build -p %%P
    cargo build -p %%P %CARGO_FLAG%
    if errorlevel 1 ( echo [ERROR] Build failed: %%P & goto :fail )
)

REM ---- 3. Stage sidecars -------------------------------------
echo.
echo == [2/4] Staging sidecar binaries ==
if not exist "src-tauri\binaries" mkdir "src-tauri\binaries"
for %%P in (gltf_convertor rpim_pc rpim_solver) do (
    set "SRC=%ROOT%\target\%PROFILE%\%%P.exe"
    if not exist "!SRC!" ( echo [ERROR] Missing !SRC! & goto :fail )
    copy /y "!SRC!" "src-tauri\binaries\%%P-%TRIPLE%.exe" >nul
    echo   staged %%P-%TRIPLE%.exe
)

REM ---- 4. Frontend dependencies ------------------------------
echo.
echo == [3/4] Frontend dependencies ==
if not exist "frontend\node_modules" (
    pushd frontend
    call npm ci
    if errorlevel 1 ( popd & echo [ERROR] npm ci failed & goto :fail )
    popd
) else (
    echo   node_modules present - skipping npm ci
)

REM ---- 5. Launch ---------------------------------------------
echo.
echo == [4/4] Launching ElectroniX ==
set "TAURI=%ROOT%\frontend\node_modules\.bin\tauri.cmd"
if not exist "%TAURI%" ( echo [ERROR] Tauri CLI missing at %TAURI% & goto :fail )

if "%MODE%"=="release" (
    call "%TAURI%" build --no-bundle
    if errorlevel 1 ( echo [ERROR] tauri build failed & goto :fail )
    set "EXE=%ROOT%\target\release\electronix.exe"
    if not exist "!EXE!" ( echo [ERROR] Missing !EXE! & goto :fail )
    echo   Starting !EXE!
    start "" "!EXE!"
) else (
    REM Runs the Vite dev server + Rust backend and opens the app window.
    call "%TAURI%" dev
)

endlocal
exit /b 0

:fail
echo.
echo   Build/launch aborted.
endlocal
pause
exit /b 1
