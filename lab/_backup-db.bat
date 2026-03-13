@echo off
setlocal

set SCRIPT_DIR=%~dp0
set REPO_ROOT=%SCRIPT_DIR%..

call "%SCRIPT_DIR%lab-config.bat"

for %%I in ("%REPO_ROOT%") do set PROJECT=%%~nxI

set DB=%REPO_ROOT%\%OTREE_DIR%\db.sqlite3

for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyyy-MM-dd_HH-mm-ss"') do set TS=%%i

set TARGET=%BACKUP_ROOT%\%PROJECT%\%TS%

mkdir "%TARGET%" >nul 2>&1

if not exist "%DB%" (
    echo ERROR: database not found:
    echo %DB%
    pause
    exit /b 1
)

copy "%DB%" "%TARGET%\db.sqlite3" >nul
if errorlevel 1 (
    echo ERROR: backup copy failed.
    echo Source: %DB%
    echo Target: %TARGET%\db.sqlite3
    pause
    exit /b 1
)

echo Project: %PROJECT% > "%TARGET%\session-log.txt"
echo Backup time: %TS% >> "%TARGET%\session-log.txt"
echo Source DB: %DB% >> "%TARGET%\session-log.txt"
echo Backup path: %TARGET%\db.sqlite3 >> "%TARGET%\session-log.txt"

echo.
echo Backup created:
echo %TARGET%\db.sqlite3
echo.