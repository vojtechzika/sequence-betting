@echo off
setlocal

set SCRIPT_DIR=%~dp0
for %%I in ("%SCRIPT_DIR%..") do set OTREE_ROOT=%%~fI

call "%SCRIPT_DIR%lab-config.bat"

for %%I in ("%OTREE_ROOT%\..") do set PROJECT=%%~nxI

set DB=%OTREE_ROOT%\db.sqlite3

for /f %%i in ('powershell -NoProfile -Command "Get-Date -Format yyyy-MM-dd_HH-mm-ss"') do set TS=%%i

set TARGET=%BACKUP_ROOT%\%PROJECT%\%TS%

if not exist "%DB%" (
    echo No database found yet, skipping backup.
    echo.
    exit /b 0
)

mkdir "%TARGET%" >nul 2>&1

copy "%DB%" "%TARGET%\db.sqlite3" >nul
if errorlevel 1 (
    echo ERROR: backup failed.
    pause
    exit /b 1
)

echo Project: %PROJECT% > "%TARGET%\session-log.txt"
echo Backup time: %TS% >> "%TARGET%\session-log.txt"
echo Source DB: %DB% >> "%TARGET%\session-log.txt"
echo Backup path: %TARGET%\db.sqlite3 >> "%TARGET%\session-log.txt"

echo Backup created:
echo %TARGET%\db.sqlite3
echo.