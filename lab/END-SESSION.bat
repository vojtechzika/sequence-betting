@echo off
setlocal

set SCRIPT_DIR=%~dp0
set REPO_ROOT=%SCRIPT_DIR%..

call "%SCRIPT_DIR%lab-config.bat"

echo Stopping oTree server...
taskkill /FI "WINDOWTITLE eq oTree Server*" /T /F >nul 2>&1

call "%SCRIPT_DIR%_backup-db.bat"
if errorlevel 1 (
    echo ERROR: Backup failed.
    echo.
    pause
    exit /b 1
)

for %%I in ("%REPO_ROOT%") do set PROJECT=%%~nxI

if not exist "%BACKUP_ROOT%\%PROJECT%" (
    echo ERROR: Backup folder not found:
    echo %BACKUP_ROOT%\%PROJECT%
    echo.
    pause
    exit /b 1
)

echo Opening backup folder...
explorer "%BACKUP_ROOT%\%PROJECT%"

pause