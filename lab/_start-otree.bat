@echo off
setlocal

set SCRIPT_DIR=%~dp0
set REPO_ROOT=%SCRIPT_DIR%..

call "%SCRIPT_DIR%lab-config.bat"

if not exist "%REPO_ROOT%\venv\Scripts\activate.bat" (
    echo ERROR: Python environment not found.
    echo Run setup-lab.bat first.
    exit /b 1
)

if not exist "%REPO_ROOT%\%OTREE_DIR%\settings.py" (
    echo ERROR: oTree project folder not found:
    echo %REPO_ROOT%\%OTREE_DIR%
    exit /b 1
)

call "%REPO_ROOT%\venv\Scripts\activate.bat"
if errorlevel 1 (
    echo ERROR: Failed to activate virtual environment.
    exit /b 1
)

cd /d "%REPO_ROOT%\%OTREE_DIR%"

set OTREE_PRODUCTION=1
set OTREE_ADMIN_PASSWORD=%OTREE_ADMIN_PASSWORD%

start "oTree Server" cmd /k "title oTree Server && otree prodserver %OTREE_PORT%"