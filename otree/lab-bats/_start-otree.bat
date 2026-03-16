@echo off
setlocal

set SCRIPT_DIR=%~dp0
for %%I in ("%SCRIPT_DIR%..") do set OTREE_ROOT=%%~fI

call "%SCRIPT_DIR%lab-config.bat"

set VENV_ACTIVATE=%OTREE_ROOT%\venv\Scripts\activate.bat
set ROOM_FILE=%OTREE_ROOT%\%ROOM_LABEL_FILE%

if not exist "%VENV_ACTIVATE%" (
    echo ERROR: Python environment not found.
    echo Run setup-lab.bat first.
    exit /b 1
)

if not exist "%OTREE_ROOT%\settings.py" (
    echo ERROR: settings.py not found in:
    echo %OTREE_ROOT%
    exit /b 1
)

if not exist "%ROOM_FILE%" (
    echo ERROR: room label file not found:
    echo %ROOM_FILE%
    exit /b 1
)

call "%VENV_ACTIVATE%"
if errorlevel 1 (
    echo ERROR: Failed to activate virtual environment.
    exit /b 1
)

cd /d "%OTREE_ROOT%"

set OTREE_PRODUCTION=1
set OTREE_ADMIN_PASSWORD=%OTREE_ADMIN_PASSWORD%

start "oTree Server" cmd /k "title oTree Server && otree prodserver %OTREE_PORT%"