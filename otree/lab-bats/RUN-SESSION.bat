@echo off
setlocal

set SCRIPT_DIR=%~dp0
for %%I in ("%SCRIPT_DIR%..") do set OTREE_ROOT=%%~fI

call "%SCRIPT_DIR%lab-config.bat"

set VENV_DIR=%OTREE_ROOT%\venv
set ROOM_FILE=%OTREE_ROOT%\%ROOM_LABEL_FILE%

where py >nul 2>&1 || (
    echo ERROR: Python launcher not found.
    echo Install Python %PYTHON_VERSION% and make sure the py launcher is available.
    echo.
    pause
    exit /b 1
)

if not exist "%OTREE_ROOT%\settings.py" (
    echo ERROR: settings.py not found in:
    echo %OTREE_ROOT%
    echo.
    pause
    exit /b 1
)

if not exist "%OTREE_ROOT%\requirements.txt" (
    echo ERROR: requirements.txt not found in:
    echo %OTREE_ROOT%
    echo.
    pause
    exit /b 1
)

if not exist "%ROOM_FILE%" (
    echo ERROR: room label file not found:
    echo %ROOM_FILE%
    echo.
    echo Check ROOM_LABEL_FILE in lab-config.bat or create the file.
    echo.
    pause
    exit /b 1
)

echo Creating virtual environment...
py -%PYTHON_VERSION% -m venv "%VENV_DIR%"
if errorlevel 1 (
    echo ERROR: Failed to create virtual environment.
    echo.
    pause
    exit /b 1
)

call "%VENV_DIR%\Scripts\activate.bat"
if errorlevel 1 (
    echo ERROR: Failed to activate virtual environment.
    echo.
    pause
    exit /b 1
)

cd /d "%OTREE_ROOT%"

python -m pip install --upgrade pip
if errorlevel 1 (
    echo ERROR: Failed to upgrade pip.
    echo.
    pause
    exit /b 1
)

pip install -r requirements.txt
if errorlevel 1 (
    echo ERROR: Failed to install dependencies from requirements.txt.
    echo.
    pause
    exit /b 1
)

echo.
echo Setup complete.
echo.
pause