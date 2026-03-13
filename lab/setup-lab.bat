@echo off
setlocal

set SCRIPT_DIR=%~dp0
set REPO_ROOT=%SCRIPT_DIR%..

call "%SCRIPT_DIR%lab-config.bat"

where py >nul 2>&1 || (
    echo ERROR: Python launcher not found.
    echo Install Python %PYTHON_VERSION% and make sure the py launcher is available.
    echo.
    pause
    exit /b 1
)

if not exist "%REPO_ROOT%\%OTREE_DIR%\settings.py" (
    echo ERROR: oTree project folder not found:
    echo %REPO_ROOT%\%OTREE_DIR%
    echo.
    pause
    exit /b 1
)

echo Creating virtual environment...
py -%PYTHON_VERSION% -m venv "%REPO_ROOT%\venv"
if errorlevel 1 (
    echo ERROR: Failed to create virtual environment.
    echo.
    pause
    exit /b 1
)

call "%REPO_ROOT%\venv\Scripts\activate.bat"
if errorlevel 1 (
    echo ERROR: Failed to activate virtual environment.
    echo.
    pause
    exit /b 1
)

cd /d "%REPO_ROOT%\%OTREE_DIR%"

if not exist requirements.txt (
    echo ERROR: requirements.txt not found in:
    echo %CD%
    echo.
    pause
    exit /b 1
)

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