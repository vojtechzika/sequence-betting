@echo off
setlocal

REM ------------------------------------------------
REM Start experiment session
REM ------------------------------------------------

set SCRIPT_DIR=%~dp0
for %%I in ("%SCRIPT_DIR%..") do set OTREE_ROOT=%%~fI

call "%SCRIPT_DIR%lab-config.bat"

set VENV_ACTIVATE=%OTREE_ROOT%\venv\Scripts\activate.bat
set ROOM_FILE=%OTREE_ROOT%\%ROOM_LABEL_FILE%

if not exist "%VENV_ACTIVATE%" (
    echo ERROR: Python environment not found.
    echo.
    echo Please run setup-lab.bat first.
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

if not exist "%ROOM_FILE%" (
    echo ERROR: room label file not found:
    echo %ROOM_FILE%
    echo.
    echo Check ROOM_LABEL_FILE in lab-config.bat or create the file.
    echo.
    pause
    exit /b 1
)

tasklist /FI "WINDOWTITLE eq oTree Server*" | find /I "cmd.exe" >nul
if %errorlevel%==0 (
    echo ERROR: oTree server already running.
    echo.
    echo Close the existing server before starting a new session.
    echo.
    pause
    exit /b 1
)

echo.
echo Creating pre-session backup...
call "%SCRIPT_DIR%_backup-db.bat" >nul 2>&1

echo.
echo Starting experiment...
echo.

call "%SCRIPT_DIR%_start-otree.bat"
if errorlevel 1 (
    echo ERROR: Failed to start oTree server.
    echo.
    pause
    exit /b 1
)

timeout /t 3 >nul

echo Opening admin page...
start http://localhost:%OTREE_PORT%

echo.
echo ======================================
echo Session running
echo Use END-SESSION.bat when finished
echo ======================================
echo.
