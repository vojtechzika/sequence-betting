@echo off
setlocal

REM ------------------------------------------------
REM Start experiment session
REM ------------------------------------------------

set SCRIPT_DIR=%~dp0
set REPO_ROOT=%SCRIPT_DIR%..

call "%SCRIPT_DIR%lab-config.bat"

REM ------------------------------------------------
REM Check that virtual environment exists
REM ------------------------------------------------

if not exist "%REPO_ROOT%\venv\Scripts\activate.bat" (
    echo ERROR: Python environment not found.
    echo.
    echo Please run setup-lab.bat first.
    echo.
    pause
    exit /b 1
)

REM ------------------------------------------------
REM Prevent starting server twice
REM ------------------------------------------------

tasklist /FI "WINDOWTITLE eq oTree Server*" | find /I "cmd.exe" >nul
if %errorlevel%==0 (
    echo ERROR: oTree server already running.
    echo.
    echo Close the existing server before starting a new session.
    echo.
    pause
    exit /b
)

REM ------------------------------------------------
REM Pre-session backup
REM ------------------------------------------------

echo.
echo Creating pre-session backup...
call "%SCRIPT_DIR%_backup-db.bat" >nul 2>&1

echo.
echo Starting experiment...
echo.

REM ------------------------------------------------
REM Start server
REM ------------------------------------------------

call "%SCRIPT_DIR%_start-otree.bat"

REM Wait for server startup
timeout /t 3 >nul

REM ------------------------------------------------
REM Open useful pages
REM ------------------------------------------------

echo Opening admin page...
start http://localhost:%OTREE_PORT%/admin

echo Opening participant room...
start http://localhost:%OTREE_PORT%/room/%ROOM_NAME%/

echo.
echo ======================================
echo Session running
echo Use END-SESSION.bat when finished
echo ======================================
echo.

pause