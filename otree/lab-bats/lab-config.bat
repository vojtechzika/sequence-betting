@echo off

REM ------------------------------------------------
REM Project configuration for lab scripts
REM ------------------------------------------------

REM Python version used for venv
set PYTHON_VERSION=3.10

REM Port for server
set OTREE_PORT=8000

REM Auth Level
set OTREE_AUTH_LEVEL=STUDY

REM Admin password
set OTREE_ADMIN_PASSWORD=1324

REM Room name used in experiment
set ROOM_NAME=lab

REM Relative path inside otree folder to participant label file
set ROOM_LABEL_FILE=_rooms\lab_labels.txt

REM Backup location
set BACKUP_ROOT=D:\otree-backups
