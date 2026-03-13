oTree Lab Scripts
=================

These scripts help run oTree experiments in the laboratory.

Operator workflow
-----------------

1. First-time setup

   Run:

   setup-lab.bat

   This will:
   - create the Python virtual environment
   - install dependencies from requirements.txt

2. Start the experiment

   Double-click:

   RUN-SESSION.bat

   This will:
   - check the environment
   - create a pre-session backup of the database
   - start the oTree server
   - open the admin page
   - open the participant room page

3. Run the experiment with participants.

4. End the experiment

   Double-click:

   END-SESSION.bat

   This will:
   - stop the oTree server
   - create a backup of the database
   - open the backup folder

Backups
-------

Backups are stored in:

D:\otree-backups\<project-name>\timestamp

Each backup contains:

- db.sqlite3
- session-log.txt

Configuration
-------------

Project settings are defined in:

lab-config.bat

Common settings include:
- OTREE_DIR
- PYTHON_VERSION
- OTREE_PORT
- ROOM_NAME
- BACKUP_ROOT

Internal scripts
----------------

The following scripts are internal and should not be run directly:

- _start-otree.bat
- _backup-db.bat