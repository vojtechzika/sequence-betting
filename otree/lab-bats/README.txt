oTree Lab Scripts
=================

These scripts help run oTree experiments in the laboratory.
The system is designed so that the operator only needs to run two scripts.

Project structure
-----------------

repo
│
└─ otree
   ├─ settings.py
   ├─ requirements.txt
   ├─ db.sqlite3
   ├─ venv
   ├─ _rooms
   │    └─ lab_labels.txt
   └─ lab-bats
        ├─ setup-lab.bat
        ├─ RUN-SESSION.bat
        ├─ END-SESSION.bat
        ├─ lab-config.bat
        ├─ _start-otree.bat
        └─ _backup-db.bat


Operator workflow
-----------------

1. First-time setup

   Run:

   setup-lab.bat

   This will:
   - verify Python installation
   - create the Python virtual environment in

     otree\venv

   - install dependencies from requirements.txt

   This step only needs to be done once on each computer.


2. Start the experiment

   Double-click:

   RUN-SESSION.bat

   This will:

   - verify the Python environment exists
   - verify settings.py exists
   - verify the participant room file exists
   - create a pre-session backup of the database
   - start the oTree server
   - open the oTree admin page


3. Run the experiment with participants.


4. End the experiment

   Double-click:

   END-SESSION.bat

   This will:

   - stop the oTree server
   - create a post-session backup of the database
   - open the backup folder


Participant access
------------------

Participants join using the room URL:

http://localhost:8000/room/lab

Participant labels are loaded from:

otree\_rooms\lab.txt


Backups
-------

Backups are stored in:

D:\otree-backups\<project-name>\timestamp

Each backup contains:

- db.sqlite3
- session-log.txt

Two backups are created automatically:

1) before the session starts
2) after the session ends


Configuration
-------------

Project settings are defined in:

lab-config.bat

Important settings include:

- PYTHON_VERSION      (Python used to create venv)
- OTREE_PORT          (server port)
- OTREE_ADMIN_PASSWORD
- ROOM_NAME
- ROOM_LABEL_FILE
- BACKUP_ROOT


Internal scripts
----------------

The following scripts are internal and should not be run directly:

- _start-otree.bat
- _backup-db.bat


Important
---------

The operator should normally only run:

RUN-SESSION.bat
END-SESSION.bat

All other scripts are part of the internal infrastructure.
