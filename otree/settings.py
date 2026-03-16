from os import environ

DEBUG = False

SESSION_CONFIGS = [
    dict(
        name='sb26',
        display_name='Sequence Betting with Optimism (2026)',
        num_demo_participants=1,
        app_sequence=[
            'intro',
            'sequences',
            'mpl',
            'lotr',
            'payout',
            'debriefing',
        ],
        treatment_weights={'m25': 0.75, 'm19': 0.25},
        treatment_multipliers={'m25': 2.5, 'm19': 1.9},
        experiment_version='sb26_v1_1_preregistered',
    ),
]

SESSION_CONFIG_DEFAULTS = dict(
    real_world_currency_per_point=1.00,
    participation_fee=150.00,
    doc="",
)

PARTICIPANT_FIELDS = ['treatment', 'multiplier', 'experiment_version']
SESSION_FIELDS = ['experiment_version']

LANGUAGE_CODE = 'cs'

USE_POINTS = True
POINTS_CUSTOM_NAME = 'ECU'

REAL_WORLD_CURRENCY_CODE = 'CZK'

ADMIN_USERNAME = 'sekvence'
ADMIN_PASSWORD = 'sekvence'

SECRET_KEY = 'lab-secret-key'

DEMO_PAGE_INTRO_HTML = ""

ROOMS = [
    dict(
        name='lab',
        display_name='Laboratory',
        participant_label_file='_rooms/lab_labels.txt',
    ),
]
