from os import environ

DEBUG = False

SESSION_CONFIGS = [
    dict(
        name='sb26',
        display_name='Sequence Betting with Optimism (2026)',
        num_demo_participants=1,
        app_sequence=[
        'intro', # General Intro + Consent
        'sequences', # Sequence Betting Task
        'mpl', # Holt-Laury
        'lotr', # LOT-R Optimism
        'payout', # Payoff screen
        'debriefing'
        ],  

        # Treatment controls
        treatment_weights={'m25': 0.75, 'm19': 0.25},      
        treatment_multipliers={'m25': 2.5, 'm19': 1.9},

        experiment_version = 'sb26_v1_preregistered',
    ),
]



# if you set a property in SESSION_CONFIG_DEFAULTS, it will be inherited by all configs
# in SESSION_CONFIGS, except those that explicitly override it.
# the session config can be accessed from methods in your apps as self.session.config,
# e.g. self.session.config['participation_fee']

SESSION_CONFIG_DEFAULTS = dict(
    real_world_currency_per_point=1.00,   # 1 ECU = 1 CZK (adjust if needed)
    participation_fee=150.00,
    doc=""
)

PARTICIPANT_FIELDS = ['treatment', 'multiplier', 'experiment_version']
SESSION_FIELDS = ['experiment_version']

LANGUAGE_CODE = 'cs'

# Experimental units
USE_POINTS = True
POINTS_CUSTOM_NAME = 'ECU'


# Real payment currency (only for final payout)
REAL_WORLD_CURRENCY_CODE = 'CZK'
#REAL_WORLD_CURRENCY_PER_POINT = 1

ADMIN_USERNAME = 'admin'
ADMIN_PASSWORD = environ.get('OTREE_ADMIN_PASSWORD')

DEMO_PAGE_INTRO_HTML = ""

SECRET_KEY = '4416572140507'

ROOMS = [
    dict(
        name='lab',
        display_name='Laboratory',
        participant_label_file='_rooms/lab_labels.txt',  # optional
    ),
]


