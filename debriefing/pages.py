from otree.api import *


class Survey(Page):
    form_model = 'player'
    form_fields = [
        'belief_independence',
        'reliance_on_sequence',
        'perceived_realism_history',
        'action_seeking',
        'enjoyment',
        'fatigue',
        'self_risk_tolerance',
        'used_strategy',
        'strategy_text',
        'comment',
    ]


class Finished(Page):
    pass


page_sequence = [Survey, Finished]