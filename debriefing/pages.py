from otree.api import *

class Debriefing(Page):
    form_model = 'player'
    form_fields = ['used_strategy', 'strategy_text', 'belief_sequences', 'confidence', 'comment']


class Finished(Page):
    pass


page_sequence = [Debriefing, Finished]