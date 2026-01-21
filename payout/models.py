from otree.api import *

class C(BaseConstants):
    NAME_IN_URL = 'payout'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1

class Subsession(BaseSubsession): pass
class Group(BaseGroup): pass

class Player(BasePlayer):
    paid_seq_round = models.IntegerField()
    paid_seq_amount = models.CurrencyField()

    paid_mpl_row = models.IntegerField()
    paid_mpl_choice = models.StringField()
    paid_mpl_amount = models.CurrencyField()

    total_paid = models.CurrencyField()