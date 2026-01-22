from otree.api import *


class C(BaseConstants):
    NAME_IN_URL = 'debriefing'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1


class Subsession(BaseSubsession):
    pass


class Group(BaseGroup):
    pass


class Player(BasePlayer):
    used_strategy = models.StringField(
        label="Použil(a) jste nějakou strategii při rozhodování?",
        choices=[['yes', 'Ano'], ['no', 'Ne']],
        widget=widgets.RadioSelect,
    )
    strategy_text = models.LongStringField(
        label="Pokud ano, popište ji stručně:",
        blank=True,
    )
    belief_sequences = models.LongStringField(
        label="Měl(a) jste pocit, že některé sekvence jsou pravděpodobnější než jiné? Proč?",
        blank=True,
    )
    confidence = models.IntegerField(
        label="Jak jste si byl(a) jistý(á) svými rozhodnutími? (0–10)",
        min=0, max=10,
    )
    comment = models.LongStringField(
        label="Další komentář (volitelné):",
        blank=True,
    )