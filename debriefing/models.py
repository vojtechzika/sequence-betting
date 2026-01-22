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
        label="Používal(a) jste během úlohy se sekvencemi nějakou strategii?",
        choices=[['yes', "Ano"], ['no', "Ne"], ['unsure', "Nevím / nejsem si jistý(á)"]],
        blank=True
    )
    strategy_text = models.LongStringField(
        label="Pokud ano, můžete ji stručně popsat?",
        blank=True
    )
    belief_sequences = models.StringField(
        label="Měl(a) jste pocit, že některé sekvence zvyšují pravděpodobnost jedné strany?",
        choices=[['yes', "Ano"], ['no', "Ne"], ['unsure', "Nevím / nejsem si jistý(á)"]],
        blank=True
    )
    confidence = models.IntegerField(
        label="Jak jistý(á) jste si byl(a) svými rozhodnutími v úloze se sekvencemi?",
        choices=list(range(1, 8)),
        blank=True
    )
    comment = models.LongStringField(
        label="Máte k experimentu nějaký komentář?",
        blank=True
    )