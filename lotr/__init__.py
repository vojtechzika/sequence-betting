from otree.api import *

doc = "LOT-R survey (raw answers only)."

class C(BaseConstants):
    NAME_IN_URL = 'lotr'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1

    # Choose your Likert scale here
    LIKERT = [
        [0, 'rozhodně nesouhlasím'],
        [1, 'nesouhlasím'],
        [2, 'ani nesouhlasím, ani souhlasím'],
        [3, 'souhlasím'],
        [4, 'rozhodně souhlasím'],
    ]

    # Paste your final Czech item texts here (10 items).
    ITEMS = [
        "V nejistých dobách obvykle očekávám to nejlepší.",
        "Je pro mě snadné se uvolnit.",
        "Pokud se v mém životě něco může pokazit, tak se to pokazí.",
        "Na svou budoucnost se vždy dívám optimisticky.",
        "Velmi si užívám společnost svých přátel.",
        "Je pro mě důležité neustále něco dělat.",
        "Téměř nikdy neočekávám, že věci půjdou tak, jak bych si přál/a.",
        "Nenechám se snadno rozhodit.",
        "Málokdy počítám s tím, že se mi přihodí něco dobrého.",
        "Celkově očekávám, že se mi v životě přihodí více dobrého než špatného.",
    ]


class Subsession(BaseSubsession):
    pass

class Group(BaseGroup):
    pass

class Player(BasePlayer):
    lotr_1  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_2  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_3  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_4  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_5  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_6  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_7  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_8  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_9  = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)
    lotr_10 = models.IntegerField(choices=C.LIKERT, widget=widgets.RadioSelect)


class LOTR(Page):
    form_model = 'player'
    form_fields = [f'lotr_{i}' for i in range(1, 11)]

    @staticmethod
    def vars_for_template(player: Player):
        return dict(items=C.ITEMS)


page_sequence = [LOTR]