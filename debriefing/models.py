from otree.api import *


class C(BaseConstants):
    NAME_IN_URL = 'survey'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1


class Subsession(BaseSubsession):
    pass


class Group(BaseGroup):
    pass


class Player(BasePlayer):

    # --- A) 50/50 transparency + hypothetical sequences ---

    belief_independence = models.IntegerField(
        label=(
            "1. Do jaké míry jste věřil(a), že každý hod má pravděpodobnost 50 % \
            na Hlavu a 50 % na Orel bez ohledu na zobrazenou historii šesti hodů? \
            (0 = vůbec ne, 10 = zcela)"
        ),
        min=0, max=10,
    )

    reliance_on_sequence = models.IntegerField(
        label=(
            "2. Do jaké míry jste při rozhodování zohledňoval(a) zobrazenou historii šesti hodů, \
            i když jste věděl(a), že mince je spravedlivá? \
            (0 = vůbec ne, 10 = velmi silně)"
        ),
        min=0, max=10,
    )

    perceived_realism_history = models.IntegerField(
        label=(
            "3. Do jaké míry jste vnímal(a) zobrazené historie šesti hodů jako „reálné historie“ \
            (tj. že skutečně předcházely dalšímu hodu), oproti tomu, že šlo jen o různé možné příklady? \
            (0 = jen možné příklady, 10 = reálná historie)"
        ),
        min=0, max=10,
    )

    # --- B) action-seeking / engagement / fatigue ---

    action_seeking = models.IntegerField(
        label=(
            "4. Jak často jste sázel(a) hlavně proto, že Vám přišlo lepší „něco udělat“ \
            než nesázet? (0 = nikdy, 10 = velmi často)"
        ),
        min=0, max=10,
    )

    enjoyment = models.IntegerField(
        label=(
            "5. Nakolik Vás úloha bavila? (0 = vůbec, 10 = velmi)"
        ),
        min=0, max=10,
    )

    fatigue = models.IntegerField(
        label=(
            "6. Nakolik jste se během úlohy nudil(a) nebo byl(a) unavený(á)? (0 = vůbec, 10 = velmi)"
        ),
        min=0, max=10,
    )

    # --- C) risk preference triangulation ---

    self_risk_tolerance = models.IntegerField(
        label=(
            "7. Jak byste celkově popsal(a) svou ochotu podstupovat finanční riziko? \
            (0 = velmi nerad(a) riskuji, 10 = velmi rád(a) riskuji)"
        ),
        min=0, max=10,
    )

    # --- D) strategy + open text ---

    used_strategy = models.StringField(
        label="8. Používal(a) jste při rozhodování nějaké pravidlo nebo strategii?",
        choices=[
            ['yes', 'Ano'],
            ['no', 'Ne'],
            ['unsure', 'Nejsem si jistý(á)'],
        ],
        widget=widgets.RadioSelect,
    )

    strategy_text = models.LongStringField(
        label="9. Pokud ano, stručně ji popište: (volitelné)",
        blank=True,
    )

    comment = models.LongStringField(
        label="10. Máte k této části studie jakýkoli další komentář? (volitelné)",
        blank=True,
    )