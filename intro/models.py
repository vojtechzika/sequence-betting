from otree.api import *
import random

class C(BaseConstants):
    NAME_IN_URL = 'intro'
    PLAYERS_PER_GROUP = None
    NUM_ROUNDS = 1


class Subsession(BaseSubsession):
    def creating_session(self):
        if self.round_number != 1:
            return

        cfg = self.session.config
        weights = cfg.get('treatment_weights', {'m25': 0.5, 'm19': 0.5})
        multipliers = cfg.get('treatment_multipliers', {'m25': 2.5, 'm19': 1.9})

        if 'treat_counts' not in self.session.vars:
            self.session.vars['treat_counts'] = {k: 0 for k in weights.keys()}
            self.session.vars['treat_total_assigned'] = 0

        counts = self.session.vars['treat_counts']

        for pl in self.get_players():
            if 'treatment' in pl.participant.vars:
                continue

            total = self.session.vars['treat_total_assigned']
            if total == 0:
                treat = random.choices(list(weights.keys()), weights=list(weights.values()), k=1)[0]
            else:
                gaps = {}
                for t, w in weights.items():
                    gaps[t] = w - (counts[t] / total)

                max_gap = max(gaps.values())
                candidates = [t for t, g in gaps.items() if g == max_gap]
                treat = random.choices(candidates, weights=[weights[c] for c in candidates], k=1)[0]

            pl.participant.vars['treatment'] = treat
            pl.participant.vars['m'] = multipliers[treat]

            counts[treat] += 1
            self.session.vars['treat_total_assigned'] += 1


class Group(BaseGroup):
    pass


class Player(BasePlayer):
    consent = models.BooleanField(
        label="Souhlasím s účastí ve studii.",
        choices=[
            [True, "Ano"],
            [False, "Ne"],
        ],
    )

    age = models.IntegerField(
        label="Váš věk:",
        min=18,
        max=99
    )

    sex = models.StringField(
        label="Pohlaví:",
        choices=[
            ["M", "Muž"],
            ["F", "Žena"],
            ["O", "Nechci uvést"],
        ]
    )