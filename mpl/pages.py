from otree.api import Currency as c
from ._builtin import Page, WaitPage
from .config import Constants



# helper to inject lottery labels everywhere we need them
def _lottery_ctx():
    return dict(
        lottery_a_lo=c(Constants.lottery_a_lo),
        lottery_a_hi=c(Constants.lottery_a_hi),
        lottery_b_lo=c(Constants.lottery_b_lo),
        lottery_b_hi=c(Constants.lottery_b_hi),
    )

class Instructions(Page):
    def is_displayed(self):
        return self.subsession.round_number == 1

    def vars_for_template(self):
        ctx = dict(num_choices=len(self.participant.vars['mpl_choices']))
        ctx.update(_lottery_ctx())
        return ctx

class Decision(Page):
    form_model = 'player'

    def get_form_fields(self):
        form_fields = [list(t) for t in zip(*self.participant.vars['mpl_choices'])][1]
        if Constants.one_choice_per_page:
            page = self.subsession.round_number
            return [form_fields[page - 1]]
        return form_fields

    def vars_for_template(self):
        total = len(self.participant.vars['mpl_choices'])
        page = self.subsession.round_number
        progress = page / total * 100

        if Constants.one_choice_per_page:
            ctx = dict(
                page=page,
                total=total,
                progress=progress,
                choices=[self.participant.vars['mpl_choices'][page - 1]],
            )
        else:
            ctx = dict(choices=self.participant.vars['mpl_choices'])

        ctx.update(_lottery_ctx())
        return ctx

    def before_next_page(self):
        round_number = self.subsession.round_number
        form_fields = [list(t) for t in zip(*self.participant.vars['mpl_choices'])][1]
        indices = [list(t) for t in zip(*self.participant.vars['mpl_choices'])][0]
        index = indices[round_number - 1]

        if Constants.one_choice_per_page:
            current_choice = getattr(self.player, form_fields[round_number - 1])
            self.participant.vars['mpl_choices_made'][index - 1] = current_choice
            #if index == self.participant.vars['mpl_index_to_pay']:
            #    self.player.set_payoffs()
            if index == self.participant.vars['mpl_index_to_pay']:
                self.player.payoff = c(0)
            if round_number == Constants.num_choices:
                self.player.set_consistency()
                self.player.set_switching_row()
        else:
            for j, choice in zip(indices, form_fields):
                choice_i = getattr(self.player, choice)
                self.participant.vars['mpl_choices_made'][j - 1] = choice_i
            #self.player.set_payoffs()
            self.player.payoff = c(0)
            self.player.set_consistency()
            self.player.set_switching_row()

class Results(Page):
    def is_displayed(self):
        return (self.subsession.round_number == Constants.num_rounds) if Constants.one_choice_per_page else True

    def vars_for_template(self):
        choices = [list(t) for t in zip(*self.participant.vars['mpl_choices'])]
        indices = choices[0]
        index_to_pay = self.participant.vars['mpl_index_to_pay']
        round_to_pay = indices.index(index_to_pay) + 1
        if Constants.one_choice_per_page:
            ctx = dict(
                choice_to_pay=[self.participant.vars['mpl_choices'][round_to_pay - 1]],
                option_to_pay=self.player.in_round(round_to_pay).option_to_pay,
                payoff=self.player.in_round(round_to_pay).payoff,
            )
        else:
            ctx = dict(
                choice_to_pay=[self.participant.vars['mpl_choices'][round_to_pay - 1]],
                option_to_pay=self.player.option_to_pay,
                payoff=self.player.payoff,
            )
        ctx.update(_lottery_ctx())
        return ctx

page_sequence = [Decision]
if Constants.instructions:
    page_sequence.insert(0, Instructions)
if Constants.results:
    page_sequence.append(Results)