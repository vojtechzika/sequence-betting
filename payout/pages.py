from otree.api import *
import random
from mpl.config import Constants as M


def fmt_czk(x):
    """
    Format float CZK with decimal comma and 2 decimals, plus 'Kč'.
    """
    s = f"{float(x):.2f}".replace('.', ',')
    return f"{s} Kč"


class Payout(Page):

    def _find_sequences_players(self, all_players):
        out = []
        for pl in all_players:
            if hasattr(pl, 'trial_id') and hasattr(pl, 'seq') and hasattr(pl, 'stake') and hasattr(pl, 'side'):
                out.append(pl)
        return out

    def _find_mpl_player(self, all_players):
        for pl in all_players:
            if hasattr(pl, 'choice_1') and hasattr(pl, 'choice_2'):
                return pl
        return None

    def _ensure_computed(self):
        me = self.player
        par = self.participant

        # prevent redraw on refresh
        if me.field_maybe_none('paid_seq_round') is not None and me.field_maybe_none('paid_mpl_row') is not None:
            return

        all_players = par.get_players()

        # ---------- SEQUENCES ----------
        seq_players = self._find_sequences_players(all_players)
        seq_players = [pl for pl in seq_players if pl.field_maybe_none('round_earnings') is not None]

        if not seq_players:
            raise Exception("No sequence rounds with round_earnings found.")

        chosen_seq = random.choice(seq_players)

        me.paid_seq_round = chosen_seq.round_number
        me.paid_seq_amount = float(chosen_seq.round_earnings)

        # admin Payments
        chosen_seq.payoff = cu(chosen_seq.round_earnings)

        # ---------- MPL ----------
        mpl_player = self._find_mpl_player(all_players)
        if mpl_player is None:
            raise Exception("MPL player not found.")

        i = par.vars.get('mpl_index_to_pay')
        choice_field = par.vars.get('mpl_choice_to_pay')

        if i is None or choice_field is None:
            i = random.randint(1, M.num_choices)
            choice_field = f'choice_{i}'
            par.vars['mpl_index_to_pay'] = i
            par.vars['mpl_choice_to_pay'] = choice_field

        option = getattr(mpl_player, choice_field, None)
        if option not in ('A', 'B'):
            option = 'A'

        r = random.randint(1, M.num_choices)

        if option == 'A':
            amount = M.lottery_a_hi if r <= i else M.lottery_a_lo
        else:
            amount = M.lottery_b_hi if r <= i else M.lottery_b_lo

        me.paid_mpl_row = i
        me.paid_mpl_choice = option
        me.paid_mpl_amount = float(amount)

        # admin Payments
        mpl_player.payoff = cu(amount)

        # total (points)
        me.total_paid = me.paid_seq_amount + me.paid_mpl_amount

        # audit trail
        par.vars['seq_paid_round'] = me.paid_seq_round
        par.vars['seq_paid_amount'] = me.paid_seq_amount
        par.vars['mpl_paid_row'] = i
        par.vars['mpl_paid_option'] = option
        par.vars['mpl_draw'] = r
        par.vars['mpl_paid_amount'] = me.paid_mpl_amount

    def vars_for_template(self):
        self._ensure_computed()

        rw_per_point = float(self.session.config.get('real_world_currency_per_point', 1))
        showup_czk = float(self.session.config.get('participation_fee', 0))

        paid_seq_czk = self.player.paid_seq_amount * rw_per_point
        paid_mpl_czk = self.player.paid_mpl_amount * rw_per_point
        variable_czk = paid_seq_czk + paid_mpl_czk
        total_czk = showup_czk + variable_czk

        return dict(
            showup_fee=fmt_czk(showup_czk),
            variable_fee=fmt_czk(variable_czk),
            total_fee=fmt_czk(total_czk),

            paid_seq_round=self.player.paid_seq_round,
            paid_seq_amount=fmt_czk(paid_seq_czk),

            paid_mpl_row=self.player.paid_mpl_row,
            paid_mpl_amount=fmt_czk(paid_mpl_czk),
        )


page_sequence = [Payout]