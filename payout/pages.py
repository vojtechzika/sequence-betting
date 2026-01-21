from otree.api import *
import random
from mpl.config import Constants as M


class Payout(Page):

    def _find_sequences_players(self, all_players):
        # sequences players have fields like trial_id/seq/stake/side
        out = []
        for pl in all_players:
            if hasattr(pl, 'trial_id') and hasattr(pl, 'seq') and hasattr(pl, 'stake') and hasattr(pl, 'side'):
                out.append(pl)
        return out

    def _find_mpl_player(self, all_players):
        # mpl player has choice_1..choice_10 fields
        for pl in all_players:
            if hasattr(pl, 'choice_1') and hasattr(pl, 'choice_2'):
                return pl
        return None

    def _ensure_computed(self):
        me = self.player
        par = self.participant

        # already computed (prevents redraw on refresh)
        if me.field_maybe_none('paid_seq_round') is not None and me.field_maybe_none('paid_mpl_row') is not None:
            return

        all_players = par.get_players()

        # -------- SEQUENCES: pick 1 round ----------
        seq_players = self._find_sequences_players(all_players)
        seq_players = [pl for pl in seq_players if pl.field_maybe_none('round_earnings') is not None]

        if not seq_players:
            raise Exception("No sequences rounds with round_earnings found for this participant.")

        chosen_seq = random.choice(seq_players)
        me.paid_seq_round = chosen_seq.round_number
        me.paid_seq_amount = cu(chosen_seq.round_earnings)
        chosen_seq.payoff = me.paid_seq_amount

        # -------- MPL: pay one row ----------
        mpl_player = self._find_mpl_player(all_players)
        if mpl_player is None:
            raise Exception("Could not locate MPL player (choice_1/choice_2 not found).")

        # If MPL already set these vars during creating_session, use them.
        # Otherwise (e.g., old session), create them now so payout is stable.
        i = par.vars.get('mpl_index_to_pay')
        choice_field = par.vars.get('mpl_choice_to_pay')

        if i is None or choice_field is None:
            i = random.randint(1, M.num_choices)
            choice_field = f'choice_{i}'
            par.vars['mpl_index_to_pay'] = i
            par.vars['mpl_choice_to_pay'] = choice_field

        option = getattr(mpl_player, choice_field, None)
        if option not in ('A', 'B'):
            option = 'A'  # policy if blank/missing

        # Correct random draw: 1..num_choices
        r = random.randint(1, M.num_choices)

        if option == 'A':
            amount = M.lottery_a_hi if r <= i else M.lottery_a_lo
        else:
            amount = M.lottery_b_hi if r <= i else M.lottery_b_lo

        me.paid_mpl_row = i
        me.paid_mpl_choice = option
        me.paid_mpl_amount = cu(amount)

        # store on mpl player for audit/export if field exists
        try:
            mpl_player.round_earnings = float(amount)
        except Exception:
            pass

        mpl_player.payoff = cu(amount)

        me.total_paid = me.paid_seq_amount + me.paid_mpl_amount

        # also store audit in participant vars (helps if you export)
        par.vars['seq_paid_round'] = me.paid_seq_round
        par.vars['seq_paid_amount'] = float(chosen_seq.round_earnings)
        par.vars['mpl_paid_row'] = i
        par.vars['mpl_paid_option'] = option
        par.vars['mpl_draw'] = r
        par.vars['mpl_paid_amount'] = float(amount)

    def vars_for_template(self):
        self._ensure_computed()
        return dict(
            paid_seq_round=self.player.paid_seq_round,
            paid_seq_amount=self.player.paid_seq_amount,
            paid_mpl_row=self.player.paid_mpl_row,
            paid_mpl_choice=self.player.paid_mpl_choice,
            paid_mpl_amount=self.player.paid_mpl_amount,
            total_paid=self.player.total_paid,
        )


page_sequence = [Payout]