from otree.api import *
from .models import C
import random


# =========================
# INTRO (was sequences_intro)
# =========================

class TaskInstructions(Page):
    def is_displayed(self):
        return self.round_number == 1

    def vars_for_template(self):
        return dict(
            n_rounds=C.NUM_ROUNDS,     
            n_blocks=C.NUM_BLOCKS,
            endowment=C.ENDOWMENT,
            multiplier=self.participant.vars['m'],
        )


class PracticeFixation(Page):
    template_name = 'sequences/Fixation.html'

    def is_displayed(self):
        return self.round_number == 1 and C.FIXATION_MS > 0

    def vars_for_template(self):
        return dict(
            fixation_ms=C.FIXATION_MS,
            is_practice=True,
        )


class PracticeTrial(Page):
    template_name = 'sequences/Trial.html'
    form_model = 'player'
    form_fields = ['aoi_boxes_px', 'viewport_w_px', 'viewport_h_px', 'dpr']

    def is_displayed(self):
        return self.round_number == 1

    def vars_for_template(self):
        N = C.NUM_ROUNDS
        B = C.NUM_BLOCKS

        # keep your original "fake" position for the progress bar
        t_sim = min(20, N)

        base, rem = divmod(N, B)
        sizes = [base + (i < rem) for i in range(B)]

        blk, acc = 1, 0
        for i, sz in enumerate(sizes, start=1):
            if t_sim <= acc + sz:
                blk, pos = i, t_sim - acc
                break
            acc += sz

        return dict(
            is_practice=True,
            t=t_sim,
            n_rounds=N,
            n_blocks=B,
            blk=blk,
            pos=pos,
            trial_id='practice',
            uid='practice',
            block=blk,
            pos_in_block=pos,
            seq="XYXYXY",
            outcomes_left='X',
            outcomes_right='Y',
            nickname_left='Hlava',
            nickname_right='Orel',
            endowment=C.ENDOWMENT,
            multiplier=self.participant.vars['m'],
        )

    def before_next_page(self):
        pv = self.participant.vars
        for k in ['aoi_boxes_px', 'viewport_w_px', 'viewport_h_px', 'dpr']:
            v = self.player.field_maybe_none(k)
            if v not in (None, '', 0):
                pv[k] = v


class Comprehension(Page):
    form_model = 'player'
    form_fields = ['cq_keep_endowment', 'cq_multiplier', 'cq_payment', 'cq_recency', 'cq_toss']


    def is_displayed(self):
        return self.round_number == 1

    def vars_for_template(self):
        E = C.ENDOWMENT
        m = self.participant.vars.get('m')
        m_txt = str(m) if m is not None else "m"

        return dict(
            cq_keep_endowment_choices=[
                (1, "0 ECU"),
                (2, f"{E} ECU"),
                (3, f"{int(E * m)} ECU"),
            ],
            cq_multiplier_choices=[
                (1, "100 ECU"),
                (2, "Nevsazená částka"),
                (3, f"{m_txt}-násobek vsazené částky"),
                (4, f"Nevsazená částka plus {m_txt}-násobek vsazené částky"),
            ],
            cq_payment_choices=[
                (1, "v jednom náhodně vybraném kole"),
                (2, "ve všech kolech"),
                (3, "v posledním kole"),
            ],
            cq_recency_choices=[
                (1, "Vpravo"),
                (2, "Vlevo"),
            ],
            cq_toss_choices=[
                (1, "je v zobrazené sekvenci nejčastější"),
                (2, "je v zobrazené sekvenci poslední"),
                (3, "padne v hodu spravedlivou mincí, který by následoval po zobrazené sekvenci"),
            ],
        )

    def error_message(self, values):
        errors = {}

        correct = {
            'cq_keep_endowment': 2,
            'cq_multiplier': 4,
            'cq_payment': 1,
            'cq_recency': 1,
            'cq_toss': 3,
        }

        missing_msg = "Odpovězte prosím na tuto otázku."
        wrong_msg = "Tato otázka je zodpovězena nesprávně."

        for field, corr in correct.items():
            v = values.get(field)
            if v is None:
                errors[field] = missing_msg
            elif v != corr:
                errors[field] = wrong_msg

        self.participant.vars['cq_failed'] = bool(errors)
        return errors if errors else None

# =========================
# REAL TASK (your original sequences/pages.py)
# =========================

class Fixation(Page):
    def is_displayed(self):
        return C.FIXATION_MS > 0 and self.round_number <= C.NUM_ROUNDS

    def vars_for_template(self):
        return dict(fixation_ms=C.FIXATION_MS)


class Trial(Page):
    form_model = 'player'
    form_fields = ['side', 'stake', 'screen_time_ms']

    def is_displayed(self):
        return self.player.field_maybe_none('side') is None

    def error_message(self, values):
        s = values.get('side')
        b = values.get('stake')
        if s == 'NB':
            if (b or 0) != 0:
                return "Při volbě „Nesázet“ musí být sázka 0."
        else:
            if not s:
                return "Zvolte stranu nebo „Nesázet“."
            if b is None or b < 1:
                return "Zadejte sázku alespoň 1 ECU."
            if b > C.ENDOWMENT:
                return f"Maximální sázka je {C.ENDOWMENT} ECU."

    def vars_for_template(self):
        entry = self.participant.vars['manifest'][self.round_number - 1]
        p = self.player
        p.t = entry['t']
        p.trial_id = entry['trial_id']
        p.seq = entry['seq']
        p.uid = entry['uid']
        p.block = entry['block']
        p.pos_in_block = entry['pos_in_block']
        p.m_used = self.participant.vars['m']
        p.treatment = self.participant.vars.get('treatment', '')

        man = self.participant.vars['manifest']
        block_sizes, cnt = [], 0
        cur_block = man[0]['block'] if man else None
        for row in man:
            if row['block'] == cur_block:
                cnt += 1
            else:
                block_sizes.append(cnt)
                cur_block, cnt = row['block'], 1
        if cnt:
            block_sizes.append(cnt)
        block_sizes_csv = ",".join(str(x) for x in block_sizes)

        if not p.field_maybe_none('button_order'):
            btns = list(dict.fromkeys(list(C.OUTCOMES)))
            random.shuffle(btns)
            p.button_order = ''.join(btns)
        btns = list(p.button_order)

        left = btns[0] if btns else None
        right = btns[1] if len(btns) > 1 else None

        return dict(
            t=p.t,
            trial_id=p.trial_id,
            seq=p.seq,
            uid=p.uid,
            block=p.block,
            pos_in_block=p.pos_in_block,
            blk=p.block,
            pos=p.pos_in_block,

            outcomes_left=left,
            outcomes_right=right,
            nickname_left=C.NICKNAMES.get(left, left) if left else None,
            nickname_right=C.NICKNAMES.get(right, right) if right else None,

            endowment=C.ENDOWMENT,
            multiplier=self.participant.vars['m'],
            n_blocks=C.NUM_BLOCKS,
            n_rounds=C.NUM_ROUNDS,
            block_sizes_csv=block_sizes_csv,
            is_practice=False,
            button_order=p.button_order,
        )

    def before_next_page(self, timeout_happened=False):
        p = self.player
        m = self.participant.vars['m']

        p.m_used = m
        p.treatment = self.participant.vars.get('treatment', '')

        outcomes = list(dict.fromkeys(list(C.OUTCOMES)))
        realized = random.choice(outcomes)
        p.realized = realized

        if p.side == 'NB':
            earnings = C.ENDOWMENT
            p.win = False
        else:
            p.win = (p.side == realized)
            if p.win:
                earnings = (C.ENDOWMENT - p.stake) + m * p.stake
            else:
                earnings = C.ENDOWMENT - p.stake

        p.round_earnings = float(earnings)
        p.payoff = cu(0)


class Break(Page):
    timer_text = ''  # hides oTree's yellow countdown box
    
    def is_displayed(self):
        if C.BREAK_SECONDS <= 0 or self.round_number >= C.NUM_ROUNDS:
            return False
        man = self.participant.vars['manifest']
        return man[self.round_number]['block'] != man[self.round_number - 1]['block']

    def get_timeout_seconds(self):
        return C.BREAK_SECONDS

    def vars_for_template(self):
        return dict(
            blk=self.player.block,
            n_blocks=C.NUM_BLOCKS,
            break_seconds=C.BREAK_SECONDS,
        )


page_sequence = [
    TaskInstructions,
    PracticeFixation,
    PracticeTrial,
    Comprehension,
    Fixation,
    Trial,
    Break
]