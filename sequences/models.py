from otree.api import *
import itertools, random, uuid

class C(BaseConstants):
    """
    Experiment constants:
    - Either use CUSTOM_SEQUENCES, or auto-generate all |OUTCOMES|^LENGTH unique sequences.
    - SEQUENCES is a list of dicts [{trial_id:int, seq:str}, ...].
    - NUM_ROUNDS = len(SEQUENCES).
    """
    NAME_IN_URL = 'sequences'
    PLAYERS_PER_GROUP = None

    # === Experiment parameters ===
    LENGTH = 6                  # length of each sequence (if generated automatically)
    OUTCOMES = "HO"             # symbols for outcomes (e.g. "HT" for coin, "123456" for dice)
    NICKNAMES = {         # what participants actually see
        "H": "Hlava",
        "O": "Orel",
    }
    ENDOWMENT = 100             # initial endowment

    # === Blocking and Breaking ===
    NUM_BLOCKS = 4              # set to 1 for no breaks
    BREAK_SECONDS = 30          # 0 = no break screen; otherwise auto-continue after N seconds

    # Fixation cross duration in ms.
    # If set to 0, the fixation page is skipped entirely.
    FIXATION_MS = 500

    # === Custom sequences ===
    # If non-empty, overrides automatic generation.
    # Example: ["HHHHHH","TTTTTT"] or any strings you want to use.
    CUSTOM_SEQUENCES = []

    # === Sequence set with fixed IDs (1..N) ===
    if CUSTOM_SEQUENCES:
        _seqs = CUSTOM_SEQUENCES
    else:
        _seqs = [''.join(bits) for bits in itertools.product(OUTCOMES, repeat=LENGTH)]

    SEQUENCES = [
        dict(trial_id=i + 1, seq=s)   # fixed mapping: trial_id -> seq
        for i, s in enumerate(_seqs)
    ]

    # === Number of rounds ===
    NUM_ROUNDS = len(SEQUENCES)


class Subsession(BaseSubsession):
    def creating_session(self):
        for p in self.get_players():
            # don't recreate if already present (reconnects)
            if 'manifest' in p.participant.vars:
                continue

            # per-participant randomized order (without mutating constants)
            seqs = random.sample(C.SEQUENCES, k=len(C.SEQUENCES))

            # balanced block sizes from NUM_BLOCKS
            nb = max(1, int(C.NUM_BLOCKS))
            base, rem = divmod(len(seqs), nb)
            block_sizes = [base + (1 if b < rem else 0) for b in range(nb)]

            # assign block + position; build manifest
            manifest, i = [], 0
            for b, sz in enumerate(block_sizes, start=1):
                for pos in range(1, sz + 1):
                    entry = seqs[i]; i += 1
                    manifest.append(dict(
                        t=i,                                 # participant-specific order (1..N)
                        trial_id=entry['trial_id'],          # fixed across participants
                        seq=entry['seq'],                    # sequence string
                        block=b,                             # 1..NUM_BLOCKS
                        pos_in_block=pos,                    # 1..block_size
                        uid=str(uuid.uuid4()),               # unique trial key
                        n_blocks=nb,
                        n_rounds=len(seqs),
                    ))
            p.participant.vars['manifest'] = manifest


class Group(BaseGroup):
    pass


class Player(BasePlayer):
    # trial metadata (saved each round)
    t = models.IntegerField()                  # participant-specific order
    trial_id = models.IntegerField()           # fixed global ID (same for all)
    seq = models.StringField()                 # sequence string (e.g., "HHTTHH")
    uid = models.StringField()                 # unique trial UUID
    block = models.IntegerField()
    pos_in_block = models.IntegerField()

    # payout audit
    m_used = models.FloatField()
    treatment = models.StringField(blank=True)

    # realized outcome + earnings for this round
    realized = models.StringField(blank=True)
    win = models.BooleanField()
    round_earnings = models.FloatField()

    # decisions
    # choices = OUTCOMES + NB; stored as single-char codes (or strings if OUTCOMES uses words)
    side = models.StringField(choices=list(set(list(C.OUTCOMES))) + ['NB'])
    stake = models.IntegerField(min=0, max=C.ENDOWMENT, initial=0)
    screen_time_ms = models.IntegerField()
    button_order = models.StringField(blank=True) 

     # AOI fields
    aoi_boxes_px = models.LongStringField(blank=True)
    viewport_w_px = models.IntegerField(blank=True)
    viewport_h_px = models.IntegerField(blank=True)
    dpr = models.FloatField(blank=True)

    # Control questions (neutral codes; correctness is enforced in the Page)
    cq_keep_endowment = models.IntegerField(
        label="1. Pokud se v kole vybraném k vyplacení rozhodnete 'Nesázet', Vaše odměna bude:",
        choices=[
            [1, "OPTION_1"],
            [2, "OPTION_2"],
            [3, "OPTION_3"],
        ],
        widget=widgets.RadioSelect,
    )

    cq_multiplier = models.IntegerField(
        label="2. Pokud si v kole vybraném k vyplacení rozhodnete vsadit a sázka bude úspěšná, Vaše odměna bude:",
        choices=[
            [1, "OPTION_1"],
            [2, "OPTION_2"],
            [3, "OPTION_3"],
            [4, "OPTION_4"],
        ],
        widget=widgets.RadioSelect,
    )

    cq_toss = models.IntegerField(
        label="3. Sázka je úspěšná, pokud Vámi vybraná strana odpovídá straně, která:",
        choices=[
            [1, "OPTION_1"],
            [2, "OPTION_2"],
            [3, "OPTION_3"],
        ],
        widget=widgets.RadioSelect,
    )

    cq_payment = models.IntegerField(
        label="4. Odměna z této části studie bude určena na základě rozhodnutí:",
        choices=[
            [1, "OPTION_1"],
            [2, "OPTION_2"],
            [3, "OPTION_3"],
        ],
        widget=widgets.RadioSelect,
    )

    cq_recency = models.IntegerField(
        label="5. Který symbol (H = Hlava, O = Orel) je v zobrazené šestici hodů nejnovější (poslední)?",
        choices=[
            [1, "OPTION_1"],
            [2, "OPTION_2"],
        ],
        widget=widgets.RadioSelect,
    )

    

