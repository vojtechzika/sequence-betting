from otree.api import *
from .models import Player

class Consent(Page):
    form_model = 'player'
    form_fields = ['consent']

    def error_message(self, values):
        if not values.get('consent'):
            return (
                "Pro pokračování ve studii je nutné udělit souhlas s podmínkami účasti. "
                "Pokud s podmínkami nesouhlasíte a nechcete se studie zúčastnit, "
                "zvedněte prosím ruku a vyčkejte příchodu administrátora."
            )


class Demographics(Page):
    form_model = 'player'
    form_fields = ['age', 'sex']


page_sequence = [Consent, Demographics]