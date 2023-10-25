from functools import cached_property
from Signal import Signal
from SpirowareSF6MBW import SpirowareSF6MBW


class SpirowareSF6MBWashin(SpirowareSF6MBW):
    """
    Class to deal specifically with SF6 MBW wash-ins, as opposed to wash-outs.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method_name = 'Spiroware SF6 Washin'

    @cached_property
    def prephase_breaths(self):
        """
        :return: simply returns the global prephase breaths, as the prephase for the washin is the same as the prephase
        for the combined measurement.
        """
        return self.prephase_breaths_raw

    @cached_property
    def tracer(self):
        """
        :return: Tracer concentration (SF6). This needs to be calculated in the level of SF6 WashIN, as it is uniquely
        different from the way other tracer gas concentrations are calculated.
        """
        return Signal(
            name='Tracer (inverted for washin)',
            unit=self.sf6.unit,
            data=self.config.sf6_tracer_gas - self.sf6.data
        )

    @cached_property
    def tracer_normalized(self):
        """"
        :return: End-tidal concentration divided by initial concentration. Used for the determination of the washout end
        criterion.
        """
        return self.tracer_endtidal / 0.04

    @cached_property
    def washout_breaths(self):
        """
        :return:
        """
        return self.washin_breaths_raw



