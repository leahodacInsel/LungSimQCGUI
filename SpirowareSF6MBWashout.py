from functools import cached_property
from SpirowareSF6MBW import SpirowareSF6MBW
import numpy as np


class SpirowareSF6MBWashout(SpirowareSF6MBW):
    """
    Class to deal specifically with SF6 MBW wash-outs, as opposed to wash-ins.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method_name = 'Spiroware SF6 Washout'

    @cached_property
    def tracer(self):
        """
        :return: Tracer concentration (SF6). This needs to be calculated in the level of SF6 WashOUT, as the WashIn
        sister class calculates a tracer gas signal completely differently.
        """
        return self.sf6

    @cached_property
    def prephase_breaths(self):
        """
        :return: simply returns the global prephase breaths, as the prephase for the washin is the same as the prephase
        for the combined measurement.
        """
        return np.arange(start=self.washin_end, stop=self.washin_breaths_raw[-1]+1)

    @cached_property
    def washout_breaths(self):
        """
        :return:
        """
        return self.washout_breaths_raw
