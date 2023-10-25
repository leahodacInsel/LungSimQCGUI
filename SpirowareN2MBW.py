from Signal import Signal
from SpirowareMBW import SpirowareMBW
from functools import cached_property
from LungSimException import LungSimException
import numpy as np


class SpirowareN2MBW(SpirowareMBW):
    """
    Sub-class of Spiroware Algorithms which uses nitrogen (N2) as its tracer gas. This class contains certain functions
    entirely specific to this type of analysis, such as the N2 calculation from O2 and CO2, phase detection specific to
    nitrogen measurements etc.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.method_name = 'Spiroware N2 Washout'

    @cached_property
    def argon(self):
        """
        Argon signal as a fraction, as a function of calculated N2
        :return: Argon fraction signal
        """
        return Signal(
            name='Argon (calculated)',
            unit='si',
            data=self.n2.data * (1/self.config.fraction_n2_inert - 1)
        )

    @cached_property
    def n2(self):
        """
        Calculate nitrogen signal from O2 and CO2 signal
        :return: Nitrogen signal as a fraction
        """
        return Signal(
            name='Nitrogen (calculated)',
            unit='si',
            data=(1 - self.o2.data - self.co2.data) * self.config.fraction_n2_inert
        )

    @cached_property
    def o2(self):
        """
        Final O2 signal, after drift correction.
        :return:
        """
        # Non-elegant way of solving an issue that passing a cached_property to this function which causes a
        # LungSimException will cause an error even though this function can very well be completed without error if one
        # of the breaths properties causes an exception. They are therefore here given a special treatment which will
        # allow the signal processing to complete, while the outcomes still cause their LSExceptions.
        # TODO: Potentially reorganise this so it is less confusing
        try:
            pp = self.prephase_breaths
        except LungSimException:
            pp = np.array([]).astype(int)
        try:
            wi = self.washin_breaths
        except LungSimException:
            wi = np.array([]).astype(int)
        try:
            wo = self.washout_breaths
        except LungSimException:
            wo = np.array([]).astype(int)

        if self.config.o2_drift_option == 'breath-by-breath':
            return Signal(
                name='O2 (drift-corrected)',
                unit=self.o2_cutout.unit,
                data=self.o2_drift_calculation(
                    o2=self.o2_cutout.data,
                    co2=self.co2.data,
                    volume=self.volume.data,
                    breaths=self.breaths,
                    config=self.config,
                    pp=pp,
                    wi=wi,
                    wo=wo,
                    target_pp=self.config.n2_target_prephase,
                    target_wi=0,
                    target_wo=self.config.n2_target_washout
                )
            )
        elif self.config.o2_drift_option == 'maximum':
            return Signal(
                name='O2 (maximum-corrected)',
                unit=self.o2_cutout.unit,
                data=self.o2_maximum_calculation(
                    o2=self.o2_cutout.data,
                    target=self.config.n2_target_washout,
                    samples=self.config.drift_sample_number,
                    window=self.config.drift_correction_window
                )
            )
        else:
            return self.o2_cutout

    @cached_property
    def prephase_breaths(self):
        """
        Function to return which breaths belong to the pre-phase portion of an N2 Spiroware measurement. Based around
        which breath is the first to have a maximal inspiratory O2 value above the threshold defined by Spiroware
        (default = 60%).
        :return: a range of breath numbers which belong to the washout portion of a measurement
        """
        max_o2 = np.array([np.max(self.o2_cutout.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])
        above_threshold = max_o2 > self.config.o2_start_wo
        if above_threshold.any():
            return np.arange(
                start=0,
                stop=self.first_nonzero(above_threshold)
            ).astype(int)
        else:
            msg = 'Prephase detection error as no inspiration has O2 above {}%.'.format(self.config.o2_start_wo * 100)
            raise LungSimException(output=msg)

    @cached_property
    def tracer(self):
        """
        Tracer gas (test gas) used for the MBW trial (here N2). Used here separately so the grandmother equations which
        reference tracer still work.
        :return: tracer gas signal in fractions
        """
        return self.n2

    @cached_property
    def washout_breaths(self):
        """
        Function to return which breaths belong to the washout portion of an N2 Spiroware measurement. Based around
        which breath is the first to have a maximal inspiratory O2 value above the threshold defined by Spiroware
        (default = 60%).
        :return: a range of breath numbers which belong to the washout portion of a measurement
        """
        max_o2 = np.array([np.max(self.o2_cutout.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])
        above_threshold = max_o2 > self.config.o2_start_wo
        if above_threshold.any():
            return np.arange(
                start=self.first_nonzero(above_threshold),
                stop=self.breaths.number
            ).astype(int)
        else:
            msg = 'Washout detection error as no inspiration has O2 above {}%.'.format(self.config.o2_start_wo*100)
            raise LungSimException(output=msg)




