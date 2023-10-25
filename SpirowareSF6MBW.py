from functools import cached_property
from LungSimException import LungSimException
from Signal import Signal
from SpirowareMBW import SpirowareMBW
import numpy as np


class SpirowareSF6MBW(SpirowareMBW):
    """
    Class for the calculation of MBW measurements performed using sulphur hexafluoride (SF6). Contains the shared
    code of the SF6 Washout and Washin algorithms. Based on "Novel methodology to perform sulfur hexafluoride (
    SF6)-based multiple-breath wash-in and washout in infants using current commercially available equipment" as
    implemented in Spiroware 3.2.1+
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @cached_property
    def ar(self):
        """
        :return: Argon concentration
        """
        # Real contributions of argon and sf6 to the composition of atmospheric air and the artificial gas mixture
        # that is used for the washout. Factor for the determination of argon fraction as a function of sf6
        # concentration. Named "s" after the convention of the Gustaffson paper referenced above.
        s = (self.config.ar_tracergas - self.config.ar_air) / (self.config.sf6_tracergas - self.config.sf6_air)

        return Signal(
            name='Argon (calculated)',
            unit='si',
            data=self.config.ar_air + s * self.sf6.data
        )

    @cached_property
    def n2(self):
        """
        :return: Nitrogen concentration
        """
        return Signal(
            name='Nitrogen (calculated)',
            unit='si',
            data=1 - self.o2.data - self.co2.data - self.sf6.data - self.ar.data
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
            pp = self.prephase_breaths_raw
        except LungSimException:
            pp = np.array([]).astype(int)
        try:
            wi = self.washin_breaths_raw
        except LungSimException:
            wi = np.array([]).astype(int)
        try:
            wo = self.washout_breaths_raw
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
                    target_pp=self.config.sf6_target_prephase,
                    target_wi=self.config.sf6_target_washin,
                    target_wo=self.config.sf6_target_washout
                )
            )
        else:
            return self.o2_cutout

    @cached_property
    def prephase_breaths_raw(self):
        """
        Function to return which breaths belong to the pre-phase portion of an SF6 Spiroware measurement. Based around
        which breath is the first to have a median inspiratory MMss value above half of the minimum-maximum difference.
        :return: a range of breath numbers which belong to the washout portion of a measurement.
        NOTE: This algorithm does not correspond to the Spiroware method, which instead reportedly checks for the
        first breath above 2% preliminary calculated SF6, to determine washin-start. The method used here
        performs better at replicating Spiroware than the algorithm described to us by Ecomedics.
        """
        sf6_max = np.array([np.max(self.sf6_estimate.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])

        max_above_threshold = sf6_max > self.config.sf6_start_wi
        if max_above_threshold.any():
            return np.arange(
                start=0,
                stop=self.first_nonzero(max_above_threshold)
            ).astype(int)
        else:
            msg = 'Prephase detection error as no inspiration has SF6 above {}%.'.format(0.02 * 100)
            raise LungSimException(output=msg)

    def calculate_sf6(self, mmside, co2, o2):
        """
        :return: calculated sf6 fraction as a function of MMside, CO2 and O2 signals. Packaged into a separate function
        so that sf6 can be calculated as any point, even with incomplete signals.
        """
        # Ratios of molar mass to kappa (for easier calculations going forward)
        f_o2 = self.config.mass['o2'] / self.config.kappa['o2']
        f_co2 = self.config.mass['co2'] / self.config.kappa['co2']
        f_sf6 = self.config.mass['sf6'] / self.config.kappa['sf6']
        f_n2 = self.config.mass['n2'] / self.config.kappa['n2']
        f_ar = self.config.mass['ar'] / self.config.kappa['ar']

        # Real contributions of argon and sf6 to the composition of atmospheric air and the artificial gas mixture
        # that is used for the washout. Factor for the determination of argon fraction as a function of sf6
        # concentration. Named "s" after the convention of the Gustaffson paper referenced above.
        s = (self.config.ar_tracer_gas - self.config.ar_air) / (self.config.sf6_tracer_gas - self.config.sf6_air)

        # Calculation of the sf6 fraction
        sf6 = (mmside / self.config.kappa['air'] -
               f_n2 * (1 - o2 - co2 - self.config.ar_air) -
               f_o2 * o2 - f_co2 * co2 - f_ar * self.config.ar_air) / \
              (f_sf6 - (1 + s) * f_n2 + s * f_ar)

        return sf6

    @cached_property
    def sf6(self):
        """
        :return: Sulphur hexafluoride (SF6) signal.
        """
        return Signal(
            name='SF6 (calculated)',
            unit='si',
            data=self.calculate_sf6(mmside=self.mmside.data, co2=self.co2.data, o2=self.o2.data)
        )

    @cached_property
    def sf6_estimate(self):
        """
        :return: Sulphur hexafluoride (SF6) signal, but calculated before drift correction. This is used by Spiroware
        to calculate the phase detection, which depends on SF6 concentration. The final O2 signal at that point is
        however not available, as the drift correction depends in turn on the phase detection.
        """
        return Signal(
            name='SF6 (calculated)',
            unit='si',
            data=self.calculate_sf6(mmside=self.mmside.data, co2=self.co2.data, o2=self.o2_cutout.data)
        )

    @cached_property
    def tracer_endtidal(self):
        """
        End-tidal (from x% to y% of expired volume of a breath) sampling of the tracer gas signal.
        :return: Concentration of each breath in an end-tidal interval
        """
        return np.array(
            [
                self.end_tidal(
                    signal=self.tracer.data[self.breaths.exp_range(i)],
                    flow=self.flow.data[self.breaths.exp_range(i)],
                    time=self.flow.time[self.breaths.exp_range(i)],
                    window=self.config.endtidal_window
                ) for i in self.breaths.numbers
            ]
        )

    @cached_property
    def washin_breaths_raw(self):
        """
        Function to return which breaths belong to the pre-phase portion of an SF6 Spiroware measurement. Based around
        which breath is the first to have a median inspiratory MMss value above half of the minimum-maximum difference.
        :return: a range of breath numbers which belong to the washout portion of a measurement.
        NOTE: This algorithm does not correspond to the Spiroware method, which instead reportedly checks for the
        first breath above 2% preliminary calculated SF6, to determine washin-start. The method used here
        performs better at replicating Spiroware than the algorithm described to us by Ecomedics.
        NOTE: This corresponds to the overall phase detection, where we detect prephase, washin, washout, hence _raw.
        For the purposes of calculating outcomes etc. we consider different ranges the prephase and washout, depending
        on whether we are analysing a washin, or washout. By using this slightly counterintuitive notation we can
        simplify a lot of code.
        """
        max_insp = np.array([np.max(self.sf6_estimate.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])
        min_insp = np.array([np.min(self.sf6_estimate.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])
        min_exp = np.array([np.min(self.sf6_estimate.data[self.breaths.exp_range(i)]) for i in self.breaths.numbers])

        exp_above_threshold = min_exp > self.config.sf6_min_wi
        if not exp_above_threshold.any():
            concentration = np.round(self.config.sf6_min_wi * 100, 2)
            msg = 'Washin expirations do not reach at least {}%.'.format(concentration)
            raise LungSimException(output=msg)

        breath_wi_complete = self.first_nonzero(exp_above_threshold)
        candidates = self.breaths.numbers >= breath_wi_complete
        insp_below_threshold = np.all([candidates, min_insp < self.config.sf6_start_wo], 0)
        insp_above_threshold = max_insp > self.config.sf6_start_wi

        if not insp_below_threshold.any():
            concentration = self.config.sf6_start_wo * 100
            msg = 'Washin detection error as no inspiration (after washin) has SF6 below {}%.'.format(concentration)
            raise LungSimException(output=msg)

        elif not insp_above_threshold.any():
            concentration = self.config.sf6_start_wi * 100
            msg = 'Washin detection error as no inspiration has SF6 above {}%.'.format(concentration)
            raise LungSimException(output=msg)

        else:
            return np.arange(
                start=self.first_nonzero(insp_above_threshold),
                stop=self.first_nonzero(insp_below_threshold)
            ).astype(int)

    @cached_property
    def washin_end(self):
        """
        After washing in past the point of 3.9% of exhaled SF6, the measurement transitions from washin to washout.
        Spiroware cuts the measurement at this point, and treats the resulting two halves as separate "washout"
        measurements. Will be used for the calculation of the properties prephase_breaths und washout_breaths for the
        daughter classes of this one.
        :return: breath number where the washin ends according to Spiroware's algorithm.
        """
        # Target concentration of SF6 tracer gas, above which the washin in considered complete.
        target = self.config.tracer_washin_end
        # Find target breath, as the first of three breaths which are above the target concentration. We
        # construct a logical array of those end-tidal concentrations which lie below the target, circ-shift the array
        # twice and add the three arrays together. The first entry to be equal to 3 corresponds to the first entry
        # which is followed by two valid end-tidal concentrations.
        original = self.tracer_endtidal[self.washin_breaths_raw] > target
        shift_one = np.roll(original, -1)
        shift_one[-1] = 0
        shift_two = np.roll(original, -2)
        shift_two[-2:] = 0
        logical_array = original & shift_one & shift_two
        # If there are array indices that fulfil the above criterium of three consecutive breaths above the target,
        # choose the first one of these for the washout end
        if logical_array.any():
            washin_end = logical_array.nonzero()[0][0] + 3
            washin_end = self.washin_breaths_raw[washin_end]

        # If three consecutive breaths aren't found, find the last dip above the target and return it instead.
        else:
            cutoff = self.tracer_endtidal[self.washin_breaths_raw] - target
            crossing = cutoff * np.roll(cutoff, 1) <= 0
            negative = cutoff > 0
            # If no dip below the critical condition is found, return nothing. Otherwise return the last dip below.
            if negative.any():
                washin_end = (crossing & negative).nonzero()[0][-1]
                washin_end = self.washin_breaths_raw[washin_end]

            else:
                msg = 'The washin end of {}% is not reached'.format(target*100)
                raise LungSimException(output=msg)

        return washin_end

    @cached_property
    def washout_breaths_raw(self):
        """
        Function to return which breaths belong to the pre-phase portion of an SF6 Spiroware measurement. Based around
        which breath is the first to have a median inspiratory MMss value above half of the minimum-maximum difference.
        :return: a range of breath numbers which belong to the washout portion of a measurement.
        NOTE: This algorithm does not correspond to the Spiroware method, which instead reportedly checks for the
        first breath above 2% preliminary calculated SF6, to determine washin-start. The method used here
        performs better at replicating Spiroware than the algorithm described to us by Ecomedics.
        NOTE: This corresponds to the overall phase detection, where we detect prephase, washin, washout, hence _raw.
        For the purposes of calculating outcomes etc. we consider different ranges the prephase and washout, depending
        on whether we are analysing a washin, or washout. By using this slightly counterintuitive notation we can
        simplify a lot of code.
        """
        max_insp = np.array([np.max(self.sf6_estimate.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])
        min_insp = np.array([np.min(self.sf6_estimate.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])
        min_exp = np.array([np.min(self.sf6_estimate.data[self.breaths.exp_range(i)]) for i in self.breaths.numbers])

        exp_above_threshold = min_exp > self.config.sf6_min_wi
        if not exp_above_threshold.any():
            concentration = np.round(self.config.sf6_min_wi * 100, 2)
            msg = 'Washin expirations do not reach at least {}%.'.format(concentration)
            raise LungSimException(output=msg)

        breath_wi_complete = self.first_nonzero(exp_above_threshold)
        candidates = self.breaths.numbers >= breath_wi_complete
        insp_below_threshold = np.all([candidates, min_insp < self.config.sf6_start_wo], 0)
        insp_above_threshold = max_insp > self.config.sf6_start_wi

        if not insp_below_threshold.any():
            concentration = self.config.sf6_start_wo * 100
            msg = 'Washout detection error as no inspiration (after washin) has SF6 below {}%.'.format(concentration)
            raise LungSimException(output=msg)

        elif not insp_above_threshold.any():
            concentration = self.config.sf6_start_wi * 100
            msg = 'Washout detection error as no inspiration has SF6 above {}%.'.format(concentration)
            raise LungSimException(output=msg)

        else:
            return np.arange(
                start=self.first_nonzero(insp_below_threshold),
                stop=self.breaths.number
            ).astype(int)
