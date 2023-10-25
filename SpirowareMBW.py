from Breaths import Breaths
from functools import cache
from functools import cached_property
from LungSimException import LungSimException
from MBW import MBW
from math import nan
from math import log
from numba import njit
import numpy as np
from scipy.integrate import cumtrapz
from scipy.signal import firwin
from scipy.signal import lfilter
from scipy.signal import lfilter_zi
from scipy.signal import bilinear
from Signal import Signal
from SpirowareConfig import SpirowareConfig
from SpirowareQualityControl import SpirowareQualityControl
from UnitConversion import UnitConversion
from copy import deepcopy


class SpirowareMBW(MBW):
    """
    The SpirowareMBW class will contain all the attributes and methods required to go from a desired Output to the raw
    data contained in the SpirowareFile. This class collects all the methods and attributes which are unique to Spiroware,
    but not unique to the individual kinds of washout that Spiroware can handle (i.e. N2WO, SF6WI, SF6WO etc.)
    """

    def __init__(self, **kwargs):
        """
        A Spiroware algorithm almost certainly requires a File, a Config, and Results to be computed. This is inherited
        from the super class, MBW, which shares these traits
        """
        super().__init__(**kwargs)

    @classmethod
    def atp_correction(cls, co2, h_insp, t_insp, h_flowhead, t_flowhead, p_ambient):
        """
        Correction applied to CO2 signal to calculate the concentration of CO2 in air if water vapor was not taken
        into account
        :param co2: co2 signal data
        :param h_insp: humidity of air at the end of the correction (0% by default)
        :param t_insp: temperature of air in room
        :param h_flowhead: estimated humidity at the flowhead when exhaling (60% by default)
        :param t_flowhead: estimated temperature at the flowhead when exhaling
        :param p_ambient: ambient pressure in room
        :return: corrected co2 signal data
        """
        pressure_inspired = cls.vapor_pressure(h_insp, t_insp)
        pressure_flowhead = cls.vapor_pressure(h_flowhead, t_flowhead)

        correction = (p_ambient - pressure_inspired) / (p_ambient - pressure_flowhead)

        return co2 * correction

    @cached_property
    def breaths(self):
        """
        The breaths property contains the final breath indices. It is created from the raw breaths object originally
        created from an uncut flow signal. After the parts of the signal between and outside breaths have been removed
        by the cutout function, the breath indices need to be readjusted.
        :return: A Breaths object containing the final breath indices
        """
        # We start by copying breaths_raw, as we do not want to affect it with our calculations
        output = deepcopy(self.breaths_raw)

        # Only perform this correction if the signal between breaths is actually removed. Otherwise it is unnecessary.
        if self.config.cutout_option:
            # The +1 is important to function like Spiroware.
            reference = np.copy(output.exp_end) + 1
            cut_away = np.zeros_like(reference)
            # We initialise the last element of cut_away as it will be the first to be called in the for loop below
            # (0-1=-1)
            cut_away[-1] = reference[-1]
            # Determine by how much each breath is set back from its original position
            for i in output.numbers:
                cut_away[i] = cut_away[i - 1] - reference[i - 1] + output.insp_start[i]
            # Reduce the original breath indices by cut_away
            output.insp_start -= cut_away
            output.insp_end -= cut_away
            output.exp_start -= cut_away
            output.exp_end -= cut_away

        return output

    @cached_property
    def breaths_raw(self):
        """
        Method to perform breath detection. Breath detection is an important step in the signal processing of Spiroware
        measurements, where different breath detection algorithms may produce different outcomes. This function controls
        which breath detection is used, depending on the configuration.
        :return: A Breaths object containing breath indices
        """
        # Breath detection algorithm used by Spiroware
        if self.config.breath_option == 'spiroware':
            return self.breath_detection_spiroware(
                time=self.flow_filtered.time,
                flow=self.flow_filtered.data,
                minvol=self.config.breath_detection_minvol
            )
        # Breath detection algorithm based on flow + co2 signal
        elif self.config.breath_option == 'co2':
            return self.breath_detection_co2(
                flow=self.flow_filtered.data,
                co2=self.co2_crosstalk.data,
                threshold=self.config.co2_threshold
            )
        # Breath detection algorithm based on minimum volume, similar to Spiroware, but working properly.
        elif self.config.breath_option == 'minvol':
            return self.breath_detection_minvol(
                time=self.flow_filtered.time,
                flow=self.flow_filtered.data,
                minvol=self.config.breath_detection_minvol
            )

    @staticmethod
    def breath_detection_spiroware(time, flow, minvol):
        """
        Complicated breath detection algorithm based primarily on a minimum-volume criterion. Once a volume between two
        zero-crossings reaches a certain volume, it is considered to be sufficient. However, there are major exceptions
        to this:
        1. If two valid expirations are separated by an invalid inspiration, the inspiration will become valid, causing
        miniature breaths to appear in the analysis.
        2. If two valid inspirations are separated by an invalid expiration, the previous valid inspiration will
        typically be deleted together with the invalid expiration.
        :param time: time array
        :param flow: flow signal in m^3/s
        :param minvol: a minimum volume in m^3 which acts as an acceptance criterion for a breath.
        :return: A Breaths object containing breath indices.
        """
        # Find inspiration start and end
        flow_positive = flow > 0
        # Create a shifted logical array of the flow signal
        flow_positive_original = flow_positive[0:-1]
        flow_positive_shift = flow_positive[1:]
        # Find the indices where original flow is negative and shifted flow is positive, i.e. sign changed from
        # expiration to inspiration. Then inspiration start is one point later.
        insp_start = (~ flow_positive_original & flow_positive_shift).nonzero()[0] + 1
        # Find the indices where original flow is positive and shifted flow is negative, i.e. sign changed from
        # inspiration to expiration. This is exactly where the inspiration ends.
        insp_end = (flow_positive_original & ~ flow_positive_shift).nonzero()[0]

        # Find expiration start and end:
        flow_negative = flow < 0
        flow_negative_original = flow_negative[0:-1]
        flow_negative_shift = flow_negative[1:]
        # Find the indices where original flow is positive and shifted flow is negative, i.e. sign changed from
        # inspiration to expiration. Then expiration start is one point later
        exp_start = (~ flow_negative_original & flow_negative_shift).nonzero()[0] + 1
        # Find the indices where original flow is negative and shifted flow is positive, i.e. sign changed from
        # expiration to inspiration. This is exactly where the expiration ends
        exp_end = (flow_negative_original & ~ flow_negative_shift).nonzero()[0]

        # If any of the above has a size of zero, raise an exception:
        if insp_start.size * insp_end.size * exp_start.size * exp_end.size < 1:
            raise LungSimException(output='Failed to detect the minimum of one complete breath (insp + exp)')

        # Now we want pairs of inspirations and expirations. Experiment must start with an inspiration start
        if insp_end[0] < insp_start[0]:
            insp_end = np.delete(insp_end, 0, axis=0)
        if exp_start[0] < insp_start[0]:
            exp_start = np.delete(exp_start, 0, axis=0)
        if exp_end[0] < insp_start[0]:
            exp_end = np.delete(exp_end, 0, axis=0)

        # After having deleted some of the relevant indices, check again whether the size of arrays has dropped to 0. If
        # yes, raise exception.
        if insp_start.size * insp_end.size * exp_start.size * exp_end.size < 1:
            raise LungSimException(output='Failed to detect the minimum of one complete breath (insp + exp)')

        # Must end with an exspiration end
        if insp_start[-1] > exp_end[-1]:
            insp_start = np.delete(insp_start, -1, axis=0)
        if insp_end[-1] > exp_end[-1]:
            insp_end = np.delete(insp_end, -1, axis=0)
        if exp_start[-1] > exp_end[-1]:
            exp_start = np.delete(exp_start, -1, axis=0)

        # After having deleted some of the relevant indices, check again whether the size of arrays has dropped to 0. If
        # yes, raise exception.
        if insp_start.size * insp_end.size * exp_start.size * exp_end.size < 1:
            raise LungSimException(output='Failed to detect the minimum of one complete breath (insp + exp)')

        # Integrate the signal once, which then avoids errors with the cumtrapz function
        volume_total = cumtrapz(y=flow, x=time, initial=0)
        vol_insp = volume_total[insp_end] - volume_total[insp_start]
        vol_exp = volume_total[exp_start] - volume_total[exp_end]

        # Logical vector of valid inspirations and expirations i.e. those breaths above volume threshold defined in
        # A-file header or default options
        valid_insp = vol_insp >= minvol
        valid_exp = vol_exp >= minvol

        # If there are either no valid inspirations or expirations, the algorithm cannot continue beyond this point.
        if np.sum(valid_insp) < 1 or np.sum(valid_exp) < 1:
            raise LungSimException(output='Failed to detect breath (insp + exp) with sufficient volume')

        # Spiroware convention: a trial always starts with a valid inspiration and a valid expiration
        first_valid_insp = valid_insp.nonzero()[0][0]

        # Initialize some parameters used in actual breath detection. Caution!: This might limit the max number of
        # breaths to 1000. This is unlikely to happen
        breath_number = -1
        last_valid_insp_number = -1
        length = 1000
        real_insp_start = np.zeros(length)
        real_insp_end = np.zeros(length)
        real_exp_start = np.zeros(length)
        real_exp_end = np.zeros(length)

        try:
            # Actual breath detection starting with the first valid inspiration (This block of code contains some odd
            # behavior, but this is intended to replicate Spiroware 3.2.1)
            for k in range(first_valid_insp, max(np.shape(exp_end))):

                # For every originally detected breath, we check whether it contains a valid expiration
                if valid_exp[k] == 1:
                    # If it does contain a valid expiration, we raise the breath number by 1.
                    breath_number += 1
                    # We accept this expiration start as a member of our final selection of expiration starts
                    real_exp_start[breath_number] = exp_start[k]
                    # We also accept the expiration end of the same index as a member of the final selection. It is pretty
                    # safe to assume that this is the next zero crossing following this expiration start.
                    real_exp_end[breath_number] = exp_end[k]

                    # If, in addition, this index contains a valid inspiration...
                    if valid_insp[k] == 1:
                        # ...then we accept both the inspiration start index and the end index
                        real_insp_start[breath_number] = insp_start[k]
                        real_insp_end[breath_number] = insp_end[k]
                    # ...however, if the inspiration is not valid...
                    else:
                        # ...and the index is the same as the one of the last valid insp number(??)... What purpose does
                        # this serve?
                        if last_valid_insp_number == breath_number:
                            # ... then the insp start and end of this index are the ones calculated below.
                            real_insp_start[breath_number] = last_valid_insp_start
                            real_insp_end[breath_number] = last_valid_insp_end
                        else:
                            real_insp_start[breath_number] = real_exp_end[breath_number - 1] + 1
                            real_insp_end[breath_number] = real_exp_start[breath_number] - 1
                else:
                    # If the index does not contain a valid expiration, but a valid inspiration...
                    if valid_insp[k] == 1:
                        # ...we store this inspiration's information in the last_valid_insp variables. We then wait to
                        # see if in the next round we have a valid expiration. If we have both a valid inspiration and
                        # expiration, this information goes unused. If there is only a valid expiration, we use this
                        # information if it immediately preceded it.
                        last_valid_insp_start = insp_start[k]
                        last_valid_insp_end = insp_end[k]
                        last_valid_insp_number = breath_number + 1

        except IndexError:
            raise LungSimException(output='Spiroware breath detection failed.')

        # Removing zero entries
        first_zero_entry = (real_insp_start == 0).nonzero()[0][0]

        if first_zero_entry:
            insp_start = real_insp_start[0:first_zero_entry]
            insp_end = real_insp_end[0:first_zero_entry]
            exp_start = real_exp_start[0:first_zero_entry]
            exp_end = real_exp_end[0:first_zero_entry]
        else:
            insp_start = real_insp_start
            insp_end = real_insp_end
            exp_start = real_exp_start
            exp_end = real_exp_end

        # Create a Breaths object
        return Breaths(
            insp_start=insp_start,
            insp_end=insp_end,
            exp_start=exp_start,
            exp_end=exp_end
        )

    @staticmethod
    def breath_detection_co2(flow, co2, threshold):
        """
        Breath detection function based on "Breath detection algorithms affect multiple-breath washout outcomes in
        pre-school and school age children" (2022?) by Oestreich et al.
        Improved breath detection which takes CO2 signal to help identify breaths.
        :param flow: Flow signal
        :param co2: CO2 signal (synchronised to flow)
        :param threshold: CO2 concentration above which zero crossings are identified as belonging to exps, and below
            which they are identified as insps.
        :return: Breaths object containing start and stop indices.
        """
        # A logical array of where the flow crosses 0
        zero_crossing = (flow * np.roll(flow, 1)) <= 0
        # Array of indexes where the flow crosses 0
        zero_indexes = zero_crossing.nonzero()[0]
        # Returns an array of the co2 values of the above zero crossings. In an ideal breath detection, these would
        # alternative between values around ~5% for beginnings of inspirations and ~0% for the beginnings of
        # expirations.
        crossing = co2[zero_indexes]
        # Instead of percentage values, we take whether or not the co2 values of the crossings lie above or below the
        # approximate half-way point between 0% and 5%. This is a semi-arbitrary choice, as you could argue that it
        # should be calculated from the signal. Potentially could be changed in the future
        crossing = crossing > threshold
        # A logical array which is non-zero for indices where the above arbitrary line is crossed.
        sign_change = crossing[0:-1] ^ crossing[1:]
        # Returns the indices of those non-zero values (i.e. the ones we want to keep). In short, if there are many
        # zero-crossings which all have a co2 signal above the threshold, we ignore all but the last one. Same for
        # those with values below threshold.
        keep = sign_change.nonzero()[0]
        if keep.size == 0:
            raise LungSimException(output='No flow reversal in flow signal. Breath detection impossible.')

        # Cuts the two arrays of zero-indexes and co2-crossings to only the indexes we want to keep
        zero_indexes = zero_indexes[keep]
        crossing = crossing[keep]
        # Finds the first and last inspiration. If
        if crossing.nonzero()[0].size == 0:
            raise LungSimException(output='Flow signal too short to detect multiple breaths.')

        first_crossing = crossing.nonzero()[0][0]
        last_crossing = crossing.nonzero()[0][-1]
        # Only keeps the zero indexes from first to last inspiration. We need the last inspiration crossing index
        # only to calculate where the last expiration ends.
        zero_indexes = zero_indexes[first_crossing:last_crossing]
        # Because of the above operation, we now know that the sequence of zero-crossings begins with an inspiration
        # and ends with an expiration.
        inspirations = np.arange(start=0, stop=np.shape(zero_indexes)[0], step=2)
        expirations = np.arange(start=1, stop=np.shape(zero_indexes)[0] - 1, step=2)
        # The inspiration starts coincide with the above indexes, but we ignore the last one because we want to end
        # with an expiration.
        insp_start = zero_indexes[inspirations[0:-1]]
        # The expiration starts coincide with the above expirations.
        exp_start = zero_indexes[expirations]
        # The ends of inspirations are exactly one index before the starts of expirations.
        insp_end = zero_indexes[expirations] - 1
        # And the expiration ends are one index before the starts of the inspirations. We ignore the first one
        # because we want to start with an inspiration.
        exp_end = zero_indexes[inspirations[1:]] - 1

        return Breaths(
            insp_start=insp_start,
            insp_end=insp_end,
            exp_start=exp_start,
            exp_end=exp_end
        )

    @staticmethod
    def btps_correction(flow, factor_insp, factor_exp):
        """
        BTPS stands for body temperature and pressure, saturated conditions. This correction multiplies any
        incoming flow signal by how much inspired air will expand in volume inside the lungs, and multiplies any
        outgoing flow by how much it has decreased in volume since being expired.
        :param flow: flow data (numpy array): expirations negative, inspirations positive.
        :param factor_insp: multiplier for inspiratory flow
        :param factor_exp: multiplier for expiratory flow
        :return: copy of input flow signal with multiplier factors applied
        """
        output = np.copy(flow)
        return output * ((output > 0) * factor_insp + (output <= 0) * factor_exp)

    @cached_property
    def btps_factor_insp(self):
        """
        Function to calculate the factor by which inspired volume needs to be multiplied in order to achieve its volume
        once it has been warmed up to body temperature and saturated with body humidity.
        :return: inspiratory btps factor
        """
        # Magnus equation estimation of the vapor pressure of water in air with a given temp and humidity.
        vapor_p_body = self.vapor_pressure(humidity=self.config.humidity_body, temp=self.config.temp_body)
        vapor_p_inspired = self.vapor_pressure(humidity=self.config.humidity_inspired, temp=self.config.temp_inspired)
        # "Partial" pressure of air without water inside the body and outside
        p_body = self.config.pressure_ambient - vapor_p_body
        p_inspired = self.config.pressure_ambient - vapor_p_inspired
        # Temperature inside the body and room in Kelvin
        t_body = self.config.temp_body - self.config.temp_kelvin_zero
        t_inspired = self.config.temp_inspired - self.config.temp_kelvin_zero
        # Volumes of the same air packet in the room vs in the body
        v_body = t_body / p_body
        v_inspired = t_inspired / p_inspired
        return v_body / v_inspired

    @cached_property
    def btps_factor_exp(self):
        """
        Function to calculate the factor by which expired volume needs to be multiplied in order to retrieve its volume
        of when it was still warmed up to body temperature and saturated with body humidity.
        :return: expiratory btps factor
        """
        # Magnus equation estimation of the vapor pressure of water in air with a given temp and humidity.
        vapor_p_body = self.vapor_pressure(humidity=self.config.humidity_body, temp=self.config.temp_body)
        vapor_p_flowhead = self.vapor_pressure(humidity=self.config.humidity_flowhead, temp=self.config.temp_flowhead)
        # "Partial" pressure of air without water inside the body and outside
        p_body = self.config.pressure_ambient - vapor_p_body
        p_flowhead = self.config.pressure_ambient - vapor_p_flowhead
        # Temperature inside the body and flowhead in Kelvin
        t_body = self.config.temp_body - self.config.temp_kelvin_zero
        t_flowhead = self.config.temp_flowhead - self.config.temp_kelvin_zero
        # Volumes of the same air packet in the flowhead vs in the body
        v_body = t_body / p_body
        v_flowhead = t_flowhead / p_flowhead
        return v_body / v_flowhead

    @cached_property
    def cev_alv(self):
        """
        The cumulative expired volume (CEV) at alveolar opening (alv) represents the total volume expired from washout
        start to the breath where expired concentration falls below 2.5% of the initial concentration. It is
        additionally adjusted for the volume of dead space which does not contribute to the washout for each breath.
        This total dead space (machine + physiological) volume is subtracted from each expiration of the washout
        breaths, regardless of their volume.
        :return: Cumulative Expired Volume at alveolar opening for each breath.
        """
        cev_alv = np.zeros_like(self.volume_expired)
        total_dead_space = self.machine_dead_space + self.physio_dead_space
        cev_alv[self.washout_breaths] = np.cumsum(self.volume_expired[self.washout_breaths] - total_dead_space)
        return cev_alv

    @cached_property
    def co2_atp(self):
        """
        CO2-signal corrected for atmospheric temperature and pressure. The CO2 sensor measures partial pressure, and in
        order to receive a measure of concentration from this measure we need to calculate the contribution of water
        vapor to total pressure, and adjust the CO2 value for it
        :return: ATP-corrected CO2 signal
        """
        if self.config.atp_option:
            return Signal(
                name='CO2 (ATP-corrected)',
                unit=self.co2_raw.unit,
                data=self.atp_correction(
                    co2=self.co2_raw.data,
                    h_insp=self.config.humidity_inspired,
                    t_insp=self.config.temp_inspired,
                    h_flowhead=self.config.humidity_flowhead,
                    t_flowhead=self.config.temp_flowhead,
                    p_ambient=self.config.pressure_ambient
                )
            )
        else:
            return self.co2_raw

    @cached_property
    def co2_cev(self):
        """
        Cumulative net expired volume of CO2
        """
        co2_cev = np.zeros_like(self.co2_expired_net)
        co2_cev[self.washout_breaths] = np.cumsum(self.co2_expired_net[self.washout_breaths])
        return co2_cev

    @cached_property
    def co2_crosstalk(self):
        """
        co2 signal after crosstalk correction. Crosstalk correction corrects the co2 concentration for the presence of
        Co2 in the gas mixture starting from Spiroware 3.3.1 and onwards. This crosstalk mode is enabled when the option
        crosstalk_option is 'quadratic' and disabled when it is 'one-way'.
        :return: crosstalk corrected co2 signal
        """
        return Signal(
            name='CO2 (crosstalk corrected)',
            unit=self.co2_filtered.unit,
            data=self.co2_crosstalk_calculation(
                config=self.config,
                co2_name=self.file.co2_name,
                o2=self.o2_spedup.data,
                co2=self.co2_filtered.data
            )
        )

    @staticmethod
    def co2_crosstalk_calculation(config, co2_name, o2, co2):
        """
        CO2 crosstalk correction. This function has some additional pit-falls to consider. Historically, the CO2 signal
        would be crosstalk corrected by a factor that corresponds to 100% O2. This corrected signal was written into the
        A-files as "CO2". Once the crosstalk correction was implemented, it was removed, and the uncorrected signal was
        called "CO2_RAW". If we therefore work with a file which has a pre-corrected signal, we need to undo that
        correction before applying the O2-concentration-dependent one.
        Additionally, up until Spiroware 3.3.0, the O2-dependent CO2-correction was different ('one-way') than after and
        including Spiroware 3.3.1 (see "Wyler et al 2021 Correction of ...."
        :param o2: Oxygen signal data as fraction
        :param co2: CO2 signal data as fraction
        :return: corrected CO2 signal data, depending on which kind of crosstalk correction should be applied.
        """
        # Correction of old signal processing in legacy files
        if co2_name == 'CO2':
            co2 /= 1.075

        # Apply crosstalk correction
        if config.crosstalk_option == 'one-way':
            co2 *= config.gain_co2 * o2 + config.offset_co2
            return co2

        elif config.crosstalk_option == 'quadratic':
            E = config.crosstalk_E
            F = config.crosstalk_F
            G = config.crosstalk_G
            return co2 + E * co2 + F * np.square(co2) + G * co2 * o2

    @cached_property
    def co2(self):
        """
        co2 signal after the breath-detection related cutting out of in-between breaths sections of the signal.
        :return: co2 signal after cutout
        """
        if self.config.cutout_option:
            return Signal(
                name='co2 (cutout)',
                unit=self.co2_crosstalk.unit,
                data=self.cutout(
                    input=self.co2_crosstalk.data,
                    breaths=self.breaths_raw
                )
            )
        else:
            return self.co2_crosstalk

    @cached_property
    def co2_delayed(self):
        """
        The raw co2 signal is delay corrected by a somewhat complicated algorithm which is under certain circumstances
        defective, and certainly not maximally precise. However, this is a centrally important function for the correct
        interpretation of raw Spiroware signals. As the Cco2 signal, the co2 signal and the co2 signal are measured at
        very different places in the experimental setup, it becomes necessary to align them in time.
        :return: delay-corrected co2 signal
        """
        return Signal(
            name='co2 (delay-corrected)',
            unit=self.co2_raw.unit,
            data=self.delay_correction(
                input=self.co2_atp.data,
                delay_times=self.delay_times['co2']
            )
        )

    @cached_property
    def co2_endtidal(self):
        """
        End-tidal (from x% to y% of expired volume of a breath) sampling of the CO2 gas signal.
        :return: Concentration of each breath in an end-tidal interval
        """
        return np.array(
            [
                self.end_tidal(
                    signal=self.co2.data[self.breaths.exp_range(i)],
                    flow=self.flow.data[self.breaths.exp_range(i)],
                    time=self.flow.time[self.breaths.exp_range(i)],
                    window=self.config.endtidal_window
                ) for i in self.breaths.numbers
            ]
        )

    @cached_property
    def co2_expired(self):
        """
        Expired CO2
        """
        return self.co2_volume.data[self.breaths.exp_start] - self.co2_volume.data[self.breaths.exp_end]

    @cached_property
    def co2_expired_net(self):
        """
        Net expired CO2
        """
        return self.co2_expired - self.co2_inspired

    @cached_property
    def co2_filtered(self):
        """
        The CO" signal is filtered with a lowpass filter to smooth it. This leads to signal delay and some loss of data.
        The Spiroware algorithm's default filters the signal with a harsh lowpass filter with a cutoff of 2Hz, which 
        leads to signal delay especially in the rise and fall moments.
        :return: Filtered CO2 signal
        """
        if self.config.co2_filter_option:
            return Signal(
                name='co2 (filtered)',
                unit=self.co2_trimmed.unit,
                data=self.filter_lowpass_spiroware(
                    input=self.co2_trimmed.data,
                    sampling_rate=self.co2_trimmed.sampling_rate,
                    cutoff=self.config.co2_filter_cutoff,
                    window=self.config.co2_filter_window
                )
            )
        else:
            return self.co2_trimmed

    @cached_property
    def co2_inspired(self):
        """
        Inspired CO2
        We made it depend on the same reinspired rule as the tracer.
        """
        if self.config.tracer_reinspired_option == 'integral':
            return self.co2_volume.data[self.breaths.insp_end] - self.co2_volume.data[self.breaths.insp_start]
        elif self.config.tracer_reinspired_option == 'non-integral':
            return self.config.dead_space_postsensor * np.roll(self.co2_endtidal, 1) * self.btps_factor_insp

    @cached_property
    def co2_raw(self):
        """
        The raw co2 signal, rescaled to be in fractions
        :return: raw co2 signal as a fraction
        """
        file_signal = self.file.signals['co2']
        output_unit = 'si'
        conversion = UnitConversion(input_unit=file_signal.unit, output_unit=output_unit)
        return Signal(
            name='co2 (raw, rescaled)',
            unit=output_unit,
            data=conversion.convert(input=file_signal.data)
        )

    @cached_property
    def co2_trimmed(self):
        """
        After delay-correction, there are parts of the signal which no longer make any sense, or are simply empty. The
        delay correction produces a Signal which indicates which parts of the signals going into the delay correction
        are considered to be valid still. This function reduces the co2 signal down to that size.
        :return: co2 Signal trimmed to size after delay correction.
        """
        return Signal(
            name='CO2 (trimmed)',
            unit=self.co2_delayed.unit,
            data=self.co2_delayed.data[self.valid_after_delay_correction.data]
        )

    @cached_property
    def co2_volume(self):
        """
        Signal array of respired CO2 volume, modified by a threshold rule.
        """
        return Signal(
            name='Tracer Volume (instantaneous)',
            unit=self.flow.unit,
            data=cumtrapz(self.flow.data * self.co2.data, self.flow.time, initial=0)
            # data=cumtrapz(self.flow.data * self.threshold(self.co2.data, 0.005), self.flow.time, initial=0)
        )

    @cached_property
    def config(self):
        """
        We want the config object-attribute of the SpirowareMBW-like classes to only be initialized the first time the
        config is actually used, as opposed to the time when it is initialized.
        :return: A SpirowareConfig object, adapted
        """
        config = SpirowareConfig()
        config.apply_dictionary(self.file.header.translated)
        config.apply_dictionary(self.options)
        config.apply_file_set(self.file.set)

        return config

    @staticmethod
    def cutout(input, breaths):
        """
        This function exists to cut signals to the same length after breath detection by removing all non-breath parts
        of the signal.
        :param input: Signal object from which all non-breath ranges are to be cut out.
        :param breaths: Breaths object containing information about breath starts and ends.
        :return: output: A trimmed copy of the input signal.
        """
        output = np.copy(input)
        # Offset is a variable which changes after each breaths "in-between-spaces" are removed sequentially.
        offset = 0
        for i in breaths.numbers:
            # Part of the original signal to be copied
            copy_range = breaths.range(i)
            # Part of the output signal this range will be pasted into
            paste_range = np.arange(start=offset, stop=breaths.length(i) + offset).astype(int)
            output[paste_range] = input[copy_range]
            # New offset, which begins after the signal just pasted
            offset += breaths.length(i)
        # After all breaths have been moved, remove the tail end of the signal
        output = output[:int(offset)]

        return output

    @classmethod
    def delay_correction(cls, input, delay_times):
        """
        Function to combine an array of signal with an array of delay times to return a new version of the input signal
        which is delay corrected by the values in delay_times.
        :param input: input numpy array to be delayed
        :param delay_times: input numpy array
        :return:
        """
        # Sampling time step
        dt = 1 / 200
        # Calculate the delays as integer multiples of the sampling time step
        delay_indices = cls.round_half_up_array_safe(delay_times / dt)
        length = input.size
        indices = np.arange(length)
        output = np.zeros(length)

        for j in range(1, -1, -1):
            # Calculate the indices where the old array entries will land, making sure the arrays don't loop around
            index_new = indices + j - delay_indices
            index_new = cls.contain_within(index_new, [0, length - 1]).astype(int)

            # Calculate the indices in the old signals which correspond to the new indices
            index_old = indices + j
            index_old = cls.contain_within(index_old, [0, length - 1]).astype(int)

            # Shifting the signals
            output[index_new] = input[index_old]

        return output

    def delay_resync(self):
        """
        This function takes as an argument the four Signals flow, co2, o2, and mmside. It attempt to sequentially align
        them to provide the best possible synchronisation values.
        For CO2: Analogous to how Spiroware calculates its own default delay values, we look at each inspiration,
        calculate the moment when enough air has been inspired for it to have traversed the post-cap deadspace, and then
        calculate by how much we would need to shift the co2 signal for its half-rise point to align with this moment.
        For O2 and MMside: We calculate which shift would produce best alignment (via x-correlation, or noise reduction)
        :return: Returns the delay "deltas", i.e. how much a signal should be shifted, for each of the co2, o2 and
        mmside Signals, for best synchronisation with the
        """
        # Sampling time step (5ms)
        dt = 0.005
        # Range (backwards and forwards) in indices that is checked for synchronisation quality if shifted there (1s).
        sync_dist = 50

        # Relevant signals
        co2 = self.co2_filtered.data
        o2 = self.o2_spedup.data
        mmside = self.mmside_spedup.data

        # CO2 delay can be computed or determined in three ways: a) 'fixed', where the CO2 delay is forced to a value
        # determined in the config, b) 'calibration', where the value from the header or failing that the default
        # delay is taken, and c) 'spiroware', where LungSim tries to compute the CO2 delay analogous to how Spiroware
        # performs its delay time calibration.
        delay_co2 = self.config.delay_co2_resync
        difference = self.config.delay_co2_resync - self.config.delay_co2
        co2 = np.roll(
            a=co2,
            shift=int(difference / dt)
        )

        # O2 delay
        # We cut the co2 signal on both sides to allow comparison of full arrays even when O2 and MMside are shifted.
        indices_co2 = np.arange(start=sync_dist, stop=co2.size-sync_dist)
        co2_short = co2[indices_co2]
        sync_range = np.arange(start=-sync_dist, stop=sync_dist+1)
        normalization = np.zeros(sync_range.size)
        noise = np.zeros(sync_range.size)

        # Calculate the "noise" of the gas sum of o2 and co2, as quantified by the sum of the absolute value of the
        # first discrete differential. When that is minimal, visual noise in the signal is minimal. (I made this up, but
        # it works, so why not use it).
        for delay in sync_range:
            indices_o2 = indices_co2 + delay
            o2_short = o2[indices_o2]
            index = delay - sync_range[0]
            evaluation = co2_short + o2_short
            normalization[index] = self.spikyness(co2_short) + self.spikyness(o2_short)
            noise[index] = self.spikyness(evaluation)

        noise = noise / normalization

        function_o2 = noise
        function_o2 = (function_o2 - np.max(function_o2)) / (np.min(function_o2) - np.max(function_o2))

        o2_delta_index = np.argmin(noise)
        delay_o2 = self.config.delay_o2 + sync_range[o2_delta_index] * dt
        o2 = np.roll(
            a=o2,
            shift=-int(sync_range[o2_delta_index])
        )

        o2_short = o2[indices_co2]
        # Calculate the "noise" of the inert portion molar mass, as quantified by the sum of the absolute value of the
        # first discrete differential. When that is minimal, visual noise in the signal is minimal.
        for delay in sync_range:
            indices_mmside = indices_co2 + delay
            mmside_short = mmside[indices_mmside]
            index = delay - sync_range[0]
            evaluation = mmside_short - self.config.mass['co2'] * co2_short #- self.config.mass['o2'] * o2_short
            normalization[index] = self.spikyness(mmside_short) + self.spikyness(self.config.mass['co2'] * co2_short)
            noise[index] = self.spikyness(evaluation)

        noise = noise / normalization

        function_mm = noise
        function_mm = (function_mm - np.max(function_mm)) / (np.min(function_mm) - np.max(function_mm))

        mmside_delta_index = np.argmin(noise)
        delay_mmside = self.config.delay_mmside + sync_range[mmside_delta_index] * dt

        return delay_co2, delay_o2, delay_mmside, function_o2, function_mm

    def delay_resync_2(self):
        """
        This function takes as an argument the four Signals flow, co2, o2, and mmside. It attempt to sequentially align
        them to provide the best possible synchronisation values.
        For CO2: Analogous to how Spiroware calculates its own default delay values, we look at each inspiration,
        calculate the moment when enough air has been inspired for it to have traversed the post-cap deadspace, and then
        calculate by how much we would need to shift the co2 signal for its half-rise point to align with this moment.
        For O2 and MMside: We calculate which shift would produce best alignment (via x-correlation, or noise reduction)
        :return: Returns the delay "deltas", i.e. how much a signal should be shifted, for each of the co2, o2 and
        mmside Signals, for best synchronisation with the
        """
        # Sampling time step (5ms)
        dt = 0.005
        # Range (backwards and forwards) in indices that is checked for synchronisation quality if shifted there (1s).
        sync_dist = 50

        # Relevant signals
        co2 = self.co2_filtered.data
        o2 = self.o2_spedup.data
        mmside = self.mmside_spedup.data

        # Shortcut
        shift1 = self.config.delay_co2_resync - self.config.delay_co2

        # Shift the co2 signal so that the correlations happen directly with the shifted signal
        co2 = np.roll(
            a=co2,
            shift=int(shift1 / dt)
        )

        shift2, _, function_o2 = self.signal_correlate(co2, o2, sync_dist)
        _, shift3, function_mm = self.signal_correlate(co2, mmside, sync_dist)

        delay_co2 = self.config.delay_co2_resync
        delay_o2 = self.config.delay_o2 - shift2 * dt
        delay_mmside = self.config.delay_mmside - shift3 * dt

        function_o2 = np.flip(function_o2, 0)
        function_o2 = (function_o2 - np.max(function_o2)) / (np.min(function_o2) - np.max(function_o2))

        function_mm = np.flip(function_mm, 0)
        function_mm = (function_mm - np.min(function_mm)) / (np.max(function_mm) - np.min(function_mm))

        return delay_co2, delay_o2, delay_mmside, function_o2, function_mm

    def delay_resync_3(self):
        """
        This function takes as an argument the four Signals flow, co2, o2, and mmside. It attempt to sequentially align
        them to provide the best possible synchronisation values.
        For CO2: Analogous to how Spiroware calculates its own default delay values, we look at each inspiration,
        calculate the moment when enough air has been inspired for it to have traversed the post-cap deadspace, and then
        calculate by how much we would need to shift the co2 signal for its half-rise point to align with this moment.
        For O2 and MMside: We calculate which shift would produce best alignment (via x-correlation, or noise reduction)
        :return: Returns the delay "deltas", i.e. how much a signal should be shifted, for each of the co2, o2 and
        mmside Signals, for best synchronisation with the
        """
        # Sampling time step (5ms)
        dt = 0.005
        # CO2 -based breath detection threshold (2%). We use this BD because it's more reliable than Spiroware's method.
        threshold = 0.02
        # Number of sampling points after breath beginning in which the maximum of the CO2 signal is determined (200ms)
        co2_beginning = 40
        sync_dist = 50

        # Relevant signals
        flow = self.flow_filtered.data
        co2 = self.co2_filtered.data
        o2 = self.o2_spedup.data
        mmside = self.mmside_spedup.data

        # Robust breath detection
        try:
            breaths = self.breath_detection_co2(flow=flow, co2=co2, threshold=threshold)
        except LungSimException:
            raise LungSimException(output='Breath detection during re-synchronisation failed.')

        # For each breath, compute the delay between the air reaching the CO2 sensor, and half-rise time
        co2_deltas = np.array([])
        o2_deltas = np.array([])
        mmside_deltas = np.array([])
        back_shift = 40
        front_shift = 200
        for i in breaths.numbers:
            try:
                start = breaths.insp_start[i] - back_shift
                # Find the closest moment when the CO2 signal dips below the halfway-point of the min-max span. The
                # max is defined as the maximum of the CO2 signal within the first 200ms seconds of the breath.
                co2_breath = co2[start:]
                co2_target = 0.5 * np.max(co2_breath[:co2_beginning])
                co2_reached = self.first_nonzero(co2_breath < co2_target)
                # Append the difference between the above two results to the co2_deltas collection

                # O2 delay
                o2_breath = o2[start:]
                o2_sample = o2[start:start + front_shift]
                top = np.max(o2_sample)
                bottom = np.min(o2_sample)

                o2_target = bottom + 0.5 * (top - bottom)
                o2_reached = self.first_nonzero(o2_breath > o2_target)

                # MMside
                mmside_breath = mmside[start:]
                mmside_sample = mmside[start:start + front_shift]
                top = np.max(mmside_sample)
                bottom = np.min(mmside_sample)
                start = mmside_sample[0]
                mmside_target = bottom + 0.5 * (top - bottom)
                if np.abs(top - start) > np.abs(bottom - start):
                    mmside_reached = self.first_nonzero(mmside_breath > mmside_target)
                else:
                    mmside_reached = self.first_nonzero(mmside_breath < mmside_target)

                co2_deltas = np.append(arr=co2_deltas, values=co2_reached)
                o2_deltas = np.append(arr=o2_deltas, values=o2_reached)
                mmside_deltas = np.append(arr=mmside_deltas, values=mmside_reached)

            except:
                # If one or the other of the two moments cannot be found, ignore breath
                continue

        if co2_deltas.size == 0 or o2_deltas.size == 0 or mmside_deltas.size == 0:
            msg = 'Resynchronisation impossible: Insufficient normal breaths for CO2-delay.'
            raise LungSimException(output=msg)

        # The delay for the CO2 signal is then the original delay minus the one subtracted above plus the one
        # calculated by shifting. The original minus the shift above is equivalent to the mod of the original delay.
        co2_diff = self.config.delay_co2_resync - self.config.delay_co2
        delay_co2 = self.config.delay_co2_resync
        co2_deltas = co2_deltas - np.round(co2_diff / dt)
        delay_o2 = self.config.delay_o2 + np.median(o2_deltas - co2_deltas) * dt
        delay_mmside = self.config.delay_mmside + np.median(mmside_deltas - co2_deltas) * dt

        # From the collection of observed delay times, create a quadratic error function centred around the peak 
        # of the histogram for O2
        f, _ = np.histogram(
            a=o2_deltas-co2_deltas,
            bins=np.arange(start=-sync_dist-0.5, stop=sync_dist+1.5)
        )
        j = np.arange(start=-sync_dist, stop=sync_dist+1)
        g = j - j[np.argmax(f)]
        h = g*(-g)
        function_o2 = (h - np.min(h)) / (np.max(h)-np.min(h))

        # Same process as above for MM
        f, _ = np.histogram(
            a=mmside_deltas-co2_deltas,
            bins=np.arange(start=-sync_dist-0.5, stop=sync_dist+1.5)
        )
        j = np.arange(start=-sync_dist, stop=sync_dist+1)
        g = j - j[np.argmax(f)]
        h = g*(-g)
        function_mm = (h - np.min(h)) / (np.max(h)-np.min(h))

        return delay_co2, delay_o2, delay_mmside, function_o2, function_mm

    @cached_property
    def delay_times(self):
        """
        The delay time property returns a DelayTimes object, which contains just the delay times for the signals which
        need to be delay corrected (CO2, O2, and MMss). When the delay_times object is created, a very heavy computation
        uses the sidestream flow signal to compute the delay times. This must happen only once.
        In some cases it is necessary to recompute better base delay values than the ones found in the header, or some
        mean provided in the options. In that case, this function will create a separate analysis of the same type with
        the default values, find ideal adjusted values for best synchronisation, and provide them here instead of the
        default. Needless to say this is also very slow and computationally heavy.
        :return: delay times dictionary with delays for each signal in [s]
        """
        if self.config.resync_option:
            # Make a copy of the analysis so far
            decoy = deepcopy(self)
            # Adjust the config to accept unchanged delay values
            decoy.config.resync_option = False
            # Extract the relevant signals for resynchronisation from the decoy analysis
            self.d_co2_1, self.d_o2_1, self.d_mmside_1, f1_o2, f1_mm = decoy.delay_resync()
            self.d_co2_2, self.d_o2_2, self.d_mmside_2, f2_o2, f2_mm = decoy.delay_resync_2()
            self.d_co2_3, self.d_o2_3, self.d_mmside_3, f3_o2, f3_mm = decoy.delay_resync_3()
            delays = np.arange(start=-50, stop=51)
            # Translate the outputs of the delay functions into single delay values
            delay_co2 = np.mean((self.d_co2_1, self.d_co2_2, self.d_co2_3))
            delay_o2 = self.config.delay_o2 + delays[np.argmax(f1_o2 * f2_o2 * f3_o2)] * 0.005
            delay_mmside = self.config.delay_mmside + delays[np.argmax(f1_mm * f2_mm * f3_mm)] * 0.005
            # If the three delays are not in logical order, delay correction becomes impossible or illogical.
            if delay_co2 < delay_o2 < delay_mmside:
                self.config.delay_co2 = delay_co2
                self.config.delay_o2 = delay_o2
                self.config.delay_mmside = delay_mmside
            elif delay_co2 < delay_o2 and self.file.method == 'n2_mbw':
                self.config.delay_co2 = delay_co2
                self.config.delay_o2 = delay_o2
                self.config.delay_mmside = delay_o2 + 0.1
            else:
                # Take the default delay values if the above fails.
                raise LungSimException(output='Re-synchronisation failed. Delays inconsistent.')
                #pass

            decoy.reset()

        if self.config.delay_option == 'dynamic':
            # Calculate delay per sampling time step
            return self.delay_time_calculation_dynamic(
                flowside=self.flowside_raw,
                static_delay_co2=self.config.delay_co2,
                static_delay_o2=self.config.delay_o2,
                static_delay_mmside=self.config.delay_mmside
            )
        elif self.config.delay_option == 'static':
            # Return delay as a static value
            return {
                'mmside': np.full_like(self.flow_sidestream_corrected.data, self.config.delay_mmside),
                'o2': np.full_like(self.flow_sidestream_corrected.data, self.config.delay_o2),
                'co2': np.full_like(self.flow_sidestream_corrected.data, self.config.delay_co2)
            }

    @classmethod
    def delay_time_calculation_dynamic(cls, flowside, static_delay_co2, static_delay_o2, static_delay_mmside):
        """
        Function to calculate the delay times of the signals CO2, O2 and MMside based on the static delay values and the
        side stream flow. This is written to mimic the Spiroware method and matches exactly. It is however quite slow,
        and therefore requires some numba compilation to run decently fast.
        This function estimates the volume between the sensors, and then attempts to find the delays.
        :param flowside: Sidestream flow signal [L/s]
        :param static_delay_co2: Calibration value for CO2 delay [s]
        :param static_delay_o2: Calibration value for O2 delay [s]
        :param static_delay_mmside: Calibration value for MMss delay [s]
        :return: A delay values dictionary containing the delays for all three signals.
        """
        # Sampling time step (set here by default to be 1/200. This function will only ever be used by Spiroware, which
        # by default has 200 Hz as its sampling frequency.
        dt = 1 / 200

        # The number of compartments in the "delay line", i.e. the number of subdivisions of volume by which the
        # Nafion tube is modelled. This is an unnecessary introduction in this calculation and leads to errors described
        # in more detail below.
        compartments = 256

        # Nominal flow, i.e. average flow during the calibration. This is the flow value by which the static delay
        # values were established. It can therefore be used to estimate the volume of tube between the sensors.
        # Volume = Time * Flow
        # The default value for this is 200ml/minute.
        flow_nominal = 200e-6 / 60

        # Volume estimate between CO2 sensor and O2 sensor
        volume_partial = flow_nominal * (static_delay_o2 - static_delay_co2)
        # Volume estimate between CO2 sensor and MMside sensor
        volume_total = flow_nominal * (static_delay_mmside - static_delay_co2)
        # Volume of one division of the sidestream tube
        volume_unit = volume_total / compartments

        # Location of O2 sensor as number of complete volume compartments away from CO2 sensor
        location_right = int(np.floor(volume_partial / volume_unit))
        # Partial index of O2 sensor location (non-integer part)
        fraction_o2 = volume_partial / volume_unit - location_right
        # Weights to interpolate to location of O2 sensor
        weight_left = 1 - fraction_o2
        weight_right = fraction_o2
        location_left = location_right - 1

        # Number of compartments pushed out at each time point
        push_array = - flowside.data * (dt / volume_unit)
        # Rounded down
        push_integer = np.floor(push_array).astype(int)
        # Fraction of one of the two array entries to be interpolated
        left_fraction = push_array - push_integer
        # Fraction of the other
        right_fraction = 1 - left_fraction

        # Initializing array used for delay line interpolation
        numbers = np.arange(compartments + 1)
        # Find unique values of the push_integer array. Changed the way these arrays are computed. Basically there
        # were only a couple possible options for the push_integer array, and every time point would compute one of
        # these, and then use it. But we can prepare those few arrays, and then select from them in the for loop.
        # The list of unique values
        push_unique = np.arange(min(push_integer), max(push_integer) + 1)
        # push_alt transforms the array of integer pushes into the indexes corresponding to the prepared arrays
        push_indices = push_integer - push_unique[0]
        # a stack of arrays that will become the options for pushing
        push_stack = np.tile(numbers, (len(push_unique), 1))
        # a stack of possible push indexes
        index_right = push_stack - np.reshape(push_unique, (len(push_unique), 1))
        # the corresponding neighboring indexes
        index_left = index_right - 1
        # Eliminating out of bounds indexes
        index_right = cls.contain_within(index_right, [0, compartments])
        index_left = cls.contain_within(index_left, [0, compartments])
        # We only modify the elements of the delay line from the second index. We therefore shorten this array for later
        index_right = index_right[:, 1:]
        index_left = index_left[:, 1:]

        # Initialized delay model. Array of delays to each location in the tube as defined by the volume compartments
        # The use of one more compartment than initially declared and used to calculate the location of the O2 sensor
        # introduces another error (#2)
        # How long would it take for the initial flow signal entry to reach each point in the tube?
        delay_model = np.zeros(compartments + 1) - np.arange(0, compartments + 1) * volume_unit / flowside.data[0]

        # We send the necessary values to a numba compiled function. We use nan_to_num to make sure the values are
        # compatible.
        delay_values = cls.delay_time_numba(
            model=np.nan_to_num(delay_model),
            push=np.nan_to_num(push_indices),
            index_r=np.nan_to_num(index_right),
            index_l=np.nan_to_num(index_left),
            fraction_r=np.nan_to_num(right_fraction),
            fraction_l=np.nan_to_num(left_fraction),
            loc_r=np.nan_to_num(location_right),
            loc_l=np.nan_to_num(location_left),
            dt=dt,
            length=flowside.length
        )

        # Extract delay values
        dmmside = delay_values[:, 0] - dt
        # Interpolate between the two O2 sensor location compartments
        do2 = delay_values[:, 1] * weight_left + delay_values[:, 2] * weight_right - dt

        return {
            'mmside': dmmside + static_delay_co2,
            'o2': do2 + static_delay_co2,
            'co2': np.ones(flowside.length) * static_delay_co2
        }

    @staticmethod
    @njit(cache=True)
    def delay_time_numba(model, push, index_r, index_l, fraction_r, fraction_l, loc_r, loc_l, dt, length):
        """
        Attempt to make very fast portion of the signal processing using compiled function. This is the heaviest part of
        the Spiroware algorithm, and improving its performance leads to great improvement of performance overall.
        We use the @njit decorator because we want to ensure optimal performance, and cache the compilation for better
        future performance.
        :param model: 257x1 array of a model of the Nafion tube going from the CO2 sensor to the side stream 
            sensors. Contains the time of delay of every time point of how long it takes the sample to reach this point.
        :param index_l: 256xn array which contains the indices that correspond to a certain push speed for the left 
            side of the interpolation. n corresponds to the n different possible integer push speeds encountered in this
            measurement.
        :param index_r: same as index_left, but the right side of the interpolation. 
        :param push: timex1 array of the index of the line of right fraction. 0 corresponds to the first line of 
            right_fraction, which corresponds to copying from the same indices for the delay line. 
        :param fraction_r: non-integer contribution of the shift of the delay line. Fraction for each time point how
            much the right element of the interpolation contributes to the new delay line.
        :param fraction_l: Fraction at each time point how much the left element of the interpolation contributes to
             the new time points delay time
        :param loc_r: index of the compartment to the left of the O2 sensor
        :param loc_l: index of the compartment to the right of the O2 sensor
        :param dt: time step, usually 1/200 s
        :param length: length of the signal to be delay-corrected
        :return: delay values the two contributing compartments of O2 and for MMss
        """
        # Initializing delay arrays instead of three separate assignments, we can assign the three points of interest
        # along the delay aline to a single array 2D array entry
        delay_values = np.zeros((length, 3))

        for i in np.arange(length):
            # Interpolate the model forward.
            model[1:] = fraction_l[i] * model[index_l[push[i]]] + fraction_r[i] * model[index_r[push[i]]] + dt
            # Extract the delay values from the model at three specific locations
            delay_values[i] = (model[-1], model[loc_l], model[loc_r])

        return delay_values

    @cached_property
    def drop_initial_expected(self):
        indices = self.washout_breaths
        args = dict()
        args['variable'] = 'frc_ao'
        args['method'] = 'critical'
        args['concentration'] = 2.5
        frc_raw = self.calculate_output(args)
        frc_raw = self.frc_ao[self.washout_breaths[0]]
        ds_phys = self.physio_dead_space
        ds_total = self.physio_dead_space + self.machine_dead_space
        frc = frc_raw - ds_phys
        vt = self.volume_tidal[indices]
        k = (vt - ds_total) / (frc + vt)
        k[k < 0] = 0
        return 1-k[self.washout_breaths[0]]

    @cached_property
    def drop_initial_realized(self):
        return self.tracer_normalized[self.washout_breaths[0]]

    @staticmethod
    def end_tidal(signal, flow, time, window):
        """
        Function to calculate a sample of signal at the end (defined by window) of the breath data given in signal and
        flow
        :param signal: CO2, O2, MMss data of one breath
        :param flow: flow data of one breath
        :param time: time series corresponding to flow data
        :param window: tuple in the range [0, 1) which defines in which percentage of volume we sample the signal data.
        :return: flow-weighted average of the data in signal between window[0] and window[1] of expired volume
        """
        # If the breath is extremely short, it may cause some bugs if these conditions are not handled. We do not want
        # to raise an exception, as one zero-breath does not mean the rest of the measurement cannot be calculated.
        if signal.size == 0:
            return None

        if signal.size == 1:
            return np.mean(signal)

        # Calculate end-tidal quantities in defined windows of each breath
        volume = cumtrapz(-flow, time, initial=0)
        # Calculate volume integral of  gas
        signal_volume = cumtrapz(-flow * signal, time, initial=0)
        # Two different interpretations are possible for this line: Last entry or maximum entry.
        total_volume = np.max(volume)
        # Interval to be considered for the averaging of values
        # It seems that this method performs the best among those tested. (<99.95% accuracy)
        interval_start = (volume > window[0] * total_volume).nonzero()[0][0]
        interval_end = (volume > window[1] * total_volume).nonzero()[0][0]
        # Condition for badly detected breaths which are only a few indices long. Prevents NaN values from appearing
        # in the concentrations, but may be prone to difference from Spiroware.
        if interval_start == interval_end:
            return signal_volume[interval_start] / volume[interval_start]

        else:
            # flow weighted concentration average in interval
            volume_window = volume[interval_end] - volume[interval_start]
            signal_window = signal_volume[interval_end] - signal_volume[interval_start]
            return signal_window / volume_window

    @staticmethod
    def filter_lowpass_spiroware(input, sampling_rate, cutoff, window, initial_condition=True):
        """
        Returns a lowpass-filtered copy of the signal contained within the data of self. The filter order here being +1
        of the sampling rate * window is simply a reproduction of Spiroware code. It is unclear why this kind of filter
        was chosen, but this works to reproduce the behavior of Spiroware.
        :param input: Signal object to be filtered
        :param sampling_rate: sampling rate of the signal in [Hz]
        :param cutoff: cutoff frequency of the lowpass filter in [Hz]
        :param window: filter window in [s]
        :param initial_condition: boolean which controls whether the filter tries to find initial conditions to match
            initial array value. Only False for flow signal.
        :return: output: A lowpass-filtered copy of the signal
        """
        if input.size == 0:
            raise LungSimException(output='Signal length too short. Reduced to 0 length when filtering.')

        output = np.copy(input)
        filter_order = int(sampling_rate * window)
        filter_design = firwin(numtaps=filter_order + 1, cutoff=cutoff, pass_zero='lowpass', fs=sampling_rate)

        if initial_condition:
            filter_zi = lfilter_zi(b=filter_design, a=1)
            output, _ = lfilter(b=filter_design, a=1, x=output, zi=filter_zi*output[0])
        else:
            output = lfilter(b=filter_design, a=1, x=output)

        return output

    @cached_property
    def flow(self):
        """
        Final flow signal of a Spiroware algorithm. The flow signal of a Spiroware Algorithm after filtering, trimming,
        BTPS correction, and compensation for the side stream flow.
        :return: Flow signal completely processed
        """
        if self.config.cutout_option:
            return Signal(
                name='Flow (final)',
                unit=self.flow_filtered.unit,
                data=self.cutout(
                    input=self.flow_filtered.data,
                    breaths=self.breaths_raw
                )
            )
        else:
            return self.flow_filtered

    @cached_property
    def flow_filtered(self):
        """
        Flow signal after filtering, trimming, BTPS correction, and compensation for the side stream flow. The filter
        option controls whether any filtering is performed. If the option is off, filtering is skipped and the previous
        signal in the processing stream is returned.
        :return: Flow signal filtered
        """
        if self.config.flow_filter_option:
            return Signal(
                name='Flow (filtered)',
                unit=self.flow_btps.unit,
                data=self.filter_lowpass_spiroware(
                    input=self.flow_btps.data,
                    sampling_rate=self.flow_btps.sampling_rate,
                    cutoff=self.config.flow_filter_cutoff,
                    window=self.config.flow_filter_window,
                    initial_condition=False
                )
            )
        else:
            return self.flow_btps

    @cached_property
    def flow_btps(self):
        """
        Flow signal after trimming, BTPS correction and compensation for the side stream flow. The BTPS option controls
        whether any filtering is performed. If off, the btps correction is skipped and the delay-trimmed signal is
        returned. BTPS stands for body temperature and pressure, saturated conditions. This correction multiplies any
        incoming flow signal by how much inspired air will expand in volume inside the lungs, and multiplies any
        outgoing flow by how much it has decreased in volume since being expired.
        :return: Flow signal BTPS-corrected.
        """
        if self.config.btps_option:
            return Signal(
                name='Flow (BTPS)',
                unit=self.flow_trimmed.unit,
                data=self.btps_correction(
                    flow=self.flow_trimmed.data,
                    factor_insp=self.btps_factor_insp,
                    factor_exp=self.btps_factor_exp
                )
            )
        else:
            return self.flow_trimmed

    @cached_property
    def flow_trimmed(self):
        """
        After delay-correction, there are parts of the signal which no longer make any sense, or are simply empty. The
        delay correction produces a Signal which indicates which parts of the signals going into the delay correction
        are considered to be valid still. This function reduces the flow signal down to that size.
        :return: Flow Signal trimmed to size after delay correction.
        """
        return Signal(
            name='Flow (trimmed)',
            unit=self.flow_sidestream_corrected.unit,
            data=self.flow_sidestream_corrected.data[self.valid_after_delay_correction.data]
        )

    @cached_property
    def flow_sidestream_corrected(self):
        """
        The flow signal of a Spiroware analysis typically has to be corrected for the amount of air that is sampled from
        the main stream into the sidestream tube. The convention of that sidestream flow is such that in order to add
        the missing volume to the inspirations of the main flow signal, the negative side stream flow must be
        subtracted.
        :return:
        """
        if self.config.sampleflow_option:
            return Signal(
                name='Flow (sidestream corrected)',
                unit=self.flow_raw.unit,
                data=self.flow_raw.data + self.flowside_raw.data
            )
        else:
            return self.flow_raw

    @cached_property
    def flow_raw(self):
        """
        The raw flow signal is the signal extracted from the File containing the raw data itself. The raw data, in order
        to be easily used in the rest of the analysis, should have its units converted.
        :return: A rescaled version of the raw flow signal
        """
        file_signal = self.file.signals['flowmain']
        output_unit = 'm^3/s'
        conversion = UnitConversion(input_unit=file_signal.unit, output_unit=output_unit)
        return Signal(
            name='Flow (raw, rescaled)',
            unit=output_unit,
            data=conversion.convert(input=file_signal.data)
        )

    @cached_property
    def flowside_raw(self):
        """
        The raw flowside signal is the rescaled signal of the flow which traverses the sidestream sampling tube. In
        order for the signal to be easily used, its units have to be converted here. Convention is that negative flow
        is flow away from the flowhead towards the Exhalyzer D box.
        :return: A rescaled version of the raw sidestream flow signal
        """
        file_signal = self.file.signals['flowside']
        output_unit = 'm^3/s'
        conversion = UnitConversion(input_unit=file_signal.unit, output_unit=output_unit)
        return Signal(
            name='Sidestream Flow (raw, rescaled)',
            unit=output_unit,
            data=conversion.convert(input=file_signal.data)
        )

    @staticmethod
    def fractional_crossing(array, threshold):
        diff = array - threshold
        if not (diff < 0).any():
            return None
        x2 = (diff < 0).nonzero()[0][0]
        x1 = x2 - 1
        x = abs(diff[x1])+abs(diff[x2])
        f1 = abs(diff[x1])/x
        f2 = abs(diff[x2])/x
        return f2*x1 + f1*x2

    @cached_property
    def lcr_crude(self):
        """
        :return: Lung clearance ratio based on global averages.
        """
        # Calculate FRCao, at 2.5% threshold. This is the standard MBW outcome."
        args = dict()
        args['variable'] = 'frc_ao'
        args['method'] = 'critical'
        args['concentration'] = 2.5
        frc_raw = self.calculate_output(args)
        # Calculate vt as the mean tidal volume (mean of insp and exp) over washout to the test-end criterion."
        args = dict()
        args['variable'] = 'volume_tidal'
        args['method'] = 'mean'
        args['indices'] = 'washout_breaths_to_critical(2.5)'
        vt = self.calculate_output(args)
        # The variable frc is here now the alveolar (non-dead space) FRC.
        frc = frc_raw - self.physio_dead_space
        # Dead space is all measured physiological and machine dead space.
        ds_total = self.physio_dead_space + self.machine_dead_space
        # Calculate breaths predicted
        k = (vt-ds_total)/(frc+vt)
        bp = log(0.025)/log(1-k)
        # Calculate fractional(!) end of test (i.e. WO can end after e.g. 7.7 breaths instead of 8)
        wo_end = self.washout_end(2.5)
        bm = self.fractional_crossing(self.tracer_normalized[wo_end-1:wo_end+1], 0.025)
        bm = wo_end + bm - self.washout_breaths[0]
        # The output is the (fractional) measured breaths over (fractional) predicted breaths.
        return bm/bp

    @cached_property
    def lcr_granular(self):
        """
        :return: Lung clearance ratio based on individual breaths
        """
        indices = self.washout_breaths
        args = dict()
        args['variable'] = 'frc_ao'
        args['method'] = 'critical'
        args['concentration'] = 2.5
        frc_raw = self.calculate_output(args)
        ds_phys = self.physio_dead_space
        ds_total = self.physio_dead_space + self.machine_dead_space
        frc = frc_raw - ds_phys
        vt = self.volume_tidal[indices]
        k = (vt - ds_total) / (frc + vt)
        f = np.cumprod(1-k)
        bp = self.fractional_crossing(f, 0.025) + 1
        if bp is None:
            LungSimException(msg='No baseline test end detected.')
        wo_end = self.washout_end(2.5)
        bm = self.fractional_crossing(self.tracer_normalized[wo_end-1:wo_end+1], 0.025)
        bm = wo_end + bm - self.washout_breaths[0]
        return bm/bp

    @cached_property
    def lcr_fractional(self):
        """
        :return: Lung clearance ratio based on individual breaths
        """
        indices = self.washout_breaths
        args = dict()
        args['variable'] = 'frc_ao'
        args['method'] = 'critical'
        args['concentration'] = 2.5
        frc_raw = self.calculate_output(args)
        ds_phys = self.physio_dead_space
        ds_total = self.physio_dead_space + self.machine_dead_space
        frc = frc_raw - ds_phys
        vt = self.volume_tidal[indices]
        k = (vt - ds_total) / (frc + vt)
        r = np.cumsum(k)
        f = np.cumprod(1-k)
        bp = self.fractional_crossing(f, 0.025)
        if bp is None:
            LungSimException(msg='No baseline test end detected.')
        wo_end = self.washout_end(2.5)
        bm = self.fractional_crossing(self.tracer_normalized[wo_end-1:wo_end+1], 0.025)
        bm = wo_end + bm - self.washout_breaths[0] - 1
        first = int(np.floor(bp))
        second = first + 1
        vtp = (second-bp)*f[first] + (bp-first)*f[second]
        first = int(np.floor(bm))
        second = first + 1
        vtm = (second-bm)*f[first] + (bm-first)*f[second]
        lcr = np.log10(vtp/vtm)
        return lcr

    @cached_property
    def lcr_granular_k(self):
        """
        :return: Lung clearance ratio based on individual breaths
        """
        indices = self.washout_breaths
        args = dict()
        args['variable'] = 'frc_ao'
        args['method'] = 'critical'
        args['concentration'] = 2.5
        frc_raw = self.calculate_output(args)
        ds_phys = self.physio_dead_space
        ds_total = self.physio_dead_space + self.machine_dead_space
        frc = frc_raw - ds_phys
        vt = self.volume_tidal[indices]
        k = (vt - ds_total) / (frc + vt)
        k[k < 0] = 0

        d = np.log(1-k)/np.log(0.025)
        dc = np.cumsum(d)

        wo_end = self.washout_end(2.5)
        bm = self.log_crossing(self.tracer_normalized[wo_end-1:wo_end+1], 0.025)
        bm = wo_end + bm - self.washout_breaths[0] - 1
        first = int(np.floor(bm))
        second = first + 1
        lcr = (second-bm)*dc[first] + (bm-first)*dc[second]
        return lcr

    @cached_property
    def lcr_x(self):
        """
        :return: Lung clearance ratio based on individual breaths
        """
        indices = self.washout_breaths
        args = dict()
        args['variable'] = 'frc_ao'
        args['method'] = 'critical'
        args['concentration'] = 25.0
        frc_raw = self.calculate_output(args)
        ds_phys = self.physio_dead_space
        ds_total = self.physio_dead_space + self.machine_dead_space
        frc = frc_raw - ds_phys
        vt = self.volume_tidal[indices]
        k = (vt - ds_total) / (frc + vt)
        k[k < 0] = 0

        d = np.log(1-k)/np.log(0.25)
        dc = np.cumsum(d)

        wo_end = self.washout_end(25.0)
        bm = self.log_crossing(self.tracer_normalized[wo_end-1:wo_end+1], 0.25)
        bm = wo_end + bm - self.washout_breaths[0] - 1
        first = int(np.floor(bm))
        second = first + 1
        lcr = (second-bm)*dc[first] + (bm-first)*dc[second]
        return lcr

    @staticmethod
    def log_crossing(array, threshold):
        # Log-transforming the input
        array[array < np.min(np.abs(array))] = np.min(np.abs(array))
        array = np.log(array)
        threshold = np.log(threshold)

        diff = array - threshold
        if not (diff < 0).any():
            return None
        x2 = (diff < 0).nonzero()[0][0]
        x1 = x2 - 1
        x = abs(diff[x1])+abs(diff[x2])
        f1 = abs(diff[x1])/x
        f2 = abs(diff[x2])/x
        return f2*x1 + f1*x2

    @cached_property
    def machine_dead_space(self):
        """
        Total machine dead space in a Spiroware measurement is defined as the dead space between the airway opening /
        mouth and the bypass flow, and is composed of a volume before and after the Spiroware sensors.
        :return:
        """
        return self.config.dead_space_postsensor + self.config.dead_space_presensor

    @cached_property
    def mmside(self):
        """
        mmside signal after the breath-detection related cutting out of in-between breaths sections of the signal.
        :return: mmside signal after cutout
        """
        if self.config.cutout_option:
            return Signal(
                name='MM side (cutout)',
                unit=self.mmside_spedup.unit,
                data=self.cutout(
                    input=self.mmside_spedup.data,
                    breaths=self.breaths_raw
                )
            )
        else:
            return self.mmside_spedup

    @cached_property
    def mmside_delayed(self):
        """
        The raw MMss signal is delay corrected by a somewhat complicated algorithm which is under certain circumstances
        defective, and certainly not maximally precise. However, this is a centrally important function for the correct
        interpretation of raw Spiroware signals. As the CO2 signal, the O2 signal and the MMss signal are measured at
        very different places in the experimental setup, it becomes necessary to align them in time.
        :return: delay-corrected mmside signal
        """
        return Signal(
            name='MMside (delay-corrected)',
            unit=self.mmside_raw.unit,
            data=self.delay_correction(
                input=self.mmside_raw.data,
                delay_times=self.delay_times['mmside']
            )
        )

    @cached_property
    def mmside_filtered(self):
        """
        The mmside signal is filtered with a lowpass filter to smooth it. This leads to signal delay and some loss of 
        data. The Spiroware algorithm's default filters the signal with a harsh lowpass filter with a cutoff of 2Hz, 
        which leads to signal delay especially in the rise and fall moments.
        :return: Filtered mmside signal
        """
        if self.config.mmside_filter_option:
            return Signal(
                name='MM side (filtered)',
                unit=self.mmside_trimmed.unit,
                data=self.filter_lowpass_spiroware(
                    input=self.mmside_trimmed.data,
                    sampling_rate=self.mmside_trimmed.sampling_rate,
                    cutoff=self.config.mmside_filter_cutoff,
                    window=self.config.mmside_filter_window
                )
            )
        else:
            return self.mmside_trimmed

    @cached_property
    def mmside_raw(self):
        """
        The raw side stream molar mass signal (MMss) in kg/mol, rescaled from the file.
        :return:
        """
        file_signal = self.file.signals['mmside']
        output_unit = 'kg/mol'
        conversion = UnitConversion(input_unit=file_signal.unit, output_unit=output_unit)
        return Signal(
            name='Sidestream Molar Mass (raw, rescaled)',
            unit=output_unit,
            data=conversion.convert(input=file_signal.data)
        )
    
    @cached_property
    def mmside_spedup(self):
        """
        The mmside signal has significant time-lag with respect to the CO2 signal, which can lead to signal artifacts if
        left uncorrected. A special filter is therefore applied to it which aims to "speed up" the parts of the
        signal where any change happens. This leads to noise in the signal, but reduces a bit the spikes caused by
        different time constants of the mmside and CO2 sensors. Speed up selected signals to compensate for the rise 
        time delay of the sensors. 
        :return: Sped-up mmside signal
        """
        if self.config.mmside_speed_option:
            return Signal(
                name='MM side (sped-up)',
                unit=self.mmside_filtered.unit,
                data=self.speeding(
                    array=self.mmside_filtered.data,
                    sampling_rate=self.mmside_filtered.sampling_rate,
                    tau=self.config.mmside_filter_tau,
                    butter_cutoff=self.config.speed_butter_cutoff
                )
            )
        else:
            return self.mmside_filtered

    @cached_property
    def mmside_trimmed(self):
        """
        After delay-correction, there are parts of the signal which no longer make any sense, or are simply empty. The
        delay correction produces a Signal which indicates which parts of the signals going into the delay correction
        are considered to be valid still. This function reduces the mmside signal down to that size.
        :return: mmside Signal trimmed to size after delay correction.
        """
        return Signal(
            name='mmside (trimmed)',
            unit=self.mmside_delayed.unit,
            data=self.mmside_delayed.data[self.valid_after_delay_correction.data]
        )

    @cached_property
    def moment_ratio_1(self):
        return self.outcome_divide(self.moment_one, self.moment_zero)

    @cached_property
    def moment_ratio_1_spiroware(self):
        """
        Moment ratio outcome evaluated at the "breath table end", i.e. the third of the three breaths below the
        critical tracer concentration of 2.5%.
        :return: single value of the "moment ratio 1" in Spiroware
        """
        try:
            breath = self.spiroware_breaths[-1] - 1
        except IndexError:
            raise LungSimException(output='Appropriate washout portion not found.')

        return self.moment_ratio_1[breath]

    @cached_property
    def moment_ratio_1_TO6(self):
        """
        Moment ratio outcome evaluated at two breaths before Turnover 6 is reached.
        :return: single value of the "moment ratio 1 TO 6" in Spiroware
        """
        return self.moment_ratio_X_TOY(self.moment_ratio_1, 6)

    @cached_property
    def moment_ratio_1_TO8(self):
        """
        Moment ratio outcome evaluated at two breaths before Turnover 8 is reached.
        :return: single value of the "moment ratio 1 TO 8" in Spiroware
        """
        return self.moment_ratio_X_TOY(self.moment_ratio_1, 8)

    @cached_property
    def moment_ratio_2(self):
        return self.outcome_divide(self.moment_two, self.moment_zero)

    @cached_property
    def moment_ratio_2_spiroware(self):
        """
        Moment ratio outcome evaluated at the "breath table end", i.e. the third of the three breaths below the
        critical tracer concentration of 2.5%.
        :return: single value of the "moment ratio 2" in Spiroware
        """
        try:
            breath = self.spiroware_breaths[-1] - 1
        except IndexError:
            raise LungSimException(output='Appropriate washout portion not found.')

        return self.moment_ratio_2[breath]

    @cached_property
    def moment_ratio_2_TO6(self):
        """
        Moment ratio outcome evaluated at two breaths before Turnover 6 is reached.
        :return: single value of the "moment ratio 2 TO 6" in Spiroware
        """
        return self.moment_ratio_X_TOY(self.moment_ratio_2, 6)

    @cached_property
    def moment_ratio_2_TO8(self):
        """
        Moment ratio outcome evaluated at two breaths before Turnover 8 is reached.
        :return: single value of the "moment ratio 2 TO 8" in Spiroware
        """
        return self.moment_ratio_X_TOY(self.moment_ratio_2, 8)

    def moment_ratio_X_TOY(self, moment_ratio_X, Y):
        """
        Function to perform the calculations shared between the various moment ratio functions. We input a moment ratio
        array and a turnover value, and the function computes the Spiroware way of outputting moment ratio outcomes.
        The weird - 1 and -2 here are there to produce Spiroware values, but they make little sense.
        :param moment_ratio_X: array of breath-wise moment ratio outcomes
        :param turnover_array: array of turnover values breathwise.
        :param Y: turnover value at which the moment ratio array is evaluated.
        :return: single moment ratio outcome
        """
        # Find the breath table end. If that proves impossible, return error
        try:
            breath_end = self.spiroware_breaths[-1] + 1
        except IndexError:
            raise LungSimException(output='Appropriate washout portion not found.')

        # Find the breaths which have a turnover value greater than value Y
        above_TOY = (self.lci_ao > Y).nonzero()[0]
        if above_TOY.size == 0:
            raise LungSimException(output='Turnover {} not reached.'.format(Y))

        # To evaluate the moment ratio, select the breath 2 breaths away from the first breath to be above the target
        selection = above_TOY[0] - 2
        if breath_end < selection:
            raise LungSimException(output='Turnover {} not reached before MBW end.'.format(Y))

        return moment_ratio_X[selection]

    @cached_property
    def moment_zero(self):
        """
        # TODO: Complete moment analysis description
        Moment analysis.
        :return: Moment zero
        """
        tracer = np.copy(self.tracer_normalized)
        tracer[self.prephase_breaths] = 1
        return cumtrapz(tracer, self.lci_ao)

    @cached_property
    def moment_one(self):
        """
        Moment analysis.
        :return: Moment one
        """
        tracer = np.copy(self.tracer_normalized)
        tracer[self.prephase_breaths] = 1
        return cumtrapz(tracer * self.lci_ao, self.lci_ao)

    @cached_property
    def moment_two(self):
        """
        Moment analysis.
        :return: Moment two
        """
        tracer = np.copy(self.tracer_normalized)
        tracer[self.prephase_breaths] = 1
        return cumtrapz(tracer * self.lci_ao * self.lci_ao, self.lci_ao)

    @cached_property
    def nlci(self):
        """
        Calculates the normalizes Lung Clearance Index. Defined as the cumulative expired volume of CO2 over FRC(ao)
        """
        return self.outcome_divide(self.co2_cev, self.frc_ao)

    @cached_property
    def o2_cutout(self):
        """
        O2 signal after the breath-detection related cutting out of in-between breaths sections of the signal.
        :return: O2 signal after cutout
        """
        if self.config.cutout_option:
            return Signal(
                name='O2 (cutout)',
                unit=self.o2_crosstalk.unit,
                data=self.cutout(
                    input=self.o2_crosstalk.data,
                    breaths=self.breaths_raw
                )
            )
        else:
            return self.o2_crosstalk

    @cached_property
    def o2_crosstalk(self):
        """
        O2 signal after crosstalk correction. Crosstalk correction corrects the O2 concentration for the presence of
        CO2 in the gas mixture starting from Spiroware 3.3.1 and onwards. This crosstalk mode is enabled when the option
        crosstalk_option is 'quadratic' and disabled when it is 'one-way'.
        :return: crosstalk corrected O2 signal
        """
        if self.config.crosstalk_option == 'quadratic':
            return Signal(
                name='O2 (crosstalk corrected)',
                unit=self.o2_spedup.unit,
                data=self.o2_crosstalk_calculation(
                    config=self.config,
                    o2=self.o2_spedup.data,
                    co2=self.co2_filtered.data
                )
            )
        elif self.config.crosstalk_option == 'one-way':
            return self.o2_spedup

    @cached_property
    def o2_delayed(self):
        """
        The raw O2 signal is delay corrected by a somewhat complicated algorithm which is under certain circumstances
        defective, and certainly not maximally precise. However, this is a centrally important function for the correct
        interpretation of raw Spiroware signals. As the CO2 signal, the O2 signal and the O2 signal are measured at
        very different places in the experimental setup, it becomes necessary to align them in time.
        :return: delay-corrected O2 signal
        """
        return Signal(
            name='O2 (delay-corrected)',
            unit=self.o2_raw.unit,
            data=self.delay_correction(
                input=self.o2_raw.data,
                delay_times=self.delay_times['o2']
            )
        )

    @staticmethod
    def o2_drift_calculation(o2, co2, volume, breaths, config, pp, wi, wo, target_pp, target_wi, target_wo):
        """
        Function to calculate a "drift-corrected" (for each breath) O2 signal. It might be more accurate to describe it
        as a correction to reach the target values of O2 and CO2 for the prephase and washout respectively.
        :param o2: o2 Signal data
        :param co2: co2 Signal data
        :param volume: volume Signal data
        :param breaths: Breaths object
        :param config: Spiroware Config object containing additional info about sampling window, and allowed ranges
        :param wo: washout Breaths object
        :param wi: washin Breaths object
        :param pp: prephase Breaths object
        :param target_pp: prephase target value of O2
        :param target_wi: washin target value of O2
        :param target_wo: washout target value of O2
        :param samples: number of datapoints sampled per breath to calculate maximum
        :return: 
        """
        # Loop through all breaths, and calculate the drift correction factor
        drift = np.zeros(breaths.number)
        for i in breaths.numbers:
            # LungSim addition to prevent too small a range from being used for calculations which require at least
            # three data points. This is a consequence of the badly functioning breath detection algorithm,
            # and may need to be updated.
            if breaths.insp_range(i).size < 2:
                drift[i] = 1
                continue

            # Get the necessary signals of breath i (flow, o2, co2, time)
            o2_insp = o2[breaths.insp_range(i)]
            co2_insp = co2[breaths.insp_range(i)]
            volume_insp = volume[breaths.insp_range(i)] - volume[breaths.insp_range(i)[0]]
            # Calculate relative volume signal
            volume_insp_rel = volume_insp / volume_insp[-1]
            # Find the indices corresponding to the lower and upper limit
            index_min = (volume_insp_rel >= config.o2_drift_window[0]).nonzero()[0][0]
            index_max = (volume_insp_rel >= config.o2_drift_window[1]).nonzero()[0][0]
            # Cut the range to be between the two indices (strange indexing is intentional)
            o2_insp_range = o2_insp[index_min:index_max + 1]
            co2_insp_range = co2_insp[index_min:index_max + 1]
            # Calculate the sum of o2+co2 in this range
            sum_range = o2_insp_range + co2_insp_range
            # Sort ascending
            sum_range = np.sort(sum_range)
            # Mean of 3 highest data points is taken, if possible:
            if sum_range.size < config.drift_sample_number:
                drift[i] = 1
                continue

            peak = np.mean(sum_range[-config.drift_sample_number:])
            # Depending on which phase of the measurement the breath is in, a different correction factor is
            # calculated
            # During prephase:
            if i in pp:
                drift[i] = target_pp / peak
            # During washin:
            elif i in wi:
                drift[i] = target_wi / peak
            # During washout:
            elif i in wo:
                drift[i] = target_wo / peak

            # If the calculated correction factor for a breath is outside the boundaries, no correction is applied
            if drift[i] < config.drift_correction_window[0] or drift[i] > config.drift_correction_window[1]:
                drift[i] = 1

        # Create a corrected version of the input signal
        o2_corrected = np.copy(o2)
        for i in breaths.numbers:
            o2_corrected[breaths.range(i)] = o2_corrected[breaths.range(i)] * drift[i]

        return o2_corrected

    @cached_property
    def o2_endtidal(self):
        """
        End-tidal (from x% to y% of expired volume of a breath) sampling of the O2 gas signal.
        :return: Concentration of each breath in an end-tidal interval
        """
        return np.array(
            [
                self.end_tidal(
                    signal=self.o2.data[self.breaths.exp_range(i)],
                    flow=self.flow.data[self.breaths.exp_range(i)],
                    time=self.flow.time[self.breaths.exp_range(i)],
                    window=self.config.endtidal_window
                ) for i in self.breaths.numbers
            ]
        )

    @cached_property
    def o2_filtered(self):
        """
        The O2 signal is filtered with a lowpass filter to smooth it. This leads to signal delay and some loss of data.
        The Spiroware algorithm's default filters the signal with a harsh lowpass filter with a cutoff of 2Hz, which 
        leads to signal delay especially in the rise and fall moments.
        :return: Filtered O2 signal
        """
        if self.config.o2_filter_option:
            return Signal(
                name='O2 (filtered)',
                unit=self.o2_trimmed.unit,
                data=self.filter_lowpass_spiroware(
                    input=self.o2_trimmed.data,
                    sampling_rate=self.o2_trimmed.sampling_rate,
                    cutoff=self.config.o2_filter_cutoff,
                    window=self.config.o2_filter_window
                )
            )
        else:
            return self.o2_trimmed

    @cached_property
    def o2_raw(self):
        """
        The raw O2 signal, rescaled to be in fractions
        :return: raw O2 signal as a fraction
        """
        file_signal = self.file.signals['o2']
        output_unit = 'si'
        conversion = UnitConversion(input_unit=file_signal.unit, output_unit=output_unit)
        return Signal(
            name='O2 (raw, rescaled)',
            unit=output_unit,
            data=conversion.convert(input=file_signal.data)
        )

    @staticmethod
    def o2_maximum_calculation(o2, target, samples, window):
        """
        :param o2: o2 data before correction
        :param target: target value for maximum O2 value
        :param samples: number of samples collected to calculate maximum
        :param window: allowed correction window
        :return: corrected o2 data
        """
        o2_sorted = np.sort(o2)
        if o2_sorted.size < samples:
            return o2

        peak = np.mean(o2_sorted[-samples:])
        factor = target / peak

        # If the calculated correction factor for a breath is outside the boundaries, no correction is applied
        if factor < window[0] or factor > window[1]:
            return o2

        return o2 * factor

    @cached_property
    def o2_spedup(self):
        """
        The O2 signal has significant time-lag with respect to the CO2 signal, which can lead to signal artifacts if
        left uncorrected. A special filter is therefore applied to it which aims to "speed up" the parts of the
        signal where any change happens. This leads to noise in the signal, but reduces a bit the spikes caused by
        different time constants of the O2 and CO2 sensors. Speed up selected signals to compensate for the rise time
        delay of the sensors. v
        :return: Sped-up O2 signal
        """
        if self.config.o2_speed_option:
            return Signal(
                name='O2 (sped-up)',
                unit=self.o2_filtered.unit,
                data=self.speeding(
                    array=self.o2_filtered.data,
                    sampling_rate=self.o2_filtered.sampling_rate,
                    tau=self.config.o2_filter_tau,
                    butter_cutoff=self.config.speed_butter_cutoff
                )
            )
        else:
            return self.o2_filtered

    @cached_property
    def o2_trimmed(self):
        """
        After delay-correction, there are parts of the signal which no longer make any sense, or are simply empty. The
        delay correction produces a Signal which indicates which parts of the signals going into the delay correction
        are considered to be valid still. This function reduces the O2 signal down to that size.
        :return: O2 Signal trimmed to size after delay correction.
        """
        return Signal(
            name='O2 (trimmed)',
            unit=self.o2_delayed.unit,
            data=self.o2_delayed.data[self.valid_after_delay_correction.data]
        )

    @staticmethod
    def o2_crosstalk_calculation(config, o2, co2):
        """
        Oxygen crosstalk correction based on Wyler et al (2021). "Correction of sensor crosstalk error in Exhalyzer D
        multiple-breath washout device significantly impacts outcomes in children with cystic fibrosis"
        :param config: must be a Spiroware Config object containing the variables below
        :param o2: O2 signal data
        :param co2: CO2 signal data
        :return: Corrected O2 signal data
        """
        A = config.crosstalk_A
        B = config.crosstalk_B
        C = config.crosstalk_C
        D = config.crosstalk_D
        calib1 = config.o2_calibration_1
        calib2 = config.o2_calibration_2
        return o2 + A * (o2 - calib1) * (o2 - calib2) + B * co2 + C * np.square(co2) + D * co2 * o2

    @cached_property
    def physio_dead_space(self):
        """
        Estimate of the physiological dead space based on the CO2 signal. The volume at which the CO2 sensor reaches
        half-maximal signal is considered an estimate of the physiological dead space between airway and alveolar
        opening.
        """
        dead_space = np.zeros(self.breaths.number)
        for i in self.breaths.numbers:
            idx = self.breaths.exp_range(i)
            volume = self.volume.data[idx[0]] - self.volume.data[idx]
            co2 = self.co2.data[idx]
            threshold = 0.5 * np.max(co2)
            #threshold = 0.02
            idx_where = np.argwhere(co2 > threshold)

            #pp.plot(1e6*volume, 1e2*co2, 'k')
            #pp.axvline(volume[idx_where[0]]*1e6, color='b')
            #threshold = 0.02
            #idx_where = np.argwhere(co2 > threshold)
            #pp.axvline(volume[idx_where[0]]*1e6, color='r')

            if np.size(idx_where) == 0:
                dead_space[i] == nan
            else:
                dead_space[i] = volume[idx_where[0]]

        if np.all(np.isnan(dead_space)):
            raise LungSimException(output='All breaths failed physiological dead space calculation.')
        else:
            return np.nanmedian(dead_space) - self.config.dead_space_presensor

    @cached_property
    def quality_control(self):
        return SpirowareQualityControl(self)

    @staticmethod
    def round_half_up_array(n, decimals=0):
        """
        Round an array, where if the number falls on x.5 it will always round up. This replicates the behavior of
        whatever platform Spiroware runs on.
        :param n:
        :param decimals: Specification to which decimal the number n should be rounded.
        :return: rounded x, with x.5 rounded up always, independent of x
        """
        multiplier = 10 ** decimals
        return np.floor(n * multiplier + 0.5) / multiplier

    @classmethod
    def round_half_up_array_safe(cls, n, decimals=0):
        """
        This is a very strange work-around for a problem encountered when the delay values in the header are a multiple
        of the sampling time 0.005s, and also happen to be 12.5, 14.5, etc. In those numbers, Spiroware rounds them the
        way MATLAB rounds them, up to the nearest integer. Numpy by default rounds to the nearest even integer. The
        function below makes use of a trick to replicate MATLAB's behavior, with an additional fix for floating point
        numbers. If the delay 14.5 is stored as 14.9999999999, then this function first rounds it to 14.500.., and then
        performs the MATLAB rounding trick.
        :param n: number to be rounded
        :param decimals: number of decimals when rounding.
        :return:
        """
        multiplier = 10 ** decimals
        return np.floor(cls.round_half_up_array(n=n * multiplier + 0.5, decimals=decimals+5)) / multiplier

    @classmethod
    def threshold(cls, input, t):
        output = np.zeros_like(input)
        above = input > t
        output[above] = 1
        return output

    @cached_property
    def tracer_initial(self):
        """
        Initial tracer concentration. This is different for every washout/washin as well as for Spiroware specifically,
        because their breath selection in the pre-phase is a little bit bonkers.
        :return: A single value for the concentration of the tracer gas before the washout start.
        """
        if self.prephase_breaths.size == 0:
            msg = 'Initial tracer concentrations cannot be computed as no prephase breaths were detected.'
            raise LungSimException(output=msg)

        # Calculate the mean of the tidal volume of the prephase
        mean_volume_tidal = np.nanmean(self.volume_tidal[self.prephase_breaths])
        # Disqualify all pre-washout-breaths from the calculation of the initial tracer concentration whose expired
        # volumes deviate more than 25% from the tidal volume average. This is an attempt by Spiroware to select only
        # valid average breaths for the calculation of initial tracer concentration. However, it leads to quite a few
        # complications if the pre-breaths have irregular tidal volumes. Note that the comparison is between EXPIRED
        # VOLUMES and the mean of the TIDAL VOLUMES!
        breaths = np.copy(self.prephase_breaths)
        breaths = breaths[self.volume_expired[breaths] >= mean_volume_tidal - (0.25 * mean_volume_tidal)]
        breaths = breaths[self.volume_expired[breaths] <= mean_volume_tidal + (0.25 * mean_volume_tidal)]
        # Take the last three valid breaths. If there are fewer than three valid pre_breaths, take as many as there
        # are. If there are none, take the last breath all by itself.
        if breaths.size > 3:
            breaths = breaths[-3:]
        elif breaths.size == 0:
            breaths = self.prephase_breaths[-1]

        return np.mean(self.tracer_endtidal[breaths])

    @cached_property
    def tracer_inspired(self):
        """
        This property is governed by Spiroware options and is therefore specified for Spiroware Algorithms. If the
        integral option is chosen, the re-inspired volume of tracer gas is calculated the same as the expired volume.
        If the non-integral option is chosen, it is instead calculated as a estimate from the post-sensor dead space
        volume times the concentration of the last breath. This puts an intrinsic cap on how much tracer gas can be
        re-inspired, and presumably makes the calculations more exact in the context of infants.
        :return: Re-inspired tracer gas volume for
        each breath.
        """
        if self.config.tracer_reinspired_option == 'integral':
            return self.tracer_volume.data[self.breaths.insp_end] - self.tracer_volume.data[self.breaths.insp_start]
        elif self.config.tracer_reinspired_option == 'non-integral':
            return self.config.dead_space_postsensor * np.roll(self.tracer_endtidal, 1) * self.btps_factor_insp

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
    def tracer_normalized(self):
        """"
        :return: End-tidal concentration divided by initial concentration. Used for the determination of the washout end
        criterion.
        """
        return self.tracer_endtidal / self.tracer_initial

    @cached_property
    def tracer_slope_3(self):
        """
        The slope 3 of the tracer gas curve refers to the slope of a linear fit to the plateau of the expiratory tracer
        gas concentration. It has different fit windows for different Set sizes. It is specific to Spiroware exactly
        because the slope windows are set-dependent.
        :return: Slope 3 (absolute) of tracer gas curve during expiratory plateau.
        """
        slopes = np.zeros(self.breaths.number)

        for i in self.breaths.numbers:
            volume = self.volume.data[self.breaths.exp_start[i]] - self.volume.data[self.breaths.exp_range(i)]
            tracer = self.tracer.data[self.breaths.exp_range(i)]

            index_start = self.first_nonzero(volume > self.config.tracer_slope3[0] * np.max(volume))
            index_stop = self.first_nonzero(volume > self.config.tracer_slope3[1] * np.max(volume))

            # If start or stop indices cannot be found for a breath, make it not a number (NaN)
            if index_start.size == 0 or index_stop == 0:
                slopes[i] = nan

            # If start or stop indices are problematically close to each other, make result not a number (NaN)
            if index_stop - index_start < 1:
                slopes[i] = nan

            fit_range = np.arange(start=index_start, stop=index_stop+1)
            fit = np.polyfit(volume[fit_range], tracer[fit_range], 1)

            # If slope is negative, set it to 0
            slopes[i] = max(0, fit[0])

        return slopes

    @cached_property
    def tracer_slope_3_normalized(self):
        """
        The slope 3 of the tracer gas curve refers to the slope of a linear fit to the plateau of the expiratory tracer
        gas concentration. It has different fit windows for different Set sizes. It is specific to Spiroware exactly
        because the slope windows are set-dependent.
        :return: Slope 3 (absolute) of tracer gas curve during expiratory plateau.
        """
        slopes = np.zeros(self.breaths.number)

        for i in self.breaths.numbers:
            volume = self.volume.data[self.breaths.exp_start[i]] - self.volume.data[self.breaths.exp_range(i)]
            tracer = self.tracer.data[self.breaths.exp_range(i)]

            index_start = self.first_nonzero(volume > self.config.tracer_slope3[0] * np.max(volume))
            index_stop = self.first_nonzero(volume > self.config.tracer_slope3[1] * np.max(volume))

            # If start or stop indices cannot be found for a breath, make it not a number (NaN)
            if index_start.size == 0 or index_stop == 0:
                slopes[i] = nan

            # If start or stop indices are problematically close to each other, make result not a number (NaN)
            if index_stop - index_start < 1:
                slopes[i] = nan

            fit_range = np.arange(start=index_start, stop=index_stop+1)
            normalization = np.mean(tracer[fit_range])
            fit = np.polyfit(volume[fit_range], tracer[fit_range] / normalization, 1)

            # If slope is negative, set it to 0
            slopes[i] = max(0, fit[0])

        return slopes

    @cached_property
    def tracer_slope_3_normalized_VT(self):
        return self.tracer_slope_3_normalized * self.volume_expired

    @cached_property
    def sacin(self):
        """
        Acinar slope. Not to be confused with sacin * VT
        :return:
        """
        n = self.washout_breaths[0]
        above_min_vol = self.volume_expired[n] > 0.75 * np.mean(self.volume_tidal[self.spiroware_breaths])
        under_max_vol = self.volume_expired[n] < 1.25 * np.mean(self.volume_tidal[self.spiroware_breaths])

        if above_min_vol and under_max_vol:
            return max(0, self.tracer_slope_3_normalized[n] - self.lci_ao[n] * self.scond)
        else:
            raise LungSimException('First breath of washout automatically rejected. Sacin cannot be computed.')

    @cached_property
    def sacin_VT(self):
        """
        Acinar slope. Normalized by tidal volume
        :return:
        """
        n = self.washout_breaths[0]
        above_min_vol = self.volume_expired[n] > 0.75 * np.mean(self.volume_tidal[self.spiroware_breaths])
        under_max_vol = self.volume_expired[n] < 1.25 * np.mean(self.volume_tidal[self.spiroware_breaths])

        if above_min_vol and under_max_vol:
            return max(0, self.tracer_slope_3_normalized_VT[n] - self.lci_ao[n] * self.scond_VT)
        else:
            raise LungSimException('First breath of washout automatically rejected. Sacin cannot be computed.')

    @cached_property
    def scond(self):
        """
        :return: Slope of the conductive airways. Not to be confused with scond * VT. Not allowed to be less than 0.
        """
        if np.sum(self.scond_breaths) == 0:
            raise LungSimException('No breaths could be included for Scond analysis.')
        elif np.sum(self.scond_breaths) == 1:
            return 0
        else:
            scond = np.polyfit(
                x=self.lci_ao[self.scond_breaths],
                y=self.tracer_slope_3_normalized[self.scond_breaths],
                deg=1
            )
            return max(0, scond[0])

    @cached_property
    def scond_VT(self):
        """
        :return: Slope of the conductive airways. Not to be confused with scond * VT. Not allowed to be less than 0.
        """
        if np.sum(self.scond_breaths) == 0:
            raise LungSimException('No breaths could be included for Scond analysis.')
        elif np.sum(self.scond_breaths) == 1:
            return 0
        else:
            scond_VT = np.polyfit(
                x=self.lci_ao[self.scond_breaths],
                y=self.tracer_slope_3_normalized_VT[self.scond_breaths],
                deg=1
            )
            return max(0, scond_VT[0])

    @cached_property
    def scond_breaths(self):
        """
        These are the breaths included for the Scond analysis. They are the breaths of the washout which are within
        25% of the tidal volume, up to 2 breaths after 2.5% critical breath, and which are between 1.5 and 6 TO. Note
        that the +/-25% allowed deviation of expired volume from tidal volume is the same criterion as during the
        initial tracer gas calculation.
        :return: breath numbers of breaths included in Scond analysis
        """
        after_wo_start = self.breaths.numbers >= self.washout_breaths[0]
        before_extended_end = self.breaths.numbers <= self.spiroware_breaths[-1]
        min_turnover = self.lci_ao >= self.config.scond_window[0]
        max_turnover = self.lci_ao <= self.config.scond_window[1]
        min_volume = self.volume_expired >= 0.75 * np.mean(self.volume_tidal[self.spiroware_breaths])
        max_volume = self.volume_expired <= 1.25 * np.mean(self.volume_tidal[self.spiroware_breaths])

        return np.all([
                after_wo_start,
                before_extended_end,
                min_turnover,
                max_turnover,
                min_volume,
                max_volume
            ], 0)

    @staticmethod
    def speeding(array, sampling_rate, tau, butter_cutoff):
        """
        Implementation of Arieli speeding inspired by Arieli (1981) Corrections for
        the response time and delay of mass spectrometers". The purpose of this function is to compensate for the
        slow rise times of the side stream sensors. Air sampled at the point of interest has to fill small
        measurement chambers, which introduces a delay equivalent to a lowpass filter. By applying the inverse of a
        lowpass filter, we can in theory return the signal to a state closer to reality. This function implements a
        second order correction, as detailed in the paper. This function is essentially a noise-amplifying machine,
        so if the O2 or MM signals contain any noise above a certain frequency this function will pick that up.
        """
        # Calculate the butterworth filter inputs for signal smoothing
        b = 1
        a = np.array([1 / (2 * np.pi * butter_cutoff), 1])
        b_sm, a_sm = bilinear(b, a, sampling_rate)
        zi = lfilter_zi(b=b_sm, a=a_sm)
        # Apply smoothing filter
        array, _ = lfilter(b=b_sm, a=a_sm, x=array, zi=zi*array[0])

        # Calculate response filter inputs for signal speeding
        tau = tau + 1 / (butter_cutoff * 2 * np.pi)
        b = 1 / tau
        a = np.array([1, b])
        b, a = bilinear(b=b, a=a, fs=sampling_rate)
        # Invert the filter to turn a lowpass filter into a speeding filter
        b_sp = a
        a_sp = b
        # Calculate initial conditions
        zi = lfilter_zi(b=b_sp, a=a_sp)
        # Apply speeding filter
        array, _ = lfilter(b=b_sp, a=a_sp, x=array, zi=zi*array[0])

        return array

    @cached_property
    def spiroware_breaths(self):
        """
        This property serves to calculate a variety of Spiroware outcomes, as typically they are calculated over the
        breath table, and the breath table is cut off after 2 breaths after the washout end for 2.5% tracer. This
        allows easy range indications in the config file, where Spiroware outcome calculations can be easily
        replicated by writing: indices: spiroware_breaths
        :return: a range of breath numbers going from washout start to +2 critical breath at 2.5% tracer
        """
        # TODO: Make it so that if this calculation does not find the 2.5%, take the whole washout?
        start = self.washout_breaths[0]
        try:
            stop = self.washout_end(2.5)
        except LungSimException:
            stop = self.breaths.number

        return self.adjust_range(
            base_range=np.arange(start=start, stop=stop),
            prepend=0,
            append=3,
            boundaries=[0, self.breaths.number-1]
        )

    @cached_property
    def valid_after_delay_correction(self):
        """
        After the dynamic delay correction, we have to determine how much of the signals is still usable. As the MMss
        signal is the most affected, we check how much of the tail-end of the MMss signal has only zeros in it, and
        create a bool Signal object, which we can use to shorten all the other signals.
        :return: Signal containing a boolean (1=still valid, 0=invalid) used for shortening signals.
        """
        valid = np.zeros_like(self.mmside_delayed.data)
        valid[0:self.mmside_delayed.last_nonzero] = 1
        return valid.astype(bool)

    @cached_property
    def washin_breaths(self):
        """
        There is no washin portion in Spiroware Algorithms, but we define it here so we can use the parent
        function's way to return the array of phase names. While SF6 washouts (one class of daughter algorithms) have a
        washin portion, once the measurement is cut in two, it results again in a washout, and a pseudo-washout
        (washin portion flipped).
        :return: empty int array
        """
        return np.array([]).astype(int)

    @cache
    def washout_end(self, concentration):
        """
        Returns the critical breath with respect to concentration.
        :param concentration: concentration of tracer gas in [%]
        :return: breath number where the normalized exhaled concentration has fallen below this value.
        """
        if not (isinstance(concentration, int) or isinstance(concentration, float)):
            msg = 'Washout end method requires threshold. "{}" not supported.'.format(concentration)
            raise LungSimException(output=msg)

        # Convert the string of percentages into a fraction
        target = float(concentration) / 100
        # Find target breath, as the first of three breaths which are below a certain target concentration. We
        # construct a logical array of those end-tidal concentrations which lie below the target, circ-shift the array
        # twice and add the three arrays together. The first entry to be equal to 3 corresponds to the first entry
        # which is followed by two valid end-tidal concentrations.
        original = self.tracer_normalized[self.washout_breaths] < target
        shift_one = np.roll(original, -1)
        shift_one[-1] = 0
        shift_two = np.roll(original, -2)
        shift_two[-2:] = 0
        logical_array = original & shift_one & shift_two
        # If there are array indices that fulfil the above criterium of three consecutive breaths below the target,
        # choose the first one of these for the washout end
        if logical_array.any():
            washout_end = logical_array.nonzero()[0][0]
            washout_end = self.washout_breaths[washout_end]

        # If three consecutive breaths aren't found, find the last dip below the target and return it instead.
        else:
            cutoff = self.tracer_normalized[self.washout_breaths] - target
            crossing = cutoff * np.roll(cutoff, 1) <= 0
            negative = cutoff < 0
            # If no dip below the critical condition is found, return nothing. Otherwise return the last dip below.
            if negative.any():
                # washout_end = (crossing & negative).nonzero()[0][-1]
                # washout_end = self.washout_breaths[washout_end]
                msg = 'No three consecutive breaths below {}%'.format(concentration)
                raise LungSimException(output=msg)
            # Potentially risky way of noting that the washout end does not exist
            else:
                msg = 'The specified concentration of {}% is not reached'.format(concentration)
                raise LungSimException(output=msg)

        return washout_end
