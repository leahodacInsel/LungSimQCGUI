from functools import cached_property
from math import nan
from math import isnan
from numpy import log
import numpy as np
from scipy.interpolate import interp1d


class SpirowareQualityControl:

    def __init__(self, spiroware_mbw):
        self.mbw = spiroware_mbw

    @cached_property
    def breath_duration_insp(self):
        return self.self_similarity(self.mbw.breath_duration_insp)

    @cached_property
    def breath_duration_exp(self):
        return self.self_similarity(self.mbw.breath_duration_exp)

    @cached_property
    def concentration_initial(self):
        """
        :return: returns a 0 if the average initial concentration is below 77% for N2 washouts, and below 3.9% for SF6
        washouts. This indicates that there was not sufficient time in between measurements for the lung to wash in
        again for N2 measurements, or not a sufficiently long washin time for SF6 measurements. returns a 1 if this
        condition is fulfilled.
        """
        if self.mbw.method_name == 'Spiroware N2 Washout':
            cutoff = 0.77
        elif self.mbw.method_name == 'Spiroware SF6 Washin' or self.mbw.method_name == 'Spiroware SF6 Washout':
            cutoff = 0.039

        return int(self.mbw.tracer_initial > cutoff)

    @cached_property
    def efficiency(self):
        """
        :return: returns the fraction of volume efficiency of breathing.
        """
        # Calculate FRCao, at 2.5% threshold. This is the standard MBW outcome."
        args = dict()
        args['variable'] = 'frc_ao'
        args['method'] = 'critical'
        args['concentration'] = 2.5
        frc_raw = self.mbw.calculate_output(args)

        # Calculate the "alveolar" part of the FRC volume
        frc_alv = frc_raw - self.mbw.physio_dead_space
        # Dead space originating purely from the device
        ds_machine = self.mbw.machine_dead_space
        # Total dead space of the setup (machine + physiological)
        ds_tot = self.mbw.machine_dead_space + self.mbw.physio_dead_space

        # Given the real alveolar FRC and dead spaces, calculate the expected number of
        # washout breaths required for a perfectly homogenous lung.
        vt_all = np.linspace(start=0, stop=10 * frc_alv, num=10000)
        k = (vt_all - ds_tot) / (frc_alv + vt_all)
        bp = log(0.025) / log(1 - k)

        # Calculate how much of the Vt of each breath would be counted towards the LCI
        vt_counted = (vt_all - ds_machine) / frc_raw
        vt_counted[vt_counted <= 0] = nan

        # Calculate the expected LCI for each Vt
        lci_est = bp * vt_counted
        lci_est[lci_est <= 0] = nan

        # Calculate the efficiency as "how much more LCI would I expect based purely on Vt, as a function of DS and FRC.
        # efficiency = 1 - (lci_est - np.nanmin(lci_est)) / np.nanmin(lci_est)
        efficiency = 1 - np.abs(1 - lci_est / np.nanmin(lci_est))
        efficiency[efficiency <= 0] = 0

        return vt_all, efficiency

    @cached_property
    def efficiency_max(self):
        vt_all, efficiency = self.efficiency
        return vt_all[np.nanargmax(efficiency)]

    @cached_property
    def efficiency_score(self):
        # Calculate vt as the mean tidal volume (mean of insp and exp) over washout to the test-end criterion."
        vt_all, efficiency = self.efficiency

        args = dict()
        args['variable'] = 'volume_tidal'
        args['method'] = 'mean'
        args['indices'] = 'washout_breaths_to_critical(2.5)'
        vt = self.mbw.calculate_output(args)

        # Evaluate the real vt in the function efficiency(vt_all)
        interpolation = interp1d(x=vt_all, y=efficiency)
        return interpolation(vt)

    @cached_property
    def efficiency_score2(self):
        """
        :return: returns the fraction of volume efficiency of breathing.
        """
        vt_all, efficiency = self.efficiency

        # Calculate vt as the mean tidal volume (mean of insp and exp) over washout to the test-end criterion."
        args = dict()
        args['variable'] = 'volume_tidal'
        args['method'] = 'mean'
        args['indices'] = 'washout_breaths_to_critical(2.5)'
        vt = self.mbw.calculate_output(args)

        # Evaluate the real vt in the function efficiency(vt_all)
        ideal = vt_all[np.nanargmax(efficiency)]
        score = 1 - np.abs((vt-ideal)/ideal)
        if score < 0:
            return 0
        else:
            return score

    @cached_property
    def end_of_test_reached(self):
        """
        :return: returns 0 if there are not three consecutive breaths with less than the target concentration of 2.5%
        within the measurement. Returns 1 if this is the case.
        """
        # Target of an MBW is still basically always 2.5%.
        target = 0.025
        # Find target breath, as the first of three breaths which are below a certain target concentration. We
        # construct a logical array of those end-tidal concentrations which lie below the target, circ-shift the array
        # twice and add the three arrays together. The first entry to be equal to 3 corresponds to the first entry
        # which is followed by two valid end-tidal concentrations.
        original = self.mbw.tracer_normalized[self.mbw.washout_breaths] < target
        shift_one = np.roll(original, -1)
        shift_one[-1] = 0
        shift_two = np.roll(original, -2)
        shift_two[-2:] = 0
        logical_array = original & shift_one & shift_two
        # If there are array indices that fulfil the above criterium of three consecutive breaths below the target,
        # choose the first one of these for the washout end
        if logical_array.any():
            return 1
        else:
            return 0

    @staticmethod
    def first_nonzero(array):
        """
        returns the index of the first nonzero element of a numpy array
        """
        nz = (array != 0).nonzero()
        if np.size(nz) == 0:
            return np.size(array)
        else:
            return (array != 0).nonzero()[0][0]

    @cached_property
    def flow_insp_bottom(self):
        return self.self_similarity(self.flow_sampling['bottom_insp'])

    @cached_property
    def flow_insp_mean(self):
        return self.self_similarity(self.flow_sampling['mean_insp'])

    @cached_property
    def flow_insp_peak(self):
        return self.self_similarity(self.flow_sampling['peak_insp'])

    @cached_property
    def flow_exp_bottom(self):
        return self.self_similarity(self.flow_sampling['bottom_exp'])

    @cached_property
    def flow_exp_mean(self):
        return self.self_similarity(self.flow_sampling['mean_exp'])

    @cached_property
    def flow_exp_peak(self):
        return self.self_similarity(self.flow_sampling['peak_exp'])

    @cached_property
    def flow_sampling(self):
        """
        Function to return several samples of flow from each breath:
        top 20% of flow speeds
        bottom 20% of flow speeds
        mean flow speeds
        for insp and exp
        :return:
        """
        sampling = dict()
        sampling['bottom_insp'] = np.zeros(self.mbw.breaths.number)
        sampling['peak_insp'] = np.zeros(self.mbw.breaths.number)
        sampling['mean_insp'] = np.zeros(self.mbw.breaths.number)
        sampling['bottom_exp'] = np.zeros(self.mbw.breaths.number)
        sampling['peak_exp'] = np.zeros(self.mbw.breaths.number)
        sampling['mean_exp'] = np.zeros(self.mbw.breaths.number)

        fraction = 0.2
        for i in self.mbw.breaths.numbers:
            # Sampling within inspiration
            indices = self.mbw.breaths.insp_range(i)
            sorted_flow = np.sort(self.mbw.flow.data[indices])
            sample_points = int(np.ceil(sorted_flow.size * fraction))
            sampling['bottom_insp'][i] = np.mean(sorted_flow[:sample_points])
            sampling['peak_insp'][i] = np.mean(sorted_flow[-sample_points:])
            sampling['mean_insp'][i] = np.mean(sorted_flow)

            # Sampling within expiration
            indices = self.mbw.breaths.exp_range(i)
            sorted_flow = np.sort(self.mbw.flow.data[indices])
            sample_points = int(np.ceil(sorted_flow.size * fraction))
            sampling['bottom_exp'][i] = np.mean(sorted_flow[:sample_points])
            sampling['peak_exp'][i] = np.mean(sorted_flow[-sample_points:])
            sampling['mean_exp'][i] = np.mean(sorted_flow)

        return sampling

    @cached_property
    def co2_drift(self):
        """
        :return: Returns a measure of the drift of the CO2 signal over the course of the washout
        """
        breaths = self.mbw.washout_breaths_to_critical(2.5)
        volume_median = np.median(self.mbw.volume_expired[breaths])
        co2_et = np.zeros(np.size(breaths))

        k = 0

        for i in breaths:

            idx = self.mbw.breaths.exp_range(i)
            volume_breath = self.mbw.volume.data[idx[0]] - self.mbw.volume.data[idx]

            if np.max(volume_breath) < volume_median:
                co2_et[k] = -1
                k += 1
                continue

            co2_breath = self.mbw.co2.data[idx]
            interpolation = interp1d(volume_breath, co2_breath)
            co2_et[k] = interpolation(volume_median)
            k += 1

        valid = co2_et > 0
        x = (breaths[valid] - np.min(breaths)) / (np.max(breaths)-np.min(breaths))
        y = co2_et[valid]

        fitting = np.polyfit(x, y, 1)
        start = np.polyval(fitting, 0)
        stop = np.polyval(fitting, 1)
        if (start > 0.06 or start < 0.04) or (stop > 0.06 or stop < 0.04):
            # return 0
            return 1 - np.abs(fitting[0]) / 0.02
            #return (fitting[0]) / 0.01
        else:
            return 1 - np.abs(fitting[0]) / 0.02
            #return (fitting[0]) / 0.01

    @staticmethod
    def self_similarity(input):
        """
        This function returns an array of the same shape as input, where each entry corresponds to how far each array
        entry is distant from the median of the array. The output is normalized so that inputs at 0 receive an output of
        0, and inputs at 2x the median receive 0 also.
        This implies that inputs should all either be strictly positive or negative.
        :param input: np array (most likely corresponding to breath-by-breath results.
        :return: np array of the normalized distance to the array median
        """
        median = np.median(input)
        distance = np.abs((input - median) / median)  # Alternatively: (input-median)**2
        score = 1 - distance
        score[score < 0] = 0
        return score

    @cached_property
    def score(self):
        return np.mean(a=(self.breath_duration_insp, self.breath_duration_exp, self.flow_exp_mean,
                          self.flow_exp_peak, self.flow_insp_mean, self.flow_insp_peak,
                          self.volume_inspired, self.volume_expired),
                       axis=0)

    @cached_property
    def score_end_of_test(self):
        """
        :return: Returns the average self-similarity score (mixed insp and exp) for the first three breaths of the
        washout phase below the target concentration.
        """
        breaths = np.arange(start=self.mbw.washout_end(2.5), stop=self.mbw.washout_end(2.5)+3)
        return np.mean(self.score[breaths])

    @cached_property
    def score_first_washout(self):
        """
        :return: Returns the average self-similarity score (mixed insp and exp) for the first three breaths of the
        washout phase.
        """
        breaths = self.mbw.washout_breaths[:3]
        return np.mean(self.score[breaths])

    @cached_property
    def score_last_prephase(self):
        """
        :return: Returns the average self-similarity score (mixed insp and exp) for the last three breaths of the
        prephase.
        """
        breaths = self.mbw.prephase_breaths[-3:]
        return np.mean(self.score[breaths])

    @cached_property
    def score_rest_of_washout(self):
        """
        :return: Returns the average self-similarity score (mixed insp and exp) for the breaths of the
        washout phase not captured by any other criteria.
        """
        breaths = self.mbw.washout_breaths[3:self.mbw.washout_end(2.5)]
        return np.mean(self.score[breaths])

    @cached_property
    def volume_inspired(self):
        return self.self_similarity(self.mbw.volume_inspired)

    @cached_property
    def volume_expired(self):
        return self.self_similarity(self.mbw.volume_expired)

    @cached_property
    def log_slope(self):

        slopes = np.zeros(self.mbw.breaths.number)

        for i in self.mbw.breaths.numbers:

            indices = self.mbw.breaths.exp_range(i)
            volume = self.mbw.volume.data[indices[0]] - self.mbw.volume.data[indices]
            relevant = volume > (self.mbw.physio_dead_space + self.mbw.machine_dead_space)

            if np.sum(relevant) < 20:
                slopes[i] = -1
                continue

            co2 = self.mbw.co2.data[indices[relevant]]
            tracer = self.mbw.tracer.data[indices[relevant]]/self.mbw.tracer_initial

            fit = np.polyfit(co2, tracer, 1)

            if fit[0] > 20 or fit[0] < 0:
                slopes[i] = -1
                continue

            slopes[i] = np.polyval(fit, 0.05)
            # slopes[i] = fit[0]

        keep = slopes > 0
        slopes[slopes < 0] = 1000
        slopes = log(slopes)
        slopes = slopes[self.mbw.washout_breaths_to_critical(2.5)]
        keep = keep[self.mbw.washout_breaths_to_critical(2.5)]

        slopes = slopes[keep]

        return slopes

    @cached_property
    def uptick(self):
        """
        :return: returns a measure of how much uptick in end-tidal concentrations is detected.
        """
        slopes = self.log_slope
        number = np.size(slopes)
        uptick = np.zeros(number)

        for i in np.arange(start=0, stop=number):

            current = slopes[i]
            behind = slopes[:i]
            ahead = slopes[i+1:]

            correct = np.sum(ahead < current) + np.sum(behind > current)
            incorrect = np.sum(ahead > current) + np.sum(behind < current)

            uptick[i] = correct / (correct + incorrect)

        return np.mean(uptick)

    @cached_property
    def uptick2(self):
        """
        :return: returns a measure of how much uptick in end-tidal concentrations is detected.
        """
        slopes = np.exp(self.log_slope)
        differences = np.diff(slopes)
        score = 1 - np.sum(differences[differences > 0]) / (np.max(slopes)-np.min(slopes))
        if score < 0:
            return 0
        else:
            return score

    @cached_property
    def uptick3(self):
        """
        :return: returns a measure of how much uptick in end-tidal concentrations is detected.
        """
        slopes = self.log_slope
        differences = np.diff(slopes)
        score = 1 - np.sum(differences[differences > 0]) / (np.max(slopes)-np.min(slopes))
        if score < 0:
            return 0
        else:
            return score
