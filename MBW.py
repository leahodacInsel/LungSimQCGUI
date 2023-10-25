from Breaths import Breaths
import numpy as np
from functools import cache
from functools import cached_property
from scipy.integrate import cumtrapz
from scipy.signal import find_peaks
from LungSimException import LungSimException
from Signal import Signal
from UnitConversion import UnitConversion


class MBW:
    """
    A (multiple-breath) washout algorithm is a set of functions and data that are typically encountered in any
    multiple-breath washout analysis. There exist a variety of methods and protocols, and a myriad of options of how to
    analyse any given multiple-breath washout experiment. Here are a few of those:
    - Ecomedics / Spiroware Nitrogen MBW
    - Ecomedics / Spiroware Sulphur Hexafluoride MBW
    - ndd / Wbreath Nitrogen MBW
    - ndd / Wbreath Sulphur Hexafluoride Mainstream Molar mass MBW
    - Custom Lungsim Sulphur Hexafluoride Mainstream Molar mass MBW

    These individual algorithms possess many functions which are only useful for their specific signal processing, and
    outcome computation. That is why their unique functions will be collected in classes which inherit from this one,
    while here we collect the functions that are useful for all algorithms to have access to.

    A washout will always contain at least a flow signal and a matching tracer gas signal. The algorithms will typically
    differ radically in how they get to the tracer gas signal, but it is useful to keep in mind that some steps are
    shared across several algorithms. Those should be made available to the daughter classes here. They include custom
    breath detection functions, some common algorithm steps such as Body-temperature/pressure-saturated (BTPS) flow
    correction, and computations of outcomes such as Lung Clearance Index (LCI) and functional residual capacity (FRC)
    according to consensus guidelines.

    As a general principle, whenever an algorithm should become available to multiple kinds of analysis, it may be
    helpful to include it here.

    Each MBW will be associated with a File, another LungSim class. For each run, additional options can be
    given to it by the MBWAnalysis class. These can include certain parameters forced to be changed across files.
    """

    def __init__(self, **kwargs):
        self.file = kwargs.pop('file')
        self.options = kwargs.pop('options', None)
        self.version = '2.04 (17.01.2022)'

    @classmethod
    def adjust_range(cls, base_range, prepend, append, boundaries):
        """
        Function to adjust base ranges by prepending or appending x consecutive breaths
        :param base_range: base range of breath numbers
        :param prepend: single number of how many breaths (without gaps) are added at the start of the range
        :param append: single number of how many breaths (without gaps) are added at the end of the range
        :param boundaries: boundaries of range
        :return: new range
        """
        new_range = np.copy(base_range)
        if prepend != 0:
            prepend_range = np.arange(start=-prepend, stop=base_range[0])
            new_range = np.insert(arr=new_range, obj=[0], values=prepend_range, axis=0)
        if append != 0:
            append_range = np.arange(start=base_range[-1], stop=base_range[-1] + append) + 1
            new_range = np.append(arr=new_range, values=append_range, axis=0)

        # Cap the range to the size of outcome.
        new_range = new_range[new_range <= boundaries[1]]
        new_range = new_range[new_range >= boundaries[0]]

        return new_range

    @cached_property
    def all_breaths(self):
        """
        Equivalent to the ranges used in the washout phases
        :return: Indices of all breaths detected
        """
        return self.breaths.numbers

    @classmethod
    def breath_detection_minvol(cls, time, flow, minvol):
        """
        Breath detection method based solely on flow which uses a properly implemented minimum volume criterion
        defined by the user. This function does not need to be a @cached_property as breaths itself is cached, and
        controls which function is called to compute the breath detection.
        :param time: time data numpy array
        :param flow: flow data numpy array
        :param minvol: minimum volume criterion for a breath to be detected
        :return: Breaths structure containing breath indices
        """
        # Find the maxima and minima of the volume signal, and sort them into one array
        volume = cumtrapz(y=flow, x=time, initial=0)
        maxima = find_peaks(volume)[0]
        minima = find_peaks(-volume)[0]
        zeros = np.sort(np.concatenate([maxima, minima]))

        # Find out whether the first breath is an inspiration or an expiration
        volumes = volume[zeros]
        volume_check = np.abs(volumes[1:] - volumes[0]) > minvol
        if not volume_check.any():
            msg = 'No breaths of sufficient size could be detected.'
            raise LungSimException(output=msg)

        first_direction = cls.first_nonzero(volume_check)
        if volumes[first_direction + 1] - volumes[0] > 0:
            current = 'insp'
        elif volumes[first_direction + 1] - volumes[0] < 0:
            current = 'exp'

        # First zero crossing
        accepted = zeros[0]

        # Iterate through the zero-crossings to find new valid zero crossings, starting from the 2nd one.
        for i in np.arange(start=1, stop=zeros.size):

            if current == 'insp':
                # If the current zero crossing is not part of the collection of signal peaks, skip it.
                if not zeros[i] in maxima:
                    continue

                # If the current zero crossing does not satisfy the volume criterion, skip it.
                if accepted.size == 1:
                    last = accepted
                else:
                    last = accepted[-1]

                if volumes[i] - volume[last] < minvol:
                    continue

                # If the above criteria are met, we still only accept this breath if the next half-breath criterion is
                # satisfied before an even higher volume peak is found.
                valid_minima = minima[minima > zeros[i]]
                # Relative volume to current index
                volume_diff = volumes[i] - volume[valid_minima]
                # First next breath below the necessary volume criterion.
                volume_check = volume_diff > minvol
                if not volume_check.any():
                    # If no such zero crossing exists, set value to maximum
                    first_below = np.max(zeros) + 1
                else:
                    first_below = valid_minima[cls.first_nonzero(volume_check)]

                # Find the next zero crossing with a higher volume
                valid_maxima = maxima[maxima > zeros[i]]
                volume_diff = volume[valid_maxima] - volumes[i]
                volume_check = volume_diff > 0
                # If no such zero crossing exists, set the value to an even higher maximum. Why? Below we compare the
                # first above with the first below value, and in the case where there is neither, we still want to
                # accept the last zero crossing as the end of the measurement. This automatically satisfies the same
                # question.
                if not volume_check.any():
                    first_above = np.max(zeros) + 2
                else:
                    first_above = valid_maxima[cls.first_nonzero(volume_check)]

                # Only if the next half - breath satisfies the criterion, before a zero - crossing of even higher volume
                # is found, do we accept the zero crossing and switch the search to the other half - breath.
                if first_below < first_above:
                    accepted = np.append(arr=accepted, values=zeros[i])
                    current = 'exp'

            elif current == 'exp':
                # 1:1 but flipped of the above
                if not zeros[i] in minima:
                    continue

                if accepted.size == 1:
                    last = accepted
                else:
                    last = accepted[-1]

                if volume[last] - volumes[i] < minvol:
                    continue

                valid_maxima = maxima[maxima > zeros[i]]
                volume_diff = volume[valid_maxima] - volumes[i]
                volume_check = volume_diff > minvol
                if not volume_check.any():
                    first_above = np.max(zeros) + 1
                else:
                    first_above = valid_maxima[cls.first_nonzero(volume_check)]

                valid_minima = minima[minima > zeros[i]]
                volume_diff = volumes[i] - volume[valid_minima]
                volume_check = volume_diff > 0
                if not volume_check.any():
                    first_below = np.max(zeros) + 2
                else:
                    first_below = valid_minima[cls.first_nonzero(volume_check)]

                if first_above < first_below:
                    accepted = np.append(arr=accepted, values=zeros[i])
                    current = 'insp'

        # Sort the breaths so they start with an inspiration and end with an expiration
        breath_volume = volume[accepted[1:]-1] - volume[accepted[:-1]]
        breath_volume = np.append(arr=breath_volume, values=volume[-1]-volume[accepted[-1]])

        first_last_check = breath_volume > 0
        if not first_last_check.any():
            msg = 'No breaths of sufficient size could be detected.'
            raise LungSimException(output=msg)

        first_inspiration = cls.first_nonzero(first_last_check)
        last_expiration = cls.last_nonzero(first_last_check)

        accepted = accepted[first_inspiration:last_expiration+1]
        # Because of the above operation, we now know that the sequence of zero-crossings begins with an inspiration
        # and ends with an expiration.
        inspirations = np.arange(start=0, stop=accepted.size, step=2)
        expirations = np.arange(start=1, stop=accepted.size - 1, step=2)
        # The inspiration starts coincide with the above indexes, but we ignore the last one because we want to end
        # with an expiration.
        insp_start = accepted[inspirations[:-1]]
        # The expiration starts coincide with the above expirations.
        exp_start = accepted[expirations]
        # The ends of inspirations are exactly one index before the starts of expirations.
        insp_end = accepted[expirations] - 1
        # And the expiration ends are one index before the starts of the inspirations. We ignore the first one
        # because we want to start with an inspiration.
        exp_end = accepted[inspirations[1:]] - 1

        return Breaths(
            insp_start=insp_start,
            insp_end=insp_end,
            exp_start=exp_start,
            exp_end=exp_end
        )

    @cached_property
    def breath_duration(self):
        """
        :return: Duration of each breath's inspiration + expiration
        """
        return self.breath_duration_insp + self.breath_duration_exp

    @cached_property
    def breath_duration_insp(self):
        """
        :return: Duration of each breath's inspiration.
        """
        return self.flow.time[self.breaths.insp_end] - self.flow.time[self.breaths.insp_start]

    @cached_property
    def breath_duration_exp(self):
        """
        :return: Duration of each breath's exppiration
        """
        return self.flow.time[self.breaths.exp_end] - self.flow.time[self.breaths.exp_start]

    def calculate_output(self, args):
        """
        Based on a MBW-like object and specs in the result_config file, return output.
        Specs include the following keywords:
        variable: Name of cached_property in analysis. This can be an array of results or a single number / string.
        method: 'critical' for specific turnover moments, 'mean', 'sd', 'cv', 'count' etc.
        indices: 'all', 'washout_breaths' etc.
        prepend/append: how many breaths are added before and after the range
        rounding:
        """
        # Variable handling. If no variable is specified in the output, raise exception.
        if 'variable' not in args:
            raise LungSimException(output='Specify variable name in output config.')
        elif not isinstance(['variable'], str):
            variable = str(args['variable'])
        else:
            variable = args['variable']

        # If specified variable does not exist, raise exception
        try:
            variable = eval('self.' + variable)
        except AttributeError:
            msg = 'Specified variable "{}" does not exist in this algorithm.'.format(variable)
            raise LungSimException(msg)

        # Method
        if 'method' not in args:
            method = 'all'
        elif not isinstance(['method'], str):
            method = str(args['method'])
        else:
            method = args['method']

        # Indices / ranges
        if method != 'single' and method != 'critical':
            if 'indices' not in args:
                indices = 'all_breaths'
            elif not isinstance(args['indices'], str):
                indices = str(args['indices'])
            else:
                indices = args['indices']

            try:
                indices = eval('self.' + indices)
            except AttributeError:
                msg = 'Specified indices range "{}" does not exist in this algorithm.'.format(indices)
                raise LungSimException(msg)

        if method == 'single':
            indices = 'all'

        # Concentration
        if method == 'critical':
            if 'concentration' in args:
                try:
                    concentration = float(args['concentration'])
                    indices = self.washout_end(concentration=concentration)
                except ValueError:
                    msg = 'Specified critical concentration "{}" formatted incorrectly'.format(args['concentration'])
                    raise LungSimException(output=msg)
            else:
                raise LungSimException(output='Please specify concentration for "critical" method')

        # Prepend / append
        if 'prepend' in args:
            try:
                prepend = int(args['prepend'])
            except ValueError:
                msg = 'Specified prepend indices "{}" not formatted correctly. Enter number'.format(args['prepend'])
                raise LungSimException(msg)
        else:
            prepend = 0

        if 'append' in args:
            try:
                append = int(args['append'])
            except ValueError:
                msg = 'Specified append indices "{}" not formatted correctly. Enter number'.format(args['append'])
                raise LungSimException(msg)
        else:
            append = 0

        # Adjusting index range. Only necessary and possible if method is not 'single'
        if method != 'single' and method != 'critical':
            indices = self.adjust_range(
                base_range=indices,
                prepend=prepend,
                append=append,
                boundaries=[0, variable.size - 1]
            )

        # Unit conversion from program-internal SI units to the desired output unit, if supported.
        if 'unit' in args:
            conversion = UnitConversion(input_unit='si', output_unit=args['unit'])
            variable = conversion.convert(variable)

        # Apply method to specified range of variable
        return self.method_manager(
            outcome=variable,
            indices=indices,
            method=method
        )

    @cached_property
    def cev_bp(self):
        """
        The cumulative expired volume (CEV) at the bypass flow (bp) represents the total volume expired from
        washout start to the breath where expired concentration falls below 2.5% of the initial concentration. No
        further adjustments to volume have to be made, as all air measured here contributes to washing out (i.e. there
        is no machine dead space beyond the bp point.
        :return: cumulative expired volume at the bypass (unadjusted CEV).
        """
        cev_bp = np.zeros_like(self.volume_expired)
        cev_bp[self.washout_breaths] = np.cumsum(self.volume_expired[self.washout_breaths])
        return cev_bp

    @cached_property
    def cev_ao(self):
        """
        The cumulative expired volume (CEV) at airway opening (ao) represents the total volume expired from washout
        start to the breath where expired concentration falls below 2.5% of the initial concentration. It is
        additionally adjusted for the volume of dead space which "does not contribute to the washout" for each breath.
        This total machine dead space volume is subtracted from each expiration of the washout breaths, regardless of
        their volume.
        :return: Cumulative Expired Volume at airway opening for each breath.
        """
        cev_ao = np.zeros_like(self.volume_expired)
        cev_ao[self.washout_breaths] = np.cumsum(self.volume_expired[self.washout_breaths] - self.machine_dead_space)
        return cev_ao

    @staticmethod
    def contain_within(values, bounds):
        """
        Function to ensure that all numbers in input array are contained within the bounds established by bounds
        :param bounds: a pair of two numbers, which limit the values that are allowed in input
        :param values: input numbers
        :return: input, but with all numbers greater than bounds[1] reduced to bounds[1] and same for the other side.
        """
        values[values < bounds[0]] = bounds[0]
        values[values > bounds[1]] = bounds[1]

        return values

    @staticmethod
    def first_nonzero(array):
        """
        returns the index of the first nonzero element of a numpy array
        """
        return (array != 0).nonzero()[0][0]

    @cached_property
    def frc_ao(self):
        """
        The functional residual capacity at airway opening (FRCao) is another main outcome of LCI measurements. It is
        calculated by measuring the expired volume of tracer gas (N2, SF6, or other), and assessing how much of the
        initial volume of FRC was occupied by said gas. FRC corresponds here to the volume after expiration of the last
        breath before the start of the washout procedure. Airway opening specifies that we subtract all the volume
        between the sensor point and the airway opening from the final result.
        :return: functional residual capacity at airway opening (FRCao)
        """
        frc_ao = self.frc_gs
        frc_ao[self.washout_breaths] -= self.config.dead_space_presensor
        return frc_ao

    @cached_property
    def frc_gs(self):
        """
        Functional residual capacity at the gas-sampling / sensor point.
        """
        frc_gs = np.zeros_like(self.volume_expired)
        wo = self.washout_breaths
        frc_gs[wo] = self.tracer_cev[wo] / (self.tracer_initial - self.tracer_endtidal[wo])
        return frc_gs

    @staticmethod
    def last_nonzero(array):
        """
        returns the index of the last nonzero element of a numpy array
        """
        return (array != 0).nonzero()[0][-1]

    @cached_property
    def lci_ao(self):
        """
        The Lung Clearance Index is the principal outcome of a MBW measurement. It can be calculated in a variety of
        ways: At the gas sampling point, at the alveolar opening, or like here: at the airway opening (ao). This
        involves subtracting from the denominator all the dead space between the airway opening and the sensors, and
        subtracting from the cumulative expired volume for each breath the total machine dead space.
        :return: Lung Clearance Index in Turnovers at airway opening for each breath[TO]
        """
        return self.outcome_divide(self.cev_ao, self.frc_ao)

    @classmethod
    def method_manager(cls, outcome, indices='all', method='all'):
        """
        Function central for the calculation of single value outcomes from an array of breath-wise outcomes. Expired
        volume is a property with an outcome for every breath, but in order to output results into a single outcome
        table, there are many ways to extract information from it. Mean, SD, specific single values related to TO, etc.
        :param outcome: Base array of outcomes, typically either
        :param indices: indices across which the method is applied
        :param method: method applied to the outcome (mean, SD, min, max, CV, all, etc.)
        :return: method applied to the outcomes in the indices specified.
        """
        # Return outcomes based on the method chosen
        if method == 'mean':
            return np.mean(outcome[indices])
        if method == 'median':
            return np.median(outcome[indices])
        elif method == 'sd':
            return np.std(outcome[indices])
        elif method == 'cv':
            return np.std(outcome[indices]) / np.mean(outcome[indices])
        elif method == 'max':
            return np.max(outcome[indices])
        elif method == 'min':
            return np.min(outcome[indices])
        elif method == 'count':
            return outcome[indices].size
        elif method == 'all' or method == 'critical':
            return outcome[indices]
        elif method == 'single':
            return outcome
        else:
            raise LungSimException(output='Specified method "{}" not supported.'.format(method))

    @staticmethod
    def spikyness(signal):
        """
        Function to quantify how "spiky" a signal is, i.e. how many ups and downs are present.
        :param signal: signal to be quantified for spikyness
        :return: single number relatively quantifying spikyness
        """
        return np.sum(np.abs(np.diff(signal)))

    @staticmethod
    def outcome_divide(num, den):
        """
        Array division which only divides the entries where the denominator is not zero, and leaves those at 0
        :param num: Numerator array
        :param den: Denominator array
        :return: Division array
        """
        result = np.zeros_like(num)
        if (den != 0).size == 0:
            return result
        else:
            result[den != 0] = num[den != 0] / den[den != 0]
            return result

    @cached_property
    def peak_flow_insp(self):
        """
        :return: Peak value of the flow signal within a breath.
        """
        return np.array([np.max(self.flow.data[self.breaths.insp_range(i)]) for i in self.breaths.numbers])

    @cached_property
    def peak_flow_exp(self):
        """
        :return: Peak value of the flow signal within a breath.
        """
        return np.array([np.max(self.flow.data[self.breaths.exp_range(i)]) for i in self.breaths.numbers])

    @cached_property
    def phases(self):
        """
        Array containing strings with the names of the phases each breath belongs to.
        :return: array containing 'prephase', 'washin', or 'washout' for each breath.
        """
        output = np.full(self.breaths.number, '').astype('object')
        output[self.prephase_breaths] = 'pre-phase'
        output[self.washout_breaths] = 'test-phase'
        return output

    def reset(self):
        """
        This function is centrally important (for now), even though it shouldn't exist. Somehow, when computing lots of
        MBW measurements, a lot of garbage accumulates progressively in the RAM. For some reason, the MBW objects are
        not properly destroyed, and it's unclear to me exactly why. The "solution" to this, for now, is after completing
        all computations relating to one analysis, it is reset, and its properties and attributes are emptied.
        This function is called in the loop going over all measurements.
        """
        dictkeys = list(self.__dict__.keys())
        for key in dictkeys:
            self.__dict__.pop(key)

    @staticmethod
    def signal_correlate(x, y, sweep):
        """
        Custom signal correlation function.
        :param x: Signal 1
        :param y: Signal 2
        :param sweep: Local distance over which the correlation is performed
        :return: minimum value within the sweep, maximum value, and the all the values
        """
        k = 0
        lags = np.arange(start=-sweep, stop=sweep + 1)
        value = np.zeros(np.shape(lags))
        for i in lags:
            value[k] = np.dot(x[sweep:-sweep], np.roll(y, i)[sweep:-sweep])
            k = k + 1

        return np.argmin(value) - sweep, np.argmax(value) - sweep, value

    @cached_property
    def time(self):
        """
        The time signal is used very sparingly, as it is a late addition. It exists solely to simplify the outputs in
        B-files. Usually the time "signal" is a feature of the signal class, as it accompanies another signal, as
        opposed to being a real signal of its own.
        """
        return Signal(
            name='Time',
            unit='s',
            data=np.arange(start=0, stop=self.flow.size)/200
        )

    @cached_property
    def tlc_ao(self):
        """
        The tidal lung capacity is the volume defined as the volume in the lung after tidal inspiration. Here at airway
        opening (ao).
        """
        tlc_ao = self.tlc_gs
        tlc_ao[self.washout_breaths] -= self.config.dead_space_presensor
        return tlc_ao

    @cached_property
    def tlc_gs(self):
        """
        The tidal lung capacity is the volume defined as the volume in the lung after tidal inspiration. Here at gas
        sampling (gs) point.
        """
        tlc_gs = np.zeros_like(self.volume_expired)
        wo = self.washout_breaths
        tlc_gs[wo] = self.tracer_cev_plus[wo] / (self.tracer_initial - self.tracer_endtidal[wo])
        return tlc_gs

    @cached_property
    def tracer_cev(self):
        """
        Cumulative expired volume of tracer gas starting from the washout start.
        :return: Cumulative expired volume (so far) for each breath, starting from washout start.
        """
        tracer_cev = np.zeros_like(self.tracer_expired_net)
        tracer_cev[self.washout_breaths] = np.cumsum(self.tracer_expired_net[self.washout_breaths])
        return tracer_cev

    @cached_property
    def tracer_cev_plus(self):
        """
        Cumulative expired volume of tracer gas starting from the half-breath before washout start. This will be used
        to calculate the tidal lung capacity rather than functional residual capacity.
        """
        tracer_cev = np.zeros_like(self.tracer_expired_net)
        tracer_cev[self.washout_breaths] = np.cumsum(self.tracer_expired_net[self.washout_breaths])
        "Critical addition: We also count the expiration of the last breath before washout start"
        tracer_cev[self.washout_breaths] += self.tracer_expired[self.washout_breaths[0]-1]
        return tracer_cev

    @cached_property
    def tracer_expired_net(self):
        """
        :return: Volume of tracer gas expired in each breath minus inspired volume
        """
        return self.tracer_expired - self.tracer_inspired

    @cached_property
    def tracer_expired(self):
        """
        :return: Expired volume of tracer gas for each breath separately.
        """
        return self.tracer_volume.data[self.breaths.exp_start] - self.tracer_volume.data[self.breaths.exp_end]

    @cached_property
    def tracer_volume(self):
        """
        :return: Volume signal for the tracer gas. Represents tracer gas moving at any given moment in and out of the
        patient, in m^3
        """
        return Signal(
            name='Tracer Volume (instantaneous)',
            unit=self.flow.unit,
            data=cumtrapz(self.flow.data * self.tracer.data, self.flow.time, initial=0)
        )

    @staticmethod
    def vapor_pressure(humidity, temp):
        """
        Calculates the water vapor pressure at a given humidity and temperature according to the Magnus formula.
        :param humidity: Humidity as a fraction
        :param temp: Temperature in degrees C
        :return: water vapor pressure in Pa
        """
        return 611.2 * humidity * np.exp(17.62 * temp / (temp + 243.12))

    @cached_property
    def volume(self):
        """
        The volume property is a Signal object, a simple integration of the corrected flow signal
        :return: The volume Signal based on recorded flow.
        """
        return Signal(
            name='Volume (continuous)',
            unit='m^3',
            data=cumtrapz(y=self.flow.data, x=self.flow.time, initial=0)
        )

    @cached_property
    def volume_abs(self):
        """
        The volume_abs property is a Signal object, which creates a continuous x-axis that mirrors volume
        :return:
        """
        return Signal(
            name='Volume (absolute)',
            unit='m^3',
            data=cumtrapz(y=np.abs(self.flow.data), x=self.flow.time, initial=0)
        )

    @cached_property
    def volume_inspired(self):
        """
        The inspired volume refers to the total inspired volume per breath
        :return: Volume inspired for each breath
        """
        return self.volume.data[self.breaths.insp_end] - self.volume.data[self.breaths.insp_start]

    @cached_property
    def volume_expired(self):
        """
        The expired volume refers to the total expired volume per breath.
        :return: Volume expired for each breath
        """
        return self.volume.data[self.breaths.exp_start] - self.volume.data[self.breaths.exp_end]

    @cached_property
    def volume_tidal(self):
        """
        Mean of the expired and inspired volume per breath. Not to be confused with expired volume, or volume normalized
        by body weight etc.
        :return: Mean of inspired and expired volume for each breath
        """
        return (self.volume_expired + self.volume_inspired) / 2

    @cache
    def washout_breaths_to_critical(self, concentration):
        """
        Range of breaths from washout start to any given washout end including only the first critical breath.
        :param concentration: concentration in %
        :return: range of breath numbers
        """
        start = self.washout_breaths[0]
        stop = self.washout_end(concentration)

        return np.arange(start=start, stop=stop + 1)
