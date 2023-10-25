from Config import Config
from LungSimException import LungSimException


class SpirowareConfig(Config):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # TODO: Remove unused parameters

        # ATP and BTPS parameters
        self.sampleflow_option = True
        self.atp_option = True
        self.btps_option = True
        self.pressure_ambient = 96000
        self.humidity_inspired = 0
        self.humidity_flowhead = 0.6
        self.humidity_body = 1
        self.temp_inspired = 26.0
        self.temp_flowhead = 30.0
        self.temp_body = 37.0
        self.temp_kelvin_zero = - 273.15

        # Delay parameters
        self.delay_option = 'dynamic'
        self.delay_co2 = 0.07
        self.delay_o2 = 0.65
        self.delay_mmside = 0.75
        self.flow_nominal = 3.3333333333333333e-06

        # Resynchronisation parameters
        self.resync_option = False
        self.resync_co2_option = 'fixed'
        self.delay_co2_resync = 0.07

        # Filter parameters
        self.flow_filter_option = True
        self.flow_filter_cutoff = 2
        self.flow_filter_window = 0.25
        self.co2_filter_option = True
        self.co2_filter_cutoff = 2
        self.co2_filter_window = 0.25
        self.o2_filter_option = True
        self.o2_filter_cutoff = 2
        self.o2_filter_window = 0.25
        self.mmside_filter_option = True
        self.mmside_filter_cutoff = 5
        self.mmside_filter_window = 0.25

        # Speeding options
        self.speed_butter_cutoff = 8
        self.o2_speed_option = True
        self.mmside_speed_option = True
        self.o2_filter_tau = 0.03
        self.mmside_filter_tau = 0.035

        # Crosstalk options
        self.crosstalk_option = 'quadratic'
        self.gain_co2 = 0.09
        self.offset_co2 = 0.9855
        self.crosstalk_A = - 0.0129935471586500
        self.crosstalk_B = 0.00177517115362386
        self.crosstalk_C = 0.117943271855761
        self.crosstalk_D = 0.156764534371431
        self.crosstalk_E = - 0.0194694850557080
        self.crosstalk_F = 0.449457902448671
        self.crosstalk_G = 0.100878454479871
        self.o2_calibration_1 = 1
        self.o2_calibration_2 = 0.2094

        # Breath detection
        self.breath_option = 'spiroware'
        self.co2_threshold = 0.02
        self.cutout_option = True

        # Phase detection
        self.o2_start_wo = 0.6
        self.sf6_start_wi = 0.02
        self.sf6_min_wi = 0.035
        self.sf6_start_wo = 0.02
        self.tracer_washin_end = 0.039

        # Drift correction
        self.o2_drift_option = 'breath-by-breath'
        self.o2_drift_window = [0.7, 0.9]
        self.n2_target_washout = 1.0
        self.n2_target_prephase = 0.2094
        self.sf6_target_washout = 0.2094
        self.sf6_target_washin = 0.21
        self.sf6_target_prephase = 0.2094
        self.drift_correction_window = [0.98, 1.02]
        self.drift_sample_number = 3

        # Physical parameters
        self.mass = {
            'ar': 0.0399480,
            'co2': 0.0440095,
            'n2': 0.0280134,
            'o2': 0.0319988,
            'sf6': 0.1460554
        }
        self.kappa = {
            'ar': 1.6696,
            'co2': 1.2941,
            'n2': 1.4013,
            'o2': 1.3967,
            'sf6': 1.0984,
            'air': 1.4018,
        }
        self.fraction_n2_inert = 0.9881799083753069
        self.ar_air = 0.0093
        self.ar_tracer_gas = 0.0
        self.sf6_air = 0.0
        self.sf6_tracer_gas = 0.04
        self.endtidal_window = [0.9, 0.95]
        self.tracer_reinspired_option = 'integral'

        # Dead spaces
        self.dead_space_set = None
        self.dead_space_postsensor_set3 = 22e-6
        self.dead_space_presensor_set3 = 30e-6
        self.dead_space_postsensor_set2 = 9.5e-6
        self.dead_space_presensor_set2 = 30e-6
        self.dead_space_postsensor_set1 = 3.5e-6
        self.dead_space_presensor_set1 = 2e-6

        # Minimum volume breath detection parameter
        self.breath_min_vol_set3 = 25e-6
        self.breath_min_vol_set2 = 15e-6
        self.breath_min_vol_set1 = 2e-6

        # Slope III parameters
        self.tracer_slope3_set1 = [0.7, 0.95]
        self.tracer_slope3_set2 = [0.65, 0.95]
        self.tracer_slope3_set3 = [0.5, 0.95]
        self.scond_window = [1.5, 6]

    def apply_file_set(self, set_name):
        """
        Function to modify the config variables which depend on the (Dead space reducer) set. This is the sequence in
        which this is determined in the overall program.
        1. Default set is None
        2. We apply the options, if option sets the set, we don't modify it any further
        3. We check the header of the file
        4. We apply the set in the variables if the dead_space_set option is still None
        :param set_name:
        """
        if self.dead_space_set is None:
            self.dead_space_set = set_name
            
    @property
    def breath_detection_minvol(self):
        """
        The behavior of some values of the SpirowareConfig are a bit unique, including the machine pre- and postcap
        dead spaces. We want to be able to set default values for the individual sets, and then let the analysis
        figure out which set it is analysing independently. We also want to be able to force a certain set on an
        analysis. This implementation of the parameters affected by set choice as properties allows them to be
        initialized as a function of other information the first time we need them, without disrupting analysis that
        is independent of set choice should the set choice be ambiguous.
        :return: Breath detection minimum volume
        """
        if self.dead_space_set.lower() == 'set1':
            return self.breath_min_vol_set1
        elif self.dead_space_set.lower() == 'set2':
            return self.breath_min_vol_set2
        elif self.dead_space_set.lower() == 'set3':
            return self.breath_min_vol_set3
        else:
            msg = 'Information about Set insufficient. Set manually in config file.'
            raise LungSimException(output=msg)

    @property
    def dead_space_postsensor(self):
        """
        The behavior of some values of the SpirowareConfig are a bit unique, including the machine pre- and postcap
        dead spaces. We want to be able to set default values for the individual sets, and then let the analysis
        figure out which set it is analysing independently. We also want to be able to force a certain set on an
        analysis. This implementation of the parameters affected by set choice as cached properties allows them to be
        initialized as a function of other information the first time we need them, without disrupting analysis that
        is independent of set choice should the set choice be ambiguous.
        :return: Post-cap machine dead space
        """
        if self.dead_space_set.lower() == 'set1':
            return self.dead_space_postsensor_set1
        elif self.dead_space_set.lower() == 'set2':
            return self.dead_space_postsensor_set2
        elif self.dead_space_set.lower() == 'set3':
            return self.dead_space_postsensor_set3
        else:
            msg = 'Information about Set insufficient. Set manually in config file.'
            raise LungSimException(output=msg)

    @property
    def dead_space_presensor(self):
        """
        The behavior of some values of the SpirowareConfig are a bit unique, including the machine pre- and postcap
        dead spaces. We want to be able to set default values for the individual sets, and then let the analysis
        figure out which set it is analysing independently. We also want to be able to force a certain set on an
        analysis. This implementation of the parameters affected by set choice as cached properties allows them to be
        initialized as a function of other information the first time we need them, without disrupting analysis that
        is independent of set choice should the set choice be ambiguous.
        :return: Pre-cap machine dead space
        """
        if self.dead_space_set.lower() == 'set1':
            return self.dead_space_presensor_set1
        elif self.dead_space_set.lower() == 'set2':
            return self.dead_space_presensor_set2
        elif self.dead_space_set.lower() == 'set3':
            return self.dead_space_presensor_set3
        else:
            msg = 'Information about Set insufficient. Set manually in config file.'
            raise LungSimException(output=msg)

    @property
    def tracer_slope3(self):
        """
        The behavior of some values of the SpirowareConfig are a bit unique, including the machine pre- and postcap
        dead spaces. We want to be able to set default values for the individual sets, and then let the analysis
        figure out which set it is analysing independently. We also want to be able to force a certain set on an
        analysis. This implementation of the parameters affected by set choice as cached properties allows them to be
        initialized as a function of other information the first time we need them, without disrupting analysis that
        is independent of set choice should the set choice be ambiguous.
        :return: Relative volume range used for the fit of slope 3 in the calculation of slope 3 for scond/sacin
        analysis.
        """
        if self.dead_space_set.lower() == 'set1':
            return self.tracer_slope3_set1
        elif self.dead_space_set.lower() == 'set2':
            return self.tracer_slope3_set2
        elif self.dead_space_set.lower() == 'set3':
            return self.tracer_slope3_set3
        else:
            msg = 'Information about Set insufficient. Set manually in config file.'
            raise LungSimException(output=msg)
