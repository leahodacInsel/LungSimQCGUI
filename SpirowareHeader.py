class SpirowareHeader:
    """
    Class to deal with converting Spiroware Header to and from a txt-file string into a dictionary which can be applied
    to a SpirowareConfig object. Two items in the Spiroware header have different behavior from the others, the Pre and
    Post-cap dead space. They affect all dead spaces upon initialization. Then if the dead spaces for either individual
    sets or all sets have to be overwritten by the options file, they can be.
    """
    def __init__(self, **kwargs):
        self.raw_string = kwargs.pop('raw_string', '')
        self.translation = {
            'DSR': {'dead_space_set': 0},
            'T': {'temp_inspired': 1},
            'P': {'pressure_ambient': 100},
            'BTPS': {'btps_option': 0},
            'H': {'humidity_inspired': 0.01},
            'FT': {'temp_flowhead': 1},
            'FH': {'humidity_flowhead': 0.01},
            'BT': {'temp_body': 1},
            'BH': {'humidity_body': 0.01},
            'ManBTPS': {'X': 0},
            'ManFIC': {'X': 0},
            'FIC': {'X': 0},
            'FEC': {'X': 0},
            'CO2C': {'atp_option': 0},
            'VS': {'breath_detection_minvol': 1e-6},
            'Pre': {'dead_space_presensor_set1': 1e-6,
                    'dead_space_presensor_set2': 1e-6,
                    'dead_space_presensor_set3': 1e-6},
            'Post': {'dead_space_postsensor_set1': 1e-6,
                     'dead_space_postsensor_set2': 1e-6,
                     'dead_space_postsensor_set3': 1e-6},
            'SFC': {'sampleflow_option': 0},
            'O2RTC': {'o2_speed_option': 0},
            'O2RT': {'o2_filter_tau': 1},
            'MMssRTC': {'mmside_speed_option': 0},
            'MMssRT': {'mmside_filter_tau': 1},
            'O2': {'delay_o2': 1},
            'CO2': {'delay_co2': 1},
            'MMss': {'delay_mmside': 1},
            'WOO2': {'sf6_target_washout': 0.01},
            'WOCO2': {'X': 0},
            'WOAR': {'ar_air': 0.01},
            'WIO2': {'sf6_target_washin': 0.01},
            'WICO2': {'X': 0},
            'WIAR': {'ar_tracer_gas': 0.01},
            'WIHE': {'X': 0},
            'WISF6': {'sf6_tracer_gas': 0.01},
            'ByH': {'X': 0}
        }
        self.raw = self.string_translation(self.raw_string)
        self.translated = self.headers_translation(self.raw)

    @classmethod
    def string_translation(cls, string):

        string = string.replace('(', ';')
        string = string.replace(')', ';')
        item_list = string.split(";")
        item_list = [i for i in item_list if i]
        raw = dict()
        for item in item_list:
            split_items = item.split("=")
            if len(split_items) > 1:
                raw[split_items[0]] = split_items[1]

        return raw

    def headers_translation(self, raw):
        """
        translates raw header into translated header and returns translated header
        :return: translated header
        """
        translated = dict()
        for key_raw, value_raw in raw.items():
            if key_raw in self.translation:
                for key_translated, value_translated in self.translation[key_raw].items():
                    # If the value in the translation dictionary is anything other than 0, then that parameter is a
                    # number to be converted by the value found there. Otherwise, it is to be treated as a boolean
                    # option type translation. If it is not boolean, it is a string and to be copied.
                    if bool(value_translated):
                        translated[key_translated] = float(value_translated) * float(value_raw)
                    elif value_raw == 'True' or value_raw == '1':
                        translated[key_translated] = True
                    elif value_raw == 'False' or value_raw == '0':
                        translated[key_translated] = False
                    else:
                        translated[key_translated] = value_raw

        return translated
