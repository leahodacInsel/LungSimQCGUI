from LungSimException import LungSimException
from copy import deepcopy


class UnitConversion:
    """
    A class to handle different kinds of unit conversions. Used to convert signals to and from SI units, to and from
    the units used as conventions in different MBW files. Additionally used during the calculation of outcomes to
    achieve the desired units for the outputs.
    The lookup table is a needlessly complicated way of collecting a variety of conversions from one format to another.
    """

    def __init__(self, input_unit, output_unit):
        self.input_unit = input_unit
        self.output_unit = output_unit
        self.lookup = {
            'Volume': {
                'si': 1.0,
                'm^3': 1.0,
                'dm^3': 1e-3,
                'l': 1e-3,
                'dl': 1e-4,
                'cl': 1e-5,
                'cm^3': 1e-6,
                'ml': 1e-6,
                'mm^3': 1e-9
            },
            'Flow': {
                'si': 1.0,
                'm^3/s': 1.0,
                'dm^3/s': 1e-3,
                'l/s': 1e-3,
                'dl/s': 1e-4,
                'cl/s': 1e-5,
                'cm^3/s': 1e-6,
                'ml/s': 1e-6,
                'mm^3/s': 1e-9,
                'm^3/min': 1.0/60,
                'dm^3/min': 1e-3/60,
                'l/min': 1e-3/60,
                'dl/min': 1e-4/60,
                'cl/min': 1e-5/60,
                'cm^3/min': 1e-6/60,
                'ml/min': 1e-6/60,
                'mm^3/min': 1e-9/60
            },
            'Fraction': {
                'si': 1.0,
                '': 1.0,
                '%': 1e-2
            },
            'MolarMass': {
                'si': 1.0,
                'kg/mol': 1.0,
                'g/mol': 1e-3
            },
            'Turnover': {
                'si': 1.0,
                'TO': 1.0
            },
            'Time': {
                'si': 1.0,
                's': 1.0,
                'min': 60,
                'ms': 1e-3,
            },
            'General': {
                '1e9': 1e-9,
                '1e8': 1e-8,
                '1e7': 1e-7,
                '1e6': 1e-6,
                '1e5': 1e-5,
                '1e4': 1e-4,
                '1e3': 1e-3,
                '1e2': 1e-2,
                '1e1': 1e-1,
                'si': 1,
                '1e-1': 1e1,
                '1e-2': 1e2,
                '1e-3': 1e3,
                '1e-4': 1e4,
                '1e-5': 1e5,
                '1e-6': 1e6,
                '1e-7': 1e7,
                '1e-8': 1e8,
                '1e-9': 1e9,
            }
        }

    @property
    def factor(self):
        success = 0
        for key in self.lookup:
            entry = self.lookup[key]
            # If both units appear within the same category, calculate the conversion factor and indicate success.
            if self.input_unit in entry and self.output_unit in entry:
                factor_value = entry[self.input_unit] / entry[self.output_unit]
                success = 1
                break
        if success == 0:
            message = 'Unit conversion from ' + self.input_unit + ' to ' + self.output_unit + ' failed.'
            raise LungSimException(output=message)

        return factor_value

    def convert(self, input):
        output = deepcopy(input)
        return output * self.factor

