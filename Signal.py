import numpy as np


class Signal:
    """
    A signal is one of the basic components of a MBW analysis. A MBW measurement will typically include many different
    types of signals (such as Flow, Tracer gas, etc.), and they all share certain characteristics. These include name,
    unit, length, maximum, minimum, a long numpy array of the signal itself etc.
    """
    def __init__(self, **kwargs):
        self.name = kwargs.pop('name', None)
        self.unit = kwargs.pop('unit', None)
        self.data = kwargs.pop('data', None)
        self.sampling_rate = kwargs.pop('sampling_rate', 200)

    @property
    def length(self):
        return self.data.shape[0]

    @property
    def size(self):
        return self.data.size

    @property
    def last_nonzero(self):
        try:
            return self.data.nonzero()[0][-1]
        except IndexError:
            return self.length - 1

    @property
    def time(self):
        """
        The time signals don't need to be separate signals. They are generally just an accompaniment of another signal.
        :return: an array of the same shape as the signal, filled with time points based on the sampling rate, in [s].
        """
        return np.arange(self.length) / self.sampling_rate
