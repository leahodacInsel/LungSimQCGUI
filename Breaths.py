import numpy as np


class Breaths:
    """
    Object to carry information about the output of a breath detection algorithm. Contains the inspiration and
    expiration starts and ends, index ranges for each breath, number of breaths etc.
    """
    def __init__(self, insp_start, insp_end, exp_start, exp_end):
        self.insp_start = insp_start.astype(int)
        self.insp_end = insp_end.astype(int)
        self.exp_start = exp_start.astype(int)
        self.exp_end = exp_end.astype(int)

    @property
    def number(self):
        """
        :return: number of breaths or length of breaths structure
        """
        return self.insp_start.shape[0]

    @property
    def numbers(self):
        """
        Index range going from 0 to the size of breaths. This allows easy for loops across all breaths.
        :return: range of breath numbers
        """
        return np.arange(start=0, stop=self.number).astype(int)

    def insp_range(self, i):
        """
        returns a numpy array of the inspiratory range of breath number i. We add +1 to insp_end in order to include it
        in the resulting range
        :param i: breath number
        :return: range from insp_start to insp_end
        """
        return np.arange(start=self.insp_start[i], stop=self.insp_end[i] + 1).astype(int)

    def insp_length(self, i):
        """
        returns the inspiratory length of breath i.
        :param i: breath number
        """
        return self.insp_end[i] - self.insp_start[i] + 1

    def exp_range(self, i):
        """
        returns a numpy array of the expiratory range of breath number i. We add +1 to exp_end in order to include it
        in the resulting range
        :param i: breath number
        :return: range from exp_start to exp_end
        """
        return np.arange(start=self.exp_start[i], stop=self.exp_end[i] + 1).astype(int)

    def exp_length(self, i):
        """
        returns the expiratory length of breath i.
        :param i: breath number
        """
        return self.exp_end[i] - self.exp_start[i] + 1

    def range(self, i):
        """
        returns a numpy array of the whole breath range of breath number i. We add +1 to exp_end in order to include it
        in the resulting range
        :param i: breath number
        :return: range from insp_start to exp_end
        """
        return np.arange(start=self.insp_start[i], stop=self.exp_end[i] + 1).astype(int)

    def length(self, i):
        """
        returns the length of breath i.
        :param i: breath number
        :return: length as the defined by the difference between exp_end and insp_start
        """
        return self.exp_end[i] - self.insp_start[i] + 1
