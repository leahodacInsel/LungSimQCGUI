import numpy as np
from easygui import fileopenbox
from LungSimException import LungSimException
from Signal import Signal
from functools import lru_cache


class File:
    """
    This class is intended to handle multiple-breath washout files. This acts as the mother class to the classes of
    specific raw data of each MBW platform. There are specific functions related to data files which are only relevant
    to Spiroware files, for example. The File class will handle all functions that are common to MBW files, and the
    device-specific File classes will inherit from this one.

    The most common way for a file to be instantiated will be via a path to the file.
    There is another use to being able to instantiate a file from existing signals, so that signal data can be written
    to files.
    """

    def __init__(self, **kwargs):
        """
        fullname: path + filename
        filename: This corresponds to the name of the file itself. .txt files are the standard for Spiroware files.
            The default value for the name is None, to avoid files being saved without a valid name given.
        path: Path to the file in question. Default is None, to avoid files being saved.
        signals: Dictionary containing the Signal objects contained within the file. Default is an empty
            dictionary.
        first_line: string of the first line of the file, which usually contains signal names and potentially
            additional information
        """
        self.fullname = kwargs.pop('fullname', None)
        self.filename = kwargs.pop('filename', None)
        self.path = kwargs.pop('path', None)
        self.signals = kwargs.pop('signals', dict())
        self.first_line = kwargs.pop('first_line', None)

    @staticmethod
    def filename_from_fullname(fullname):
        """
        Function to extract from a full name (path + filename) just the filename, by splitting the fullname by the
        file/folder separators and selecting the last item
        :param fullname: path + filename
        :return: filename
        """
        return fullname.split('\\')[-1]

    @staticmethod
    @lru_cache
    def first_line_from_fullname(fullname):
        """
        Function to centralise reading of first line of file in one place. Function is static so it can be accessed by
        classmethods used as constructors, and cached because reading a .txt file can be costly in processing.
        :param fullname: path + filename of .txt file
        :return: First line of file
        """
        file = open(fullname)
        first_line = file.readline()
        file.close()
        return first_line

    @classmethod
    def from_fullname(cls, fullname):
        """
        Function to redirect the program to the various different data-format opening functions implemented. For now:
        .txt
        :param fullname:
        :return: class instance initiated via .txt file
        """
        if '.txt' in fullname:
            return cls.from_txt(fullname=fullname)
        else:
            raise LungSimException(output="Non-standard Spiroware input file (not .txt)")

    @classmethod
    def from_selection_box(cls):
        fullname = fileopenbox(
            msg="Select single MBW file",
            multiple=False
        )
        return cls.from_fullname(fullname=fullname)

    @classmethod
    def from_txt(cls, fullname):
        """
        Instantiates a File from a full name (path + filename), by using the methods of file to extract the
        necessary information.
        :param fullname:
        :return:
        """
        file = cls(
            signals=cls.signals_from_txt(fullname=fullname),
            filename=cls.filename_from_fullname(fullname=fullname),
            path=cls.path_from_fullname(fullname=fullname),
            fullname=fullname,
            first_line=cls.first_line_from_fullname(fullname=fullname)
        )
        return file

    @classmethod
    def path_from_fullname(cls, fullname):
        """
        Function to extract from a full name (path + filename) just the name of the path. We do this by extracting the
        filename and then replacing the filename in the original fullname with an empty string.
        :param fullname: path + filename
        :return: path
        """
        remove = cls.filename_from_fullname(fullname=fullname)
        return fullname.replace(remove, '')

    @staticmethod
    def signal_names():
        """
        Aliases of signal names which can be encountered in MBW files. The signal names dictionary
        contains for each signal the various spellings which could be encountered in file headers, but which belong
        to the same signal. The list below is not exhaustive, nor is it contained to the smallest possible number
        of names. The final name as it is interpreted by LungSim is the dictionary key, for example for Helium it
        would be "he".
        """
        return {
            'he': ['he', 'He'],
            'sf6': ['sf6', 'SF6'],
            'mmside': ['mmside', 'MMss'],
            'o2': ['o2', 'O2'],
            'flowside': ['flowside', 'SampleFlow'],
            'n2': ['n2', 'N2'],
            'res': ['res', 'Res'],
            'delay_o2': ['delay_o2', 'DelayO2'],
            'delay_mmside': ['delay_mmside', 'DelayMMss'],
            'time': ['time', 'Time', 'Tid'],
            'mmmain': ['mmmain', 'MMms', 'MM'],
            'ar': ['ar', 'Ar', 'AR'],
            'co2': ['co2', 'CO2', 'CO2_', 'CO2_RAW'],
            'co2star': ['co2star', 'CO2*'],
            'flowmain': ['flowmain', 'Flow', 'Flow_BTPS', 'flow_btps']
        }

    @staticmethod
    def signal_units():
        """
        The signal units dictionary contains information about the units of the signal found within
        MBW files. Should any type of file have different units a daughter class can be created.
        """
        return {
            'he': '%',
            'sf6': '%',
            'mmside': 'g/mol',
            'o2': '%',
            'flowside': 'l/s',
            'n2': '%',
            'res': '',
            'delay_o2': 's',
            'delay_mmside': 's',
            'time': 'ms',
            'mmmain': 'g/mol',
            'ar': '%',
            'co2': '%',
            'co2star': '%',
            'flowmain': 'l/s'
        }

    @classmethod
    def signals_from_txt(cls, fullname):
        """
        Function to extract a dictionary of signals from a txt file. Aliases of the same kind of signal are converted
        to a standardized name. The aliases of signal names, as well as the units corresponding to individual signals
        are stored in a DictionaryAliases object.
        aliases: To read signals we need a dictionary of aliases of signal names. From there the necessary information
        can be imported to interpret which kinds of signals are contained within a file.
        :param fullname: path + filename of file from which signals are to be extracted
        :return: dictionary of signals recognized within the file.
        """
        # Attempt to read signal data, if unsuccessful raise exception
        try:
            signal_data = np.loadtxt(fullname, skiprows=1)
        except ValueError:
            raise LungSimException(output='.txt file did not contain appropriate signal data.')
        # Extract signals names from file
        header = cls.first_line_from_fullname(fullname)
        header_items = header.split()
        number_of_signals = np.shape(signal_data)[1]
        # Go through the list of signals inherent to the File class, and if a list contains a tag that can be found in
        # the header of the file to be read, enter it into the signal dictionary
        signals = dict()
        signal_names = cls.signal_names()
        signal_units = cls.signal_units()

        for signal_type in signal_names:
            for spelling in signal_names[signal_type]:
                # If the header name appears in any of the possible tags, and said tag actually corresponds to data from
                # the input_data...:
                if spelling in header_items:
                    if header_items.index(spelling) < number_of_signals:
                        # Create a signal object from the correct column of the file data
                        signals[signal_type] = Signal(
                            name=signal_type,
                            unit=signal_units[signal_type],
                            data=np.array(signal_data[:, header_items.index(spelling)]),
                        )
                    break

        return signals


