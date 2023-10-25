from datetime import datetime
from easygui import fileopenbox
from easygui import diropenbox
from LungSimException import LungSimException
from math import nan
import numpy as np
# Is needed even if shaded out, for excel writing of pandas
import openpyxl
from os import listdir
from pandas import DataFrame
from LungSimProgressBar import LungSimProgressBar
from SpirowareFile import SpirowareFile
from SpirowareN2MBW import SpirowareN2MBW
from SpirowareSF6MBWashout import SpirowareSF6MBWashout
from SpirowareSF6MBWashin import SpirowareSF6MBWashin
from UnitConversion import UnitConversion
from YAMLTools import dict_from_yaml


class MBWAnalysis:
    """
    Whenever a proper washout analysis is performed, we need a bunch of information. This information, and how the files
    are analysed and their outputs exported, is governed by this class. A MBW analysis has the following components:
    - A filename or a list of filenames to be analysed
    - A results configuration which determines what outcomes should be calculated
    - An options configuration which governs options set by the user (optional)

    It then has additional attributes which it makes based on the above:
    - An output trial table which contains the desired trial table outcomes for all files.

    It also has the methods required to interpret the results configuration file and execute it.
    """

    def __init__(self, **kwargs):
        self.filenames = kwargs.pop('filenames', None)

        self.options_config = kwargs.pop('options_config', None)
        if self.options_config is not None:
            self.options = dict_from_yaml(self.options_config)
        else:
            self.options = dict()

        self.results_path = kwargs.pop('results_path', None)
        self.current_index = 0
        self.current_file = None
        self.current_analysis = None
        self.trial_table = self.initial_trial_table

    @classmethod
    def from_fileopenbox(cls, **kwargs):
        """
        Function to initialize a Washout Analysis instance from the easy gui interface. Some parameters may be pre-
        declared using the kwargs, those left empty will cause an easygui window to pop up asking for input.
        :param kwargs: pre-declared attributes, and optional choice to determine which method is used for the selection
        of filenames ('folder' or 'files').
        :return: An instance of a MBWAnalysis
        """
        filenames = kwargs.pop('filenames', None)
        filenames_method = kwargs.pop('filenames_method', None)
        options_config = kwargs.pop('options_config', None)
        results_path = kwargs.pop('results_path', None)

        # Parameter config
        if options_config is None:
            options_config = fileopenbox(msg='Select a option config file.', multiple=False)

        # Filenames
        if filenames is None:
            if filenames_method == 'files' or filenames_method is None:
                filenames = fileopenbox(msg='Select a list of files to be analysed.', multiple=True)
            elif filenames_method == 'folder':
                folder = diropenbox(msg='Select a folder containing the files to be analysed.')
                filenames = listdir(folder)

        # Output folder
        if results_path is None:
            results_path = diropenbox(msg='Select folder where outputs are saved.')

        return cls(
            filenames=filenames,
            options_config=options_config,
            results_path=results_path
        )

    @property
    def file_number(self):
        return len(self.filenames)

    @property
    def file_numbers(self):
        return np.arange(start=0, stop=self.file_number)

    def calculate_all(self):
        """
        Function to calculate everything, going backwards from the result_config
        :return: None
        """
        # Initializing progress bar
        bar = LungSimProgressBar('Analyzing', max=self.file_number)

        for i in self.file_numbers:
            self.current_index = i
            # Attempt to parse the file and return one of the supported file types. If parsing fails, carry on with the
            # analysis while skipping the file.
            try:
                self.current_file = self.parse_file(self.filenames[i])
            except LungSimException:
                continue

            for self.current_analysis in self.parse_analysis(self.current_file, self.options):
                self.calculate_single()
                self.write_signal_file()
                self.write_breath_table()
                self.current_analysis.reset()
            bar.next()

        if 'write_trial_table' in self.options:
            if self.options['write_trial_table']:
                self.write_trial_table()

        bar.finish()

    def calculate_single(self):
        """
        Calculate the outcomes for a single file, identified by current_analysis. Always include the filename as a
        separate output. Otherwise either export the result, or the relevant LungSimException error message as an
        output.
        :return: None
        """
        # Add filename to trial table
        self.trial_table['Filename'] = np.append(
            arr=self.trial_table['Filename'],
            values=self.filenames[self.current_index].split('\\')[-1]
        )

        # Add all other values to trial table
        if 'trial_table' in self.options:
            for trial_output in self.options['trial_table']:
                try:
                    values = self.current_analysis.calculate_output(
                        args=self.options['trial_table'][trial_output]
                    )
                except LungSimException as error:
                    values = error.output

                # It is useful to comment these three lines for the purposes of debugging, but for release versions it
                # provides excellent robustness to long reloads to simply catch all unexpected errors in one place.
                except:
                    values = 'Unexpected error in file {}.'.format(self.current_file.filename)
                    print('Unexpected error in file {}.'.format(self.current_file.filename))

                self.trial_table[trial_output] = np.append(
                    arr=self.trial_table[trial_output],
                    values=values
                )

    @property
    def initial_trial_table(self):
        """
        Initializes the trial table based on the number of files in filenames, and the number of desired trial table
        outputs in result_config. We initialize each array as an object array as the outputs can be numbers or error
        messages interchangeably.
        :return: Dictionary of empty outputs.
        """
        trial_table = dict()
        trial_table['Filename'] = np.array([]).astype('object')
        if 'trial_table' in self.options:
            for output in self.options['trial_table']:
                trial_table[output] = np.array([]).astype('object')

        return trial_table

    def load_single(self):
        """
        Function to load / prepare the MBW measurement which is in position of current_index in the filelist
        """
        # TODO: Very clumsy and mostly redundant with the "calculate all" funtion above. Rework this class.
        self.current_file = self.parse_file(self.filenames[self.current_index])
        self.current_analysis = self.parse_analysis(self.current_file, self.options)[-1]

    @staticmethod
    def parse_file(filename):
        """
        :param filename: full name of file (path + filename)
        :return: File object based on file characteristics
        """
        # TODO: Parse nature of input file and automatically select which file type
        if '.txt' in filename:
            return SpirowareFile.from_fullname(filename)
        else:
            raise LungSimException(output='Only .txt files currently supported.')

    @staticmethod
    def parse_analysis(file, options):
        """
        Function to parse which kind of analysis should be returned depending on the attributes of the input file
        :param file: File-like object
        :param options:
        :return: A list of MBW-like objects
        """
        if isinstance(file, SpirowareFile):
            if file.method == 'n2_mbw':
                return [
                    SpirowareN2MBW(file=file, options=options)
                ]
            elif file.method == 'sf6_mbw':
                return [
                   SpirowareSF6MBWashin(file=file, options=options),
                   SpirowareSF6MBWashout(file=file, options=options)
                ]
            else:
                raise LungSimException(output='Only N2 Washout and SF6 Washout/Washin currently supported')
        else:
            raise LungSimException(output='Only Spiroware Washout currently supported')

    def write_signal_file(self):
        """
        If desired, format and write a .txt file containing the signals specified in the config
        :return: None
        """
        try:
            if 'write_signal_files' in self.options:
                if self.options['write_signal_files']:
                    if 'signal_files' in self.options:
                        # Determine maximum length
                        signal_lengths = np.array([])
                        # Eliminate elements that cannot be calculated
                        dictkeys = list(self.options['signal_files'])

                        for name in dictkeys:
                            try:
                                details = self.options['signal_files'][name]
                                signal = eval('self.current_analysis.'+details['signal'])
                                signal_lengths = np.append(arr=signal_lengths, values=signal.length)
                            except AttributeError:
                                self.options['signal_files'].pop(name)
                                signal_lengths = np.append(arr=signal_lengths, values=1)
                            except LungSimException:
                                signal_lengths = np.append(arr=signal_lengths, values=1)
                                pass
                                # self.options['signal_files'].pop(name)

                        columns = len(self.options['signal_files'])
                        max_length = int(np.max(signal_lengths))
                        output = np.full((max_length, columns), nan).astype('object')
                        k = 0
                        delimiter = '\t'
                        header = ''
                        fmt = list()
                        for name in self.options['signal_files']:
                            try:
                                details = self.options['signal_files'][name]
                                signal = eval('self.current_analysis.'+details['signal'])
                                if 'unit' in details:
                                    try:
                                        conversion = UnitConversion(input_unit='si', output_unit=details['unit'])
                                        output[:signal.length, k] = signal.data * conversion.factor
                                        fmt.append('%4f')
                                    except LungSimException as error:
                                        output[0, k] = error.output
                                        output[1:max_length, k] = ''
                                        fmt.append('%s')

                                else:
                                    output[:signal.length, k] = signal.data
                                    fmt.append('%4f')

                            except LungSimException as error:
                                output[0, k] = error.output
                                output[1:max_length, k] = ''
                                fmt.append('%s')

                            header = header + name + delimiter
                            k += 1

                        now = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
                        # prefix = 'B-'+self.current_analysis.method_name+'-'+now+'-'
                        prefix = 'B-'
                        if 'A-' in self.current_file.filename:
                            filename_new = self.current_file.filename.replace('A-', prefix)
                        else:
                            filename_new = prefix+self.current_file.filename

                        save_name = self.results_path+'\\'+filename_new
                        np.savetxt(fname=save_name,
                                   X=output,
                                   fmt=fmt,
                                   delimiter=delimiter,
                                   header=header,
                                   comments='')

        except:
            # This is terrible terrible code practise. It simply serves to make sure that analyses of large numbers of
            # files encounter no possibility of failing.
            # TODO: Make this obsolete by writing a better B-file writing function.
            pass

    def write_breath_table(self):
        """
        Function to write breath-tables (one file per input file, containing breath-by-breath information
        :return:
        """
        try:
            if 'write_breath_table' in self.options:
                if self.options['write_breath_table']:
                    if 'breath_table' in self.options:
                        # Determine maximum length
                        output_lengths = np.array([])
                        dictkeys = list(self.options['breath_table'])

                        for name in dictkeys:
                            try:
                                details = self.options['breath_table'][name]
                                output = eval('self.current_analysis.'+details['variable'])
                                output_lengths = np.append(arr=output_lengths, values=output.size)
                            except AttributeError:
                                self.options['breath_table'].pop(name)
                            except LungSimException:
                                pass

                        max_length = int(np.max(output_lengths))
                        table = dict()
                        for name in self.options['breath_table']:
                            try:
                                table_variable = np.full(max_length, nan).astype('object')
                                details = self.options['breath_table'][name]
                                variable = eval('self.current_analysis.'+details['variable'])
                                if 'unit' in details:
                                    try:
                                        conversion = UnitConversion(input_unit='si', output_unit=details['unit'])
                                        table_variable[:variable.size] = variable * conversion.factor
                                    except LungSimException as error:
                                        table_variable[0] = error.output
                                else:
                                    table_variable[:variable.size] = variable

                            except LungSimException as error:
                                table_variable[0] = error.output

                            table[name] = table_variable

                        now = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
                        prefix = 'BreathTable-'+self.current_analysis.method_name+'-'+now+'-'
                        if 'A-' in self.current_file.filename:
                            filename_new = self.current_file.filename.replace('A-', prefix)
                        else:
                            filename_new = prefix+self.current_file.filename

                        if '.txt' in self.current_file.filename:
                            filename_new = filename_new.replace('.txt', '.xlsx')
                        else:
                            filename_new = filename_new+'.xlsx'

                        save_name = self.results_path+'\\'+filename_new
                        table = DataFrame(table)
                        table.to_excel(excel_writer=save_name, index=False)

        except:
            # TODO: Improve this. Just here to make sure that the run never crashes due to B-table writing.
            pass

    def write_trial_table(self):
        """
        Write trial table to current output folder
        :return: None
        """
        now = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        output_name = self.results_path + "\\" "Trial_table_" + now + ".xlsx"
        output = DataFrame(self.trial_table)
        output.to_excel(excel_writer=output_name, index=False)
