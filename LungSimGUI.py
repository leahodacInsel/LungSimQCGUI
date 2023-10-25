from easygui import fileopenbox
from easygui import diropenbox
import tkinter as tk
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from MBWAnalysis import MBWAnalysis
from LungSimException import LungSimException


class FileNavigation(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg='white')
        self.parent = parent

        # Listbox
        self.listbox_frame = None
        self.listbox = None
        self.y_scrollbar = None
        self.x_scrollbar = None

        # Buttons
        self.load_button = None
        self.next_button = None
        self.previous_button = None

        self.init_listbox()
        self.init_buttons()

    def init_listbox(self):
        # Create a frame for the listbox and the scrollbar
        self.listbox_frame = tk.Frame(self)
        self.listbox_frame.pack(side=tk.TOP, fill=tk.X, expand=1)
        self.listbox_frame.rowconfigure(0, weight=1)
        self.listbox_frame.rowconfigure(1, weight=0)
        self.listbox_frame.columnconfigure(0, weight=1)
        self.listbox_frame.columnconfigure(1, weight=0)
        # Creating a Scrollbar and
        # attaching it to root window
        self.y_scrollbar = tk.Scrollbar(self.listbox_frame, orient=tk.VERTICAL)
        self.x_scrollbar = tk.Scrollbar(self.listbox_frame, orient=tk.HORIZONTAL)
        # Adding Scrollbar to the right side of root window
        self.y_scrollbar.grid(column=1, row=0, sticky='ns')
        self.x_scrollbar.grid(column=0, row=1, sticky='ew')
        # Initialize listbox containing names of A-files
        self.listbox = tk.Listbox(self.listbox_frame, height=20)
        # Adding Listbox to the left side of frame
        # self.listbox.pack(side = tk.LEFT, fill = tk.BOTH, expand = 1)
        self.listbox.grid(column=0, row=0, sticky='nsew')
        # Insert elements into the listbox
        for i in self.parent.filenames:
            insert = i.split('\\')[-1]
            self.listbox.insert(tk.END, insert)

        # Attaching Listbox to Scrollbar. Since we need to have a vertical scroll we use yscrollcommand
        self.listbox.config(yscrollcommand=self.y_scrollbar.set)
        self.listbox.config(xscrollcommand=self.x_scrollbar.set)
        # setting scrollbar command parameter to listbox.yview method its yview because we need to have a vertical view
        self.y_scrollbar.config(command=self.listbox.yview)
        self.x_scrollbar.config(command=self.listbox.xview)

    def init_buttons(self):
        # Load button
        self.load_button = tk.Button(master=self, text='Load', bg='light gray', command=self.load_file)
        self.load_button.pack(side=tk.TOP, fill=tk.X, expand=True)
        # Previous button
        self.previous_button = tk.Button(master=self, text='Previous', bg='light gray', command=self.previous_file)
        self.previous_button.pack(side=tk.TOP, fill=tk.X, expand=True)
        # Next button
        self.next_button = tk.Button(master=self, text='Next', bg='light gray', command=self.next_file)
        self.next_button.pack(side=tk.TOP, fill=tk.X, expand=True)

    def load_file(self):
        self.parent.mbw_analysis.current_index = self.listbox.curselection()[0]
        self.parent.refresh_all()

    def previous_file(self):
        self.parent.mbw_analysis.current_index = max(0, self.parent.index - 1)
        self.parent.refresh_all()

    def next_file(self):
        max_index = len(self.parent.filenames)
        self.parent.mbw_analysis.current_index = min(max_index - 1, self.parent.index + 1)
        self.parent.refresh_all()


class IntroWindow(tk.Frame):
    # Intro screen consists of three buttons:
    # 1. Load config button (opens a window for the user to select config file)
    # 2. Load A-file-list button (opens a window for the user to select a list of A-files)
    # 3. Start button (activates only when a config file and a A-file-list have been selected.
    #    Activates the main window with the config file and A-file-list, then destroys the intro screen)
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        # Parent of the intro screen is root
        self.parent = parent
        self.pack(side="top", fill="both", expand=True)

        # Main information of intro window
        self.filenames = None
        self.options_config = None
        self.results_path = None

        # Title of the window
        self.parent.wm_title("LungSim MBW Quality Control")

        # Lifts the window and makes it zoomed
        self.parent.lift()
        self.parent.attributes('-topmost', True)
        self.parent.after_idle(self.parent.attributes, '-topmost', False)

        # Window resolution
        x_resolution = 400
        y_resolution = 400

        # Turn resolution parameters into a string
        self.parent.geometry(str(x_resolution) + 'x' + str(y_resolution))

        # Define button states
        self.files_ready = False
        self.options_config_ready = False
        self.results_path_ready = False

        # Button attributes
        button_height = 3
        button_width = 30
        button_x_padding = 5
        button_y_padding = 15

        # Create the buttons
        # Button for loading config file
        self.config_button = tk.Button(master=self,
                                       text='Load config file',
                                       bg='light gray',
                                       command=self.load_config,
                                       height=button_height,
                                       width=button_width)

        # Button for loading files
        self.files_button = tk.Button(master=self,
                                      text='Load A-files',
                                      bg='light gray',
                                      command=self.load_files,
                                      height=button_height,
                                      width=button_width)

        # Button for defining results directory
        self.results_button = tk.Button(master=self,
                                        text='Select results directory',
                                        bg='light gray',
                                        command=self.load_results_path,
                                        height=button_height,
                                        width=button_width)

        # Button for starting a batch run
        self.start_button = tk.Button(master=self,
                                      text='Go!',
                                      bg='light gray',
                                      command=self.start_main,
                                      height=button_height,
                                      width=button_width)

        # Pack the buttons
        self.config_button.pack(side=tk.TOP, padx=button_x_padding, pady=button_y_padding)
        self.files_button.pack(side=tk.TOP, padx=button_x_padding, pady=button_y_padding)
        self.results_button.pack(side=tk.TOP, padx=button_x_padding, pady=button_y_padding)
        self.start_button.pack(side=tk.TOP, padx=button_x_padding, pady=button_y_padding)
        # Disable the go button by default
        self.start_button['state'] = 'disabled'

    def load_config(self):
        self.options_config = fileopenbox(msg='Select config file', multiple=False, filetypes="*.yaml")
        # Check if a file was selected by the user
        if self.options_config is None:
            self.options_config_ready = False
            self.config_button.configure(bg='light gray')
        else:
            self.options_config_ready = True
            self.config_button.configure(bg='green')
        # Updates the check on whether program can be launched
        self.start_button_update()

    def load_files(self):
        self.filenames = fileopenbox(msg='Select A-files', multiple=True, filetypes="*txt")
        if self.filenames is None:
            self.files_ready = False
            self.files_button.configure(bg='light gray')
        else:
            self.files_ready = True
            self.files_button.configure(bg='green')
        # Updates the check on whether program can be launched
        self.start_button_update()

    def load_results_path(self):
        self.results_path = diropenbox(msg='Select folder where outputs are saved.')
        if self.results_path is None:
            self.results_path_ready = False
            self.results_button.configure(bg='light gray')
        else:
            self.results_path_ready = True
            self.results_button.configure(bg='green')
        # Updates the check on whether program can be launched
        self.start_button_update()

    def start_button_update(self):
        if self.files_ready and self.options_config_ready and self.results_path_ready:
            self.start_button['state'] = 'normal'
        else:
            self.start_button['state'] = 'disabled'

    def start_main(self):
        # Removes current frame
        self.pack_forget()
        # Starts the QC run
        MainWindow(self.parent, self.options_config, self.filenames, self.results_path)


class ResultsWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg='#5A5A5A')
        self.parent = parent
        self.text = ''
        self.label = None

    def refresh(self):
        """
        Refresh display text based on state of MBW analysis
        """
        # Remove previous results if they exist
        for widget in self.winfo_children():
            widget.destroy()

        text = ''
        text += 'Results:\n'
        # LCI 2.5%
        args = {
            'variable': 'lci_ao',
            'method': 'critical',
            'concentration': 2.5,
            'unit': 'TO'}
        text += 'LCI 2.5% [TO]: {:.3f}\n'.format(self.mbw.calculate_output(args))

        # FRC
        args = {
            'variable': 'frc_ao',
            'method': 'critical',
            'concentration': 2.5,
            'unit': 'l'}
        text += 'FRC [L]: {:.3f}\n'.format(self.mbw.calculate_output(args))

        self.label = tk.Label(master=self, text=text, justify=tk.LEFT, anchor='nw', wraplength=200)
        self.label.pack(side=tk.TOP, fill=tk.BOTH, padx=5, pady=5, expand=1)

    @property
    def mbw(self):
        """
        Shortcut property for simpler access to current signals
        :return: currently activated MBW object
        """
        return self.parent.mbw_analysis.current_analysis


class MainWindow(tk.Frame):
    def __init__(self, parent, options_config, filenames, results_path):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        # Fullscreen
        self.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Initialize a mbw analysis
        # self.results_path = r"C:\Users\I0326962\Local Files\Former Desktop"
        self.results_path = results_path
        # self.files = [r"C:\Users\I0326962\Local Files\Former Desktop\Ecomedics reloads\Set 3
        # raw\A-20180608-132334-HIRLEA02071999-N2MultiBreathWashoutTest-Set3.txt"]
        self.filenames = filenames
        # self.config = r"C:\Github\LungSim config files\TestConfig_GUI.yaml"
        self.options_config = options_config
        self.mbw_analysis = MBWAnalysis.from_fileopenbox(
            filenames=self.filenames,
            options_config=self.options_config,
            results_path=self.results_path
        )
        self.mbw_analysis.load_single()

        # Divide the window into grids (The main item is the graphs at 1,1, the rest organizes around it)
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=5)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        # Assemble elements of main window
        self.plots = SignalWindow(self)  # .grid(column=0, row=0, sticky='nsew')
        self.plots.grid(column=0, row=0, sticky='nsew')
        self.plots.refresh()

        self.menu_bar = MenuBar(self)

        self.results = ResultsWindow(self)
        self.results.grid(column=1, row=0, sticky='nsew')
        self.results.refresh()

        self.file_navigation = FileNavigation(self)
        self.file_navigation.grid(column=2, row=0, sticky='nsew')

        # Change icon to LS icon
        self.parent.iconbitmap('LungSimIcon.ico')

    def refresh_all(self):
        """
        Function to refresh all components of the UI (opon file change)
        """
        self.mbw_analysis.load_single()
        self.plots.refresh()
        self.results.refresh()


class MenuBar:
    def __init__(self, parent):
        self.parent = parent
        self.menu_bar = tk.Menu(parent)

        # Define the file menu
        self.file_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.file_menu.add_command(label='Load config file', command=self.config_switch)
        self.file_menu.add_command(label='Load A-files', command=self.afile_switch)
        self.file_menu.add_command(label='Select results folder', command=self.results_folder_switch)
        self.menu_bar.add_cascade(label='File', menu=self.file_menu)

        # Define signal menu
        self.signal_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.signal_menu.add_command(label='Tracer', command=self.tracer_switch)
        self.signal_menu.add_command(label='Flow', command=self.flow_switch)
        self.signal_menu.add_command(label='O2', command=self.o2_switch)
        self.signal_menu.add_command(label='CO2', command=self.co2_switch)
        self.signal_menu.add_command(label='MM ss', command=self.mmside_switch)
        self.signal_menu.add_separator()
        self.signal_menu.add_command(label='Volume x-axis', command=self.x_axis_switch)
        self.signal_menu.add_command(label='End-tidal concentrations', command=self.endtidal_switch)
        self.menu_bar.add_cascade(label="Signals", menu=self.signal_menu)

        # Add the menu bar to the root window
        self.parent.parent.config(menu=self.menu_bar)

        self.bold_items()

    def config_switch(self):
        options_config = fileopenbox(msg='Select config file', multiple=False, filetypes="*.yaml")
        # Check if a file was selected by the user
        if options_config is not None:
            self.parent.options_config = options_config

        self.parent.mbw_analysis = MBWAnalysis.from_fileopenbox(
            filenames=self.parent.filenames,
            options_config=self.parent.options_config,
            results_path=self.parent.results_path
        )

        # Recalculate outcomes
        self.parent.refresh_all()

    def afile_switch(self):
        return

    def results_folder_switch(self):
        return

    def tracer_switch(self):
        self.parent.plots.signal_info["tracer"]["plot"] = not self.parent.plots.signal_info["tracer"]["plot"]
        self.parent.plots.refresh()
        self.bold_items()

    def co2_switch(self):
        self.parent.plots.signal_info["co2"]["plot"] = not self.parent.plots.signal_info["co2"]["plot"]
        self.parent.plots.refresh()
        self.bold_items()

    def o2_switch(self):
        self.parent.plots.signal_info["o2"]["plot"] = not self.parent.plots.signal_info["o2"]["plot"]
        self.parent.plots.refresh()
        self.bold_items()

    def flow_switch(self):
        self.parent.plots.signal_info["flow"]["plot"] = not self.parent.plots.signal_info["flow"]["plot"]
        self.parent.plots.refresh()
        self.bold_items()

    def mmside_switch(self):
        self.parent.plots.signal_info["mmside"]["plot"] = not self.parent.plots.signal_info["mmside"]["plot"]
        self.parent.plots.refresh()
        self.bold_items()

    def bold_items(self):
        # Configure items in signal menu
        item_list = ['tracer', 'flow', 'o2', 'co2', 'mmside']
        for item in item_list:
            if self.parent.plots.signal_info[item]['plot']:
                self.signal_menu.entryconfigure(item_list.index(item), font="SegoeUI 9 bold")
            else:
                self.signal_menu.entryconfigure(item_list.index(item), font="SegoeUI 9")

    def x_axis_switch(self):
        if self.parent.plots.x_axis_mode == 'time':
            self.parent.plots.x_axis_mode = 'volume'
            self.signal_menu.entryconfigure(6, label="Time x-axis")

        elif self.parent.plots.x_axis_mode == 'volume':
            self.parent.plots.x_axis_mode = 'time'
            self.signal_menu.entryconfigure(6, label="Volume x-axis")

        self.parent.plots.refresh()

    def endtidal_switch(self):

        self.parent.plots.plot_endtidal = not self.parent.plots.plot_endtidal

        if self.parent.plots.plot_endtidal:
            self.signal_menu.entryconfigure(7, font='SegoeUI 9 bold')
        elif self.parent.plots.plot_endtidal:
            self.signal_menu.entryconfigure(7, font='SegoeUI 9')

        self.parent.plots.refresh()


class SignalWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg='#5A5A5A')
        self.parent = parent
        self.fig = None
        self.axes = None
        self.graphs = None
        self.toolbar = None
        self.current_subplot = 0
        self.x_axis_mode = 'time'
        self.plot_endtidal = False
        """
        Plotting behavior: For each main signal, this dictionary contains:
        plot: whether the signal should be plotted
        nr: which subplot number this graph is on
        name: the LungSim-internal name of said signal
        label: the label on the y-axis
        scale_y: the scaling factor in the y-direction
        """
        self.signal_info = {
            'tracer': {
                'plot': True,
                'nr': 0,
                'name': 'tracer',
                'label': 'Tracer [%]',
                'scale_y': 100},
            'flow': {
                'plot': True,
                'nr': 0,
                'name': 'flow',
                'label': 'Flow [L/s]',
                'scale_y': 1000},
            'o2': {
                'plot': False,
                'nr': 0,
                'name': 'o2',
                'label': 'O2 [%]',
                'scale_y': 100},
            'co2': {
                'plot': True,
                'nr': 0,
                'name': 'co2',
                'label': 'CO2 [%]',
                'scale_y': 100},
            'mmside': {
                'plot': False,
                'nr': 0,
                'name': 'mmside',
                'label': 'MMss [g/mol]',
                'scale_y': 1000},
            'mmmain': {
                'plot': False,
                'nr': 0,
                'name': 'mmmain',
                'label': 'MMms [g/mol]',
                'scale_y': 1000}
        }

    def refresh(self):
        """
        Destroy previous graphs, create pyplot graph and replace it, according to the state of signal_info
        """
        for widget in self.winfo_children():
            widget.destroy()

        self.fig, self.axes = plt.subplots(self.subplot_nr, 1, sharex=True, facecolor="#5A5A5A")

        # If there is a single axis, we need to turn it into a list of length 1 for the code below to work properly
        if self.subplot_nr == 1:
            self.axes = ([self.axes])

        # Draw the measurement on a figure
        self.draw_signals()

        # Draw the end-tidal concentrations
        self.draw_endtidal()

        # Commit the graph created to the GUI
        self.insert_figure()

    def draw_endtidal(self):
        # Draw the end-tidal signal concentrations for Tracer, CO2, and O2, if set by the user
        if not self.plot_endtidal:
            return

        #
        if self.signal_info['tracer']['plot']:
            endtidal_indices = self.mbw.breaths.exp_end
            x_data = self.x_data[endtidal_indices]
            y_data = self.mbw.tracer_endtidal * self.signal_info['tracer']['scale_y']

            plot_nr = self.signal_info['tracer']['subplot_nr']

            self.axes[plot_nr].plot(x_data, y_data, 'or', markersize=2)

        if self.signal_info['co2']['plot']:
            endtidal_indices = self.mbw.breaths.exp_end
            x_data = self.x_data[endtidal_indices]
            y_data = self.mbw.co2_endtidal * self.signal_info['co2']['scale_y']

            plot_nr = self.signal_info['co2']['subplot_nr']

            self.axes[plot_nr].plot(x_data, y_data, 'or', markersize=2)

        if self.signal_info['o2']['plot']:
            endtidal_indices = self.mbw.breaths.exp_end
            x_data = self.x_data[endtidal_indices]
            y_data = self.mbw.o2_endtidal * self.signal_info['o2']['scale_y']

            plot_nr = self.signal_info['o2']['subplot_nr']

            self.axes[plot_nr].plot(x_data, y_data, 'or', markersize=2)

    @property
    def x_data(self):
        if self.x_axis_mode == 'time':
            # Time signal in [s]
            x_data = self.mbw.co2.time
        elif self.x_axis_mode == 'volume':
            # Continuous volume signal in [L]
            x_data = self.mbw.volume_abs.data * 1000

        return x_data

    def draw_signals(self):
        # Iterate through the signals and draw the graphs
        self.current_subplot = 0
        for signal in self.signal_info:
            # Unpack information for this signal
            info = self.signal_info[signal]
            # If plotting is off, skip
            if not info['plot']:
                continue

            # Add to the signal info on which subplot they are located
            self.signal_info[signal]['subplot_nr'] = self.current_subplot

            # Plot signal as a function of time
            color = 'black'
            linewidth = 0.5
            label = info['label']
            try:
                signal_object = eval('self.mbw.' + info['name'])
                # Choose the correct x-axis
                x_data = self.x_data
                y_data = signal_object.data * info['scale_y']
            except LungSimException:
                x_data = 0
                y_data = 0

            # Plot
            self.axes[self.current_subplot].plot(x_data, y_data, color=color, linewidth=linewidth, label=label)
            # Background color
            # self.axes[self.current_subplot].set_facecolor("#5A5A5A")
            # Axis Labels
            plt.axes(self.axes[self.current_subplot])
            plt.ylabel(info['label'])
            # Grid and legend
            self.axes[self.current_subplot].grid()
            self.axes[self.current_subplot].legend(loc='upper right')
            # Update the signal_info structure as to the location of this plot
            self.signal_info[signal]['nr'] = self.current_subplot
            self.current_subplot += 1

        # Label the x-axis on the last plot only to save space
        plt.axes(self.axes[-1])
        if self.x_axis_mode == 'time':
            plt.xlabel('Time [s]')
        elif self.x_axis_mode == 'volume':
            plt.xlabel('Continuous volume [L]')

    def insert_figure(self):
        self.graphs = FigureCanvasTkAgg(self.fig, master=self)
        self.graphs.draw()
        self.graphs.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        # Add the toolbar
        self.toolbar = NavigationToolbar2Tk(self.graphs, self)
        self.toolbar.update()
        self.graphs._tkcanvas.pack()

    @property
    def mbw(self):
        """
        Shortcut property for simpler access to current signals
        :return: currently activated MBW object
        """
        return self.parent.mbw_analysis.current_analysis

    @property
    def subplot_nr(self):
        """
        :return: number of "activated" plots, i.e. how many plots are supposed to be shown in the SignalWindow frame at
        this time.
        """
        nr = 0
        for signal in self.signal_info:
            nr += self.signal_info[signal]['plot']

        return nr


def run_gui():
    root = tk.Tk()
    IntroWindow(root)
    root.mainloop()


def run_test():
    root = tk.Tk()
    config = r"\\filer300\USERS3007\I0337516\Desktop\data\Config_noncutout_331_precap_24ml.yaml"
    files = [r"\\filer300\USERS3007\I0337516\Desktop\data\A-20220620-155737-WAETRI30052013-N2MultiBreathWashoutTest-Set2.txt", r"\\filer300\USERS3007\I0337516\Desktop\data\A-20220620-160354-WAETRI30052013-N2MultiBreathWashoutTest-Set2.txt"]
    results = r"\\filer300\USERS3007\I0337516\Desktop\data"
    MainWindow(root, config, files, results)
    root.mainloop()


run_test()

