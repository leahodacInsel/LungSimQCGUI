from tkinter import *
from tkinter import ttk
import matplotlib.patches
import pandas as pd
import sys
import numpy as np
import os
import time
from ttkthemes import ThemedTk
from PIL import Image, ImageTk
from easygui import fileopenbox
from easygui import diropenbox
import tkinter as tk
from matplotlib import pyplot as plt
from matplotlib.widgets import RectangleSelector
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from MBWAnalysis import MBWAnalysis
from LungSimException import LungSimException


# 25.10.23, 11:23 - begin to implement new intro window on this test_lea.py version

# WARNING: no protection against loading previous_session but with different files
# add: if select artefact, no zoom anymore --> or button to activate zoom

# Set of static variables for themes, colors and font (easy to change and try different combination)
green_insel = '#009670'
pink_insel = '#e5cbd6'
grey_insel = '#677078'
salmon_insel = '#ee8f7a'
blue_insel = '#6ca5da'
light_green_insel = '#d1e2bc'
beige_insel = '#e4dbcf'
purple_insel = '#8c7482'
brown_insel = '#9b6051'
dark_blue_insel = '#4a7094'
light_grey = '#BFBFBF'

color_background = 'white'
color_summary_column = green_insel  # color background summary column
color_curve = blue_insel  # color curve plot
color_selected_row_tabs = light_grey  # selected row
color_buttons = light_green_insel  # color other buttons
color_export_button = light_green_insel  # export button
color_curve_in_rectangle = salmon_insel  # select part of the curve color
color_deactivated_button = light_green_insel  # deactivated button
color_tables = light_green_insel  # tables
color_intro_window = light_green_insel

name_theme = "adapta"
main_font = 'Helvetica 10'
font_matplotlib = 'Calibri'  # font for matplotlib
font_small = 10  # global
font_medium = 12  # title like "select type of artefacts
font_tall = 14  # "Summary"


def style_config():
    # watch out: need to used ttk.Button and not Button in order for the style to be applied
    # however, push-push button hard to do with ttk.Button --> button in the class "gradeWindow" are not ttk and so the style is directly written there
    style = ttk.Style()
    style.theme_use(name_theme)
    style.map('TRadiobutton', background=[('selected', color_background)])
    style.configure('B2.TRadiobutton', font=(main_font, font_medium), background=color_background)



    style.configure("Treeview.Heading", font=(main_font, font_small))
    style.configure("Treeview",
                    background=color_tables,
                    foreground='black',
                    rowheight=25,
                    fieldbackground=color_tables)
    style.map('Treeview', background=[('selected', color_selected_row_tabs)])


def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)


class IntroWindow(tk.Frame):

    def __init__(self, parent, options_config, filenames):
        tk.Frame.__init__(self, parent, background=color_intro_window)
        # Parent of the intro screen is root
        self.parent = parent
        self.pack(side="top", fill="both", expand=True)

        # Main information of intro window
        self.options_config = options_config
        self.filenames = filenames
        self.prev_session_filename = None
        # self.results_path = None

        # Title of the window
        self.parent.title('MBW quality grading')
        self.parent.iconbitmap(resource_path('logo_insel.ico'))

        # Lifts the window and makes it zoomed
        self.parent.lift()
        self.parent.attributes('-topmost', True)
        self.parent.after_idle(self.parent.attributes, '-topmost', False)

        # Window resolution
        x_resolution = 800
        y_resolution = 400

        # Turn resolution parameters into a string
        self.parent.geometry(str(x_resolution) + 'x' + str(y_resolution))

        # Define button states
        self.files_ready = False
        self.prev_session_file_ready = False
        self.options_config_ready = False

        # Button attributes
        button_height = 3
        button_width = 30
        button_x_padding = 5
        button_y_padding = 15

        # Create the buttons
        # Button for loading file
        # self.config_button = tk.Button(master=self,
        #                                text='Load config file',
        #                                bg=color_intro_window,
        #                                command=self.load_config,
        #                                height=button_height,
        #                                width=button_width)
        #
        # self.files_button = tk.Button(master=self,
        #                               text='Load file',
        #                               bg=color_intro_window,
        #                               command=self.load_data,
        #                               height=button_height,
        #                               width=button_width)

        # # Button for defining results directory
        # self.results_button = tk.Button(master=self,
        #                                 text='Select results directory',
        #                                 bg=color_intro_window,
        #                                 command=self.load_results_path,
        #                                 height=button_height,
        #                                 width=button_width)

        # Button to load a previous session
        self.previous_session_button = tk.Button(master=self,
                                                 text='Load previous grading session',
                                                 bg=color_intro_window,
                                                 command=self.load_previous_session,
                                                 height=button_height,
                                                 width=button_width)

        # Button for starting app
        self.start_button = tk.Button(master=self,
                                      text='Go!',
                                      bg=green_insel,
                                      fg='white',
                                      command=self.start_main,
                                      height=button_height,
                                      width=button_width - 15)
        # Disable the go button by default
        # self.start_button['state'] = 'disabled'

        # logo

        img_original = Image.open(resource_path("logo_insel_kk_big.png"))
        self.img = ImageTk.PhotoImage(img_original.resize((300, 62)))  # original is 1523 x 315
        self.logo = Label(self, image=self.img, borderwidth=0, background=color_intro_window)

        # Create the labels
        self.welcome_label = ttk.Label(self, text='Welcome to the MBW quality control app!',
                                       background=color_intro_window, font=(main_font, font_tall), relief=FLAT)
        # self.config_label = ttk.Label(self, text='Select config file.', background=color_intro_window,
        #                               font=(main_font, font_medium), justify=LEFT, relief=FLAT)
        # self.files_label = ttk.Label(self, text='Select A-files. ', background=color_intro_window,
        #                              font=(main_font, font_medium), justify=LEFT, relief=FLAT)
        # self.results_path_label = ttk.Label(self, text='Select results path. ', background=color_intro_window,
        #                                     font=(main_font, font_medium), justify=LEFT, relief=FLAT)
        self.previous_session_label = ttk.Label(self,
                                                text='Open a previous session to modify \nor terminate the grading. [OPTIONAL]',
                                                background=color_intro_window, font=(main_font, font_medium), relief=FLAT)
        self.start_label = ttk.Label(self,
                                  text='Begin a new quality control session.',
                                                background=color_intro_window, font=(main_font, font_medium), relief=FLAT)

        # Divide the window into grids (The main item is the graphs at 1,1, the rest organizes around it)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=1)
        # self.rowconfigure(4, weight=1)
        # self.rowconfigure(5, weight=1)
        # self.rowconfigure(6, weight=1)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        self.logo.grid(row=0, column=0, columnspan=3)
        self.welcome_label.grid(row=1, column=0, columnspan=3)
        # self.files_label.grid(row=2, column=0, columnspan=2)
        # self.files_button.grid(row=2, column=2)
        # self.config_label.grid(row=3, column=0, columnspan=2)
        # self.config_button.grid(row=3, column=2)
        # self.results_path_label.grid(row=4, column=0, columnspan=2)
        # self.results_button.grid(row=4, column=2)
        self.previous_session_label.grid(row=2, column=0)
        self.previous_session_button.grid(row=3, column=0)
        self.start_label.grid(row=2, column=2)
        self.start_button.grid(row=3, column=2)

    # def load_config(self):
    #     self.options_config = fileopenbox(msg='Select config file', multiple=False, filetypes="*.yaml")
    #     # Check if a file was selected by the user
    #     if self.options_config is None:
    #         self.options_config_ready = False
    #     else:
    #         self.options_config_ready = True
    #         self.config_button.configure(bg=green_insel, foreground='white')
    #     # Updates the check on whether program can be launched
    #     self.start_button_update()

    # def load_data(self):
    #     self.filenames = fileopenbox(msg='Select A-files', multiple=True, filetypes="*txt")
    #     if self.filenames is None:
    #         self.files_ready = False
    #     else:
    #         self.files_ready = True
    #         self.files_button.configure(bg=green_insel, foreground='white')
    #     # Updates the check on whether program can be launched
    #     self.start_button_update()

    # def load_results_path(self):
    #     self.results_path = diropenbox(msg='Select folder where outputs are saved.')
    #     if self.results_path is None:
    #         self.results_path_ready = False
    #         self.results_button.configure(bg='light gray')
    #     else:
    #         self.results_path_ready = True
    #         self.results_button.configure(bg='green')
    #     # Updates the check on whether program can be launched
    #     self.start_button_update()

    def load_previous_session(self):
        self.prev_session_filename = fileopenbox(msg='Select files', multiple=False, filetypes="*xlsx")
        if self.prev_session_filename is None:
            self.prev_session_file_ready = False
            self.previous_session_button.configure(bg=color_intro_window)
        else:
            self.prev_session_file_ready = True
            self.previous_session_button.configure(bg=green_insel, foreground='white')

    def start_button_update(self):
        if self.files_ready:
            self.start_button['state'] = 'normal'
        else:
            self.start_button['state'] = 'disabled'

    def start_main(self):
        # Removes current frame
        self.pack_forget()
        # Starts the QC run
        MainWindow(self.parent, self.options_config, self.filenames, self.prev_session_filename)


class MainWindow(tk.Frame):
    def __init__(self, parent, options_config, filenames, prev_session_filename, results_path=""):
        tk.Frame.__init__(self, parent, background=color_background)
        self.parent = parent

        # -------------------------Style
        style_config()
        self.parent.title('MBW quality grading')

        # self.parent.iconbitmap(resource_path('logo_insel.ico'))

        # Fullscreen
        self.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.filenames = filenames
        self.options_config = options_config
        self.results_path = results_path
        self.previous_session_filename = prev_session_filename
        self.current_trial = 0
        self.nb_trials = len(self.filenames)
        self.abno_table = None

        # MBW
        self.mbw_analysis = MBWAnalysis.from_fileopenbox(
            filenames=self.filenames,
            options_config=self.options_config,
            results_path=self.results_path
        )
        self.parent.current_index = self.current_trial  # way to link implementation of app (spirometry) and MBW existing variables)
        self.mbw_analysis.load_single()

        # Divide the window into grids (The main item is the graphs at 1,1, the rest organizes around it)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=1)
        self.rowconfigure(4, weight=1)
        self.rowconfigure(5, weight=1)
        self.rowconfigure(6, weight=1)

        self.columnconfigure(0, weight=5)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        # Assemble elements of main window
        # -------------------------------------------------------- Column 0

        #  Signal window
        self.plots = SignalWindow(self)
        self.plots.grid(column=0, row=0, sticky='nsew', rowspan=5)
        self.plots.refresh()

        self.num_vals = NumericalValuesWindow(self)
        self.num_vals.grid(column=0, row=5, sticky='nsew')
        self.num_vals.refresh()

        # Explanation
        self.explanation = ExplanationWindow(self)
        self.explanation.grid(column=0, row=6, sticky='nsew')


        # -------------------------------------------------------- Column 1
        self.labelArtefact = ttk.Label(self, text='Select artefact: \n', font=(main_font, font_medium),
                                       background=color_background,
                                       justify=tk.LEFT)
        self.labelArtefact.grid(column=1, row=0, sticky='sw')

        #  Radio button with abnormalities choices
        self.artefact_choice = SelectArtefactsWindow(self)
        self.artefact_choice.grid(column=1, row=1, sticky='nw')
        self.artefact_choice.refresh()

        self.labelGrade = ttk.Label(self, text='Grade the trial: \n', font=(main_font, font_medium),
                                    background=color_background,
                                    justify=tk.LEFT)
        self.labelGrade.grid(column=1, row=2, sticky='sw')
        #  Grades
        self.grade = SelectGradeWindow(self)
        self.grade.grid(column=1, row=3, sticky='nw')


        self.labelRemark = ttk.Label(self, text='Add a comment:\n', font=(main_font, font_medium),
                                     background=color_background, justify=tk.LEFT)
        self.labelRemark.grid(column=1, row=4, sticky='sw')
        #  Remark
        self.remark = RemarkWindow(self)
        self.remark.grid(column=1, row=5, sticky='nw')
        self.remark.refresh()

        self.logo = LogoWindow(self)
        self.logo.grid(column=1, row=6, sticky='se')

        # -------------------------------------------------------- Column 2
        #  Recorded abnormalities
        self.abno_table = AbnormalityTableWindow(self)
        self.abno_table.grid(column=2, row=0, sticky='nsew', rowspan=2)
        self.abno_table.refresh()

        # Recorded graded trial progression
        self.graded_trials = TrialsProgressionWindow(self)
        self.graded_trials.grid(column=2, row=2, sticky='nsew', rowspan=3)

        # File navigation
        self.file_navigation = FileNavigation(self)
        self.file_navigation.grid(column=2, row=5, sticky='nsew')

        # Export button
        self.export = ExportWindow(self)
        self.export.grid(column=2, row=6, sticky='nsew')

        # --------------------------------------------------------
        # Load previous data
        if not self.previous_session_filename is None:
            self.abno_table.savedAbnoTable, self.grade.savedGrades, self.remark.savedRemark = self.load_previous_session()
            self.abno_table.refresh()
            self.graded_trials.refresh()

        self.bindings()

    def bindings(self):
        # bind arrows
        self.parent.bind('<Left>', lambda event: self.file_navigation.previous_file())
        self.parent.bind('<Right>', lambda event: self.file_navigation.next_file())

    def refresh_all(self):
        """
        Function to refresh all components of the UI (upon file change)
        """
        self.mbw_analysis.current_index = self.current_trial
        self.mbw_analysis.load_single()
        self.file_navigation.refresh()
        self.plots.refresh()
        self.artefact_choice.refresh()
        self.grade.refresh()
        self.remark.refresh()
        self.abno_table.refresh()
        self.explanation.refresh()
        self.num_vals.refresh()
        self.graded_trials.refresh()


    def load_previous_session(self):
        prev_abno = pd.read_excel(self.previous_session_filename, sheet_name='Abnormalities')
        prev_grade = pd.read_excel(self.previous_session_filename, sheet_name='Grades')
        prev_remark = pd.read_excel(self.previous_session_filename, sheet_name='Remarks')

        prev_remark['remark'] = prev_remark['remark'].fillna('')

        return prev_abno, prev_grade, prev_remark


class LogoWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent

        img_original = Image.open(resource_path("logo_insel_kk_big.png"))
        self.img = ImageTk.PhotoImage(img_original.resize((180, 37)))  # original is 1523 x 315
        self.logo = Label(self, image=self.img, borderwidth=0, background=color_background)
        self.logo.pack(fill=BOTH)


class SignalWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent
        self.fig = None
        self.axes = None
        self.graphs = None
        self.toolbar = None
        self.annotation = None

        self.b_save_abno = None

        self.current_subplot = 0
        self.x_axis_mode = 'time'
        self.plot_endtidal = False
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

        self.drawn_signals = {
            'tracer': [],
            'flow': [],
            'co2': [],
            'x_data': []
        }

    def refresh(self):
        """
        Destroy previous graphs, create pyplot graph and replace it, according to the state of signal_info
        """
        for widget in self.winfo_children():
            widget.destroy()
        plt.close(self.fig)

        self.fig, self.axes = plt.subplots(self.subplot_nr, 1, sharex=True, facecolor=color_background)

        # If there is a single axis, we need to turn it into a list of length 1 for the code below to work properly
        if self.subplot_nr == 1:
            self.axes = ([self.axes])

        # update vector and values attribute (tracer, flow, co2, time)
        # self.update_vals()

        # Draw the measurement on a figure
        self.draw_signals()

        # Draw the end-tidal concentrations
        self.draw_endtidal()

        # Commit the graph created to the GUI
        self.insert_figure()

        # rectangle annotation
        self.annotation = FigureAnnotation(self)



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


    def translate_plot_name_idx(self, val):
        if type(val) == str:
            if val == 'Tracer':
                plot_idx = 0
            elif val == 'Flow':
                plot_idx = 1
            elif val == 'CO2':
                plot_idx = 2
            else:
                plot_idx = np.nan
            return plot_idx

        else:
            if val == 0 : plot_name = 'Tracer'
            elif val == 1: plot_name = 'Flow'
            elif val == 2: plot_name = 'CO2'
            else:
                plot_name = np.nan
            return plot_name



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

                self.drawn_signals['x_data'] = x_data
                self.drawn_signals[signal] = y_data

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
            #self.axes[self.current_subplot].legend(loc='upper right')
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
        self.graphs.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        # Add the toolbar
        self.toolbar = NavigationToolbar2Tk(self.graphs, self)
        self.toolbar.update()
        self.toolbar.config(background=color_background)
        self.toolbar.config(highlightcolor=color_background)

        self.toolbar.winfo_children()[-2].config(background=color_background)
        self.toolbar.winfo_children()[-1].config(background=color_background)

        self.graphs._tkcanvas.pack()

        self.graphs.mpl_connect('button_press_event', self.on_click)


    def on_click(self, event):
        state = self.toolbar.mode
        print(state)
        if state == 'zoom rect':
            self.parent.explanation.update_label("turn_off_zoom")







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


class NumericalValuesWindow(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent
        self.LCI = None
        self.label = None

    def refresh(self):
        for widget in self.winfo_children():
            widget.destroy()

        self.LCI = self.calculate_LCI()

        # self.label.destroy()
        self.label = Label(self, text='', bg='White', font=(main_font, font_small), background=color_background,
                           justify=tk.LEFT)
        self.label.config(text='LCI:  ' + "{:.2f}".format(self.LCI) + ' [TO]')
        self.label.pack()

    def calculate_LCI(self):
        # calculate LCI here
        # Specify output (here LCI)
        args = dict()
        args['variable'] = 'lci_ao'
        args['method'] = 'critical'
        args['concentration'] = 2.5
        lci = self.parent.mbw_analysis.current_analysis.calculate_output(args)

        return lci


class ExplanationWindow(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent
        self.text = None
        self.label = None
        self.init_label()

    def init_label(self):
        self.text = 'Grade the quality of the displayed trial and/or select ' \
                    'an abnormality that you see appearing on the curve. ' \
                    '\nYou can navigate through the \"Previous\" and \"Next\" ' \
                    'button or in the Progression table.'
        self.label = Label(self, text=self.text, background=color_background, font=(main_font, font_small),
                           justify=LEFT)
        self.label.pack(fill='x')

    def update_label(self, text2display):
        if text2display == 'general':
            self.text = 'Grade the quality of the displayed trial and/or select ' \
                        'an abnormality that you see appearing on the curve. ' \
                        '\nYou can navigate through the \"Previous\" and \"Next\" ' \
                        'button or in the Progression table.'
        elif text2display == 'draw':
            self.text = 'Draw a rectangle on the flow-volume plot'
        elif text2display == 'how2draw_save':
            self.text = 'Drag and drop the mouse to create a rectangle around ' \
                        'the located abnormality.\nOnce it is done, click the ' \
                        '\"Save selected abnormality button\".'
        elif text2display == 'save_comm':
            self.text = 'Click on the \"Save\" button to save your comment.'

        elif text2display == 'turn_off_zoom':
            self.text = '"Zoom" button has to be deactivated to be able to select an artefact.'
        self.label.config(text=self.text)

    def refresh(self):
        self.update_label('general')


class SelectArtefactsWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent
        self.artefact = IntVar()
        self.artefactLabel = 'None'
        self.Labels = ['Cough', 'Suboptimal PEF', 'Irregularity', 'Glottic Closure', 'Obstruction', 'Leak',
                       'Early termination', 'Undefined']

        self.save_coord_button = None

        self.init_buttons()

    def init_buttons(self):
        for idx in range(len(self.Labels)):
            ttk.Radiobutton(self, text=self.Labels[idx], variable=self.artefact, value=idx + 1,
                            command=self.artefact_selection, style="TRadiobutton").pack(side=tk.TOP, fill='y',
                                                                                        anchor='w')

        self.save_coord_button = Button(self, text='Save selected abnormality', command=self.save_coord_table,
                                        background=color_buttons, relief=FLAT)
        self.save_coord_button.pack(side=RIGHT)

    def save_coord_table(self):

        # add entry to df with all info
        d = [self.parent.current_trial,
             self.parent.artefact_choice.artefactLabel,
             self.parent.plots.translate_plot_name_idx(self.parent.plots.annotation.static_rectangle_plot_idx),
             self.parent.plots.annotation.x1,
             self.parent.plots.annotation.y1,
             self.parent.plots.annotation.x2,
             self.parent.plots.annotation.y2,
             self.parent.filenames[self.parent.current_trial]
             ]

        self.parent.abno_table.savedAbnoTable.loc[len(self.parent.abno_table.savedAbnoTable)] = d
        self.parent.abno_table.refresh()
        self.parent.artefact_choice.deactivate_save_coord_button()
        self.parent.plots.annotation.deactivate_annotation()
        self.parent.explanation.refresh()
        self.refresh()

    def artefact_selection(self):
        # self.parent.abno_table.stored_rectangle.reset_rectangle()
        self.parent.plots.annotation.activate_annotation()
        self.artefactLabel = self.Labels[self.artefact.get() - 1]
        self.parent.explanation.update_label('draw')

    def refresh(self):
        self.artefact.set(0)
        self.artefactLabel = 'None'

    def activate_save_coord_button(self):
        if not self.save_coord_button is None:
            self.save_coord_button['state'] = 'normal'
            self.save_coord_button.configure(bg=color_buttons)
            self.parent.explanation.update_label('how2draw_save')

    def deactivate_save_coord_button(self):
        if not self.save_coord_button is None:
            self.save_coord_button['state'] = 'disabled'
            self.save_coord_button.configure(bg=color_deactivated_button)


class FigureAnnotation(object):  # TODO: adapt the self.plot_idx depending on which plot is hovered on

    def __init__(self, parent):
        self.parent = parent
        self.flag_annotation_activated = False
        self.x1 = None
        self.y1 = None
        self.x2 = None
        self.y2 = None
        self.plot_idx = 0  # which plot is clicked on
        self.selector_rectangle = None
        self.static_rectangle = None
        self.static_rectangle_plot_idx = None

    # activate rectangle only when an artefact name is selected in the radiobutton
    def activate_annotation(self):
        self.parent.graphs.mpl_connect('axes_enter_event', self.on_enter_axis)
        self.flag_annotation_activated = True

    # deactivate when abnormality saved or plot changes, before selection of new artefact typ
    def deactivate_annotation(self):
        self.parent.graphs.mpl_disconnect(self.on_enter_axis)
        self.reset_rectangles()
        self.flag_annotation_activated = False

    # update the subplot number when the mouse hovers on it
    def on_enter_axis(self, event):

        if event.inaxes == self.parent.axes[0]:
            self.plot_idx = 0
        if event.inaxes == self.parent.axes[1]:
            self.plot_idx = 1
        if event.inaxes == self.parent.axes[2]:
            self.plot_idx = 2

        self.init_selector()



    def init_selector(self):  # need to reinit a RectangleSelector every time we change subplots
        if self.flag_annotation_activated:
            self.selector_rectangle = RectangleSelector(self.parent.axes[self.plot_idx], self.selector_callback,
                                                        button=[1, 3],  # don't use middle button
                                                        minspanx=5, minspany=5,
                                                        spancoords='pixels',
                                                        interactive=True,
                                                        use_data_coordinates=True,
                                                        props=dict(color=salmon_insel)
                                                        )



    def selector_callback(self, eclick, erelease):
        """eclick and erelease are the press and release events"""
        self.remove_static_rectangle()
        self.x1, self.y1 = eclick.xdata, eclick.ydata
        self.x2, self.y2 = erelease.xdata, erelease.ydata


        if self.selector_rectangle.artists[0].get_width() == 0:  # self.toggle_selector.artists[0] is the Rectangle object
            self.parent.parent.artefact_choice.deactivate_save_coord_button()
        else:
            self.parent.parent.artefact_choice.activate_save_coord_button()

        self.selector_rectangle.set_visible(False)
        self.draw_static_rectangle_from_coordinates(self.plot_idx)

    def reset_rectangles(self):
        if self.selector_rectangle is not None:
            for artist in self.selector_rectangle.artists:
                artist.set_visible(False)
            self.selector_rectangle.set_active(False)
            self.selector_rectangle.update()
        self.remove_static_rectangle()

    def draw_static_rectangle_from_coordinates(self, plot_selected):
        self.static_rectangle = matplotlib.patches.Rectangle(xy=(self.x1, self.y1), width=(self.x2 - self.x1),
                                                             height=(self.y2 - self.y1), color=salmon_insel)
        self.parent.axes[plot_selected].add_patch(self.static_rectangle)
        self.static_rectangle_plot_idx = plot_selected
        self.parent.axes[plot_selected].figure.canvas.draw()


    def remove_static_rectangle(self):
        if self.static_rectangle is not None:
            self.static_rectangle.set_visible(False)
            # self.stored_rectangle = None
            self.parent.axes[self.static_rectangle_plot_idx].figure.canvas.draw()


class SelectGradeWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent
        self.grade = IntVar()
        self.gradeLabel = ''
        self.Labels = ['Unacceptable', 'Unclear', 'Acceptable']
        self.savedGrades = pd.DataFrame(columns=['curve_index', 'grade', 'gradeLabel', 'filename'],
                                        index=np.arange(self.parent.nb_trials))
        self.savedGrades['grade'] = 0
        self.savedGrades['curve_index'] = np.arange(self.parent.nb_trials)

        self.init_buttons()

    def init_buttons(self):

        for idx in range(len(self.Labels)):
            ttk.Radiobutton(self, text=self.Labels[idx], variable=self.grade, value=idx + 1,
                            command=self.grade_selection, style='B2.TRadiobutton').pack(anchor='w')



    def grade_selection(self):

        if self.grade.get() == self.savedGrades.at[self.parent.current_trial, 'grade']:
            self.grade.set(0)
            self.gradeLabel = ''
        else:
            self.gradeLabel = self.Labels[self.grade.get() - 1]

        self.savedGrades.at[self.parent.current_trial, 'curve_index'] = self.parent.current_trial
        self.savedGrades.at[self.parent.current_trial, 'gradeLabel'] = self.gradeLabel
        self.savedGrades.at[self.parent.current_trial, 'grade'] = self.grade.get()
        self.savedGrades.at[self.parent.current_trial, 'filename'] = self.parent.filenames[self.parent.current_trial]

        self.parent.graded_trials.refresh()

    def refresh(self):
        self.grade.set(self.savedGrades.at[self.parent.current_trial, 'grade'])
        self.gradeLabel = self.savedGrades.at[self.parent.current_trial, 'gradeLabel']
        self.grade.set(0)

        for idx, row in self.savedGrades.iterrows():
            if row['curve_index'] == self.parent.current_trial:
                self.grade.set(row['grade'])
                self.gradeLabel = row['gradeLabel']


class RemarkWindow(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent
        self.remark = None

        self.entry_var = StringVar()
        self.entry_var.trace("w", self.on_entry_change)
        self.entry = Entry(self, textvariable=self.entry_var)
        self.entry.config(fg='blue')
        self.entry.pack(side=tk.TOP, fill='x')

        self.clearButton = tk.Button(self, text='Clear', command=self.clear_remark, width=6, background=color_buttons,
                                     relief=FLAT)
        self.clearButton.pack(side=tk.RIGHT)

        self.saveButton = tk.Button(self, text="Save", command=self.saveRemark, width=6, background=color_buttons,
                                    relief=FLAT)
        self.saveButton.pack(side=tk.LEFT)

        self.savedRemark = pd.DataFrame(columns=['curve_index', 'remark', 'filename'], index=np.arange(self.parent.nb_trials))
        self.savedRemark['remark'] = ''
        self.savedRemark['curve_index'] = np.arange(self.parent.nb_trials)

    def saveRemark(self):
        self.remark = self.entry.get()
        # empty and disable saved button
        self.savedRemark.at[
            self.parent.current_trial, 'curve_index'] = self.parent.current_trial
        self.savedRemark.at[self.parent.current_trial, 'remark'] = self.remark
        self.savedRemark.at[self.parent.current_trial, 'filename'] = self.parent.filenames[self.parent.current_trial]
        self.entry.config(fg=green_insel)
        self.parent.graded_trials.refresh()
        self.parent.explanation.update_label('general')

    def clear_remark(self):
        self.entry.delete(0, END)
        self.remark = ''
        self.savedRemark.at[self.parent.current_trial, 'remark'] = ''
        self.entry.config(fg='black')
        self.parent.graded_trials.refresh()

    def on_entry_change(self, *args):
        if self.entry_var.get() == self.savedRemark.at[self.parent.current_trial, 'remark']:
            self.entry.config(fg=green_insel)

        else:
            self.entry.config(fg='black')

        self.parent.explanation.update_label('save_comm')

    def refresh(self):

        self.entry.delete(0, END)
        self.remark = self.savedRemark.at[self.parent.current_trial, 'remark']

        if not self.savedRemark.at[self.parent.current_trial, 'remark'] == '':
            self.entry.insert(0, self.savedRemark.at[self.parent.current_trial, 'remark'])
            self.entry.config(fg=green_insel)

        else:
            self.entry.config(fg='black')


class AbnormalityTableWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_summary_column)
        self.parent = parent

        self.abno_table = None
        self.init_abno_table()

        self.selected_row = None
        self.selection = []
        self.delete_row_button = None

        self.savedAbnoTable = pd.DataFrame(
            columns=['curve_index', 'abnormality', 'plot_selected', 'x1', 'y1', 'x2', 'y2', 'filename'])
        self.init_buttons()

    def init_abno_table(self):

        self.abno_table = ttk.Treeview(self, height=5, style='Treeview')
        self.abno_table['columns'] = ('curve_index', 'abnormality', 'plot_selected', 'x1', 'y1', 'x2', 'y2')

        # column
        self.abno_table.column("#0", width=0, stretch=NO)
        self.abno_table.column("curve_index")
        self.abno_table.column("abnormality", width=80)
        self.abno_table.column("plot_selected", width=40)
        self.abno_table.column("x1", width=20)
        self.abno_table.column("y1")
        self.abno_table.column("x2")
        self.abno_table.column("y2")
        # title

        self.abno_table.heading("#0", text="")
        self.abno_table.heading("curve_index", text="n°")
        self.abno_table.heading("abnormality", text="Abnormalites # " + str(self.parent.current_trial))
        self.abno_table.heading("plot_selected", text="Plot")
        self.abno_table.heading("x1", text="x1")
        self.abno_table.heading("y1", text="y1")
        self.abno_table.heading("x2", text="x2")
        self.abno_table.heading("y2", text="y2")

        self.abno_table["displaycolumns"] = ("abnormality", "plot_selected", 'x1')

        self.abno_table.pack(fill=X)
        self.abno_table.bind('<ButtonRelease>', self.on_click_release)
        self.abno_table.bind('<Button>', self.on_click)

        return self.abno_table

    def refresh(self):
        for widget in self.winfo_children():
            widget.destroy()

        self.init_abno_table()
        self.init_buttons()

        idxs = []
        for idx, row in self.savedAbnoTable.iterrows():
            if row['curve_index'] == self.parent.current_trial:
                idxs.append(idx)
        for i in idxs:
            vals = self.savedAbnoTable.loc[
                i, ['curve_index', 'abnormality', 'plot_selected', 'x1', 'y1', 'x2', 'y2', 'filename']].tolist()
            self.abno_table.insert(parent='', index='end', text='', values=vals)

        # self.deactivate_save_coord_button()
        self.parent.artefact_choice.deactivate_save_coord_button()

        self.selected_row = None

        if self.selected_row is None:
            self.delete_row_button['state'] = 'disabled'
            self.delete_row_button.configure(bg=color_deactivated_button)

    def init_buttons(self):

        self.delete_row_button = Button(self, text="Delete", command=self.delete_row, background=color_buttons,
                                        relief=FLAT)
        self.delete_row_button.pack(side=RIGHT)

    def on_click(self, event):
        if len(self.abno_table.selection()) > 0:
            item = self.abno_table.selection()[0]
            self.abno_table.selection_remove(item)
            return "break"

    def on_click_release(self, event):
        self.selected_row = self.abno_table.focus()
        if len(self.abno_table.selection()) > 0:
            self.draw_annotation()
            self.delete_row_button['state'] = 'normal'  # button do not activate
            self.delete_row_button.configure(bg=color_buttons)

        else:
            self.delete_row_button['state'] = 'disabled'
            self.delete_row_button.configure(bg=color_deactivated_button)
            self.parent.plots.annotation.remove_static_rectangle()

    def draw_annotation(self):
        vals = self.abno_table.item(self.selected_row)['values']
        if len(vals) > 0:
            for idx, row in self.savedAbnoTable.iterrows():
                if row['curve_index'] == vals[0] and row['abnormality'] == vals[1] and row['x1'] == float(vals[3]):
                    self.parent.plots.annotation.x1 = row['x1']
                    self.parent.plots.annotation.y1 = row['y1']
                    self.parent.plots.annotation.x2 = row['x2']
                    self.parent.plots.annotation.y2 = row['y2']
                    plot_idx = self.parent.plots.translate_plot_name_idx(row['plot_selected'])
                    self.parent.plots.annotation.draw_static_rectangle_from_coordinates(plot_idx)
                    break

    def delete_row(self):
        if self.selected_row is not None:
            # remove entry from global saved data
            vals = self.abno_table.item(self.selected_row)['values']

            if len(vals) > 0:
                for idx, row in self.savedAbnoTable.iterrows():
                    if row['curve_index'] == vals[0] and row['abnormality'] == vals[1] and row['x1'] == float(vals[3]):
                        self.savedAbnoTable = self.savedAbnoTable.drop([idx])
                        self.savedAbnoTable.reset_index(drop=True, inplace=True)
                        break

        # remove entry from widget
        if len(self.abno_table.selection()) > 0:
            selected_item = self.abno_table.selection()[0]
            self.parent.plots.annotation.reset_rectangles()
            self.abno_table.delete(selected_item)
            self.delete_row_button['state'] = 'disabled'


class TrialsProgressionWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_summary_column)

        self.parent = parent
        self.progression_tab = None
        self.selected_row = None
        self.init_table()
        self.refresh()

    def init_table(self):
        self.progression_tab = ttk.Treeview(self, height=5, style='Treeview')

        self.progression_tab['columns'] = ('curve_index', 'grade', 'remark')

        # width of last column of mainWindowGrid is defined with this table's column widths
        self.progression_tab.column("#0", width=0, stretch=NO)
        self.progression_tab.column("curve_index", width=30)
        self.progression_tab.column("grade", width=80)
        self.progression_tab.column("remark", width=80)

        self.progression_tab.heading("#0", text="")
        self.progression_tab.heading("curve_index", text="n°")
        self.progression_tab.heading("grade", text="Grade")
        self.progression_tab.heading("remark", text="Remark")

        for i in range(0, self.parent.nb_trials):
            vals = [str(i), '']
            self.progression_tab.insert(parent='', index='end', text='', values=vals, iid=str(i))

        # Constructing vertical scrollbar with treeview
        vertscrlbar = ttk.Scrollbar(self, orient="vertical", command=self.progression_tab.yview)
        self.progression_tab.configure(yscrollcommand=vertscrlbar.set)

        vertscrlbar.pack(side=RIGHT, fill=Y)
        # vertscrlbar.place(x=200, y=50, height=80)
        self.progression_tab.pack(fill=BOTH)

        self.progression_tab.bind('<ButtonRelease>', self.on_click_release)

    def refresh(self):
        for i in range(0, self.parent.nb_trials ):
            gra = self.parent.grade.savedGrades.at[i, 'gradeLabel']
            rem = self.parent.remark.savedRemark.at[i, 'remark']
            if type(gra) == str:
                self.progression_tab.set(i, 'grade', gra)
            if type(rem) == str:
                self.progression_tab.set(i, 'remark', rem)

    def on_click_release(self, event):
        self.selected_row = self.progression_tab.focus()
        if len(self.progression_tab.selection()) > 0:
            self.parent.current_trial = int(self.selected_row)
        self.parent.refresh_all()
        # self.parent.grade.refresh()

    def refresh_highlight_row(self):
        self.progression_tab.selection_set(self.parent.current_trial)
        self.progression_tab.yview_moveto(self.parent.current_trial / len(self.progression_tab.get_children()))


class FileNavigation(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_summary_column)
        self.parent = parent

        # Buttons
        self.next_button = None
        self.trial_label = None
        self.previous_button = None

        self.init_buttons()

    def init_buttons(self):
        # Next button
        self.next_button = tk.Button(master=self, text='>', bg=color_buttons, command=self.next_file, width=2,
                                     relief=FLAT)
        self.next_button.pack(side=tk.RIGHT, fill=tk.X)

        # trial_label
        self.trial_label = Label(self, text=str(self.parent.current_trial) + "/" + str(self.parent.nb_trials - 1),
                                 bg=color_summary_column)
        self.trial_label.pack(side=tk.RIGHT, fill=tk.X)

        # Previous button
        self.previous_button = tk.Button(master=self, text='<', bg=color_buttons, command=self.previous_file, width=2,
                                         relief=FLAT)
        self.previous_button.pack(side=tk.RIGHT, fill=tk.X)

    def update_label(self):
        self.trial_label.config(text=str(self.parent.current_trial) + "/" + str(self.parent.nb_trials - 1))

    def previous_file(self):
        self.parent.current_trial = max(0, self.parent.current_trial - 1)
        self.parent.mbw_analysis.current_index = self.parent.current_trial
        self.parent.refresh_all()

    def next_file(self):
        max_index = len(self.parent.filenames)
        self.parent.current_trial = min(max_index - 1, self.parent.current_trial + 1)
        self.parent.mbw_analysis.current_index = self.parent.current_trial
        self.parent.refresh_all()

    def refresh(self):
        self.update_label()
        # update numerical values
        # self.parent.num_vals.refresh()
        self.parent.graded_trials.refresh_highlight_row()


class ExportWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_summary_column)

        self.parent = parent
        self.abno_tab = None
        self.grade_tab = None
        self.remark_tab = None
        self.fileName = 'res_' + time.strftime("%Y%m%d-%H%M%S") + '.xlsx'
        self.export_path = None
        self.confirm_export_window = None
        self.export_button = tk.Button(master=self, text='EXPORT', bg=color_export_button, command=self.confirm_export,
                                       width=15, height=3, relief=FLAT)
        self.export_button.pack(side=RIGHT)

    def confirm_export(self):
        self.confirm_export_window = Toplevel(self)
        self.confirm_export_window.geometry("500x150")
        self.confirm_export_window.title("Confirm export")
        Label(self.confirm_export_window,
              text="\n\nDo you want to export " + self.fileName + " ?").pack()

        B1 = Button(self.confirm_export_window, text="Yes", command=self.export, width=8)
        B1.place(x=160, y=80)

        B2 = Button(self.confirm_export_window, text="No", command=self.confirm_export_window.destroy, width=8)
        B2.place(x=260, y=80)

    def export(self):
        self.export_path = diropenbox(msg='Select folder where outputs are saved.')
        if self.export_path is None:
            print('No export path')
        else:
            self.abno_tab = self.parent.abno_table.savedAbnoTable
            self.grade_tab = self.parent.grade.savedGrades
            self.remark_tab = self.parent.remark.savedRemark

            self.abno_tab.sort_values(by='curve_index')
            self.grade_tab.sort_values(by='curve_index')
            self.remark_tab.sort_values(by='curve_index')

            with pd.ExcelWriter((self.export_path + '\\' + self.fileName)) as writer:
                # use to_excel function and specify the sheet_name and index
                # to store the dataframe in specified sheet
                self.grade_tab.to_excel(writer, sheet_name="Grades", index=False)
                self.abno_tab.to_excel(writer, sheet_name="Abnormalities", index=False)
                self.remark_tab.to_excel(writer, sheet_name="Remarks", index=False)
            self.confirm_export_window.destroy()
            self.success_window()

    def success_window(self):
        success_export_window = Toplevel(self)
        success_export_window.geometry("500x150")
        success_export_window.title("Success export")
        Label(success_export_window,
              text="\n\nExport of result file " + self.fileName + "\n in the folder " + self.export_path + "\nwas a success !").pack()
        B1 = Button(success_export_window, text="Ok", command=success_export_window.destroy)
        B1.pack(pady=20)


def run_gui():

    options_config = r"\\filer300\USERS3007\I0337516\Desktop\data\Config_noncutout_331_precap_24ml.yaml"
    path_files = r'\\filer300\USERS3007\I0337516\Desktop\data\test2'
    filenames = os.listdir(path_files)
    for idx, f in enumerate(filenames): filenames[idx] = os.path.join(path_files, f)

    root = ThemedTk()
    IntroWindow(root, options_config, filenames)
    root.mainloop()


def run_test():
    root = ThemedTk()
    root.title('Spirometry QC grading App')
    root.iconbitmap(resource_path("logo_insel.ico"))
    config = r"\\filer300\USERS3007\I0337516\Desktop\data\Config_noncutout_331_precap_24ml.yaml"
    path_files = r'\\filer300\USERS3007\I0337516\Desktop\data\test2'
    files = os.listdir(path_files)
    for idx, f in enumerate(files): files[idx] = os.path.join(path_files,f)
    previous_session = None
    MainWindow(root, config, files, previous_session)
    root.mainloop()


run_gui()
#run_test()
