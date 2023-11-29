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
font_small = 8  # global
font_medium = 10  # title like "select type of artefacts
font_tall = 12  # "Summary"


def style_config():
    # watch out: need to used ttk.Button and not Button in order for the style to be applied
    # however, push-push button hard to do with ttk.Button --> button in the class "gradeWindow" are not ttk and so the style is directly written there
    style = ttk.Style()
    style.theme_use(name_theme)

    style.map('TRadiobutton', background=[('selected', color_background)])
    style.map('B2.TRadiobutton', background=[('selected', color_background)])
    style.map('TCheckbutton', background=[('selected', color_background)])


    style.configure('TRadiobutton', font=(main_font, font_medium), background=color_background)
    style.configure('B2.TRadiobutton', font=(main_font, font_medium), background=color_background)
    style.configure('TCheckbutton', font=(main_font, font_medium), background=color_background)


    style.configure("Treeview.Heading", font=(main_font, font_small))
    style.configure("Treeview",
                    background=color_tables,
                    foreground='black',
                    rowheight=25,
                    fieldbackground=color_tables)
    style.map('Treeview', background=[('selected', color_selected_row_tabs)])

    plt.style.use('seaborn-v0_8-white')


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
        self.experience_info = None

        # Title of the window
        self.parent.title('MBW quality grading')
        self.parent.iconbitmap(resource_path('imgs\logo_insel.ico'))

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


        # Create the buttons
        # Button to load a previous session
        self.previous_session_button = tk.Button(master=self,
                                                 text='Load previous grading session',
                                                 bg=color_intro_window,
                                                 command=self.start_main_previous_session,
                                                 height=button_height,
                                                 width=button_width)

        # Button for starting app
        self.start_button = tk.Button(master=self,
                                      text='Go!',
                                      bg=green_insel,
                                      fg='white',
                                      command=self.start_experience_questions,
                                      height=button_height,
                                      width=button_width)
        # logo
        img_original = Image.open(resource_path("imgs\logo_insel_kk_big.png"))
        self.img = ImageTk.PhotoImage(img_original.resize((300, 62)))  # original is 1523 x 315
        self.logo = Label(self, image=self.img, borderwidth=0, background=color_intro_window)

        # Create the labels
        self.welcome_label = Label(self, text='Welcome to the MBW quality control app!',
                                       background=color_intro_window, font=(main_font, font_tall), relief=FLAT)

        self.previous_session_label = Label(self,
                                                text='open a previous session to modify \nor terminate the grading',
                                                background=color_intro_window, font=(main_font, font_medium), relief=FLAT)
        self.start_label = Label(self,
                                     text='Begin a new quality control session',
                                     background=color_intro_window, font=(main_font, font_medium), relief=FLAT)

        self.or_label = Label(self,
                                  text='or',
                                  background=color_intro_window, font=(main_font, font_medium), relief=FLAT)

        # Divide the window into grids (The main item is the graphs at 1,1, the rest organizes around it)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=1)


        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        self.logo.grid(row=0, column=0, columnspan=3)
        self.welcome_label.grid(row=1, column=0, columnspan=3)

        self.start_label.grid(row=2, column=0)
        self.start_button.grid(row=3, column=0)

        self.or_label.grid(row=2, column=1)
        self.previous_session_label.grid(row=2, column=2)
        self.previous_session_button.grid(row=3, column=2)

    def set_experience_info(self, infos):
        self.experience_info = infos

    def start_main_previous_session(self):
        self.prev_session_filename = fileopenbox(msg='Select files', multiple=False, filetypes="*xlsx")
        if self.prev_session_filename is None:
            self.prev_session_file_ready = False
            self.previous_session_button.configure(bg=color_intro_window)
        else:
            self.prev_session_file_ready = True
            self.previous_session_button.configure(bg=green_insel, foreground='white')
            # Removes current frame
            self.pack_forget()
            self.start_main()

    def start_experience_questions(self):
        # Removes current frame
        self.pack_forget()
        # Starts the QC run
        ExperienceQuestionsWindow(self.parent, self)

    def start_main(self):
        # Starts the QC run
        MainWindow(self.parent, self.options_config, self.filenames, self.prev_session_filename, self.experience_info)


class ExperienceQuestionsWindow(tk.Frame):
    def __init__(self, parent, intro_window):
        tk.Frame.__init__(self, parent, background=color_intro_window)

        self.parent = parent # parent is the root, not the intro_window
        self.intro_window = intro_window
        self.pack(side="top", fill="both", expand=True)

        self.parent.title('Experience level questions')
        self.parent.iconbitmap(resource_path('imgs\logo_insel.ico'))

        # Lifts the window and makes it zoomed
        self.parent.lift()
        self.parent.attributes('-topmost', True)
        self.parent.after_idle(self.parent.attributes, '-topmost', False)

        # Window resolution
        x_resolution = 800
        y_resolution = 400

        # Turn resolution parameters into a string
        self.parent.geometry(str(x_resolution) + 'x' + str(y_resolution))


        self.explanation_label = ttk.Label(self, text='Grader experience level', background=color_intro_window, font=(main_font, font_tall))

        self.ID_label = ttk.Label(self, text="ID:", background=color_intro_window, justify=tk.LEFT, font=(main_font, font_medium))
        self.ID_entry = ttk.Entry(self)

        self.yearsExperience_label = ttk.Label(self, text="How many years' experience do you have in quality control?", background=color_intro_window, justify=tk.LEFT, font=(main_font, font_medium))
        self.yearsExperience_entry = Entry(self)

        self.curvesNumber_label = ttk.Label(self, text="How many curves are you grading per week", background=color_intro_window, justify=tk.LEFT, font=(main_font, font_medium))
        self.curvesNumber_entry = ttk.Entry(self)

        self.comment_label = ttk.Label(self, text="Comment:", background=color_intro_window, justify=tk.LEFT, font=(main_font, font_medium))
        self.comment_entry = ttk.Entry(self)

        self.register_button = Button(self, text="Register", command=self.register, background=color_intro_window)

        # Divide the window into grids (The main item is the graphs at 1,1, the rest organizes around it)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=1)
        self.rowconfigure(4, weight=1)
        self.rowconfigure(5, weight=1)
        self.rowconfigure(6, weight=1)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        self.explanation_label.grid(column=0, row=0, columnspan=3)
        self.ID_label.grid(column=0, row=1, sticky='w', padx=100)
        self.ID_entry.grid(column=1, row=1)
        self.yearsExperience_label.grid(column=0, row=2, sticky='w', padx=100)
        self.yearsExperience_entry.grid(column=1, row=2)
        self.curvesNumber_label.grid(column=0, row=3, sticky='w', padx=100)
        self.curvesNumber_entry.grid(column=1, row=3)
        self.comment_label.grid(column=0, row=4, sticky='w', padx=100)
        self.comment_entry.grid(column=1, row=4)
        self.register_button.grid(column=2, row=6)

    def register(self):
        # Get the user input from the form
        infos = pd.DataFrame(columns=['ID', 'yearsExperiences', 'curvesNumber', 'comment'])
        infos.loc[-1] = [self.ID_entry.get(), self.yearsExperience_entry.get(), self.curvesNumber_entry.get(),
               self.comment_entry.get()]
        self.intro_window.set_experience_info(infos)
        self.pack_forget()
        self.intro_window.start_main()


class MainWindow(tk.Frame):
    def __init__(self, parent, options_config, filenames, prev_session_filename, experience_info):
        tk.Frame.__init__(self, parent, background=color_background)
        self.parent = parent

        # -------------------------Style
        style_config()
        self.parent.title('MBW quality grading')
        self.parent.iconbitmap(resource_path('imgs\logo_insel.ico'))

        # Fullscreen
        self.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.filenames = filenames
        self.options_config = options_config
        self.previous_session_filename = prev_session_filename
        self.experience_info = experience_info
        self.current_trial = 0
        self.nb_trials = len(self.filenames)
        self.abno_table = None

        # MBW
        self.mbw_analysis = MBWAnalysis.from_fileopenbox(
            filenames=self.filenames,
            options_config=self.options_config,
            results_path=''
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

        # self.num_vals = NumericalValuesWindow(self)
        # self.num_vals.grid(column=0, row=5, sticky='nsew')
        # self.num_vals.refresh()

        self.toolbar = PlotToolbarWindow(self)
        self.toolbar.grid(column=0, row=5)

        # Explanation
        self.explanation = ExplanationWindow(self)
        self.explanation.grid(column=0, row=6, sticky='nsew')


        # -------------------------------------------------------- Column 1
        self.labelArtefact = ttk.Label(self, text='Select artefact: \n', font=(main_font, font_medium),
                                       background=color_background,
                                       justify=tk.LEFT, borderwidth=0)
        self.labelArtefact.grid(column=1, row=0, sticky='nsew')

        #  Radio button with abnormalities choices
        self.artefact_choice = SelectArtefactsWindow(self)
        self.artefact_choice.grid(column=1, row=1, sticky='nsew')
        self.artefact_choice.refresh()

        self.labelGrade = ttk.Label(self, text='Grade the trial: \n', font=(main_font, font_medium),
                                    background=color_background,
                                    justify=tk.LEFT, borderwidth=0)
        self.labelGrade.grid(column=1, row=2, sticky='nsew')
        #  Grades
        self.grade = SelectGradeWindow(self)
        self.grade.grid(column=1, row=3, sticky='nsew')


        self.labelRemark = ttk.Label(self, text='Add a comment:\n', font=(main_font, font_medium),
                                     background=color_background, justify=tk.LEFT, borderwidth=0)
        self.labelRemark.grid(column=1, row=4, sticky='nsew')
        #  Remark
        self.remark = RemarkWindow(self)
        self.remark.grid(column=1, row=5, sticky='nsew')
        self.remark.refresh()

        self.logo = LogoWindow(self)
        self.logo.grid(column=1, row=6)

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
            self.experience_info, self.abno_table.savedAbnoTable, self.grade.savedGrades, self.remark.savedRemark = self.load_previous_session()
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
        # self.num_vals.refresh()
        self.graded_trials.refresh()
        self.plots.graphs.draw()

    def load_previous_session(self):
        prev_experience = pd.read_excel(self.previous_session_filename, sheet_name='Experience Info')
        prev_abno = pd.read_excel(self.previous_session_filename, sheet_name='Abnormalities')
        prev_grade = pd.read_excel(self.previous_session_filename, sheet_name='Grades')
        prev_remark = pd.read_excel(self.previous_session_filename, sheet_name='Remarks')

        prev_remark['remark'] = prev_remark['remark'].fillna('')

        return prev_experience, prev_abno, prev_grade, prev_remark


class LogoWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=0, relief='flat', bg=color_background)
        self.parent = parent

        img_original = Image.open(resource_path("imgs\logo_insel_kk_big.png"))
        self.img = ImageTk.PhotoImage(img_original.resize((150, 31)))  # original is 1523 x 315
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
        self.plot_idx_while_selecting = None

        self.current_subplot = 0
        self.x_axis_mode = 'time'
        self.signal_info = {
            'vol/flow': {
                'plot': True,
                'nr': 0,
                'name': ['volume', 'flow'],
                'label': 'Vol [mL], Flow [L/s]',
                'scale_y': 1000,
                'color': ['red', 'k']},
            'tracer': {
                'plot': True,
                'nr': 0,
                'name': 'tracer',
                'label': 'Tracer [%]',
                'scale_y': 100,
                'color': 'firebrick'},
            'co2': {
                'plot': True,
                'nr': 0,
                'name': 'co2',
                'label': 'CO2 [%]',
                'scale_y': 100,
                'color': 'green'}
        }

        self.plot_endtidal = None

    def refresh(self):
        """
        Destroy previous graphs, create pyplot graph and replace it, according to the state of signal_info
        """
        for widget in self.winfo_children():
            widget.destroy()
        plt.close(self.fig)

        self.fig, self.axes = plt.subplots(self.subplot_nr, 1, sharex=True, facecolor=color_background)
        self.fig.subplots_adjust(left=0.08, right=0.95, hspace=0.14)

        # If there is a single axis, we need to turn it into a list of length 1 for the code below to work properly
        if self.subplot_nr == 1:
            self.axes = ([self.axes])

        # Draw the measurement on a figure
        self.draw_signals()

        # Commit the graph created to the GUI
        self.insert_figure()

        # rectangle annotation
        self.annotation = FigureAnnotation(self)

    def disable_zoom(self):
        self.toolbar.zoom()

    def zoom_on(self):
        return self.toolbar.mode == 'zoom rect'

    def restore_home_view(self):
        self.toolbar.home()


    def set_plot_entidal(self, bool):
        self.plot_endtidal = bool


    # TODO: modify this function with signal info
    def translate_plot_name_idx(self, val):
        if type(val) == str:
            if val == 'Vol/Flow':
                plot_idx = 0
            elif val == 'Tracer':
                plot_idx = 1
            elif val == 'CO2':
                plot_idx = 2
            else:
                plot_idx = np.nan
            return plot_idx

        else:
            if val == 0 : plot_name = 'Vol/Flow'
            elif val == 1: plot_name = 'Tracer'
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
            linewidth = 0.8

            # Axis Labels
            plt.axes(self.axes[self.current_subplot])
            plt.ylabel(info['label'])
            # Grid and legend
            self.axes[self.current_subplot].grid()
            self.axes[self.current_subplot].spines[['right', 'top', 'bottom', 'left']].set_visible(False)
            self.axes[self.current_subplot].axhline(y=0, color='grey', linewidth=linewidth)
            self.axes[self.current_subplot].axvline(x=0, color='grey', linewidth=linewidth)

            # Plot signal as a function of time

            try:
                # Choose the correct x-axis
                x_data = self.x_data

                if type(info['name']) == str:
                    signal_object = eval('self.mbw.' + info['name'])
                    y_data = signal_object.data * info['scale_y']
                    self.axes[self.current_subplot].plot(x_data, y_data, color=info['color'], linewidth=linewidth, label=info['label'])

                elif type(info['name']) == list:
                    signal_object1 = eval('self.mbw.' + info['name'][0])
                    y_data1 = signal_object1.data * info['scale_y']
                    signal_object2 = eval('self.mbw.' + info['name'][1])
                    y_data2 = signal_object2.data * info['scale_y']
                    y_data = [y_data1, y_data2]
                    color = ['red', 'blue']
                    legends = info['label'].split(',')
                    for idx, y in enumerate(y_data):
                        self.axes[self.current_subplot].plot(x_data, y, color=info['color'][idx], linewidth=linewidth, label=legends[idx].strip())
                    self.axes[self.current_subplot].legend(ncols=len(y_data), bbox_to_anchor=(0, 1),loc='lower left', fontsize='small')

            except LungSimException:
                x_data = 0
                y_data = 0


            # Update the signal_info structure as to the location of this plot
            self.signal_info[signal]['nr'] = self.current_subplot
            self.current_subplot += 1

        # Label the x-axis on the last plot only to save space
        self.fig.align_ylabels(self.axes[:])
        plt.axes(self.axes[-1])
        if self.x_axis_mode == 'time':
            plt.xlabel('Time [s]')
        elif self.x_axis_mode == 'volume':
            plt.xlabel('Continuous volume [L]')

        # plt.subplots_adjust(left=0.08, right=0.95, hspace=0.14)

        self.draw_endOfTest_vline()
        self.draw_concentrationThreshold_hline()
        self.draw_endtidal()

    def draw_endtidal(self):
        if not self.plot_endtidal:
            return

        if self.signal_info['tracer']['plot']:

            endtidal_indices = self.mbw.breaths.exp_end[self.mbw.washout_breaths]
            x_data = self.x_data[endtidal_indices]
            y_data = self.mbw.tracer_endtidal * self.signal_info['tracer']['scale_y']
            y_data = y_data[self.mbw.washout_breaths]

            plot_nr = self.signal_info['tracer']['subplot_nr']

            self.axes[plot_nr].plot(x_data, y_data, 'or', markersize=2)

        if self.signal_info['co2']['plot']:
            endtidal_indices = self.mbw.breaths.exp_end[self.mbw.washout_breaths]
            x_data = self.x_data[endtidal_indices]
            y_data = self.mbw.co2_endtidal * self.signal_info['co2']['scale_y']
            y_data = y_data[self.mbw.washout_breaths]

            plot_nr = self.signal_info['co2']['subplot_nr']

            self.axes[plot_nr].plot(x_data, y_data, 'or', markersize=2)

    def draw_endOfTest_vline(self):
        plot_nr = self.signal_info['tracer']['subplot_nr']
        ymin, ymax = self.axes[plot_nr].get_ylim()
        try:
            x_idx = self.mbw.breaths.exp_end[self.mbw.washout_end(2.5)]
            self.axes[plot_nr].vlines(self.x_data[x_idx], ymin=ymin, ymax=ymax, linestyles="dashed", colors='black', linewidth=0.5)
            self.axes[plot_nr].text(self.x_data[x_idx]+1, 15, 'End of test', rotation=90)
        except LungSimException:
            x_idx = 0

    def draw_concentrationThreshold_hline(self):
        plot_nr = self.signal_info['tracer']['subplot_nr']
        xmin, xmax = self.axes[plot_nr].get_xlim()
        try:
            y = self.mbw.tracer_initial * (1/40) * self.signal_info['tracer']['scale_y']
            self.axes[plot_nr].hlines(y, xmin=xmin, xmax=xmax, linestyles="dashed", colors='black', linewidth=0.8)
        except LungSimException:
            y_idx = 0

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

        #self.graphs._tkcanvas.pack()

        self.graphs.mpl_connect('button_press_event', self.on_click)
        self.graphs.mpl_connect('button_release_event', self.on_release)

    def on_click(self, event):
        state = self.toolbar.mode
        if state == 'zoom rect':
            self.parent.artefact_choice.refresh()
            self.parent.explanation.update_label("turn_off_zoom")

        self.annotation.remove_static_rectangle()

        self.plot_idx_while_selecting = self.annotation.plot_idx
        self.annotation.while_selecting = True

    def on_release(self, event):
        self.annotation.while_selecting = False

        if self.annotation.plot_idx != self.plot_idx_while_selecting:
            self.annotation.reset_rectangles(self.plot_idx_while_selecting)

        self.plot_idx_while_selecting = None


class PlotToolbarWindow(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent, borderwidth=10, relief='flat', bg=color_background)
        self.parent = parent
        self.plot_endtidal = tk.BooleanVar()
        self.endtidal_checkbox = ttk.Checkbutton(self,
                                                 text='Show end tidal concentrations',
                                                 command=self.is_checked,
                                                 variable=self.plot_endtidal,
                                                 onvalue=True,
                                                 offvalue=False)
        self.endtidal_checkbox.pack()

    def is_checked(self):
        self.parent.plots.set_plot_entidal(self.plot_endtidal.get())
        self.parent.plots.refresh()


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
        self.label = Label(self, text=self.text, background=color_background, font=(main_font, font_medium),
                           justify=LEFT)
        self.update_label('general')

        self.label.pack(fill='x')

    def update_label(self, text2display):
        if text2display == 'general':
            self.text = '(1) Select and draw an abnormality you see on the trial\n' \
                        '(2) Grade the trial (3) Add a comment\n' \
                        'Navigate in trials using "<" and ">", the keyboard arrows or the progression table.'
        elif text2display == 'draw':
            self.text = 'Draw a rectangle on any plot'
        elif text2display == 'how2draw_save':
            self.text = 'Drag and drop the mouse to create a rectangle' \
                        '\naround the located abnormality. ' \
                        ' \nClick \"save abnormality\" to add it to the abnormality table '
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
        self.Labels = ['leak', 'trapped gas', 'irregular breathing', 'cough/ sigh/ swallow', 'EELV step change',
                       'Early termination', 'Hypoventilation/ Hyperventilation', 'Undefined']
        self.save_coord_button = None
        self.init_buttons()

    def init_buttons(self):
        for idx in range(len(self.Labels)):
            ttk.Radiobutton(self, text=self.Labels[idx], variable=self.artefact, value=idx + 1,
                            command=self.artefact_selection, style="TRadiobutton").pack(side=tk.TOP, fill='y',
                                                                                        anchor='w')

        self.save_coord_button = Button(self, text='Save abnormality', command=self.save_coord_table,
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
        self.parent.plots.restore_home_view()
        self.refresh()

    def artefact_selection(self):
        # self.parent.abno_table.stored_rectangle.reset_rectangle()
        self.parent.plots.annotation.activate_annotation()
        self.artefactLabel = self.Labels[self.artefact.get() - 1]
        self.parent.explanation.update_label('draw')
        if self.parent.plots.zoom_on():
            self.parent.plots.disable_zoom()

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
        self.while_selecting = False

    # activate rectangle only when an artefact name is selected in the radiobutton
    def activate_annotation(self):
        self.parent.graphs.mpl_connect('axes_enter_event', self.on_enter_axis)
        self.flag_annotation_activated = True

    # deactivate when abnormality saved or plot changes, before selection of new artefact typ
    def deactivate_annotation(self):
        self.flag_annotation_activated = False
        self.parent.graphs.mpl_disconnect(self.on_enter_axis)
        self.reset_rectangles()

    # update the subplot number when the mouse hovers on it
    def on_enter_axis(self, event):

        if event.inaxes == self.parent.axes[0]:
            self.plot_idx = 0
        if event.inaxes == self.parent.axes[1]:
            self.plot_idx = 1
        if event.inaxes == self.parent.axes[2]:
            self.plot_idx = 2

        self.init_selector()

    def init_selector(self):# need to reinit a RectangleSelector every time we change subplots
        if self.flag_annotation_activated and not self.while_selecting:  # if user change axis between press and release of selector rectangle, we should not create a new rectangle selection because we won't be able to erase the previous one
            self.selector_rectangle = RectangleSelector(self.parent.axes[self.plot_idx], self.selector_callback,
                                                        button=[1, 3],  # don't use middle button
                                                        minspanx=5, minspany=5,
                                                        spancoords='pixels',
                                                        interactive=False,
                                                        use_data_coordinates=True,
                                                        props=dict(color=salmon_insel)
                                                        )

    def selector_callback(self, eclick, erelease):
        """eclick and erelease are the press and release events"""
        #self.remove_static_rectangle()
        self.x1, self.y1 = eclick.xdata, eclick.ydata
        self.x2, self.y2 = erelease.xdata, erelease.ydata


        if self.selector_rectangle.artists[0].get_width() == 0:  # self.toggle_selector.artists[0] is the Rectangle object
            self.parent.parent.artefact_choice.deactivate_save_coord_button()
        else:
            self.parent.parent.artefact_choice.activate_save_coord_button()

        self.selector_rectangle.set_visible(False)
        self.draw_static_rectangle_from_coordinates(self.plot_idx)

    def reset_rectangles(self, idx=None):
        if idx is None:
            idx = self.plot_idx
        self.remove_selector_rectangle(idx)
        self.remove_static_rectangle()
        self.init_selector()

    def draw_static_rectangle_from_coordinates(self, plot_selected):
        self.static_rectangle = matplotlib.patches.Rectangle(xy=(self.x1, self.y1), width=(self.x2 - self.x1),
                                                             height=(self.y2 - self.y1), color=salmon_insel)
        self.parent.axes[plot_selected].add_patch(self.static_rectangle)
        self.static_rectangle_plot_idx = plot_selected
        self.parent.axes[plot_selected].figure.canvas.draw()

    def remove_static_rectangle(self):
        if self.static_rectangle is not None:
            self.static_rectangle.set_visible(False)
            self.parent.axes[self.static_rectangle_plot_idx].figure.canvas.draw()

    def remove_selector_rectangle(self, idx):
        if self.selector_rectangle is not None:
            # for artist in self.selector_rectangle.artists:
            #    artist.set_visible(False)
            self.selector_rectangle.set_active(False)
            self.selector_rectangle.set_visible(False)
            self.parent.axes[idx].figure.canvas.draw()


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
        self.abno_table.column("abnormality", width=50)
        self.abno_table.column("plot_selected", width=20)
        self.abno_table.column("x1")
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

        self.abno_table["displaycolumns"] = ("abnormality", "plot_selected")

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
            self.parent.plots.annotation.remove_static_rectangle()
            self.parent.artefact_choice.refresh()
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
        self.progression_tab.column("curve_index", width=10)
        self.progression_tab.column("grade", width=40)
        self.progression_tab.column("remark", width=40)

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
        self.experience_tab = None
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
        self.confirm_export_window.destroy()
        self.export_path = diropenbox(msg='Select folder where outputs are saved.')
        if self.export_path is None:
            print('No export path')
        else:
            self.abno_tab = self.parent.abno_table.savedAbnoTable
            self.grade_tab = self.parent.grade.savedGrades
            self.remark_tab = self.parent.remark.savedRemark
            self.experience_tab = self.parent.experience_info

            self.abno_tab.sort_values(by='curve_index')
            self.grade_tab.sort_values(by='curve_index')
            self.remark_tab.sort_values(by='curve_index')

            with pd.ExcelWriter((self.export_path + '\\' + self.fileName)) as writer:
                # use to_excel function and specify the sheet_name and index
                # to store the dataframe in specified sheet
                self.experience_tab.to_excel(writer, sheet_name="Experience Info", index=False)
                self.grade_tab.to_excel(writer, sheet_name="Grades", index=False)
                self.abno_tab.to_excel(writer, sheet_name="Abnormalities", index=False)
                self.remark_tab.to_excel(writer, sheet_name="Remarks", index=False)
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
    # bad solution
    # options_config_filename = r"\\filer300\USERS3007\I0337516\Desktop\data\Config_noncutout_331_precap_24ml.yaml"
    # path_files = r'L:\KKM_LuFu\OfficeData\01. Documentation\QC Project\11. Datasets\MBW_last100_DS\MBW_last100_A-Files'



    # # temporary solution
    # current_dir = diropenbox(msg='Select folder of the project')
    #
    # path_files = os.path.join(current_dir, "data\MBW_last100_A-Files")
    #
    # path_config = os.path.join(current_dir, "config")
    # options_config_filename = os.path.join(path_config, "Config_noncutout_331_precap_24ml.yaml")


    # Future solution
    current_dir = resource_path("")  # current dir

    path_config = os.path.join(current_dir, "config")
    options_config_filename = os.path.join(path_config, "Config_noncutout_331_precap_24ml.yaml")

    path_files = os.path.join(current_dir, "data\MBW_last100_A-Files")
    filenames = os.listdir(path_files)
    for idx, f in enumerate(filenames): filenames[idx] = os.path.join(path_files, f)

    root = ThemedTk()
    IntroWindow(root, options_config_filename, filenames)
    root.mainloop()


def run_test():
    root = ThemedTk()
    root.title('Spirometry QC grading App')
    root.iconbitmap(resource_path("imgs\logo_insel.ico"))
    config = r"\\filer300\USERS3007\I0337516\Desktop\data\Config_noncutout_331_precap_24ml.yaml"
    path_files = r'\\filer300\USERS3007\I0337516\Desktop\data\test2'
    files = os.listdir(path_files)
    for idx, f in enumerate(files): files[idx] = os.path.join(path_files,f)
    previous_session = None
    MainWindow(root, config, files, previous_session)
    root.mainloop()


run_gui()
# run_test()
