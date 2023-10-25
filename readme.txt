Author: Florian Wyler
Date: 08.07.2021

Introduction:
This is the main script to launch LungSim, a small software to analyse multiple-breath washout (MBW) data. Due to the
variety of methods and analysis algorithms that exist to analyse an equally varied array of raw data, LungSim will not
be doing one thing, but exists / will exist as:
1.  A code base for analysing MBW data via single scripts
2.  (eventually...) A UI that aims to collect the most important aspects of MBW data (the most standard forms of
    analysis and quality control) in an easy-to-access/use user interface.

In order to preserve a certain degree of simplicity that is necessitated by the small scope of the project, many
features which will be accessible to the experienced (code-literate) user will be usable via specialized scripts.
The most important features, such as visualization of signals, access to quality control, and customization of outputs,
will be controlled via a user interface.

Design principles:
1.  Things should not grow in complexity beyond what is strictly necessary. If something can be simple, it should be.
    For example, the most standard way to analyse Ecomedics files should not be bundled with other features beyond what
    is necessary, or require much user input. There should be an easy-to-use way to manage the outputs desired, such as
    via pre-made .yaml files.
    This may for example involve one script for basic analysis, another script for quality control with/without
    metadata. In these examples, a UI is not necessary for analysis, but it is valuable for quality control, and so its
    design should reflect this.
2.  As much as possible, the inputs / outputs should not change much, so standardized configurations for analysis should
    cover a wide range of functionality for the user.
3.  Those parts of the code which aim to emulate other systems should be eminently testable, and easy to check for
    continued exactness of outcomes.

Code principles:
Code should be maximally human-readable, and well documented. Speed is desired, but not at the expense of the previous
two points.
Comments should accompany the code with two basic aims:
1.  Explain as exactly as possible what a block or line of code is _supposed_ to do
2.  Explain why the code is supposed to do this.

For example:
Rationale:
The ATP correction for the flow signal is supposed to calculate the flow in and out of the lungs based on
the flow in and out of the flow-meter. We measure the volume of air after it has cooled down and lost some humidity
(the precise cooling and loss of humidity follows an assumption down to x degrees and y%. We can then calculate what the
flow in and out of the lungs would have been, assuming body temperature and body humidity.
Implementation:
This function was written to mirror Spiroware 3.2.1 processing, based on code provided by Ecomedics.


