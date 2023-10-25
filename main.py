from MBWAnalysis import MBWAnalysis
from LungSimGUI import run_gui


def run_batch():
    a = MBWAnalysis.from_fileopenbox()
    a.calculate_all()


if __name__ == '__main__':
    run_batch()
    # run_gui()
