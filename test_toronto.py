from MBWAnalysis import MBWAnalysis
from LungSimGUI import run_gui


def run_batch():
    filenames = ['C:\\Users\\I0326962\\Local Files\\Former Desktop\\Completely fresh Wbreath take\\In vivo\\Used in vivo\\Spiroware\\A-20170227-112446-MUSMAD26122016-SF6MultiBreathWashoutTest-Set1.txt']
    config = 'C:\\Users\\I0326962\\Local Files\\Former Desktop\\CLCI\\LS config file\LCR config SF6 resync.yaml'
    results = 'C:\\Users\\I0326962\\Local Files\\Former Desktop\\CLCI'
    a = MBWAnalysis(filenames=filenames, options_config=config, results_path=results)
    a.calculate_all()


if __name__ == '__main__':
    run_batch()
    # run_gui()