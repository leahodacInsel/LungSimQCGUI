from MBWAnalysis import MBWAnalysis
from LungSimGUI import run_gui
from LungSimGUI import run_test


def run_batch():
    filenames = ['C:\\Users\\I0326962\\Local Files\\Former Desktop\\Quality control features\\A-20131031-145011-FU_BILD_281_JP220707-N2MultiBreathWashoutTest-Set2.txt', 'C:\\Users\\I0326962\\Local Files\\Former Desktop\\Quality control features\\A-20160706-125931-CL180516-SF6MultiBreathWashoutTest-Set1.txt']
    config = 'C:\\Users\\I0326962\\Local Files\\Former Desktop\\Quality control features\\QC config.yaml'
    results = 'C:\\Users\\I0326962\\Local Files\\Former Desktop\\Quality control features'
    a = MBWAnalysis(filenames=filenames, options_config=config, results_path=results)
    a.calculate_all()


if __name__ == '__main__':
    #run_batch()
    #run_gui()
    run_test()
