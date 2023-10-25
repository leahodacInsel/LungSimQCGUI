from SpirowareN2MBW import SpirowareN2MBW
from SpirowareFile import SpirowareFile
from YAMLTools import dict_from_yaml


def run_test():
    # Define raw file and config file
    filename = r"\\filer300\USERS3007\I0337516\Desktop\data\A-20220620-155737-WAETRI30052013-N2MultiBreathWashoutTest-Set2.txt"
    optionsname = r"\\filer300\USERS3007\I0337516\Desktop\data\Config_noncutout_331_precap_24ml.yaml"

    # Read the files
    file = SpirowareFile.from_fullname(filename)
    options = dict_from_yaml(optionsname)

    # Declare the analysis
    mbw_analysis = SpirowareN2MBW(file=file, options=options)

    # Specify output (here LCI)
    args = dict()
    args['variable'] = 'lci_ao'
    args['method'] = 'critical'
    args['concentration'] = 2.5
    lci = mbw_analysis.calculate_output(args)
    print(lci)


if __name__ == '__main__':
    run_test()
