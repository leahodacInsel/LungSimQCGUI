from YAMLTools import dict_from_yaml


class Config:
    """
    Parent class of the Config family. Config objects contain parameters and options for the processing of MBW
    algorithms. They are stored in its attributes, and can be modified / altered en bloc using .yaml files or
    dictionaries.
    """
    def __init__(self, **kwargs):
        """
        Here is where default parameters can be stored which affect all types of Config.
        """

    def __setitem__(self, key, value):
        # TODO: Only allow compatible variable types
        if key in self.__dict__:
            self.__dict__.pop(key)
        self.__dict__[key] = value

    def apply_yaml(self, path=None):
        """
        Function to apply the contents of a .yaml options file to a config object.
        :param path: path to the location of the .yaml file to be applied to a config object.
        """
        self.apply_dictionary(dict_from_yaml(path))

    def apply_dictionary(self, dictionary=None):
        """
        Function to apply the contents of a dictionary to a config object.
        :param dictionary: dictionary with contents to be applied to a config object
        """
        if dictionary is not None:
            for key, value in dictionary.items():
                self.__setitem__(key=key, value=value)