from ruamel.yaml import YAML

def dict_from_yaml(fullname=None):
    """
    Function to apply the contents of a .yaml options file to a config object.
    """
    with open(fullname) as file:
        yaml = YAML(typ='safe')
        return yaml.load(file)
