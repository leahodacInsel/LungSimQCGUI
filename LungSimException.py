class LungSimException(Exception):
    """
    A LungSim exception is a type of exception that can occur during the attempt to calculate any given outcome.
    It is used for common and predictable failures to calculate due to missing information, missing pieces of signals or
    measurements, and will be the output of a calculation instead of the queried outcome.
    """
    def __init__(self, output="Unspecified error"):
        """
        The main content of a LungSim error is the output to the MBW object, which will be inserted into
        the results instead of the desired calculation outcome.
        :param output:
        """
        self.output = output
