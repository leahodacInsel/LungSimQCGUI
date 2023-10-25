from progress.bar import ChargingBar


class LungSimProgressBar(ChargingBar):
    """
    Ehm. The ChargingBar, and in fact all the Bars in progress library don't display in PyCharm. Hoping that this does
    not horribly break everything, this is a replacement class which simply has a different line writing function.
    """

    def writeln(self, line):
        width = len(line)
        if width < self._max_width:
            # Add padding to cover previous contents
            line += ' ' * (self._max_width - width)
        else:
            self._max_width = width
        print('\r' + line, end='')
        self.file.flush()
