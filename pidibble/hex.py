# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Function for detecting a switch from decimal to hexadecimal in integer parsing.
This is used to determine if a string should be parsed as a hexadecimal number or a plain decimal integer.
"""

class AtomSerialParser:
    """
    Callable parser for atom serial numbers that tracks the decimal-to-hex transition.
    Each instance maintains independent state, avoiding the global-variable pitfall.
    """
    def __init__(self):
        self._hex_tripped = False

    def __call__(self, arg: str) -> int:
        """
        Convert a string atom serial number to int, auto-detecting hex format.

        Parameters
        ----------
        arg : str
            The string representation of the atom serial number.

        Returns
        -------
        int
            The integer representation of the atom serial number.
        """
        assert isinstance(arg, str)
        if arg == 'nan':
            return_object = 0
        elif self._hex_tripped or any(x in arg for x in 'abcdefABCDEF'):
            return_object = int(arg, 16)
        elif '*' in arg:
            return_object = 0
        else:
            return_object = int(arg)
        if return_object > 99999 and not self._hex_tripped:
            self._hex_tripped = True
        return return_object

    def reset(self):
        """Reset the hex-tripped flag for a fresh parse."""
        self._hex_tripped = False


# Module-level instance kept for backward compatibility
_default_parser = AtomSerialParser()

def str2atomSerial(arg):
    """
    Convert a string representation of an atom serial number to an integer.  Should be used in cases were an integer series changes format from decimal representation to hexadecimal representation beyond 99999.  The transition is signaled by the presence of hexadecimal characters in the string, which sets a global flag to indicate that subsequent strings should be parsed as hexadecimal numbers.

    Parameters
    ----------
    arg : str
        The string representation of the atom serial number.

    Returns
    -------
    int
        The integer representation of the atom serial number.
    """
    return _default_parser(arg)

def hex_reset():
    """Reset the hexadecimal parsing flag."""
    _default_parser.reset()
