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


class HexSerialEncoder:
    """
    Inverse of :class:`AtomSerialParser`: render an atom serial number for a
    fixed-width field, switching from decimal to hexadecimal once any serial
    exceeds 99999 and staying hexadecimal thereafter.

    The switch is stateful and permanent, mirroring the parser's ``_hex_tripped``
    flag exactly. This matters because the parser, once tripped, reads *every*
    subsequent serial field as hex — including small back-references in
    ``CONECT`` — so those must be encoded as hex too (e.g. serial ``10`` becomes
    ``"A"``) for the file to round-trip. One encoder instance must therefore be
    threaded through a whole document's serial fields in emission order.
    """
    def __init__(self):
        self._hex_tripped = False

    def __call__(self, serial: int) -> str:
        """
        Encode one serial as a digit string (the caller pads it to width).

        Parameters
        ----------
        serial : int
            The atom serial number.

        Returns
        -------
        str
            Decimal digits before the decimal-to-hex trip, uppercase hex after.
        """
        iv = int(serial)
        if iv > 99999:
            self._hex_tripped = True
        return format(iv, 'X') if self._hex_tripped else str(iv)

    def reset(self):
        """Reset to decimal mode for a fresh document."""
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
