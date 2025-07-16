# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Function for detecting a switch from decimal to hexadecimal in integer parsing.
This is used to determine if a string should be parsed as a hexadecimal number or a plain decimal integer.
"""
__hex_tripped__=False
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
    global __hex_tripped__
    assert type(arg)==str
    if arg=='nan':
        return_object=0
    elif __hex_tripped__ or any([(x in arg) for x in 'abcdefABCDEF']):
        return_object=int(arg,16)
    elif '*' in arg:
        return_object=0
    else:
        return_object=int(arg)
    if return_object>99999 and not __hex_tripped__:
        __hex_tripped__=True
    return return_object

def hex_reset():
    """
    Reset the hexadecimal parsing flag.
    """
    global __hex_tripped__
    __hex_tripped__=False
