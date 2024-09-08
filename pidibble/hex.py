__hex_tripped__=False
def str2atomSerial(arg):
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
    global __hex_tripped__
    __hex_tripped__=False
