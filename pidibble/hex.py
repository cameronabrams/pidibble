def int_or_hex(arg,force=False):
    if type(arg)==int and not force: return arg
    assert type(arg)==str
    if force or any([(x in arg) for x in 'abcdef']):
        return int(arg,16)
    else:
        return int(arg)