"""
Unit tests for the scalar field coercers in :mod:`pidibble.baseparsers`.
These are offline tests -- no network.
"""
import math

from pidibble.baseparsers import StringParser, safe_float, str2int_sig


def test_safe_float_maps_nan_literal_to_zero():
    """A literal 'nan' in a float column becomes 0.0, not a float NaN.

    ``float('nan')`` succeeds, so without the sentinel guard the value would
    sail through the coercion untouched and poison every downstream arithmetic
    result silently -- no exception, no nonconformance recorded.
    """
    for s in ('nan', 'NaN', '     nan'):   # float columns are right-justified
        v = safe_float(s)
        assert v == 0.0, f'{s!r} -> {v!r}'
        assert not math.isnan(v), f'{s!r} leaked a NaN'
    # ordinary values are untouched
    assert safe_float('-1.25') == -1.25
    assert safe_float('  3.5 ') == 3.5


def test_safe_float_nan_column_in_a_fixed_width_record():
    """The same guard, exercised through the fixed-column path a real file takes."""
    fmt = {'x': ['Float', [31, 38]], 'y': ['Float', [39, 46]], 'z': ['Float', [47, 54]]}
    sp = StringParser(fmt, {'Float': safe_float, 'String': str, 'Integer': str2int_sig})
    # columns 1-30 are the atom/residue identifiers; 31-54 are the coordinates,
    # with 'nan' written into y as some coordinate-producing tools do
    record = 'ATOM      1  N   LEU A  34'.ljust(30) + '  12.345     nan  -7.100'
    parsed = sp.parse(record)
    assert parsed['x'] == 12.345
    assert parsed['y'] == 0.0
    assert not math.isnan(parsed['y'])
    assert parsed['z'] == -7.100
    # 'nan' is a coercible value, so it is not flagged as a nonconformance
    assert sp.nonconformances == []


def test_str2int_sig_sentinel_and_sign():
    assert str2int_sig('  42 ') == 42
    assert str2int_sig(' -42 ') == -42
    assert str2int_sig('') == -1          # blank -> sentinel
    assert str2int_sig('1A') == -1        # non-numeric -> sentinel
