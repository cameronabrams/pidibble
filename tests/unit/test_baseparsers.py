"""
Unit tests for the scalar field coercers in :mod:`pidibble.baseparsers`.
These are offline tests -- no network.
"""
import copy
import json
import math
import pickle

import pytest
import yaml

from pidibble.baseparsers import EmptyField, StringParser, safe_float, str2int_sig

TYPEMAP = {'Float': safe_float, 'String': str, 'Integer': str2int_sig}


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
    sp = StringParser(fmt, TYPEMAP)
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


# --- EmptyField: absent typed fields must not impersonate numbers -----------

# a truncated SSBOND: the residue columns are canonical, but the record stops
# before sym1/sym2/length -- the shape a minimal writer emits
SSBOND_FMT = {'serNum': ['Integer', [8, 10]], 'sym1': ['String', [60, 65]],
              'length': ['Float', [74, 78]]}
SSBOND_SHORT = 'SSBOND   1 CYS A   54    CYS A   74'


def _short_ssbond():
    sp = StringParser(SSBOND_FMT, TYPEMAP)
    return sp, sp.parse(SSBOND_SHORT)


def test_absent_typed_field_still_reads_as_empty():
    """EmptyField must satisfy every emptiness test the library and its callers
    already make, or this becomes a breaking change instead of a safety net."""
    _, parsed = _short_ssbond()
    length = parsed['length']
    assert length == ''
    assert not length
    assert isinstance(length, str)
    assert f'[{length}]' == '[]'
    assert length.strip() == ''
    # an absent *String* field keeps a plain '': there the empty string is a
    # legitimate value, not a missing one
    assert type(parsed['sym1']) is str


def test_absent_float_field_refuses_to_act_like_a_number():
    """The silent bite: a plain '' in a Float field makes ``length * 2``
    evaluate to '' by string repetition instead of failing."""
    _, parsed = _short_ssbond()
    length = parsed['length']
    for op in (lambda: length * 2, lambda: 2 * length, lambda: float(length),
               lambda: int(length), lambda: length + 1, lambda: length - 1,
               lambda: length / 2, lambda: round(length, 2), lambda: -length,
               lambda: abs(length), lambda: length % 2):
        with pytest.raises(TypeError) as exc:
            op()
        # the message must identify the field, not just the type
        assert "'length'" in str(exc.value)
        assert 'Float' in str(exc.value)
    # and it says why it is empty, and how to test for it
    assert 'absent from the source record' in repr(length)
    assert length.field == 'length' and length.typestring == 'Float'


def test_empty_field_survives_the_usual_round_trips():
    """Anything that worked on a plain '' must keep working: PyYAML in
    particular dispatches on exact type and would otherwise refuse to
    represent parsed records."""
    _, parsed = _short_ssbond()
    length = parsed['length']
    assert yaml.safe_dump({'length': length}) == yaml.safe_dump({'length': ''})
    assert yaml.dump({'length': length}) == yaml.dump({'length': ''})
    assert json.dumps({'length': length}) == '{"length": ""}'
    assert copy.deepcopy(length) == '' and pickle.loads(pickle.dumps(length)) == ''
    assert 'x' + length == 'x'          # concatenation is still meaningful
    assert ''.join(['a', length]) == 'a'


def test_uncoercible_field_reports_the_offending_value():
    """A field that *had* text but could not be coerced gets the same guard,
    with the offending value in the reason, and is still counted as a
    nonconformance."""
    sp = StringParser({'charge': ['Float', [79, 80]]}, TYPEMAP)
    parsed = sp.parse('ATOM      1  N   LEU A  34'.ljust(78) + 'O-')
    charge = parsed['charge']
    assert charge == ''                      # empty, as before
    assert isinstance(charge, EmptyField)
    assert "'O-'" in charge.reason
    with pytest.raises(TypeError, match="'charge'"):
        float(charge)
    assert [nc['kind'] for nc in sp.nonconformances] == ['not coercible to Float']


def test_str2int_sig_sentinel_and_sign():
    assert str2int_sig('  42 ') == 42
    assert str2int_sig(' -42 ') == -42
    assert str2int_sig('') == -1          # blank -> sentinel
    assert str2int_sig('1A') == -1        # non-numeric -> sentinel
