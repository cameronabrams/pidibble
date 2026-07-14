"""
Unit tests for nonconformance tracking: NonconformanceRegistry aggregation and
StringParser's per-instance collection. These are offline tests — no network.
"""
import logging
from pidibble.baseparsers import NonconformanceRegistry, StringParser


def test_registry_empty_by_default():
    reg = NonconformanceRegistry()
    assert len(reg) == 0
    assert not reg
    assert reg.types() == []


def test_registry_aggregates_by_type():
    reg = NonconformanceRegistry()
    for _ in range(5):
        reg.add('ATOM', {'field': 'charge', 'kind': 'not coercible to Float',
                         'byte_range': [79, 80], 'value': '1-'})
    reg.add('ANISOU', {'field': 'charge', 'kind': 'not coercible to Float',
                       'byte_range': [79, 80], 'value': '2+'})
    # two distinct types, not six instances
    assert len(reg) == 2
    assert bool(reg)
    types = {(rk, field): (count, br, example)
             for (rk, field, kind, count, br, example) in reg.types()}
    assert types[('ATOM', 'charge')] == (5, [79, 80], '1-')
    assert types[('ANISOU', 'charge')] == (1, [79, 80], '2+')


def test_registry_keeps_first_exemplar():
    reg = NonconformanceRegistry()
    reg.add('ATOM', {'field': 'charge', 'kind': 'k', 'byte_range': [79, 80], 'value': 'first'})
    reg.add('ATOM', {'field': 'charge', 'kind': 'k', 'byte_range': [79, 80], 'value': 'second'})
    (rk, field, kind, count, br, example), = reg.types()
    assert count == 2
    assert example == 'first'


def test_registry_report_emits_one_info_line_per_type(caplog):
    reg = NonconformanceRegistry()
    for _ in range(3):
        reg.add('ATOM', {'field': 'charge', 'kind': 'not coercible to Float',
                         'byte_range': [79, 80], 'value': '1-'})
    with caplog.at_level(logging.INFO):
        reg.report(logging.getLogger('pidibble.test'), label='6cm3')
    infos = [r.getMessage() for r in caplog.records if r.levelno == logging.INFO]
    # one header line + one per-type line
    assert len(infos) == 2
    assert '1 PDB-format nonconformance type(s)' in infos[0]
    assert '6cm3' in infos[0]
    assert 'ATOM.charge' in infos[1]
    assert '3 value(s)' in infos[1]
    assert "'1-'" in infos[1]
    assert 'cols 79-80' in infos[1]


def test_registry_report_silent_when_empty(caplog):
    reg = NonconformanceRegistry()
    with caplog.at_level(logging.INFO):
        reg.report(logging.getLogger('pidibble.test'))
    assert [r for r in caplog.records if r.levelno == logging.INFO] == []


def test_stringparser_collects_coercion_failure():
    # a single Float field at columns 1-3 fed a non-numeric value
    parser = StringParser({'q': ['Float', [1, 3]]}, {'Float': float, 'String': str})
    parser.parse('ab ')
    assert len(parser.nonconformances) == 1
    nc = parser.nonconformances[0]
    assert nc['field'] == 'q'
    assert nc['kind'] == 'not coercible to Float'
    assert nc['byte_range'] == [1, 3]
    assert nc['value'] == 'ab'


def test_stringparser_no_nonconformance_on_blank_or_valid():
    parser = StringParser({'q': ['Float', [1, 3]]}, {'Float': float, 'String': str})
    parser.parse('   ')            # blank -> '' , no coercion attempted
    parser.parse('1.5')            # valid float
    assert parser.nonconformances == []


def test_stringparser_drains_into_registry_by_record_key():
    reg = NonconformanceRegistry()
    parser = StringParser({'charge': ['Float', [79, 80]]}, {'Float': float, 'String': str})
    rec = ' ' * 78 + '1-'
    parser.parse(rec)
    for nc in parser.nonconformances:
        reg.add('ATOM', nc)
    (rk, field, kind, count, br, example), = reg.types()
    assert (rk, field, count, example) == ('ATOM', 'charge', 1, '1-')
