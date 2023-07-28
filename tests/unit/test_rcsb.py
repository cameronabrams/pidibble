import unittest
import pytest

from pidibble.rcsb import PDBParser

def test_pdbformat():
    p=PDBParser()
    expected_sections=['record_types', 'delimiters', 'record_formats']
    assert(all([x in p.pdb_format_dict.keys() for x in expected_sections]))

def test_custom_formats():
    p=PDBParser(PDBcode='test',pdb_format_file='test_pdb_format.yaml')
    p.fetch()
    p.read()
    p.parse()
    assert 'MYREC1' in p.parsed
    assert 'MYREC2' in p.parsed
    assert p.parsed['MYREC1'][0].cussword=='FUCK'
    assert p.parsed['MYREC1'][0].residue.resName=='RRR'
    assert p.parsed['MYREC1'][0].residue.chainID=='C'
    assert p.parsed['MYREC1'][0].residue.seqNum==1111
    assert p.parsed['MYREC1'][0].residue.iCode=='I'
    assert p.parsed['MYREC2'][0].cussword=='SHIT'
    assert p.parsed['MYREC2'][0].residue.resName=='XXX'
    assert p.parsed['MYREC2'][0].residue.chainID=='D'
    assert p.parsed['MYREC2'][0].residue.seqNum==2222
    assert p.parsed['MYREC2'][0].residue.iCode=='J'
    for s in p.parsed['SITE']:
        print(s.__dict__)
    assert len(p.parsed['SITE'])==4
    assert p.parsed['SITE'][0].siteID=='AC1'
    assert p.parsed['SITE'][0].residue1.resName=='HIS'
    assert len(p.parsed['SITE'][0].residues)==3
    assert len(p.parsed['SITE'][1].residues)==5
    assert len(p.parsed['SITE'][2].residues)==5
    assert len(p.parsed['SITE'][3].residues)==11
    for s in p.parsed['SITE']:
        assert s.numRes==len(s.residues)
        print(s.continuation)
    s=p.parsed['SITE'][3]

    expected_resnames=['HIS','HIS','HIS','HIS','LEU','THR','THR','TRP','HOH','HOH','HOH']
    assert expected_resnames==[r.resName for r in s.residues]

def test_parse():
    q=PDBParser(PDBcode='4zmj')
    print(f'file {q.pdb_format_file}')
    q.fetch()
    q.read()
    q.parse()
    assert 'HEADER' in q.parsed
    assert 'VIRAL PROTEIN'==q.parsed['HEADER'].classification
    assert '04-MAY-15'==q.parsed['HEADER'].depDate
    assert '4ZMJ'==q.parsed['HEADER'].idCode
    assert type(q.parsed['COMPND'].compound)==list
    assert len(q.parsed['COMPND'].compound)==11
    assert q.parsed['COMPND'].compound[0] =='MOL_ID: 1'
    assert q.parsed['COMPND'].compound[1] =='MOLECULE: ENVELOPE GLYCOPROTEIN GP160'
    assert q.parsed['COMPND'].compound[2]=='A_FAKE_TOKEN: A_FAKE_VALUE'
    assert q.parsed['COMPND'].compound[3]=='CHAIN: G'
    assert q.parsed['TITLE'].title=='CRYSTAL STRUCTURE OF LIGAND-FREE BG505 SOSIP.664 HIV-1 ENV TRIMER THIS IS A FAKE EXTRA LINE THIS IS ANOTHER FAKE EXTRA LINE'
    # print(q.parsed['COMPND'].__dict__)
    assert type(q.parsed['COMPND'].tokengroups)==dict
    assert 'MOL_ID.1' in q.parsed['COMPND'].tokengroups['compound']
    assert 'MOL_ID.2' in q.parsed['COMPND'].tokengroups['compound']
    m1=q.parsed['COMPND'].tokengroups['compound']['MOL_ID.1']
    assert m1.MOLECULE=='ENVELOPE GLYCOPROTEIN GP160'
    assert m1.MUTATION=='YES'
    assert type(q.parsed['SOURCE'].tokengroups)==dict
    assert 'MOL_ID.1' in q.parsed['SOURCE'].tokengroups['srcName']
    assert 'MOL_ID.2' in q.parsed['SOURCE'].tokengroups['srcName']
    m1=q.parsed['SOURCE'].tokengroups['srcName']['MOL_ID.1']
    assert m1.ORGANISM_SCIENTIFIC=='HUMAN IMMUNODEFICIENCY VIRUS 1'
    # print(q.parsed['COMPND'].tokens)
    # assert 'MOL_ID.2' in q.parsed['COMPND'].tokens['compound'].keys()
    # assert len(q.parsed['COMPND'].tokens['compound'])==2
    # assert type(q.parsed['COMPND'].tokens['compound']['MOL_ID.1'])==dict
    # assert type(q.parsed['COMPND'].tokens['compound']['MOL_ID.2'])==dict
    # assert len(q.parsed['COMPND'].tokens['compound']['MOL_ID.1'])==5
    # assert q.parsed['COMPND'].tokens['compound']['MOL_ID.1']['MOLECULE']=='ENVELOPE GLYCOPROTEIN GP160'
    # assert q.parsed['COMPND'].tokens['compound']['MOL_ID.1']['MUTATION']=='YES'
    # assert type(q.parsed['SOURCE'].tokens)==dict
    # assert 'MOL_ID.1' in q.parsed['SOURCE'].tokens['srcName'].keys()
    # assert 'MOL_ID.2' in q.parsed['SOURCE'].tokens['srcName'].keys()
    # assert len(q.parsed['SOURCE'].tokens['srcName'])==2
    # assert q.parsed['SOURCE'].tokens['srcName']['MOL_ID.1']['ORGANISM_SCIENTIFIC']=='HUMAN IMMUNODEFICIENCY VIRUS 1'
    assert len(q.parsed['ATOM'])==4518
    assert len(q.parsed['ANISOU'])==4518
    assert len(q.parsed['HETATM'])==338
    assert len(q.parsed['HET'])==25
    assert len(q.parsed['LINK'])==25
    assert not 'MODRES' in q.parsed
    assert len(q.parsed['REVDAT'])==6
    assert q.parsed['REVDAT'][0].records==['COMPND', 'REMARK', 'HETNAM', 'LINK', 'SITE', 'ATOM']
    assert len(q.parsed['SEQADV'])==10
    assert q.parsed['SEQADV'][0].residue.resName=='ASN'
    assert q.parsed['SEQADV'][0].residue.chainID=='G'
    assert q.parsed['SEQADV'][0].residue.seqNum==332
    assert q.parsed['SEQADV'][0].residue.iCode==''
    assert q.parsed['SEQADV'][0].database=='UNP'
    assert q.parsed['SEQADV'][0].dbAccession=='Q2N0S6'
    assert q.parsed['SEQADV'][0].dbRes=='THR'
    assert q.parsed['SEQADV'][0].dbSeq==330
    assert q.parsed['SEQADV'][0].conflict=='ENGINEERED MUTATION'
    assert len(q.parsed['SSBOND'])==11
    assert q.parsed['SSBOND'][2].residue1.chainID=='G'
    assert q.parsed['SSBOND'][2].residue1.seqNum==126
    assert q.parsed['SSBOND'][2].residue2.chainID=='G'
    assert q.parsed['SSBOND'][2].residue2.seqNum==196
    assert len(q.parsed['SEQRES'][0].resNames)==481
    assert len(q.parsed['SEQRES'][1].resNames)==153
    expected_seq='ALA GLU ASN LEU TRP VAL THR VAL TYR TYR GLY'.split()
    assert q.parsed['SEQRES'][0].resNames[:len(expected_seq)]==expected_seq
    "ATOM   4519  OD2 ASP B 664     -15.056 125.079  66.899  1.00142.18           O  "
    assert q.parsed['ATOM'][0].residue.resName=='LEU'
    assert q.parsed['ATOM'][-1].name=='OD2'
    assert q.parsed['ATOM'][0].serial==1
    assert q.parsed['ATOM'][-1].serial==4519
    assert q.parsed['ATOM'][-1].residue.resName=='ASP'
    assert len(q.parsed['TER'])==2
    assert q.parsed['TER'][0].residue.resName=='VAL'
    assert q.parsed['TER'][0].residue.chainID=='G'
    assert q.parsed['TER'][0].residue.seqNum==505
    assert q.parsed['TER'][0].residue.iCode==''

    assert q.parsed['TER'][1].residue.resName=='ASP'
    assert q.parsed['TER'][1].residue.chainID=='B'
    assert q.parsed['TER'][1].residue.seqNum==664
    assert q.parsed['TER'][1].residue.iCode==''

    # assert len(q.parsed['JRNL'])==16

    # assert len(q.parsed['REMARK'])==649 # unparsed remarks


