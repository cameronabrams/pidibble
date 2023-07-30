import unittest
import pytest
import numpy as np
from pidibble.rcsb import PDBParser

def test_pdbformat():
    p=PDBParser()
    expected_sections=['record_types', 'delimiters', 'record_formats']
    assert(all([x in p.pdb_format_dict.keys() for x in expected_sections]))

def test_custom_formats():
    p=PDBParser(PDBcode='test',pdb_format_file='test_pdb_format.yaml').parse()
    # p.fetch()
    # p.read()
    # p.parse()
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
    # for s in p.parsed['SITE']:
    #     print(s.__dict__)
    assert len(p.parsed['SITE'])==4
    assert p.parsed['SITE'][0].siteID=='AC1'
    assert p.parsed['SITE'][0].residue1.resName=='HIS'
    assert len(p.parsed['SITE'][0].residues)==3
    assert len(p.parsed['SITE'][1].residues)==5
    assert len(p.parsed['SITE'][2].residues)==5
    assert len(p.parsed['SITE'][3].residues)==11
    for s in p.parsed['SITE']:
        assert s.numRes==len(s.residues)
        # print(s.continuation)
    s=p.parsed['SITE'][3]

    expected_resnames=['HIS','HIS','HIS','HIS','LEU','THR','THR','TRP','HOH','HOH','HOH']
    assert expected_resnames==[r.resName for r in s.residues]

def test_parse():
    q=PDBParser(PDBcode='4zmj').parse()
    # print(f'file {q.pdb_format_file}')
    # q.fetch()
    # q.read()
    # q.parse()
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

    assert len(q.parsed['JRNL.AUTH'].authorList)==53
    assert q.parsed['JRNL.AUTH'].authorList[0]=='Y.DO KWON'
    assert q.parsed['JRNL.AUTH'].authorList[-1]=='P.D.KWONG'
    assert q.parsed['JRNL.TITL'].title=='CRYSTAL STRUCTURE, CONFORMATIONAL FIXATION AND ENTRY-RELATED INTERACTIONS OF MATURE LIGAND-FREE HIV-1 ENV.'
    assert q.parsed['JRNL.REF'].pubName=='NAT.STRUCT.MOL.BIOL.'
    assert q.parsed['JRNL.REF'].volume=='22'
    assert q.parsed['JRNL.REF'].page=='522'
    assert q.parsed['JRNL.REF'].year==2015
    assert q.parsed['JRNL.REFN'].issnORessn=='ESSN'
    assert q.parsed['JRNL.REFN'].issn=='1545-9985'
    assert q.parsed['JRNL.PMID'].pmid==26098315
    assert q.parsed['JRNL.DOI'].doi=='10.1038/NSMB.3051'

    assert q.parsed['REMARK.2'].resolutionmsg==['', '3.31 ANGSTROMS.']
    assert len(q.parsed['REMARK.3'].refinementdetails)==331
    assert q.parsed['REMARK.3'].refinementdetails[0]==''
    assert q.parsed['REMARK.3'].refinementdetails[3]=='AUTHORS     : PAUL ADAMS,PAVEL AFONINE,VINCENT CHEN,IAN'
    assert q.parsed['REMARK.3'].refinementdetails[-1]=='OTHER REFINEMENT REMARKS: NULL'

    assert q.parsed['REMARK.4'].formatmsg[1]=='4ZMJ COMPLIES WITH FORMAT V. 3.30, 13-JUL-11'
    assert q.parsed['REMARK.100'].processmsg[-1]=='THE DEPOSITION ID IS D_1000209535.'
    assert len(q.parsed['REMARK.200'].crystallizationmsg)==49
    assert q.parsed['REMARK.200'].crystallizationmsg[-1]=='REMARK: NULL'
    assert len(q.parsed['REMARK.290'].crystallographicautodetails)==40
    assert q.parsed['REMARK.290'].crystallographicautodetails[-1]=='REMARK: NULL'

    # print(q.parsed['REMARK.300'].tokengroups['biomoleculedeclarations'])

    assert q.parsed['REMARK.300'].tokengroups['biomoleculedeclarations']['BIOMOLECULE'].BIOMOLECULE==['1']
    assert q.parsed['REMARK.350'].tokengroups['biomoleculespecifications']['BIOMOLECULE.1'].AUTH_BIO_UNIT=='HEXAMERIC'
    assert q.parsed['REMARK.350'].tokengroups['biomoleculespecifications']['BIOMOLECULE.1'].SOFT_QUAT_STRUCT=='HEXAMERIC'
    assert q.parsed['REMARK.350'].tokengroups['biomoleculespecifications']['BIOMOLECULE.1'].CHAIN_BIOMT==['G', 'B', 'A', 'C', 'D']
    # for k in q.parsed.keys():
    #     print(k)
    # print(q.parsed['REMARK.350'])
    # print(q.parsed['REMARK.350.BIOMT'])
    assert q.parsed['REMARK.350.BIOMT'].rowNum==[1,2,3,1,2,3,1,2,3]
    assert q.parsed['REMARK.350.BIOMT'].transNum==[1,1,1,2,2,2,3,3,3]
    assert q.parsed['REMARK.350.BIOMT'].M1==[1.0,0.0,0.0,-0.5,0.866025,0.0,-0.5,-0.866025,0.0]
    rec=q.parsed['REMARK.350.BIOMT']
    M=[]
    for i in range(3):
        M.append(np.zeros((3,3)))
    for t,r,m1,m2,m3 in zip(rec.transNum,rec.rowNum,rec.M1,rec.M2,rec.M3):
        M[t-1][r-1,0]=m1
        M[t-1][r-1,1]=m2
        M[t-1][r-1,2]=m3
    
    expM=[
        np.array([
            [1.0,0.0,0.0],
            [0.0,1.0,0.0],
            [0.0,0.0,1.0]
        ]),
        np.array([
            [-0.5,-0.866025,0.0],
            [0.866025,-0.5,0.0],
            [0.0,0.0,1.0]
        ]),
        np.array([
            [-0.5,0.866025,0.0],
            [-0.866025,-0.5,0.0],
            [0.0,0.0,1.0]
        ])       
    ]
    assert all(np.array_equal(m,em) for m,em in zip(M,expM))

    assert q.parsed['REMARK.465'].freetext[1].strip()=='REMARK 465 MISSING RESIDUES'
    assert q.parsed['REMARK.465'].tables['MISSING'][0].resname=='ALA'
    assert q.parsed['REMARK.465'].tables['MISSING'][0].modelNum==''
    assert q.parsed['REMARK.465'].tables['MISSING'][0].chainID=='G'
    assert q.parsed['REMARK.465'].tables['MISSING'][0].resseqnum==31
    assert len(q.parsed['REMARK.465'].tables['MISSING'])==61

    assert q.parsed['REMARK.500'].tables['RAMA_OUTLIERS'][0].residue.resName=='PRO'
    assert q.parsed['REMARK.500'].tables['NONCISTRANS'][0].residueC.resName=='ALA'
    assert q.parsed['REMARK.500'].tables['NONCISTRANS'][0].omega==146.93
    assert q.parsed['REMARK.500'].tables['NONCISTRANS'][-1].omega==149.23

#REMARK 500 LEU G  494     GLY G  495                  149.23                    
