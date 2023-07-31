import numpy as np
from pidibble.rcsb import PDBParser, PDBRecord, get_symm_ops
import unittest
import pytest

def test_pdbformat():
    p=PDBParser()
    expected_sections=['record_formats']
    assert(all([x in p.pdb_format_dict.keys() for x in expected_sections]))

def test_custom_formats():
    p=PDBParser(PDBcode='test',pdb_format_file='test_pdb_format.yaml').parse()
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

class TestParse(unittest.TestCase):

    # def __init__(self):

    def setUp(self):
        self.P=PDBParser(PDBcode='4zmj').parse()
    
    def test_header(self):
        q=self.P
        self.assertTrue('HEADER' in q.parsed)
        self.assertEqual('VIRAL PROTEIN',q.parsed['HEADER'].classification)
        self.assertEqual('04-MAY-15',q.parsed['HEADER'].depDate)
        self.assertEqual('4ZMJ',q.parsed['HEADER'].idCode)
    
    def test_obslte(self):
        q=self.P
        self.assertTrue('OBSLTE' not in q.parsed)

    def test_title(self):
        rec=self.P.parsed['TITLE']
        self.assertEqual(rec.title,'CRYSTAL STRUCTURE OF LIGAND-FREE BG505 SOSIP.664 HIV-1 ENV TRIMER THIS IS A FAKE EXTRA LINE THIS IS ANOTHER FAKE EXTRA LINE')

    def test_compnd(self):
        q=self.P
        rec=q.parsed['COMPND']
        self.assertEqual(type(rec.compound),list)
        self.assertEqual(11,len(rec.compound))
        self.assertTrue('MOL_ID.1' in rec.tokengroups['compound'])
        m1=rec.tokengroups['compound']['MOL_ID.1']
        self.assertEqual('ENVELOPE GLYCOPROTEIN GP160',m1.MOLECULE)
        self.assertEqual('YES',m1.MUTATION)
        self.assertEqual(['G'],m1.CHAIN)

    def test_source(self):
        q=self.P
        rec=q.parsed['SOURCE']
        m1=q.parsed['SOURCE'].tokengroups['srcName']['MOL_ID.1']
        self.assertEqual(m1.ORGANISM_SCIENTIFIC,'HUMAN IMMUNODEFICIENCY VIRUS 1')

    def test_split(self):
        self.assertTrue('SPLIT' not in self.P.parsed)

    def test_sprsde(self):
        self.assertTrue('SPRSDE' not in self.P.parsed)

    def test_cispep(self):
        rec=self.P.parsed['CISPEP']
        self.assertEqual(len(rec),3)
        arec=rec[0]
        self.assertEqual(arec.residue1.resName,'ILE')
        self.assertEqual(arec.residue1.chainID,'G')
        self.assertEqual(arec.residue1.seqNum,138)
        self.assertEqual(arec.residue2.resName,'THR')
        self.assertEqual(arec.residue2.chainID,'G')
        self.assertEqual(arec.residue2.seqNum,139)
        self.assertEqual(arec.modNum,0)
        self.assertEqual(arec.measure,-24.5)
        arec=rec[-1]
        self.assertEqual(arec.residue1.resName,'SER')
        self.assertEqual(arec.residue1.chainID,'G')
        self.assertEqual(arec.residue1.seqNum,460)
        self.assertEqual(arec.residue2.resName,'THR')
        self.assertEqual(arec.residue2.chainID,'G')
        self.assertEqual(arec.residue2.seqNum,461)
        self.assertEqual(arec.modNum,0)
        self.assertEqual(arec.measure,6.44)

    def test_conect(self):
        rec=self.P.parsed['CONECT']
        self.assertEqual(len(rec),378)
        arec=rec[-8]
        self.assertEqual(arec.serial,4851)
        self.assertEqual(arec.partners,[4852,4853,4858])

# DBREF  4ZMJ G   31   507  UNP    Q2N0S6   Q2N0S6_9HIV1    30    504             
# DBREF  4ZMJ B  512   664  UNP    Q2N0S6   Q2N0S6_9HIV1   509    661 
    def test_dbref(self):
        rec=self.P.parsed['DBREF']
        self.assertEqual(len(rec),2)
        a=rec[0]
        self.assertEqual([a.idCode,a.chainID,a.seqBegin,a.insertBegin,a.seqEnd,a.insertEnd],['4ZMJ','G',31,'',507,''])
        self.assertEqual([a.database,a.dbAccession,a.dbIdCode,a.dbseqBegin,a.dbinsBegin,a.dbseqEnd,a.dbinsEnd],['UNP','Q2N0S6','Q2N0S6_9HIV1',30,'',504,''])

    def test_dbref12(self):
        self.assertTrue('DBREF1' not in self.P.parsed)
        self.assertTrue('DBREF2' not in self.P.parsed)

    def test_helix(self):
        rec=self.P.parsed['HELIX']
        self.assertEqual(len(rec),13)
        h=rec[0]
        self.assertEqual([h.serNum,h.helixID,h.helixClass,h.length],[1,'AA1',1,6])
        self.assertEqual([h.initRes.resName,h.initRes.seqNum,h.initRes.iCode,h.initRes.chainID],['ALA',58,'','G'])
        self.assertEqual([h.endRes.resName,h.endRes.seqNum,h.endRes.iCode,h.endRes.chainID],['THR',63,'','G'])
        h=rec[5]
        self.assertEqual([h.serNum,h.helixID,h.helixClass,h.length],[6,'AA6',1,10])
        self.assertEqual([h.initRes.resName,h.initRes.seqNum,h.initRes.iCode,h.initRes.chainID],['MET',475,'','G'])
        self.assertEqual([h.endRes.resName,h.endRes.seqNum,h.endRes.iCode,h.endRes.chainID],['TYR',484,'','G'])

    def test_modres(self):
        self.assertTrue('MODRES' not in self.P.parsed)

    def test_mtrix123(self):
        self.assertTrue('MTRIX1' not in self.P.parsed)
        self.assertTrue('MTRIX2' not in self.P.parsed)
        self.assertTrue('MTRIX3' not in self.P.parsed)

    def test_keywds(self):
        rec=self.P.parsed['KEYWDS']
        self.assertEqual(len(rec.keywds),5)
        exp_keywds=[x.strip() for x in 'HIV-1, ENV TRIMER, UNLIGANDED, BG505 SOSIP, VIRAL PROTEIN'.split(',')]
        self.assertEqual(rec.keywds,exp_keywds)

    def test_expdta(self):
        rec=self.P.parsed['EXPDTA']
        self.assertEqual(rec.technique,'X-RAY DIFFRACTION')

    def test_author(self):
        rec=self.P.parsed['AUTHOR']
        self.assertEqual(rec.authorList,['Y.D.KWON','P.D.KWONG'])

    def test_mdltyp(self):
        self.assertTrue('MDLTYP' not in self.P.parsed)

    def test_atoms(self):
        atoms=self.P.parsed['ATOM']
        self.assertEqual(len(atoms),4518)
        anisou=self.P.parsed['ANISOU']
        self.assertEqual(len(anisou),4518)
        hetatm=self.P.parsed['HETATM']
        self.assertEqual(len(hetatm),338)
        anatom=atoms[0]
        self.assertEqual(anatom.residue.resName,'LEU')
        self.assertEqual(anatom.residue.seqNum,34)
        self.assertEqual(anatom.residue.chainID,'G')
        self.assertEqual(anatom.name,'N')
        self.assertEqual(anatom.serial,1)
        self.assertEqual(anatom.altLoc,'')
        self.assertEqual(anatom.occupancy,1.0)
        self.assertEqual(anatom.tempFactor,137.71)
        self.assertEqual(anatom.element,'N')
        self.assertEqual(anatom.charge,'')
        anatom=atoms[456]
        self.assertEqual(anatom.serial,457)
        self.assertEqual(anatom.residue.resName,'THR')
        self.assertEqual(anatom.residue.seqNum,90)
        self.assertEqual(anatom.residue.chainID,'G')
        self.assertEqual(anatom.name,'OG1')
        self.assertEqual(anatom.altLoc,'')
        self.assertEqual(anatom.occupancy,1.0)
        self.assertEqual(anatom.tempFactor,117.68)
        self.assertEqual(anatom.element,'O')
        self.assertEqual(anatom.charge,'')

    def test_het(self):
        het=self.P.parsed['HET']
        self.assertEqual(len(het),25)
        ahet=het[-1]
        self.assertEqual(ahet.residue.resName,'NAG')
        ahet=het[2]
        self.assertEqual(ahet.residue.resName,'BMA')
    def test_links(self):
        links=self.P.parsed['LINK']
        self.assertEqual(len(links),25)
        alink=links[0]
        self.assertEqual(alink.name1,'ND2')
        self.assertEqual(alink.altLoc1,'')
        self.assertEqual(alink.residue1.resName,'ASN')
        self.assertEqual(alink.residue2.resName,'NAG')
        self.assertEqual(alink.residue1.seqNum,156)
        self.assertEqual(alink.residue2.seqNum,615)
        alink=links[-1]
        self.assertEqual(alink.name1,'O4')
        self.assertEqual(alink.altLoc1,'')
        self.assertEqual(alink.residue1.resName,'NAG')
        self.assertEqual(alink.residue2.resName,'NAG')
        self.assertEqual(alink.residue1.seqNum,1)
        self.assertEqual(alink.residue2.seqNum,2)
        self.assertEqual(alink.residue1.chainID,'D')
        self.assertEqual(alink.residue2.chainID,'D')


    def test_revdat(self):
        q=self.P
        rec=q.parsed['REVDAT']
        self.assertEqual(len(rec),6)
        self.assertEqual(rec[0].modNum,6)
        self.assertEqual(rec[0].modDate,'29-JUL-20')
        self.assertEqual(rec[0].modId,'4ZMJ')
        self.assertEqual(rec[0].modType,1)
        self.assertEqual(rec[0].records,['COMPND', 'REMARK', 'HETNAM', 'LINK', 'SITE', 'ATOM'])
        self.assertEqual(rec[-1].modType,0)

    def test_seqadv(self):
        rec=self.P.parsed['SEQADV']
        self.assertEqual(len(rec),10)
        self.assertEqual(rec[0].residue.resName,'ASN')
        self.assertEqual(rec[0].residue.chainID,'G')
        self.assertEqual(rec[0].residue.seqNum,332)
        self.assertEqual(rec[0].residue.iCode,'')
        self.assertEqual(rec[0].database,'UNP')
        self.assertEqual(rec[0].dbAccession,'Q2N0S6')
        self.assertEqual(rec[0].dbRes,'THR')
        self.assertEqual(rec[0].dbSeq,330)
        self.assertEqual(rec[0].conflict,'ENGINEERED MUTATION')

    def test_sheet(self):
        rec=self.P.parsed['SHEET']
        self.assertTrue(len(rec),45)
        sid='AA1'
        s_aa1=[x for x in rec if x.sheetID==sid]
        self.assertTrue(len(s_aa1),3)
        sid='AA7'
        s_aa7=[x for x in rec if x.sheetID==sid]
        self.assertTrue(len(s_aa7),7)
        self.assertTrue([x.initRes.resName for x in s_aa7],['LEU','ILE','ILE','HIS','SER','GLU','HIS'])

    def test_ssbond(self):
        rec=self.P.parsed['SSBOND']
        self.assertEqual(len(rec),11)
        anssbond=rec[0]
        self.assertEqual(anssbond.residue1.chainID,'G')
        self.assertEqual(anssbond.residue1.seqNum,54)
        self.assertEqual(anssbond.residue2.chainID,'G')
        self.assertEqual(anssbond.residue2.seqNum,74)
        anssbond=rec[9]
        self.assertEqual(anssbond.residue1.chainID,'G')
        self.assertEqual(anssbond.residue1.seqNum,501)
        self.assertEqual(anssbond.residue2.chainID,'B')
        self.assertEqual(anssbond.residue2.seqNum,605)

    def test_cryst1(self):
        rec=self.P.parsed['CRYST1']
        """ 
        CRYST1  107.180  107.180  103.060  90.00  90.00 120.00 P 63          6   
        """
        self.assertEqual([rec.a,rec.b,rec.c],[107.180,107.180,103.060])
        self.assertEqual([rec.alpha,rec.beta,rec.gamma],[90,90,120])
        self.assertEqual(rec.sGroup,'P 63')
        self.assertEqual(rec.z,6)

    def test_origx123(self):
        """ 
        ORIGX1      1.000000  0.000000  0.000000        0.00000                         
        ORIGX2      0.000000  1.000000  0.000000        0.00000                         
        ORIGX3      0.000000  0.000000  1.000000        0.00000  
        """
        o=[self.P.parsed[x] for x in ['ORIGX1','ORIGX2','ORIGX3']]
        self.assertEqual([o[0].O11,o[0].O12,o[0].O13,o[0].T1],[1.0,0.0,0.0,0.0])
        self.assertEqual([o[1].O21,o[1].O22,o[1].O23,o[1].T2],[0.0,1.0,0.0,0.0])
        self.assertEqual([o[2].O31,o[2].O32,o[2].O33,o[2].T3],[0.0,0.0,1.0,0.0])

    def test_scale123(self):
        """ 
        SCALE1      0.009330  0.005387  0.000000        0.00000                         
        SCALE2      0.000000  0.010773  0.000000        0.00000                         
        SCALE3      0.000000  0.000000  0.009703        0.00000  
        """
        s=[self.P.parsed[x] for x in ['SCALE1','SCALE2','SCALE3']]
        self.assertEqual([s[0].S11,s[0].S12,s[0].S13,s[0].U1],[0.00933,0.005387,0.0,0.0])
        self.assertEqual([s[1].S21,s[1].S22,s[1].S23,s[1].U2],[0.0,0.010773,0.0,0.0])
        self.assertEqual([s[2].S31,s[2].S32,s[2].S33,s[2].U3],[0.0,0.0,0.009703,0.0])

    def test_caveat(self):
        self.assertTrue('CAVEAT' not in self.P.parsed)

    def test_formul(self):
        rec=self.P.parsed['FORMUL']
        self.assertEqual(len(rec),3)
        """ 
        FORMUL   3  NAG    21(C8 H15 N O6)                                              
        FORMUL   3  BMA    C6 H12 O6                                                    
        FORMUL   3  MAN    3(C6 H12 O6)     
        """
        a=rec[0]
        self.assertEqual([a.compNum,a.hetID,a.asterisk,a.chemicalformula],[3,'NAG','','21(C8 H15 N O6)'])

    def test_hetnam(self):
        """ 
        HETNAM     NAG 2-ACETAMIDO-2-DEOXY-BETA-D-GLUCOPYRANOSE                         
        HETNAM     BMA BETA-D-MANNOPYRANOSE                                             
        HETNAM     MAN ALPHA-D-MANNOPYRANOSE   
        """
        rec=self.P.parsed['HETNAM']
        self.assertEqual(len(rec),3)
        a=rec[1]
        self.assertEqual([a.hetID,a.chemicalname],['BMA','BETA-D-MANNOPYRANOSE'])

    def test_hetsyn(self):
        self.assertTrue('HETSYN' not in self.P.parsed)

    def test_seqres(self):
        rec=self.P.parsed['SEQRES']
        self.assertEqual(len(rec),2)
        arec=rec[0]
        self.assertEqual(arec.chainID,'G')
        self.assertEqual(len(arec.resNames),arec.numRes)
        expected_seq='ALA GLU ASN LEU TRP VAL THR VAL TYR TYR GLY'.split()
        self.assertEqual(arec.resNames[:len(expected_seq)],expected_seq)
        expected_seq='CYS LYS ARG ARG VAL VAL GLY ARG ARG ARG ARG ARG ARG'.split()
        self.assertEqual(arec.resNames[-len(expected_seq):],expected_seq)
        arec=rec[1]
        self.assertEqual(arec.chainID,'B')
        self.assertEqual(len(arec.resNames),arec.numRes)
        expected_seq='ALA VAL GLY ILE GLY ALA VAL PHE LEU GLY PHE LEU GLY ALA ALA GLY SER THR MET GLY ALA ALA SER MET THR LEU THR VAL GLN ALA ARG ASN LEU LEU SER GLY ILE VAL GLN'.split()
        self.assertEqual(arec.resNames[:len(expected_seq)],expected_seq)

    def test_site(self):
        self.assertTrue('SITE' not in self.P.parsed)

    def test_endmdl(self):
        self.assertTrue('ENDMDL' not in self.P.parsed)

    def test_nummdl(self):
        self.assertTrue('NUMMDL' not in self.P.parsed)

    def test_model(self):
        self.assertTrue('MODEL' not in self.P.parsed)

    def test_end(self):
        rec=self.P.parsed['END']
        self.assertTrue(type(rec)==PDBRecord)

    def test_master(self):
        rec=self.P.parsed['MASTER']
        """ 
        MASTER      649    0   25   13   35    0    0    6 4856    2  378   49   
        """
        expres=[649,25,13,35,0,0,6,4856,2,378,49]
        res=[rec.numRemark,rec.numHet,rec.numHelix,rec.numSheet,rec.numTurn,rec.numSite,rec.numXform,rec.numCoord,rec.numTer,rec.numConect,rec.numSeq]
        self.assertEqual(expres,res)

    def test_ter(self):
        rec=self.P.parsed['TER']
        self.assertEqual(len(rec),2)
        ater=rec[0]
        self.assertEqual(ater.serial,3544)
        self.assertEqual(ater.residue.resName,'VAL')
        self.assertEqual(ater.residue.chainID,'G')
        self.assertEqual(ater.residue.seqNum,505)
        ater=rec[1]
        self.assertEqual(ater.serial,4520)
        self.assertEqual(ater.residue.resName,'ASP')
        self.assertEqual(ater.residue.chainID,'B')
        self.assertEqual(ater.residue.seqNum,664)

    def test_jrnl(self):
        rec=self.P.parsed['JRNL.AUTH']
        self.assertEqual(len(rec.authorList),53)
        self.assertEqual(rec.authorList[0],'Y.DO KWON')
        self.assertEqual(rec.authorList[1],'M.PANCERA')
        self.assertEqual(rec.authorList[-2],'J.R.MASCOLA')
        self.assertEqual(rec.authorList[-1],'P.D.KWONG')
        rec=self.P.parsed['JRNL.TITL']
        self.assertEqual(rec.title,'CRYSTAL STRUCTURE, CONFORMATIONAL FIXATION AND ENTRY-RELATED INTERACTIONS OF MATURE LIGAND-FREE HIV-1 ENV.')
        rec=self.P.parsed['JRNL.REF']
        self.assertEqual(rec.pubName,'NAT.STRUCT.MOL.BIOL.')
        self.assertEqual(rec.volume,'22')
        self.assertEqual(rec.page,'522')
        self.assertEqual(rec.year,2015)
        rec=self.P.parsed['JRNL.REFN']
        self.assertEqual(rec.issnORessn,'ESSN')
        self.assertEqual(rec.issn,'1545-9985')
        rec=self.P.parsed['JRNL.PMID']
        self.assertEqual(rec.pmid,26098315)
        rec=self.P.parsed['JRNL.DOI']
        self.assertEqual(rec.doi,'10.1038/NSMB.3051')

    def test_remark_0(self):
        self.assertTrue('REMARK.0' not in self.P.parsed)

    def test_remark_1(self):
        self.assertTrue('REMARK.1' not in self.P.parsed)

    def test_remark_290(self):
        self.assertTrue('REMARK.290' in self.P.parsed)
        self.assertTrue('REMARK.290.CRYSTSYMMTRANS' in self.P.parsed)
        rec=self.P.parsed['REMARK.290.CRYSTSYMMTRANS']
        print(rec.__dict__)
        Mlist,Tlist=get_symm_ops(rec)
        self.assertEqual(len(Mlist),6)
        self.assertTrue(np.array_equal(Mlist[0],np.identity(3)))
        self.assertTrue(np.array_equal(Mlist[1],
                                       np.array(
                                        [
                                            [-0.500000,-0.866025,0.000000],
                                            [ 0.866025,-0.500000,0.000000],
                                            [ 0.00000,  0.000000,1.000000]
                                        ])))
        self.assertTrue(np.array_equal(Mlist[-1],
                                       np.array(
                                        [
                                            [0.500000,-0.866025,0.00000],
                                            [0.866025,0.500000,0.000000],
                                            [ 0.00000,  0.000000,1.000000]
                                        ])))
        self.assertTrue(np.array_equal(Tlist[-1],np.array([0.0,0.0,51.53])))
    
    def test_remark_350(self):
        self.assertTrue('REMARK.350' in self.P.parsed)
        self.assertTrue(hasattr(self.P.parsed['REMARK.350'],'tokens'))
        rec=self.P.parsed['REMARK.350']
        self.assertEqual(rec.tokens['APPLY THE FOLLOWING TO CHAINS'],' G, B, A, C, D')
        self.assertTrue('REMARK.350.BIOMOLECULE.1' in self.P.parsed)
        rec=self.P.parsed['REMARK.350.BIOMOLECULE.1']
        Mlist,Tlist=get_symm_ops(rec)
        self.assertEqual(len(Mlist),3)
        self.assertTrue(np.array_equal(Mlist[0],np.identity(3)))
        self.assertTrue(np.array_equal(Mlist[1],
                                       np.array(
                                        [
                                            [-0.500000,-0.866025,0.000000],
                                            [ 0.866025,-0.500000,0.000000],
                                            [ 0.00000,  0.000000,1.000000]
                                        ])))
        self.assertTrue(np.array_equal(Tlist[1],np.array([107.18,185.64121,0.0])))

        """test_remark_290 
                          - rowName
                  - replNum
                  - m1
                  - m2
                  - m3
                  - t
        """

    # assert q.parsed['REMARK.2'].resolutionmsg==['', '3.31 ANGSTROMS.']
    # assert len(q.parsed['REMARK.3'].refinementdetails)==331
    # assert q.parsed['REMARK.3'].refinementdetails[0]==''
    # assert q.parsed['REMARK.3'].refinementdetails[3]=='AUTHORS     : PAUL ADAMS,PAVEL AFONINE,VINCENT CHEN,IAN'
    # assert q.parsed['REMARK.3'].refinementdetails[-1]=='OTHER REFINEMENT REMARKS: NULL'

    # assert q.parsed['REMARK.4'].formatmsg[1]=='4ZMJ COMPLIES WITH FORMAT V. 3.30, 13-JUL-11'
    # assert q.parsed['REMARK.100'].processmsg[-1]=='THE DEPOSITION ID IS D_1000209535.'
    # assert len(q.parsed['REMARK.200'].crystallizationmsg)==49
    # assert q.parsed['REMARK.200'].crystallizationmsg[-1]=='REMARK: NULL'
    # assert len(q.parsed['REMARK.290'].crystallographicautodetails)==40
    # assert q.parsed['REMARK.290'].crystallographicautodetails[-1]=='REMARK: NULL'

    # # print(q.parsed['REMARK.300'].tokengroups['biomoleculedeclarations'])

    # assert q.parsed['REMARK.300'].tokengroups['biomoleculedeclarations']['BIOMOLECULE'].BIOMOLECULE==['1']
    # assert q.parsed['REMARK.350'].tokengroups['biomoleculespecifications']['BIOMOLECULE.1'].AUTH_BIO_UNIT=='HEXAMERIC'
    # assert q.parsed['REMARK.350'].tokengroups['biomoleculespecifications']['BIOMOLECULE.1'].SOFT_QUAT_STRUCT=='HEXAMERIC'
    # assert q.parsed['REMARK.350'].tokengroups['biomoleculespecifications']['BIOMOLECULE.1'].CHAIN_BIOMT==['G', 'B', 'A', 'C', 'D']
    # # for k in q.parsed.keys():
    # #     print(k)
    # # print(q.parsed['REMARK.350'])
    # # print(q.parsed['REMARK.350.BIOMT'])
    # assert q.parsed['REMARK.350.BIOMT'].rowNum==[1,2,3,1,2,3,1,2,3]
    # assert q.parsed['REMARK.350.BIOMT'].transNum==[1,1,1,2,2,2,3,3,3]
    # assert q.parsed['REMARK.350.BIOMT'].M1==[1.0,0.0,0.0,-0.5,0.866025,0.0,-0.5,-0.866025,0.0]
    # rec=q.parsed['REMARK.350.BIOMT']
    # M=[]
    # for i in range(3):
    #     M.append(np.zeros((3,3)))
    # for t,r,m1,m2,m3 in zip(rec.transNum,rec.rowNum,rec.M1,rec.M2,rec.M3):
    #     M[t-1][r-1,0]=m1
    #     M[t-1][r-1,1]=m2
    #     M[t-1][r-1,2]=m3
    
    # expM=[
    #     np.array([
    #         [1.0,0.0,0.0],
    #         [0.0,1.0,0.0],
    #         [0.0,0.0,1.0]
    #     ]),
    #     np.array([
    #         [-0.5,-0.866025,0.0],
    #         [0.866025,-0.5,0.0],
    #         [0.0,0.0,1.0]
    #     ]),
    #     np.array([
    #         [-0.5,0.866025,0.0],
    #         [-0.866025,-0.5,0.0],
    #         [0.0,0.0,1.0]
    #     ])       
    # ]
    # assert all(np.array_equal(m,em) for m,em in zip(M,expM))

    # assert q.parsed['REMARK.465'].freetext[1].strip()=='REMARK 465 MISSING RESIDUES'
    # assert q.parsed['REMARK.465'].tables['MISSING'][0].resname=='ALA'
    # assert q.parsed['REMARK.465'].tables['MISSING'][0].modelNum==''
    # assert q.parsed['REMARK.465'].tables['MISSING'][0].chainID=='G'
    # assert q.parsed['REMARK.465'].tables['MISSING'][0].resseqnum==31
    # assert len(q.parsed['REMARK.465'].tables['MISSING'])==61

    # assert q.parsed['REMARK.500'].tables['RAMA_OUTLIERS'][0].residue.resName=='PRO'
    # assert q.parsed['REMARK.500'].tables['NONCISTRANS'][0].residueC.resName=='ALA'
    # assert q.parsed['REMARK.500'].tables['NONCISTRANS'][0].omega==146.93
    # assert q.parsed['REMARK.500'].tables['NONCISTRANS'][-1].omega==149.23

#REMARK 500 LEU G  494     GLY G  495                  149.23                    

# def test_rem1_ref():
#     p=PDBParser(PDBcode='test2').parse()

#     assert p.parsed['JRNL.AUTH'].authorList[-1]=='P.D.KWONG'
#     assert p.parsed['REMARK.1.1.AUTH'].authorList==['J.N.BREG', 'J.H.J.VAN  OPHEUSDEN', 'M.J.M.BURGERING','R.BOELENS','R.KAPTEIN']
#     assert p.parsed['REMARK.1.1.TITL'].title=='STRUCTURE OF ARC REPRESSOR  IN SOLUTION: EVIDENCE FOR A FAMILY OF B-SHEET DNA-BINDING PROTEIN'
#     assert p.parsed['REMARK.1.1.REF'].pubName=='NATURE'
#     assert p.parsed['REMARK.1.2.AUTH'].authorList==['J.N.BREG', 'R.BOELENS','A.V.E.GEORGE','R.KAPTEIN']
#     assert p.parsed['REMARK.1.2.REF'].pubName=='BIOCHEMISTRY'
#     assert p.parsed['REMARK.0.1.AUTH'].authorList==['I.P.FREELY', 'R.U.SEERIUS']
#     assert p.parsed['REMARK.0.1'].tokengroups['tokencheck']['PDB_ID'].PDB_ID=='POOP'
