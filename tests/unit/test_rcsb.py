import numpy as np
from pidibble.pdbparse import PDBParser, PDBRecord, get_symm_ops
import unittest
from itertools import product
import logging
logger=logging.getLogger(__name__)
def test_pdbformat():
    p=PDBParser()
    expected_sections=['record_formats','BASE_URL']
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

class Test4zmj(unittest.TestCase):

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
        self.assertEqual({'MOL_ID.1': ['G'], 'MOL_ID.2': ['B']},rec.get_token('CHAIN'))

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
        for i in range(1,7):
            self.assertTrue(f'REMARK.290.CRYSTSYMMTRANS.{i}' in self.P.parsed)
        rec=self.P.parsed['REMARK.290.CRYSTSYMMTRANS.1']
        logger.debug(rec.pstr())
        
        # print(rec.__dict__)
        # Mlist,Tlist=get_symm_ops(rec)
        # self.assertEqual(len(Mlist),6)
        # self.assertTrue(np.array_equal(Mlist[0],np.identity(3)))
        # self.assertTrue(np.array_equal(Mlist[1],
        #                                np.array(
        #                                 [
        #                                     [-0.500000,-0.866025,0.000000],
        #                                     [ 0.866025,-0.500000,0.000000],
        #                                     [ 0.00000,  0.000000,1.000000]
        #                                 ])))
        # self.assertTrue(np.array_equal(Mlist[-1],
        #                                np.array(
        #                                 [
        #                                     [0.500000,-0.866025,0.00000],
        #                                     [0.866025,0.500000,0.000000],
        #                                     [ 0.00000,  0.000000,1.000000]
        #                                 ])))
        # self.assertTrue(np.array_equal(Tlist[-1],np.array([0.0,0.0,51.53])))
    
    def test_remark_350(self):
        self.assertTrue('REMARK.350' in self.P.parsed)
        # self.assertTrue(hasattr(self.P.parsed['REMARK.350'],'tokens'))
        rec=self.P.parsed['REMARK.350']
        # self.assertEqual(rec.tokens['APPLY THE FOLLOWING TO CHAINS'],' G, B, A, C, D')
        # self.assertTrue('REMARK.350.BIOMOLECULE1' in self.P.parsed)
        self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM1' in self.P.parsed)
        self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM2' in self.P.parsed)
        self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM3' in self.P.parsed)
        rec=self.P.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM1']
        M,T=get_symm_ops(rec)
        self.assertTrue(np.array_equal(M,np.identity(3)))
        rec=self.P.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
        M,T=get_symm_ops(rec)
        self.assertTrue(np.array_equal(M,
                                       np.array(
                                        [
                                            [-0.500000,-0.866025,0.000000],
                                            [ 0.866025,-0.500000,0.000000],
                                            [ 0.00000,  0.000000,1.000000]
                                        ])))
        self.assertTrue(np.array_equal(T,np.array([107.18,185.64121,0.0])))

    def test_remark_465(self):
        rec=self.P.parsed['REMARK.465']
        self.assertTrue(hasattr(rec,'tables'))
        self.assertTrue('MISSING' in rec.tables)
        missing=rec.tables['MISSING']

        self.assertEqual(len(missing),61)
        m=missing[0]
        self.assertEqual(m.modelNum,'')
        self.assertEqual(m.resName,'ALA')
        self.assertEqual(m.chainID,'G')
        self.assertEqual(m.seqNum,31)
        m=missing[-1]
        self.assertEqual(m.modelNum,'')
        self.assertEqual(m.resName,'LEU')
        self.assertEqual(m.chainID,'B')
        self.assertEqual(m.seqNum,568)
        
    def test_remark_500(self):
        rec=self.P.parsed['REMARK.500']
        self.assertTrue(hasattr(rec,'tables'))
        tables=rec.tables
        self.assertTrue('CLOSE_CONTACTS_ASYMM_UNIT' in tables)
        self.assertTrue('CLOSE_CONTACTS' in tables)
        self.assertTrue('COVALENT_BOND_ANGLES' in tables)
        self.assertTrue('RAMA_OUTLIERS' in tables)
        self.assertTrue('NONCISTRANS' in tables)

        t=tables['CLOSE_CONTACTS_ASYMM_UNIT']
        self.assertEqual(len(t),2)
        r=t[0]
        self.assertEqual(r.atom1,'O')
        self.assertEqual(r.residue1.resName,'VAL')
        self.assertEqual(r.residue1.chainID,'G')
        self.assertEqual(r.residue1.seqNum,36)
        self.assertEqual(r.residue1.iCode,'')
        self.assertEqual(r.atom2,'OG1')
        self.assertEqual(r.residue2.resName,'THR')
        self.assertEqual(r.residue2.chainID,'B')
        self.assertEqual(r.residue2.seqNum,606)
        self.assertEqual(r.residue2.iCode,'')
        self.assertEqual(r.distance,2.15)
        t=tables['CLOSE_CONTACTS']
        self.assertEqual(len(t),2)
        r=t[0]
        self.assertEqual(r.atom1,'ND2')
        self.assertEqual(r.residue1.resName,'ASN')
        self.assertEqual(r.residue1.chainID,'G')
        self.assertEqual(r.residue1.seqNum,136)
        self.assertEqual(r.residue1.iCode,'')
        self.assertEqual(r.atom2,'O3')
        self.assertEqual(r.residue2.resName,'NAG')
        self.assertEqual(r.residue2.chainID,'G')
        self.assertEqual(r.residue2.seqNum,610)
        self.assertEqual(r.residue2.iCode,'')
        self.assertEqual(r.ssymop,'5564')
        self.assertEqual(r.distance,2.16)
        t=tables['COVALENT_BOND_ANGLES']
        self.assertEqual(len(t),2)
        r=t[0]
        self.assertEqual(r.modelNum,'')
        self.assertEqual(r.residue.resName,'PRO')
        self.assertEqual(r.residue.chainID,'G')
        self.assertEqual(r.residue.seqNum,79)
        self.assertEqual(r.residue.iCode,'')
        self.assertEqual(r.atom1,'C')
        self.assertEqual(r.atom2,'N')
        self.assertEqual(r.atom3,'CA')
        self.assertEqual(r.deviation,9.7)
        self.assertEqual(r.units,'DEGREES')
        r=t[1]
        self.assertEqual(r.modelNum,'')
        self.assertEqual(r.residue.resName,'PRO')
        self.assertEqual(r.residue.chainID,'G')
        self.assertEqual(r.residue.seqNum,214)
        self.assertEqual(r.residue.iCode,'')
        self.assertEqual(r.atom1,'C')
        self.assertEqual(r.atom2,'N')
        self.assertEqual(r.atom3,'CA')
        self.assertEqual(r.deviation,-9.1)
        self.assertEqual(r.units,'DEGREES')

        
        t=tables['RAMA_OUTLIERS']
        self.assertEqual(len(t),33)
        r=t[0]
        self.assertEqual(r.modelNum,'')
        self.assertEqual(r.residue.resName,'PRO')
        self.assertEqual(r.residue.chainID,'G')
        self.assertEqual(r.residue.seqNum,43)
        self.assertEqual(r.residue.iCode,'')
        self.assertEqual(r.phi,66.5)
        self.assertEqual(r.psi,-66.28)
        r=t[-1]
        self.assertEqual(r.modelNum,'')
        self.assertEqual(r.residue.resName,'GLN')
        self.assertEqual(r.residue.chainID,'B')
        self.assertEqual(r.residue.seqNum,650)
        self.assertEqual(r.residue.iCode,'')
        self.assertEqual(r.phi,-74.39)
        self.assertEqual(r.psi,-111.72)

        t=tables['NONCISTRANS']
        self.assertEqual(len(t),4)
        r=t[0]
        self.assertEqual(r.residueN.resName,'ASP')
        self.assertEqual(r.residueN.chainID,'G')
        self.assertEqual(r.residueN.seqNum,57)
        self.assertEqual(r.residueN.iCode,'')
        self.assertEqual(r.residueC.resName,'ALA')
        self.assertEqual(r.residueC.chainID,'G')
        self.assertEqual(r.residueC.seqNum,58)
        self.assertEqual(r.residueC.iCode,'')
        self.assertEqual(r.modelNum,'')
        self.assertEqual(r.omega,146.93)        
        r=t[1]
        self.assertEqual(r.residueN.resName,'GLY')
        self.assertEqual(r.residueN.chainID,'G')
        self.assertEqual(r.residueN.seqNum,354)
        self.assertEqual(r.residueN.iCode,'')
        self.assertEqual(r.residueC.resName,'ASN')
        self.assertEqual(r.residueC.chainID,'G')
        self.assertEqual(r.residueC.seqNum,355)
        self.assertEqual(r.residueC.iCode,'')
        self.assertEqual(r.modelNum,'')
        self.assertEqual(r.omega,-134.45)        

class Test_mmCIF(unittest.TestCase):
    def setUp(self):
        self.pdb=PDBParser(PDBcode='4tvp',input_format='PDB').parse().parsed
        self.mmCIF=PDBParser(input_format='mmCIF',PDBcode='4tvp').parse().parsed

    def test_cif_pdb_correspondence_atoms(self):
        # this test verifies correct parsing of mmCIF in the case where there is an exact correspondence
        # between the mmCIF file and the PDB file
        mmCIF=self.mmCIF
        pdb=self.pdb
        self.assertEqual(len(mmCIF['ATOM']),len(pdb['ATOM']))

        for ats in ['ATOM','HETATM']:
            # print(ats)
            mmCIF_atoms=mmCIF[ats]
            pdb_atoms=pdb[ats]
            self.assertEqual(len(mmCIF_atoms),len(pdb_atoms))
            for ma,pa in zip(mmCIF_atoms,pdb_atoms):
                # serial numbers in PDB files may not match those in CIF files because
                # PDB files may have TER records that have a serial number just like an atom
                # self.assertEqual(ma.serial,pa.serial)
                self.assertEqual(ma.name,pa.name)
                self.assertEqual(ma.altLoc,pa.altLoc)
                logger.debug(f'{ma.residue.chainID} {ma.residue.resName} {ma.residue.seqNum} {ma.name} == {pa.residue.chainID} {pa.residue.resName} {pa.residue.seqNum} {pa.name} ')
                self.assertEqual(ma.residue.chainID,pa.residue.chainID)
                self.assertEqual(ma.residue.resName,pa.residue.resName)
                self.assertEqual(ma.residue.seqNum,pa.residue.seqNum)
                self.assertEqual(ma.residue.iCode,pa.residue.iCode)
                self.assertEqual(ma.x,pa.x)
                self.assertEqual(ma.y,pa.y)
                self.assertEqual(ma.z,pa.z)
                self.assertEqual(ma.occupancy,pa.occupancy)
                self.assertEqual(ma.tempFactor,pa.tempFactor)
                self.assertEqual(ma.element,pa.element)
                self.assertEqual(ma.charge,pa.charge)
        # plinks=p.parsed['LINK']
        # olinks=o.parsed['LINK']
        # self.assertEqual(len(plinks),len(olinks))
        # # do they have the same order though?
        # for i in range(len(plinks)):
        #     pn_link=plinks[i]
        #     on_link=olinks[i]
        #     self.assertEqual(pn_link.name1,on_link.name1)
        #     self.assertEqual(pn_link.altLoc1,on_link.altLoc1)
        #     self.assertEqual(pn_link.residue1.chainID,on_link.residue1.chainID)
        #     self.assertEqual(pn_link.residue1.resName,on_link.residue1.resName)
        #     self.assertEqual(pn_link.residue1.seqNum,on_link.residue1.seqNum)
        #     self.assertEqual(pn_link.residue1.iCode,on_link.residue1.iCode)
        #     self.assertEqual(pn_link.name2,on_link.name2)
        #     self.assertEqual(pn_link.altLoc2,on_link.altLoc2)
        #     self.assertEqual(pn_link.residue2.chainID,on_link.residue2.chainID)
        #     self.assertEqual(pn_link.residue2.resName,on_link.residue2.resName)
        #     self.assertEqual(pn_link.residue2.seqNum,on_link.residue2.seqNum)
        #     self.assertEqual(pn_link.residue2.iCode,on_link.residue2.iCode)

        # pssbonds=p.parsed['SSBOND']
        # ossbonds=o.parsed['SSBOND']
        # self.assertEqual(len(pssbonds),len(ossbonds))
        # for i in range(len(pssbonds)):
        #     p_ss=pssbonds[i]
        #     o_ss=ossbonds[i]
        #     self.assertEqual(p_ss.residue1.chainID,o_ss.residue1.chainID)
        #     self.assertEqual(p_ss.residue1.resName,o_ss.residue1.resName)
        #     self.assertEqual(p_ss.residue1.seqNum,o_ss.residue1.seqNum)
        #     self.assertEqual(p_ss.residue1.iCode,o_ss.residue1.iCode)

        # p_missing=p.parsed['REMARK.465'].tables['MISSING']
        # o_missing=o.parsed['REMARK.465'].tables['MISSING']
        # self.assertEqual(len(p_missing),len(o_missing))
        # for i in range(len(p_missing)):
        #     pn_missing=p_missing[i]
        #     on_missing=o_missing[i]
        #     self.assertEqual(pn_missing.modelNum,on_missing.modelNum)
        #     self.assertEqual(pn_missing.resName,on_missing.resName)
        #     self.assertEqual(pn_missing.chainID,on_missing.chainID)
        #     self.assertEqual(pn_missing.seqNum,on_missing.seqNum)
        #     self.assertEqual(pn_missing.iCode,on_missing.iCode)

        # # for k in p.parsed.keys():
        # #     print(k)
        # self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM1' in p.parsed)
        # self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM1' in o.parsed)
        # prec=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM1']
        # orec=o.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM1']
        # pM,pT=get_symm_ops(prec)
        # oM,oT=get_symm_ops(orec)
        # self.assertTrue(np.allclose(pM,oM))
        # self.assertTrue(np.allclose(pT,oT))

        # self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM2' in p.parsed)
        # self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM2' in o.parsed)
        # prec=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
        # orec=o.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
        # pM,pT=get_symm_ops(prec)
        # oM,oT=get_symm_ops(orec)
        # # print(f'cif {pM} {pT}')
        # # print(f'pdb {oM} {oT}')
        # self.assertTrue(np.allclose(pM,oM))
        # self.assertTrue(np.allclose(pT,oT))
        
        # self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM3' in p.parsed)
        # self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM3' in o.parsed)
        # prec=p.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM3']
        # orec=o.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM3']
        # pM,pT=get_symm_ops(prec)
        # oM,oT=get_symm_ops(orec)
        # self.assertTrue(np.allclose(pM,oM))
        # self.assertTrue(np.allclose(pT,oT))

        # self.assertFalse('REMARK.350.BIOMOLECULE1.TRANSFORM4' in p.parsed)
        
        # prec=p.parsed['SEQADV']
        # orec=o.parsed['SEQADV']
        # self.assertEqual(len(prec),len(orec))
        # for i in range(len(prec)):
        #     psa=prec[i]
        #     osa=orec[i]
        #     self.assertEqual(psa.idCode,osa.idCode)
        #     self.assertEqual(psa.residue.chainID,osa.residue.chainID)
        #     self.assertEqual(psa.residue.resName,osa.residue.resName)
        #     self.assertEqual(psa.residue.seqNum,osa.residue.seqNum)
        #     self.assertEqual(psa.residue.iCode,osa.residue.iCode)
        #     self.assertEqual(psa.database,osa.database)
        #     self.assertEqual(psa.dbAccession,osa.dbAccession)
        #     self.assertEqual(psa.dbRes,osa.dbRes)
        #     self.assertEqual(psa.dbSeq,osa.dbSeq)
        #     self.assertEqual(psa.conflict,osa.conflict)

        # self.assertTrue(False)

    def test_cif_fetch(self):
        p=PDBParser(input_format='mmCIF',PDBcode='8fae').parse().parsed
        self.assertEqual(len(p['ATOM']),14466)
        self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM1' in p)
        rec=p['REMARK.350.BIOMOLECULE1.TRANSFORM1']
        digits=[chr(ord('0')+x) for x in range(10)]
        L=[chr(ord('A')+x) for x in range(26)]
        l=[chr(ord('a')+x) for x in range(26)]
        DL=[''.join(x[::-1]) for x in product(L,L)]
        nc=len(rec.header)
        ML=L+l+digits+DL
        TL=ML[:nc]
        bad_idx=TL.index('q')
        TL[bad_idx]='OA'
        TL.sort()
        self.assertEqual(rec.header,TL)

class TestFourLetterResNames(unittest.TestCase):
    def test_four(self):
        P=PDBParser(PDBcode='4zmj-newresnames').parse()
        atoms=P.parsed['ATOM']
        an_atom=atoms[9590]
        self.assertEqual(an_atom.residue.resName,'BGNA')
    def test_6m0j(self):
        P=PDBParser(PDBcode='6m0j').parse()
        hets=P.parsed['HETATM']
        a=hets[0]
        self.assertEqual(a.name,'ZN')
        self.assertEqual(a.residue.resName,'ZN')

class Test4tvp(unittest.TestCase):
    def setUp(self) -> None:
        self.P=PDBParser(PDBcode='4tvp').parse()
        return super().setUp()
    
    def test_title(self):
        p=self.P.parsed
        self.assertEqual(p['TITLE'].title,'CRYSTAL STRUCTURE OF THE HIV-1 BG505 SOSIP.664 ENV TRIMER ECTODOMAIN, COMPRISING ATOMIC-LEVEL DEFINITION OF PRE-FUSION GP120 AND GP41, IN COMPLEX WITH HUMAN ANTIBODIES PGT122 AND 35O22')

    def test_remark_350(self):
        self.assertTrue('REMARK.350' in self.P.parsed)
        rec=self.P.parsed['REMARK.350']
        self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM1' in self.P.parsed)
        rec=self.P.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM1']
        self.assertEqual(rec.header,['G', 'B', 'L', 'H', 'D', 'E', 'A', 'C', 'F', 'I', 'J', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T'])
        M,T=get_symm_ops(rec)
        self.assertTrue(np.array_equal(T,np.array([0.0,0.0,0.0])))
        self.assertTrue(np.array_equal(M,np.identity(3)))
        self.assertTrue('REMARK.350.BIOMOLECULE1.TRANSFORM2' in self.P.parsed)
        rec=self.P.parsed['REMARK.350.BIOMOLECULE1.TRANSFORM2']
        self.assertEqual(rec.header,['G', 'B', 'L', 'H', 'D', 'E', 'A', 'C', 'F', 'I', 'J', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T'])
        M,T=get_symm_ops(rec)
        self.assertTrue(np.array_equal(M,
                                       np.array(
                                        [
                                            [-0.500000,-0.866025,0.000000],
                                            [ 0.866025,-0.500000,0.000000],
                                            [ 0.00000,  0.000000,1.000000]
                                        ])))
        self.assertTrue(np.array_equal(T,np.array([-515.56,0.0,0.0])))

    def test_remark_375(self):
        self.assertTrue('REMARK.375' in self.P.parsed)
        rec=self.P.parsed['REMARK.375']
        self.assertTrue(hasattr(rec,'tables'))
        tables=rec.tables
        self.assertTrue('SPECIAL_POSITIONS' in tables)
        t=tables['SPECIAL_POSITIONS']
        self.assertEqual(len(t),4)
        r=t[0]
        self.assertEqual(r.atomname,'S')
        self.assertEqual(r.residue.resName,'SO4')
        self.assertEqual(r.residue.chainID,'G')
        self.assertEqual(r.residue.seqNum,606)
        self.assertEqual(r.residue.iCode,'')
    
    def test_remark_650(self):
        self.assertTrue('REMARK.650' in self.P.parsed)
        rec=self.P.parsed['REMARK.650']
        self.assertTrue(hasattr(rec,'tokengroups'))
        tg=rec.tokengroups
        self.assertTrue('freetext' in tg)
        t=tg['freetext']
        d=rec.tokengroups['freetext']['HELIX'].HELIX
        self.assertEqual(d,'AUTHOR DETERMINED')

    def test_remark_700(self):
        self.assertTrue('REMARK.700' in self.P.parsed)
        rec=self.P.parsed['REMARK.700']
        self.assertTrue(hasattr(rec,'tokengroups'))
        tg=rec.tokengroups
        self.assertTrue('freetext' in tg)
        t=tg['freetext']
        d=rec.tokengroups['freetext']['SHEET'].SHEET
        self.assertEqual(d,'AUTHOR DETERMINED')
        self.assertEqual(rec.get_token('SHEET'),'AUTHOR DETERMINED')
    def test_remark_280(self):
        self.assertTrue('REMARK.280' in self.P.parsed)
        rec=self.P.parsed['REMARK.280']
        self.assertTrue(hasattr(rec,'tokengroups'))
        tg=rec.tokengroups
        self.assertTrue('freetext' in tg)
        t=tg['freetext']
        d=rec.tokengroups['freetext']['SOLV_CONT'].SOLV_CONT
        self.assertEqual(d,72.34)
        d=rec.tokengroups['freetext']['MATT_COEF'].MATT_COEF
        self.assertEqual(d,4.45)
        v=rec.get_token('MATT_COEF')
        self.assertEqual(v,4.45)
        v=rec.get_token('CRYS_COND')
        self.assertEqual(v,'16% ISOPROPANOL, 5.32% PEG1500, 0.2M LISO4, 0.1M SODIUM ACETATE PH 5.5, VAPOR DIFFUSION, HANGING DROP, TEMPERATURE 293K')