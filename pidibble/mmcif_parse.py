from .pdbrecord import PDBRecord
import logging
logger=logging.getLogger(__name__)
def split_ri(ri):
    if ri[-1].isdigit(): # there is no insertion code
        r=int(ri)
        i=''
    else:
        r=int(ri[:-1])
        i=ri[-1]
    return r,i

def rectify(val):
    if not val:
        return ''
    if val in '.?':
        return ''
    if val.isdigit():
        return int(val)
    try:
        val=float(val)
    except:
        pass
    return val

def mmCIF_parser(data,pdb_formats,mmcif_formats):
    recdict={}
    for rectype,mapspec in mmcif_formats.items():
        rectypeparts=rectype.split('.')
        baserectype=rectypeparts[0]
        pdb_format=pdb_formats[baserectype]
        if len(rectypeparts)>1:
            subrectype=rectypeparts[1]
            if subrectype.isdigit():
                subrectype=int(subrectype)
            subrecfmt=pdb_format['subrecords']['formats'][subrectype]
            if len(rectypeparts)>2:
                embedtype=rectypeparts[2]
                embedfmt=subrecfmt['embedded_records'][embedtype]
        cifrec=data.getObj(mapspec['data_obj'])
        cifattr=cifrec.getAttributeList()
        sigattr=mapspec.get('signal_attr',None)
        sigval=mapspec.get('signal_value',None)
        use_signal=(sigattr!=None)
        attr_map=mapspec.get('attr_map',{})
        cifndata=len(cifrec)
        if cifndata==1:
            print(f'not parsing {rectype}')
            pass
        else:
            if not rectype in recdict:
                recdict[rectype]=[]
            tables=mapspec.get('tables',[])
            if not tables:
                for i in range(cifndata):
                    if not use_signal or (cifrec.getValue(sigattr,i)==sigval):
                        idict={}
                        for k,v in attr_map.items():
                            if type(v)!=dict:
                                idict[k]=rectify(cifrec.getValue(v,i))
                            else:
                                # print(f'val is dict: {v}')
                                udict={}
                                for kk,vv in v.items():
                                    if kk=='resseqnumi':
                                        ri=cifrec.getValue(vv,i)
                                        seqNum,iCode=split_ri(ri)
                                        udict['seqNum']=int(seqNum)
                                        udict['iCode']=iCode
                                    else:
                                        udict[kk]=rectify(cifrec.getValue(vv,i))
                                idict[k]=PDBRecord(udict)
                        recdict[rectype].append(PDBRecord(idict))
            else:
                for T in tables:
                    pass
                # TODO parse table!
    return recdict
