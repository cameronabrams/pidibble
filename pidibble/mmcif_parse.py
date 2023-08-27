from .pdbrecord import PDBRecord
from .baserecord import BaseRecord
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

class MMCIF_Parser:
    def __init__(self,mmcif_formats,pdb_formats):
        self.formats=mmcif_formats
        self.pdb_formats=pdb_formats
        self.global_maps={}
        self.global_ids={}

    def update_maps(self,maps,cifrec,idx):
        for mapname,mapspec in maps.items():
            if not mapname in self.global_maps:
                self.global_maps[mapname]={}
            k=mapspec['key']
            v=mapspec['value']
            key=rectify(cifrec.getValue(k,idx))
            val=rectify(cifrec.getValue(v,idx))
            if not key in self.global_maps[mapname]:
                self.global_maps[mapname][key]=val

    def update_ids(self,idmaps,cifrec,idx):
        for idname,idspec in idmaps.items():
            if not idname in self.global_ids:
                self.global_ids[idname]=[]
            thisid=rectify(cifrec.getValue(idspec,idx))
            if not thisid in self.global_ids[idname]:
                self.global_ids[idname].append(thisid)

    def gen_dict(self,cifrec,idx,**kwargs):
        idicts=[]
        attr_map=kwargs.get('attr_map',{})
        splits=kwargs.get('splits',[])
        indexes=kwargs.get('indexes',None)
        spawns_on=kwargs.get('spawns_on',None)
        map_values=kwargs.get('map_values',{})
        spawn_data=kwargs.get('spawn_data',{})
        tables=kwargs.get('tables',{})
        list_attr=kwargs.get('list_attr',{})
        idict={}
        for k,v in attr_map.items():
            if type(v)==dict:
                val=PDBRecord(self.gen_dict(v,cifrec,idx))
            elif type(v)==list:
                val=[rectify(cifrec.getValue(x,idx)) for x in v]
            else:
                val=rectify(cifrec.getValue(v,idx))
                if k=='resseqnumi':
                    seqNum,iCode=split_ri(val)
                    val={'seqNum':int(seqNum),'iCode':iCode}
                else:
                    if k in splits and ',' in val:
                        val=[rectify(x) for x in val.split(',')]
                    if k==spawns_on and ',' in val:
                        val=[rectify(x) for x in val.split(',')]
                    if k in map_values:
                        mapper=self.global_maps[map_values[k]]
                        if type(val)==list:
                            val=[mapper(x) for x in val]
                        else:
                            val=mapper(val)
            idict[k]=val
        return [idict]

    def parse(self,data):
        recdict={}
        for rectype,mapspec in self.formats.items():
            rectypeparts=rectype.split('.')
            baserectype=rectypeparts[0]
            pdb_format=self.pdb_formats[baserectype]
            if len(rectypeparts)>1:
                subrectype=rectypeparts[1]
                if subrectype.isdigit():
                    subrectype=int(subrectype)
                subrecfmt=pdb_format['subrecords']['formats'][subrectype]
                if len(rectypeparts)>2:
                    embedtype=rectypeparts[2]
                    embedfmt=subrecfmt['embedded_records'][embedtype]
            cifrec=data.getObj(mapspec['data_obj'])
            sigattr=mapspec.get('signal_attr',None)
            sigval=mapspec.get('signal_value',None)
            global_maps=mapspec.get('global_maps',{})
            global_ids=mapspec.get('global_ids',{})
            splits=mapspec.get('splits',[])
            spawns_on=mapspec.get('spawns_on',None)
            indexes=mapspec.get('indexes',None)
            map_values=mapspec.get('map_values',{})
            tables=mapspec.get('tables',{})
            use_signal=(sigattr!=None)
            attr_map=mapspec.get('attr_map',{})
            for i in range(len(cifrec)):
                if self.global_maps:
                    self.update_maps(global_maps,cifrec,i)
                if self.global_ids:
                    self.update_ids(global_ids,cifrec,i)
                if not use_signal or (cifrec.getValue(sigattr,i)==sigval):
                    idicts=self.gen_dict(cifrec,i,attr_map=attr_map,splits=splits,spawns_on=spawns_on,map_values=map_values,indexes=indexes,tables=tables)
                    for idict in idicts:
                        recdict[rectype].append(PDBRecord(idict))

            cifndata=len(cifrec)
            if cifndata==1:
                idict={}
                for k,v in attr_map.items():
                    val=rectify(cifrec.getValue(v,0))
                    if k in splits and ',' in val:
                        val=[rectify(x) for x in val.split(',')]
                    if k==spawns_on and ',' in k:
                        val=[rectify(x) for x in val.split(',')]
                    if k in map_values:
                        mapper=self.global_maps[map_values[k]]
                        if type(val)==list:
                            val=[mapper(x) for x in val]
                        else:
                            val=mapper(val)
                    idict[k]=val
                if spawns_on and type(idict[spawns_on])==list:
                    if not rectype in recdict:
                        recdict[rectype]=[]
                    for val in idict[spawns_on]:
                        jdict=idict.copy()
                        jdict[spawns_on]=val


                if indexes:
                    rectype=f'{rectype}{idict[indexes]}'
                recdict[rectype]=PDBRecord(idict)                    
            else:
                if not rectype in recdict:
                    recdict[rectype]=[]
                if not tables:
                else:
                    tabledict={}
                    for tname,tspec in tables.items():
                        tabledict[tname]=[]
                        attr_map=tspec['row_attr_map']
                        bisv=tspec.get('blank_if_single_valued',[])
                        for i in range(cifndata):
                            idict={}
                            for k,v in attr_map.items():
                                idict[k]=rectify(cifrec.getValue(v,i))
                                if k in bisv:
                                    if len(self.global_ids[k])<2:
                                        idict[k]=''
                            tabledict[tname].append(BaseRecord(idict))
                    udict={'tables':tabledict}
                    recdict[rectype]=PDBRecord(udict)
        return recdict
