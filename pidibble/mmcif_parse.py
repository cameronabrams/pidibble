from collections import UserDict
from .pdbrecord import PDBRecord
from .baserecord import BaseRecord
import logging
logger=logging.getLogger(__name__)

def split_ri(ri):
    if type(ri)==int: # this is no insertion code
        r=ri
        i=''
    elif ri[-1].isdigit(): # there is no insertion code
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

def resolve(key,aDict):
    pass

class MMCIFDict(UserDict):
    def __init__(self,data,linkers={},blankers=[' ','','?']):
        self.data=data
        self.linkers=linkers
        self.blankers=blankers
    
    def get(self,key):
        val=self[key]
        if val in self.blankers:
            return ''
        
        key_link=self.linkers.get(val,None)
        if key_link:
            if key_link in self.keys():
                val=self[key_link]
        return val


class MMCIF_Parser:
    def __init__(self,mmcif_formats,pdb_formats,cif_data):
        self.formats=mmcif_formats
        self.pdb_formats=pdb_formats
        self.global_maps={}
        self.global_ids={}
        self.cif_data=cif_data

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

    def gen_dict(self,mapspec):
        idicts=[]
        attr_map=mapspec.get('attr_map',{})
        splits=mapspec.get('splits',[])
        spawns_on=mapspec.get('spawns_on',None)
        indexes=mapspec.get('indexes',None)
        map_values=mapspec.get('map_values',{})
        tables=mapspec.get('tables',{})
        spawn_data=mapspec.get('spawn_data',{})
        tables=mapspec.get('tables',{})
        list_attr=mapspec.get('list_attr',{})
        sigattr=mapspec.get('signal_attr',None)
        sigval=mapspec.get('signal_value',None)
        use_signal=(sigattr!=None)
        global_maps=mapspec.get('global_maps',{})
        global_ids=mapspec.get('global_ids',{})
        spawns_on=mapspec.get('spawns_on',None)
        allcaps=mapspec.get('allcaps',[])
        if_dot_replace_with=mapspec.get('if_dot_replace_with',{})
        cifrec=self.cif_data.getObj(mapspec['data_obj'])
        if not tables:
            for idx in range(len(cifrec)):
                if not use_signal or (cifrec.getValue(sigattr,idx)==sigval):
                    if global_maps:
                        self.update_maps(global_maps,cifrec,idx)
                    if global_ids:
                        self.update_ids(global_ids,cifrec,idx)
                    idict={}
                    for k,v in attr_map.items():
                        if type(v)==dict:
                            resdict={kk:rectify(cifrec.getValue(o,idx)) for kk,o in v.items()}
                            if 'resseqnumi' in resdict:
                                resdict['seqNum'],resdict['iCode']=split_ri(resdict['resseqnumi'])
                            val=PDBRecord(resdict)
                        else:
                            val=rectify(cifrec.getValue(v,idx))
                            if k=='resseqnumi':
                                idict['seqNum'],idict['iCode']=split_ri(val)
                            else:
                                if k in splits and ',' in val:
                                    val=[rectify(x) for x in val.split(',')]
                                if k==spawns_on:
                                    if type(val)==str and ',' in val:
                                        val=[rectify(x) for x in val.split(',')]
                                if k in map_values:
                                    mapper=self.global_maps[map_values[k]]
                                    if type(val)==list:
                                        logger.debug(f'mapper {mapper}')
                                        logger.debug(f'list before mapping {val}')
                                        mapped_val=list(set([str(mapper[x]) for x in val]))
                                        logger.debug(f'list after mapping {mapped_val}')
                                        try:
                                            mapped_val.sort()
                                            val=mapped_val
                                        except:
                                            raise TypeError(f'could not sort list {mapped_val} at key {k}')
                                    else:
                                        val=mapper[val]
                        idict[k]=val
                        if k==indexes:
                            idict['tmp_label']=f'{k}{val}'
                    for la,vn in list_attr.items():
                        from_existing=all([x in idict for x in vn])
                        if from_existing:
                            idict[la]=[idict[x] for x in vn]
                        else:
                            idict[la]=vn
                    if spawns_on:
                        spdicts=self.gen_dict(mapspec['spawn_data'])
                        if type(idict[spawns_on])==list:
                            spawned_dicts=[]
                            for v in idict[spawns_on]:
                                sd=idict.copy()
                                sd[spawns_on]=v
                                for sp in spdicts:
                                    if sp['spawn_idx']==v:
                                        break
                                else:
                                    raise Exception(f'(list) cannot find spawn index for {spawns_on} = {v}; spdicts: {spdicts}')
                                spc=sp.copy()
                                del spc['spawn_idx']
                                spclabel=spc.get('tmp_label','')
                                if 'tmp_label' in spc:
                                    del spc['tmp_label']
                                sd.update(spc)
                                if 'tmp_label' in sd and spclabel!='':
                                    sd['tmp_label']=f'{sd["tmp_label"]}.{spclabel}'
                                spawned_dicts.append(sd)
                            idicts.extend(spawned_dicts)
                        else:
                            spawned_dicts=[]
                            v=idict[spawns_on]
                            for sp in spdicts:
                                if sp['spawn_idx']==v:
                                    break
                            else:
                                raise Exception(f'cannot find spawn index for {spawns_on} = {v}')
                            spc=sp.copy()
                            del spc['spawn_idx']
                            spclabel=spc.get('tmp_label','')
                            if 'tmp_label' in spc:
                                del spc['tmp_label']
                            idict.update(spc)
                            if 'tmp_label' in idict and spclabel!='':
                                idict['tmp_label']=f'{idict["tmp_label"]}.{spclabel}'
                            idicts.append(idict)
                    else:
                        idicts.append(idict)
        else:
            tabledict={}
            for tname,tspec in tables.items():
                tabledict[tname]=[]
                attr_map=tspec['row_attr_map']
                bisv=tspec.get('blank_if_single_valued',[])
                for i in range(len(cifrec)):
                    tdict={}
                    for k,v in attr_map.items():
                        tdict[k]=rectify(cifrec.getValue(v,i))
                        if k in bisv:
                            if len(self.global_ids[k])<2:
                                tdict[k]=''
                    tabledict[tname].append(BaseRecord(tdict))
            udict={'tables':tabledict}
            idicts.append(udict)

        if allcaps:
            for idict in idicts:
                for k,v in idict.items():
                    if k in allcaps:
                        idict[k]=v.upper()
        return idicts

    def parse(self):
        recdict={}
        for rectype,mapspec in self.formats.items():
            idicts=self.gen_dict(mapspec)
            for idict in idicts:
                this_key=idict.get('tmp_label','')
                reckey=rectype if not this_key else f'{rectype}.{this_key}'
                if reckey in recdict:
                    if not type(recdict[reckey])==list:
                        recdict[reckey]=[recdict[reckey]]
                    idict['key']=reckey
                    recdict[reckey].append(PDBRecord(idict))
                else:
                    idict['key']=reckey
                    recdict[reckey]=PDBRecord(idict)
        return recdict
