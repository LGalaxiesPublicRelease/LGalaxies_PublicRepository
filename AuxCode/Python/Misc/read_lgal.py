import numpy as np
import os


###############################################
# Code to read the INPUT trees for L-Galaxies #
###############################################
# NOT tested

struct_lgalinput = np.dtype([
('Descendant',np.int32,1),
('FirstProgenitor',np.int32,1),
('NextProgenitor',np.int32,1),
('FirstHaloInFOFgroup',np.int32,1),
('NextHaloInFOFgroup',np.int32,1),
('Len',np.int32,1),
('M_Mean200',np.float32,1),
('M_Crit200',np.float32,1),
('M_TopHat',np.float32,1),
('Pos',np.float32,3),
('Vel',np.float32,3),
('VelDisp',np.float32,1),
('Vmax',np.float32,1),
('Spin',np.float32,3),
('MostBoundID',np.int64,1),
('SnapNum',np.int32,1),
('FileNr',np.int32,1),
('SubhaloIndex',np.int32,1),
('SubHalfMass',np.int32,1)
])

def read_input_tree(folder,firstfile,lastfile,lastsnap):
    """ Reads the tree files that are used for input into L-Galaxies. """ 
    nHalos = 0
    nTrees = 0
    ngalstree = np.array([],dtype=np.int32)
    output_trees = np.array([],dtype=struct_lgalinput)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+"/trees_%03d.%d"%(lastsnap,ifile)
        f = open(filename,"rb")
        this_nTrees = np.fromfile(f,np.int32,1)[0]
        this_nHalos = np.fromfile(f,np.int32,1)[0]
        this_ngalstree = np.fromfile(f,np.int32,this_nTrees)
        this_trees = np.fromfile(f,struct_lgalinput,this_nHalos)
        nHalos += this_nHalos
        nTrees += this_nTrees
        ngalstree = np.append(ngalstree,this_ngalstree)
        output_trees = np.append(output_trees,this_trees)
    return (nHalos,nTrees,ngalstree,output_trees)

########################################
# Code to read the OUPUT snapshot data #
########################################

def read_tree(folder,file_prefix,firstfile,lastfile,props,template):
    """ Reads L-Galaxies output data files that are written in tree format.
    Returns: (nHalos,gals)
    Inputs: (folder,file_prefix,firstfile,lastfile,props,template)
    props - list of properties to return
    template - structure dtype definition for database
    NOT TESTED """
    nHalos = 0
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)
    gals = np.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"galtree_"+"%d"%(ifile)
        f = open(filename,"rb")
        one =  np.fromfile(f,np.int32,1)
        structsize = np.fromfile(f,np.int32,1)
        if(structsize != template.itemsize):
            print "size mismatch:",structsize,template.itemsize
        this_nHalos = np.fromfile(f,np.int32,1)
        nHalos += this_nHalos
        f.seek(structsize, os.SEEK_SET) 
        print "File ", ifile," nGals = ",this_nHalos
        this_addedGalaxy = np.fromfile(f,template,this_nHalos)
        addedGalaxy = np.zeros(this_nHalos,dtype=filter_dtype)
        for prop in template.names:
            if props[prop]:
                addedGalaxy[prop] = this_addedGalaxy[prop]
        gals = np.append(gals,addedGalaxy)      
        f.close()
    return (nHalos,gals)

def read_snap(folder,file_prefix,firstfile,lastfile,props,template):
    """ Reads L-Galaxy output files.
    Returns: (nTrees,nHalos,nTreeHalos,gals)
    Inputs: (folder,file_prefix,firstfile,lastfile,props,template)
    props - list of properties to return
    template - structure dtype definition for database """
    nTrees = 0
    nHalos = 0
    nTreeHalos = np.array([],dtype=np.int32)
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)
    gals = np.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        this_nTrees =  np.fromfile(f,np.int32,1)
        nTrees += this_nTrees
        this_nHalos = np.fromfile(f,np.int32,1)
        nHalos += this_nHalos
        print "File ", ifile," nGals = ",this_nHalos
        addednTreeHalos = np.fromfile(f,np.int32,this_nTrees)
        nTreeHalos = np.append(nTreeHalos,addednTreeHalos)
        this_addedGalaxy = np.fromfile(f,template,this_nHalos) # all properties
        addedGalaxy = np.zeros(this_nHalos,dtype=filter_dtype) # selected props
        for prop in template.names:
            if props[prop]:
                addedGalaxy[prop] = this_addedGalaxy[prop]
        gals = np.append(gals,addedGalaxy)      
        f.close()
    return (nTrees,nHalos,nTreeHalos,gals)

