# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np

def read_snap(folder,FirstFile,LastFile,
              props,template,RedshiftsToRead,RedshiftList):    
    """ Reads L-Galaxy output files.
    Returns: (nTrees,nHalos,nTreeHalos,gals)
    Inputs: (folder,file_prefix,FirstFile,LastFile,props,template)
    props - list of properties to return
    template - structure dtype definition from database """     
    nTrees = 0
    nHalos = 0    
    nTreeHalos = np.array([],dtype=np.int32)
    
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)    
    gals = np.array([],dtype=filter_dtype)      
    
    SnapshotList=np.array([],dtype=np.int32)
    
    for iredshift in range(0,len(RedshiftList)-1):
        if RedshiftsToRead[iredshift]:  
            print ("\n\nReading redshift: ", RedshiftList[iredshift], "\n")
            for ifile in range(FirstFile,LastFile+1):
                char_redshift="%0.2f" % RedshiftList[iredshift]
                filename = folder+'/'+'SA_z'+char_redshift+"_"+"%d"%(ifile)
                f = open(filename,"rb")
                
                this_nTrees =  np.fromfile(f,np.int32,1)
                nTrees += this_nTrees
                this_nHalos = np.fromfile(f,np.int32,1)
                nHalos += this_nHalos
                print ("File ", ifile," nGals = ",this_nHalos)  
                
                addednTreeHalos = np.fromfile(f,np.int32,this_nTrees)
                nTreeHalos = np.append(nTreeHalos,addednTreeHalos)
                this_addedGalaxy = np.fromfile(f,template,this_nHalos) # all properties
                addedGalaxy = np.zeros(this_nHalos,dtype=filter_dtype) # selected props
                
                for prop in template.names:
                    if props[prop]:
                        addedGalaxy[prop] = this_addedGalaxy[prop]
                gals = np.append(gals,addedGalaxy)      
                f.close()           
            #endfor
        #endif    
        SnapshotList=np.append(SnapshotList,gals['SnapNum'][len(gals)-1])
    #endfor
    
    return (gals, SnapshotList)

# <codecell>


# <codecell>


