
# coding: utf-8

# In[ ]:

RedshiftsToRead = [True,True,True,True,True,False]

PLANCK=1
if PLANCK: 
    RedshiftList=[0.00,0.40,1.04,2.07,3.11,3.95]   
    BoxSize_MR    = 500.* 0.960558 #full MR 
    BoxSize_MRII  = 100.* 0.960558 #full MRII      
    Hubble_h      = 0.673
    Omega_M       = 0.315 
    Omega_Lambda  = 0.683
    MaxTreeFiles  = 512
    
Datadir = '/net/bootes/export/data1/data/'

#DirName_MR = '/net/bootes/scratch2/SAM/test1/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/test1/MRII/'

DirName_MR = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MR/'
DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MRII/'



