
# coding: utf-8

# In[ ]:


#PLOT OPTIONS
opt_stellar_mass_function=1
opt_metals_vs_stellarmass=0
opt_BHBM=0
opt_SFRF=0
opt_gas_fraction=0
opt_HI_MF=0
opt_sfr_vs_stellar_mass=0
opt_ur_vs_r=0
opt_UVJ_colour=0
opt_redfraction_color_cut=0
    
opt_plot_MCMC_sample=0

Datadir = '/Users/BrunoHenriques/Desktop/OneDrive/Workspace/GitHub_PR_Hen15/AuxCode/Python/data/'
MCMCdir = '/Users/BrunoHenriques/Desktop/OneDrive/Workspace/GitHub_PR_Hen15/MCMC/'
#DirName_MR = '/net/bootes/scratch2/SAM/test1/MR/'
#DirName_MRII = '/net/bootes/scratch2/SAM/test1/MRII/'

DirName_MR = '/Users/BrunoHenriques/Desktop/OneDrive/Workspace/GitHub_PR_Hen15/output/'
DirName_MRII = '/net/bootes/scratch2/SAM/Henriques2015a/snaps/MRII/'

prefix_this_model='This Work - PLANCK1'
file_this_model='ThisWork'

do_previous_model1=1
file_previous_model1=Datadir+'Guo2013a_m05'
prefix_previous_model1='Guo2013a - WMAP7'
linestyle_previous_model1=':'

#do_previous_model2=1
#file_previous_model2=Datadir+'Henriques2013a'
#prefix_previous_model2='Henriques2013a - WMAP7'
#linestyle_previous_model2='--'

do_previous_model2=1
file_previous_model2=Datadir+'/Henriques2015a'
prefix_previous_model2='Henriques2015 - PLANCK1'
linestyle_previous_model2='--'


slope_red_fraction=[0.075,0.275, 0.3, 0.32,0.38]
offset_red_fraction=[1.85,1.213, 1.18,0.99,0.79]
minimum_y_red_fraction=[0.0,1.3,1.3,1.3,1.3]

RedshiftsToRead = [True,True,True,True,True,True,False]

CatalogType='snap'
#CatalogType='tree'

Hubble_h_WMAP1 = 0.732
Hubble_h_WMAP7 = 0.704


PLANCK=1
if PLANCK: 
    FullRedshiftList=[0.00,0.11,0.40,1.04,2.07,3.11,3.95] 
    FullSnapshotList=[58,54,47,38,30,25,22]  
    BoxSize_MR    = 500.* 0.960558 #full MR 
    BoxSize_MRII  = 100.* 0.960558 #full MRII      
    Hubble_h      = 0.673
    Omega_M       = 0.315 
    Omega_Lambda  = 0.683
    MaxTreeFiles  = 512
    



