import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys
import os.path
import cPickle
from scipy import stats

plt.rc('text', usetex=True)


#All the SIM data locations
SIMpath = ['/lustre/scratch/astro/bf77/NM_coldMvir.pkl','/lustre/scratch/astro/bf77/NM_coldMvir_mcmc.pkl','/lustre/scratch/astro/bf77/NM_coldBower.pkl','/lustre/scratch/astro/bf77/NM_coldBower_mcmc.pkl','/lustre/scratch/astro/bf77/NM_cooledMvir.pkl','/lustre/scratch/astro/bf77/NM_cooledMvir_mcmc.pkl','/lustre/scratch/astro/bf77/NM_cooledBower.pkl','/lustre/scratch/astro/bf77/NM_cooledBower_mcmc.pkl','/lustre/scratch/astro/bf77/NM_REF.pkl','/lustre/scratch/astro/bf77/NM_TESTMOREBULGES.pkl']


#Initializing parameters and mass bins
L_Path = len(SIMpath)

scale = 1.#512./100.#512./100.
redshift = 0.
Hubble_h = 0.673
xlimB=[8.0,11.5]
ylimB=[0., 1.2]
binB=[0.1,0.05]
NbinsB=[int((xlimB[1]-xlimB[0])/binB[0]),int((ylimB[1]-ylimB[0])/binB[1])]
binB=0.25


#Loading the 3 different bulge fractions obs files
data_B = np.loadtxt("/research/astro/virgo/Benoit/L-GALAXY/mcmc_NewModel/ObsConstraints/BulgeFraction_z0.10.txt",skiprows=1)
data_D = np.loadtxt("/research/astro/virgo/SAM/Obsdata/conselice2006_disk_fract.txt",skiprows=1)
data_I = np.loadtxt("/research/astro/virgo/SAM/Obsdata/conselice2006_irr_fract.txt",skiprows=1)



fig6 = plt.figure(6)#, ax6 = plt.subplots()

for i in range(L_Path):

    #Loading the data and extracting Mass components
    i_ind = i + 1
    print "\n Run : ", i
    f1 = open(SIMpath[i], 'rb')
    data = cPickle.load(f1)
    f1.close()

    BULGE = data["BulgeMass"] * 1.0e10
    STARS = data["DiskMass"] * 1.0e10
    STELLAR_MASS = BULGE + STARS

    #Units need to correspond to those of the obs files
    BulgeMassRF = np.log10(BULGE*Hubble_h)
    StellarMassRF= np.log10(STELLAR_MASS*Hubble_h)+ np.random.randn(len(STELLAR_MASS))*0.08*(1+redshift)

    #Creating the subplots. Need to change it accordingly to the number of SIM data used/wanted
    ax6 = fig6.add_subplot(5,2,i_ind)

    #Calculating/plotting the bulge fraction for bulge dominated galaxies
    Mass_arrB=np.arange(xlimB[0],xlimB[1],binB)
    BulgeFraction=np.zeros(len(Mass_arrB),dtype=np.float32)
    sel_bulge = None
    sel_allB = None
    for ll in range(0,len(Mass_arrB)):
	    sel_bulge=data[((BulgeMassRF - StellarMassRF)> -0.154902) & (StellarMassRF>Mass_arrB[ll]-binB/2.) & (StellarMassRF<Mass_arrB[ll]+binB/2.)]
	    sel_allB=data[(StellarMassRF>Mass_arrB[ll]-binB/2.) & (StellarMassRF<Mass_arrB[ll]+binB/2.)]
	    if (float(len(sel_allB)) > 0.):
		    BulgeFraction[ll]=float(len(sel_bulge))/float(len(sel_allB))
	    else:
		    BulgeFraction[ll]=0.

    ax6.plot(Mass_arrB,BulgeFraction,linestyle='-',color='r')

    #Calculating/plotting the bulge fraction for disk dominated galaxies
    sel_bulge = None #Just in case
    sel_allB = None #Just in case
    for ll in range(0,len(Mass_arrB)):
	    sel_bulge=data[((BulgeMassRF - StellarMassRF) <= -0.154902) & ((BulgeMassRF - StellarMassRF) >= -2.0) & (StellarMassRF>Mass_arrB[ll]-binB/2.) & (StellarMassRF<Mass_arrB[ll]+binB/2.)]
	    sel_allB=data[(StellarMassRF>Mass_arrB[ll]-binB/2.) & (StellarMassRF<Mass_arrB[ll]+binB/2.)]
	    if (float(len(sel_allB)) > 0.):
		    BulgeFraction[ll]=float(len(sel_bulge))/float(len(sel_allB))
	    else:
		    BulgeFraction[ll]=0.

    ax6.plot(Mass_arrB,BulgeFraction,linestyle='-',color='b')

    #Calculating/plotting the bulge fraction for irregular/no bulge galaxies
    sel_bulge = None
    sel_allB = None
    for ll in range(0,len(Mass_arrB)):
	    sel_bulge=data[((BulgeMassRF - StellarMassRF)< -2.0) & (StellarMassRF>Mass_arrB[ll]-binB/2.) & (StellarMassRF<Mass_arrB[ll]+binB/2.)]
	    sel_allB=data[(StellarMassRF>Mass_arrB[ll]-binB/2.) & (StellarMassRF<Mass_arrB[ll]+binB/2.)]
	    if (float(len(sel_allB)) > 0.):
		    BulgeFraction[ll]=float(len(sel_bulge))/float(len(sel_allB))
	    else:
		    BulgeFraction[ll]=0.

    ax6.plot(Mass_arrB,BulgeFraction,linestyle='-',color='g')

    #Plotting the observational data
    ax6.errorbar(data_B[:,0]+np.log10(0.673)+np.log10(0.673),data_B[:,1],yerr=data_B[:,2],linestyle='None',color='r',marker='+',label="Observations (MCMC)")
    ax6.errorbar(data_D[:,0]+np.log10(0.673)+np.log10(0.673),data_D[:,1],yerr=data_D[:,2],linestyle='None',color='b',marker='+',label="Observations (MCMC)")
    ax6.errorbar(data_I[:,0]+np.log10(0.673)+np.log10(0.673),data_I[:,1],yerr=data_I[:,2],linestyle='None',color='g',marker='+',label="Observations (MCMC)")

    if (i == 4):
        ax6.set_ylabel(r"BulgeFraction")
    if ((i == 8 ) | (i == 9)):
        ax6.set_xlabel(r"$M_* [M_{\odot}h^{-2}]$")


plt.savefig("PLOT_allBF.pdf")
