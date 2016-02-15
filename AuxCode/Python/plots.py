''' 
stellar_mass_function()
bluck_red_fractions()
metals_vs_stellarmass()
BHBM

gas_metallicity_gradients()
SFR_gradients()

test_plots()
'''

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.table import Table
from importlib import reload
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import sys
from scipy.ndimage import zoom

import procedures
reload (procedures)
from procedures import *
import plots_input
reload (plots_input)
from plots_input import *


def stellar_mass_vs_halo_mass(G_MR, ThisRedshiftList, pdf):

    ylim=[4.0,12.5]
    xlim=[7.5, 15.]
  
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                         'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})    
    fig = plt.figure(figsize=(10,10))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim) 

    ylab='$\mathrm{log_{10}}(M_*[h^{-2}M_{\odot}])$'       
    xlab='$\mathrm{log_{10}}(M_{200c}[h^{-2}M_{\odot}])$'
    subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)
                  
    for ii in range (0,len(ThisRedshiftList)):
  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        HaloMass=np.log10(G0_MR['Mvir']*1.e10*Hubble_h)
        subplot.scatter(HaloMass,StellarMass,s=5, color='black')
    #endfor
    
    
    plt.tight_layout()
    plt.savefig('./fig/plots_smhm.pdf')
    pdf.savefig()
    plt.close()
    
#end stellar_mass_vs_halo_mass

def stellar_mass_function(G_MR, Volume_MR, ThisRedshiftList, pdf):
           
    xlim=[7.0,12.5]
    ylim=[-6.5, 0.5]
    bin=0.1


    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                         'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

    fig = plt.figure(figsize=(8,7))
    grid = gridspec.GridSpec(2, 2)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])

        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if ii==2 or ii == 3: 
            xlab='$\mathrm{log_{10}}(M_*[h^{-2}M_{\odot}])$'
        else:
            xlab=''
        if ii==0 or ii == 2:
            ylab='$\mathrm{log_{10}}(\phi [h^3 \mathrm{Mpc^{-3}} \mathrm{log_{10}}(M^{-1})])$'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)
                
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
        if ii==1 or ii == 3:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
           
    
        #OBSERVATIONS             
        file = MCMCdir + '/ObsConstraints/StellarMassFunction_z'+char_redshift+'.txt'        
        f = open(file, 'r')     
        line = int(f.readline())     
        obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
        
        obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.
        asy_yerror = [np.log10(obs['col3']/(obs['col3']-obs['col4'])), 
                      np.log10((obs['col3']+obs['col4'])/obs['col3'])]
        subplot.errorbar(obs_xbin, np.log10(obs['col3']),yerr=asy_yerror,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
        #sub = plt.subplot(111)
    
        #PREVIOUS MODELS
        RedshiftList_OldModels=[0.1,1.0,2.0,3.0]
        if do_previous_model1==1:           
            char_old_redshift="%0.2f" % RedshiftList_OldModels[ii]
            file = file_previous_model1+'_smf_z'+char_old_redshift+'.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model1, linewidth=2)
      
        if do_previous_model2==1:           
            file = file_previous_model2+'_smf_z'+char_redshift+'.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model2, linewidth=2)
      
    
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist=np.histogram(StellarMass, bins=bin_arr, range=(xlim[0],xlim[1]))   
        subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                 color='red', linewidth=2)


        #MCMC sample
        if opt_plot_MCMC_sample==1:
            file = MCMCSampledir + 'mcmc_plus_obs0_z'+char_redshift+'.txt' 
            obs = Table.read(file, format='ascii')      
            subplot.plot(obs['col1'],np.log10(obs['col4']), color='black', linewidth=2)
            
            if ii==len(ThisRedshiftList)-1:
                plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.55, y_percentage=0.85, color='black', xlog=0, ylog=0, 
                    label='MCMC sample', fontsize=13, fontweight='normal') 
                plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.44, y_percentage=0.87, color='black', x2_percentage=0.53, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        #LABELS
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Observations used in MCMC', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
            
        if ii==len(ThisRedshiftList)-1:
            plot_label_three_models (subplot, xlim, ylim, position='bottom_left')
                                  
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.35, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal')      
        if ii==2:
             plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.075, 
                         color='black', xlog=0, ylog=0, label='SMF evo', 
                         fontsize=16, fontweight='normal') 
           
    #endfor


    plt.tight_layout()
    plt.savefig('./fig/plots_smf_evo.pdf')
    pdf.savefig()
    plt.close()
#endif stellar_mass_function==1:



def metals_vs_stellarmass(G_MR, ThisRedshiftList, pdf):
           
    xlim=[9.0,12.0]
    ylim=[-1.5, 0.5]   
    bin=0.2
        
    plot_color=['red','purple']        
    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
    fig = plt.figure(figsize=(5,4))
    subplot=plt.subplot()
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)
    
    #format axis
    majorFormatter = FormatStrFormatter('%d')
    subplot.xaxis.set_major_locator(MultipleLocator(1))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
    
    xlab='$\mathrm{log_{10}}(M_*[h^{-2}M_{\odot}])$'       
    ylab='$\mathrm{log_{10}}(Z/Z_{\odot})$'         
    subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)    
        
    for ii in range (0,len(ThisRedshiftList)):
                
        #observations from GALLAZI    
        xmass=np.arange(16)*bin+8.8
        obsp50=[-0.60,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13,0.14,0.15]
        obsp16=[-1.11,-1.07,-1.10,-1.03,-0.97,-0.90,-0.80,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04,-0.03,-0.03]
        obsp84=[-0.00,-0.00,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28,0.29,0.30]  
        subplot.errorbar(xmass, obsp50, color='blue', fmt='o')       
        subplot.plot(xmass, obsp16,color='blue', linestyle='--')
        subplot.plot(xmass, obsp84,color='blue', linestyle='--')
    
        #PREVIOUS MODELS       
        if ThisRedshiftList[ii]==0.1:
            char_redshift="%0.2f" % ThisRedshiftList[ii]
            if do_previous_model1==1: 
                file = file_previous_model1+'_metals_median_z'+char_redshift+'.txt' 
                model = Table.read(file, format='ascii')
                subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model1, linewidth=2)
      
            if do_previous_model2==1: 
                file = file_previous_model2+'_metals_z'+char_redshift+'.txt' 
                model = Table.read(file, format='ascii')
                subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model2, linewidth=2)
      
    
        #Model
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)                 
        G0_MR=G_MR[sel]   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])                                     
        Metallicity=np.log10((G0_MR['MetalsStellarMass'])/((G0_MR['StellarMass'])*0.02))    
              
        (x_binned, median, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Metallicity)    
        subplot.plot(x_binned, median,color=plot_color[ii], linewidth=2)
        if (ii==0):
            subplot.plot(x_binned, pc16,color=plot_color[ii], linewidth=2, linestyle='-')
            subplot.plot(x_binned, pc84,color=plot_color[ii], linewidth=2, linestyle='-')
            
            
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                    color='black', xlog=0, ylog=0, label='Gallazzi 2005', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.05) 
                               
        plot_label_three_models (subplot, xlim, ylim, position='bottom_left')
        
       
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_metals_vs_stellarmass.pdf')
    pdf.savefig()
    plt.close()

#end metals_vs_stellarmass

def BHBM(G_MR, ThisRedshiftList, pdf):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[8.5,12.5]
        ylim=[5.0, 10.5]
        bin=[0.1,0.05]
        Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    
        plot_color=['red','purple']        
        plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
        fig = plt.figure(figsize=(5,4))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
            
        xlab='$\mathrm{log_{10}}(M_{\mathrm{Bulge}}[h^{-2}M_{\odot}])$'       
        ylab='$\mathrm{log_{10}}(M_{\mathrm{BH}}[h^{-1}M_{\odot}])$' 
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)        
        G0_MR_unsel=G_MR[sel]   
        G0_MR=G0_MR_unsel[(G0_MR_unsel['BulgeMass'] > 0.) & (G0_MR_unsel['BlackHoleMass'] > 0.)]
        Ngals=len(G0_MR) 
       
        BulgeMass=(np.log10(G0_MR['BulgeMass']*1.e10*Hubble_h)) 
        BHMass=(np.log10(G0_MR['BlackHoleMass']*1.e10)) 
         
        #plt.scatter(BulgeMass, BHMass, s=5, color='black')          
        H, xedges, yedges = np.histogram2d(BulgeMass, BHMass, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/0.7)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
        plt.colorbar(format='%d') 
        #(ax, cmap=None, norm=None, alpha=None, values=None, boundaries=None, orientation='vertical', 
        #ticklocation='auto', extend='neither', spacing='uniform', ticks=None, format=None, 
        #drawedges=False, filled=True, extendfrac=None, extendrect=False, label='')
       
        file = Datadir + 'mcconnel2012.dat'
        obs = Table.read(file, format='ascii', data_start=20)     
        obs_x = np.log10(obs['col14']*Hubble_h_WMAP7**2)
        obs_y = np.log10(obs['col3']*Hubble_h_WMAP7**2)        
        obs_x_err=np.zeros(len(obs_x),dtype=np.float64)+0.24 
        obs_y_err = [np.log10(obs['col3']/obs['col4']),np.log10(obs['col5']/obs['col3'])]
       
        subplot.errorbar(obs_x, obs_y,xerr= obs_x_err, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
        
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='McConnel 2012', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.925, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.7, y_percentage=0.1, 
                    color='black', xlog=0, ylog=0, label='BHBM', 
                    fontsize=13, fontweight='normal') 
        
    plt.tight_layout()
    plt.savefig('./fig/plots_bhbm.pdf')
    pdf.savefig()
    plt.close()

#end BHBM

def SFRF(G_MR, Volume_MR, ThisRedshiftList, pdf):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[-1.0,3.0]
        ylim=[-6.0,-1.0]
        bin=0.1
    
        plot_color=['red','purple']        
        plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
        fig = plt.figure(figsize=(5,4))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25)) 
            
        xlab='$\mathrm{log_{10}}(\mathrm{SFR}[h^{-2}M_{\odot}yr^{-1}])$'   
        ylab='$\mathrm{log_{10}}(\phi [h^3 \mathrm{Mpc^{-3}} \mathrm{log_{10}}(SFR^{-1})])$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)         
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['Sfr']>0.]
        SFR=(np.log10(G0_MR['Sfr']*Hubble_h**2))
                
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist=np.histogram(SFR, bins=bin_arr, range=(xlim[0],xlim[1]))   
        subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                 color='red', linewidth=2)
      
        file = Datadir + 'gruppioni2015.txt'
        obs = Table.read(file, format='ascii', data_start=0, data_end=7)       
        obs_x = (obs['col1']+obs['col2'])/2.
        obs_y = np.log10(obs['col3'])       
        obs_y_err = np.log10((obs['col3']+obs['col4']))-np.log10((obs['col3']-obs['col4']))
       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
                
        #LABELS   
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label='Gruppioni 2015', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.225, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
           
        plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.1, color='black', xlog=0, ylog=0, 
                label=prefix_this_model, fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.12, color='red', x2_percentage=0.13, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.775, y_percentage=0.85, 
                    color='black', xlog=0, ylog=0, label='SFRF', 
                    fontsize=13, fontweight='normal') 
            
    plt.tight_layout()
    plt.savefig('./fig/plots_sfrf.pdf')
    pdf.savefig()
    plt.close()

#end SFRF



def gas_fraction(G_MR, ThisRedshiftList, pdf):
   
    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[8.5,12.0]
        ylim=[-2.5,1.0]
       
        bin=0.1
        plot_color=['red','purple']        
        plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
        fig = plt.figure(figsize=(5,4))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$log_{10}(M_*[h^{-2}M_{\odot}])$'
        ylab='$M_{\mathrm{Cold}}/M_*$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)        
        G0_MR=G_MR[sel]          
        G0_MR=G0_MR[G0_MR['StellarMass']>0.]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        Fraction=np.log10(G0_MR['ColdGas']*1.e10*Hubble_h)-StellarMass 
   
        (x_binned, median, pc16, pc84)=median_and_percentiles (bin, xlim[0], xlim[1], StellarMass, Fraction)   
        sel=(median!=0)        
        subplot.plot(x_binned[sel], median[sel],color=plot_color[ii], linewidth=2)     
        subplot.plot(x_binned[sel], pc16[sel],color=plot_color[ii], linewidth=2, linestyle='--')
        subplot.plot(x_binned[sel], pc84[sel],color=plot_color[ii], linewidth=2, linestyle='--')
    
    
        file = Datadir + 'peeples_2015.txt'
        obs = Table.read(file, format='ascii', data_start=0)       
        obs_x = (obs['col1']+obs['col2'])/2.
        obs_y = np.log10(obs['col3'])       
        obs_y_err = [np.log10(obs['col3']/(obs['col3']-obs['col4'])),np.log10((obs['col3']+obs['col4'])/obs['col3'])]
       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
                
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.2, 
                    color='black', xlog=0, ylog=0, label='Peeples 2015', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.225, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
           
        plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=0.15, y_percentage=0.1, color='black', xlog=0, ylog=0, 
                label=prefix_this_model, fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                x_percentage=0.04, y_percentage=0.12, color='red', x2_percentage=0.13, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.55, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label='Gas Fraction', 
                    fontsize=13, fontweight='normal')   
            
    plt.tight_layout()
    plt.savefig('./fig/plots_GMF.pdf')
    pdf.savefig()
    plt.close()

#end gas fraction

def HI_MF(G_MR, Volume_MR, ThisRedshiftList, pdf):
    for ii in range(0,len(ThisRedshiftList)):        
        
        xlim=[8.0,11.5]
        ylim=[-6.0,0.0]
        bin=0.25
    
        plot_color=['red','purple']        
        plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                             'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})
        fig = plt.figure(figsize=(5,4))
        subplot=plt.subplot()
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)    
            
        #format axis
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))
            
        xlab='$\mathrm{log_{10}}(\mathrm{M_{\mathrm{HI}}}[h^{-2}M_{\odot}])$'   
        ylab='$\mathrm{log_{10}}(\phi [h^3 \mathrm{Mpc^{-3}} \mathrm{log_{10}}(M_{\mathrm{HI}}^{-1})])$'     
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)   
            
        #PREVIOUS MODELS       
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        if do_previous_model1==1: 
            file = file_previous_model1+'_coldgas_MF.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model1, linewidth=2)
      
        if do_previous_model2==1: 
            file = file_previous_model2+'_coldgas_MF.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model2, linewidth=2) 
            
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)         
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[G0_MR['ColdGas']>0.]
        HI=(np.log10(G0_MR['ColdGas']*0.54*1.e10*Hubble_h))
        
        bin_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        hist=np.histogram(HI, bins=bin_arr, range=(xlim[0],xlim[1]))   
        subplot.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume_MR*bin)),
                 color='red', linewidth=2)
      
        #OBSERVATIONS
        h=0.75
        file = Datadir + 'zwaan2005.txt'       
        obs = Table.read(file, format='ascii')      
        obs_x = obs['col1']
        obs_y = obs['col2']-3.*np.log10(h)      
        obs_y_err = [-obs['col3'],obs['col4']]
       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='blue', color='blue')
        
        
        file = Datadir + 'haynes2011_gmf.txt'       
        obs = Table.read(file, format='ascii')      
        obs_x = obs['col1']+2.*np.log10(Hubble_h_WMAP7)
        obs_y = obs['col2']-3.*np.log10(Hubble_h_WMAP7)      
        obs_y_err = [-obs['col3'],obs['col4']]
       
        subplot.errorbar(obs_x, obs_y, yerr=obs_y_err,
                 fmt='o', markersize=5, ecolor='limegreen', color='limegreen')
         
        #LABELS        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.5, 
                    color='black', xlog=0, ylog=0, label='Zwaan 2005', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.525, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.4, 
                    color='black', xlog=0, ylog=0, label='Haynes 2011', 
                    fontsize=13, fontweight='normal') 
        plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.425, 
                    color='limegreen', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label_three_models (subplot, xlim, ylim, position='bottom_left')
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.75, y_percentage=0.8, 
                    color='black', xlog=0, ylog=0, label='HI MF', 
                    fontsize=13, fontweight='normal')   
        
    plt.tight_layout()
    plt.savefig('./fig/plots_HI_MF.pdf')
    pdf.savefig()
    plt.close()
   

#end HI_MF

def sfr_vs_stellar_mass(G_MR, ThisRedshiftList, pdf):
           
    xlim=[8.5,11.5]
    ylim=[-2.5, 3]   
    bin=[0.1,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                         'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

    fig = plt.figure(figsize=(15,4))
    grid = gridspec.GridSpec(1, 5)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.1f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*[h^{-2}M_{\odot}])$'      
        if ii==0:
            ylab='$\mathrm{log_{10}}(\mathrm{SFR}[h^{-2}M_{\odot}yr^{-1}])$'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)
        
        subplot.text(xlim[0]+2.,ylim[0]+.3,'z='+char_redshift, fontsize=16, fontweight='normal')
        
        if ii==2:
            subplot.text(xlim[0]+0.3,ylim[0]+0.5,'MS evo', fontsize=16, fontweight='normal')
        
        
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25))
    
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
           
        
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['StellarMass']>0.) & (G0_MR['Sfr']>0.)]
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        SFR=(np.log10(G0_MR['Sfr']*Hubble_h**2))
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(StellarMass, SFR, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/0.7)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)   
        
        if ii==len(ThisRedshiftList)-1:
            plt.colorbar(format='%d') 
            
        #OBSERVATIONS
        #values at all_z
        
        #ELBAZ2007
        obs_slope_elbaz2007 =[0.77, -99.0, 0.9, -99.0, -99.0]
        obs_offset_elbaz2007=[np.log10(8.7)-(0.77*11.), -99.0, np.log10(7.2)-9, -99.0, -99.0]
        obs_offset_low_elbaz2007=[np.log10(5.0)-(0.77*11.), -99.0, np.log10(3.6)-9, -99.0, -99.0]
        obs_offset_high_elbaz2007=[np.log10(16.1)-(0.77*11.), -99.0, np.log10(14.4)-9, -99.0, -99.0]
        
        #KARIM2011
        file = Datadir + 'karim2011_sfr_mass_sf.txt'       
        karim2011 = Table.read(file, format='ascii') 
        karim_low_z_limit        = karim2011['col4']
        karim_medium_mass        = karim2011['col3']
        karim_sfr                = karim2011['col19']
        karim_sfr_error_up   = karim2011['col20']
        karim_sfr_error_down = karim2011['col21']
        log_karim_sfr_error_up=np.log10((karim_sfr+karim_sfr_error_up)/karim_sfr)
        log_karim_sfr_error_down=np.log10(karim_sfr/(karim_sfr-karim_sfr_error_down))
        
        obs_x=np.arange(xlim[0], xlim[1], 0.01)
        
        if ThisRedshiftList[ii]==0.0:  
            #ELBAZ2007
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_elbaz2007[ii] + 2.*np.log10(Hubble_h_WMAP7)
            subplot.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2)            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_low_elbaz2007[ii] + 2.*np.log10(Hubble_h_WMAP7)
            subplot.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_high_elbaz2007[ii] + 2.*np.log10(Hubble_h_WMAP7)
            subplot.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')
            
        if ThisRedshiftList[ii]==0.4: 
            #KARIM2011
            sel=(karim_low_z_limit==0.2) & (karim_medium_mass>8.8)                
            subplot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                             np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]], 
                             mfc='white', markeredgecolor='limegreen', color='limegreen', fmt='o', markersize=5)
            sel=(karim_low_z_limit==0.4) & (karim_medium_mass>8.9)                
            subplot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                             np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]], 
                             color='limegreen', fmt='o', markersize=5)
            
        if ThisRedshiftList[ii]==1.0:  
             #ELBAZ2007
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_elbaz2007[ii] + 2.*np.log10(Hubble_h_WMAP7)
            subplot.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2)            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_low_elbaz2007[ii] + 2.*np.log10(Hubble_h_WMAP7)
            subplot.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')            
            obs_y=obs_x*obs_slope_elbaz2007[ii] + obs_offset_high_elbaz2007[ii] + 2.*np.log10(Hubble_h_WMAP7)
            subplot.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')
            
            #KARIM2011
            sel=(karim_low_z_limit==0.8) & (karim_medium_mass>9.1)                
            subplot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                             np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                             mfc='white', markeredgecolor='limegreen', color='limegreen', fmt='o', markersize=5)
            sel=(karim_low_z_limit==1.0) & (karim_medium_mass>9.3)                
            subplot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                             np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],  
                             color='limegreen', fmt='o', markersize=5)
            
            #Whitaker2013
            file = Datadir + 'whitaker2013_mass_vs_sfr_z0.5_1.0.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                             obs['col3'], mfc='white', markeredgecolor='blue', color='blue', fmt='o', markersize=5)
            file = Datadir + 'whitaker2013_mass_vs_sfr_z1.0_1.5.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                             obs['col3'], color='blue', fmt='o', markersize=5)
            
        if ThisRedshiftList[ii]==2.0:  
            #KARIM2011
            sel=(karim_low_z_limit==1.6) & (karim_medium_mass>9.6)                
            subplot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                             np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                             mfc='white', markeredgecolor='limegreen', color='limegreen', fmt='o', markersize=5)
            sel=(karim_low_z_limit==2.0) & (karim_medium_mass>9.8)                
            subplot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                             np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],  
                             color='limegreen', fmt='o', markersize=5)
            
            #Whitaker2013
            file = Datadir + 'whitaker2013_mass_vs_sfr_z1.5_2.0.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                             obs['col3'], mfc='white', markeredgecolor='blue', color='blue', fmt='o', markersize=5)
            file = Datadir + 'whitaker2013_mass_vs_sfr_z2.0_2.5.txt'       
            obs = Table.read(file, format='ascii') 
            subplot.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                             obs['col3'], color='blue', fmt='o', markersize=5)
            
        if ThisRedshiftList[ii]==3.0:
            #KARIM2011
            sel=(karim_low_z_limit==2.5) & (karim_medium_mass>10.0)                
            subplot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                             np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                             [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                             color='limegreen', fmt='o', markersize=5)
         
        #labels
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.82, 
                        color='black', xlog=0, ylog=0, label='Karim 2011', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.845, 
                        color='limegreen', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.9, 
                        color='black', xlog=0, ylog=0, label='Whitaker 2013', 
                        fontsize=13, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.925, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
            plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.74, color='black', xlog=0, ylog=0, 
                    label='Elbaz 2007', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.76, color='firebrick', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_sfr_vs_stellar_mass.pdf')
    pdf.savefig()
    plt.close()
#endif stellar_mass_vs_sfr



def ur_vs_r(G_MR, ThisRedshiftList, pdf):
           
    xlim=[-23.,-14.]
    ylim=[1.0, 3.2]   
    bin=[0.5,0.025]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                         'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

    fig = plt.figure(figsize=(5,4))
    subplot=plt.subplot()
    
    subplot.set_ylim(ylim), subplot.set_xlim(xlim)      
    subplot.set_xlabel('r', fontsize=16), subplot.set_ylabel('u-r', fontsize=16) 
    
    majorFormatter = FormatStrFormatter('%2d')
    subplot.xaxis.set_major_locator(MultipleLocator(2))    
    subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      

    for ii in range(0,len(ThisRedshiftList)):                     
             
        subplot.text(xlim[0]+0.2,ylim[0]+0.25,'z=0', fontsize=16, fontweight='normal')     
        subplot.text(xlim[0]+0.2,ylim[0]+0.05,'r vs u-r cut', fontsize=16, fontweight='normal')
                  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['MagDust'][:,15]<99.) & (G0_MR['MagDust'][:,17]<99.)]        
        color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]  
        Magr=G0_MR['MagDust'][:,17]         
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(Magr, color_ur, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0]*4., Nbins[0]*4.)*Ngals/(Nbins[0]**2/0.2)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
        plt.colorbar(format='%d') 
        
        #BestFit Cut   
        bin=0.01
        x_arr=np.arange(xlim[0],xlim[1]+bin,bin)        
        y_arr=(offset_red_fraction[0]-slope_red_fraction[0]*np.tanh((x_arr+18.07)/1.09))
        subplot.plot(x_arr,y_arr,color='red', linestyle='-', linewidth=2)  
        
        #OBSERVATIONAL CUT 
        Nbin=0.01
        x_arr=np.arange(xlim[0],xlim[1]+Nbin,Nbin)       
        y_arr=2.06-0.244*np.tanh((x_arr+20.07)/1.09)        
        subplot.plot(x_arr,y_arr,color='blue', linestyle='--', linewidth=2) 
        
        #LABELS
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.9, color='black', xlog=0, ylog=0, 
                    label='Baldry 2004', fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.92, color='blue', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='--', linewidth=2)
        
        plot_label (subplot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.82, color='black', xlog=0, ylog=0, 
                    label='Best Fit cut', fontsize=13, fontweight='normal') 
        plot_label (subplot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.84, color='red', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)
        
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_ur_vs_r.pdf')
    pdf.savefig()
    plt.close()
#endif ur_vs_r







def UVJ_colour(G_MR, ThisRedshiftList, pdf):
           
    xlim=[-0.5,2.5]
    ylim=[-0.5, 2.5]   
    bin=[0.05,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                         'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

    fig = plt.figure(figsize=(12,4))
    grid = gridspec.GridSpec(1, 4)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.1f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        if ii==0:
            ylabel='U-V'
        else:
             ylabel=''
        subplot.set_xlabel('V-J', fontsize=16), subplot.set_ylabel(ylabel, fontsize=16)
        
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.25)) 
        
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
        
        subplot.text(xlim[0]+2.,ylim[0]+.2,'z='+char_redshift, fontsize=16, fontweight='normal')
        
        if ii==2:
            subplot.text(xlim[0]+0.2,ylim[0]+0.2,'UVJ cut', fontsize=16, fontweight='normal')
                  
        #MODEL
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)
        
        G0_MR=G_MR[sel]   
        G0_MR=G0_MR[(G0_MR['MagDust'][:,0]<99.) & (G0_MR['MagDust'][:,2]<99.) & (G0_MR['MagDust'][:,7]<99.)]        
        color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
        color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
        Ngals=len(G0_MR)
      
        H, xedges, yedges = np.histogram2d(color_VJ, color_UV, bins=Nbins)            
        extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
        plt.subplots_adjust(bottom=0.15, left=0.15)        
        mylevels = np.linspace(1., Nbins[0], Nbins[0])*Ngals/(Nbins[0]**2/2.5)        
        H = zoom(H, 20)        
        cont=plt.contourf(H.transpose()[::], origin='lower', cmap='Greys_r', levels=mylevels, extent=extent)        
        
        if ii==len(ThisRedshiftList)-1:
            plt.colorbar(format='%d') 
        
        #BestFit Cut   
        bin=0.01
        slope=slope_red_fraction[ii+1]
        offset=offset_red_fraction[ii+1]
        minimum_y=minimum_y_red_fraction[ii+1]
        
        x_arr=np.arange(xlim[0],xlim[1]+bin,bin)
        cut1=np.zeros(len(x_arr),dtype=np.float32)+minimum_y      
        cut2=x_arr*slope+offset
  
        sel1=x_arr<((minimum_y-offset)/slope)
        subplot.plot(x_arr[sel1],cut1[sel1],color='red', linestyle='-', linewidth=2)  
         
        sel2=x_arr > ((minimum_y-offset)/slope)
        subplot.plot(x_arr[sel2],cut2[sel2],color='red', linestyle='-', linewidth=2)  
        
        #OBSERVATIONAL CUT    
        Nbin=0.01
        slope=0.88
        if(ThisRedshiftList[ii]<1.):
            offset=0.69 
        else: 
            offset=0.59

        x_arr=np.arange(xlim[0],xlim[1]+Nbin,Nbin)
        cut1=np.zeros(len(x_arr),dtype=np.float32)+1.3   
        cut2=x_arr*slope+offset
  
        sel1=x_arr<((1.3-offset)/slope)
        subplot.plot(x_arr[sel1],cut1[sel1],color='blue', linestyle='--', linewidth=2)  
  
        sel2=x_arr>((1.3-offset)/slope) 
        subplot.plot(x_arr[sel2],cut2[sel2],color='blue', linestyle='--', linewidth=2)  
        
        
        #LABELS
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.15, y_percentage=0.15, color='black', xlog=0, ylog=0, 
                        label='Muzzin 2013', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.04, y_percentage=0.17, color='blue', x2_percentage=0.13, 
                        xlog=0, ylog=0, linestyle='--', linewidth=2)
        
            plot_label (subplot, 'label', xlim, ylim, 
                        x_percentage=0.15, y_percentage=0.07, color='black', xlog=0, ylog=0, 
                        label='Best Fit cut', fontsize=13, fontweight='normal') 
            plot_label (subplot, 'line', xlim, ylim,
                        x_percentage=0.04, y_percentage=0.09, color='red', x2_percentage=0.13, 
                        xlog=0, ylog=0, linestyle='-', linewidth=2)
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_UVJ_colour.pdf')
    pdf.savefig()
    plt.close()
#endif UVJ_colour




def redfraction_color_cut(G_MR, ThisRedshiftList, pdf):
           
    xlim=[8.5,11.5]
    ylim=[0., 1.2]   
    bin=[0.1,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]

    plt.rcParams.update({'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 
                         'xtick.minor.width': 1.0, 'ytick.minor.width': 1.0})

    fig = plt.figure(figsize=(15,4))
    grid = gridspec.GridSpec(1, 5)
    grid.update(wspace=0.0, hspace=0.0)

    for ii in range(0,len(ThisRedshiftList)):
        
        char_redshift="%0.2f" % ThisRedshiftList[ii]
        
        subplot=plt.subplot(grid[ii])
        subplot.set_ylim(ylim), subplot.set_xlim(xlim)
        
        xlab='$\mathrm{log_{10}}(M_*[h^{-2}M_{\odot}])$'      
        if ii==0:
            ylab='Red Fraction'
        else:
            ylab=''      
        subplot.set_xlabel(xlab, fontsize=16), subplot.set_ylabel(ylab, fontsize=16)
        
        #subplot.text(xlim[0]+0.1,ylim[0]+.875,'z='+char_redshift, fontsize=16, fontweight='normal')
                        
    
        majorFormatter = FormatStrFormatter('%d')
        subplot.xaxis.set_major_locator(MultipleLocator(1))    
        subplot.xaxis.set_minor_locator(MultipleLocator(0.25))      
        subplot.yaxis.set_minor_locator(MultipleLocator(0.1))
    
        if ii>0:
            plt.tick_params(axis='y', which='both', left='on', labelleft='off')
                   
        #OBSERVATIONS             
        file = MCMCdir + '/ObsConstraints/RedFraction_z'+char_redshift+'.txt'        
        f = open(file, 'r')     
        line = int(f.readline())     
        obs = Table.read(file, format='ascii', data_start=1, data_end=line+1)
        
        obs_xbin=obs['col1']+(obs['col2']-obs['col1'])/2.       
        subplot.errorbar(obs_xbin, obs['col3'],obs['col4'],
                 fmt='o', markersize=5, ecolor='blue', color='blue')
        #sub = plt.subplot(111)
        
        #PREVIOUS MODELS 
        RedshiftList_OldModels=[0.1,0.4,1.,2.,3.0]
        old_char_redshift="%0.2f" % RedshiftList_OldModels[ii]
        if do_previous_model1==1: 
            file = file_previous_model1+'_redfrac_colorcut_z'+old_char_redshift+'.txt' 
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model1, linewidth=2)           
        if do_previous_model2==1: 
            file = file_previous_model2+'_redfrac_colorcut_z'+char_redshift+'.txt'  
            model = Table.read(file, format='ascii')
            subplot.plot(model['col1'],model['col2'],color='red',linestyle=linestyle_previous_model2, linewidth=2)
        
        #MODEL
        bin=0.25        
        Mass_arr=np.arange(xlim[0],xlim[1],bin)
        RedFraction=np.zeros(len(Mass_arr),dtype=np.float32)
        
        (sel)=select_current_redshift(G_MR, ThisRedshiftList, ii)        
        G0_MR=G_MR[sel]   
        StellarMass=stellar_mass_with_err(G0_MR, Hubble_h, ThisRedshiftList[ii])
        
        #z=0.
        if ThisRedshiftList[ii]==0:
            color_ur=G0_MR['MagDust'][:,15]-G0_MR['MagDust'][:,17]  
            Magr=G0_MR['MagDust'][:,17]
            
            for ll in range(0,len(Mass_arr)):
                sel_red=G0_MR[(color_ur>(offset_red_fraction[ii]-slope_red_fraction[ii]*np.tanh((Magr+18.07)/1.09))) &
                    (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                sel_all=G0_MR[(StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                
                RedFraction[ll]=float(len(sel_red))/float(len(sel_all))
                
        #z>0.
        else:
            color_UV=G0_MR['MagDust'][:,0]-G0_MR['MagDust'][:,2]  
            color_VJ=G0_MR['MagDust'][:,2]-G0_MR['MagDust'][:,7]       
                  
            for ll in range(0,len(Mass_arr)):
                sel_red=G0_MR[(((color_VJ < (minimum_y_red_fraction[ii]-offset_red_fraction[ii])/slope_red_fraction[ii]) &
                               (color_UV > minimum_y_red_fraction[ii])) |                      
                              ((color_VJ > (minimum_y_red_fraction[ii]-offset_red_fraction[ii])/slope_red_fraction[ii]) &
                               (color_UV > (color_VJ*slope_red_fraction[ii] + offset_red_fraction[ii])))) &
                              (StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                
                sel_all=G0_MR[(StellarMass>Mass_arr[ll]-bin/2.) & (StellarMass<Mass_arr[ll]+bin/2.)]
                if len(sel_all)>0.:
                    RedFraction[ll]=float(len(sel_red))/float(len(sel_all))
                else:
                    RedFraction[ll]=0.
                    
        subplot.plot(Mass_arr, RedFraction, color='red', linestyle='-', linewidth=2) 
      
    
    
        #LABELS    
        if ii==0:
            plot_label (subplot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.91, 
                        color='black', xlog=0, ylog=0, label='Observations used in MCMC', 
                        fontsize=12, fontweight='normal') 
            plot_label (subplot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.935, 
                        color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.03) 
        
        plot_label (subplot, 'label', xlim, ylim, x_percentage=0.05, y_percentage=0.65, 
                    color='black', xlog=0, ylog=0, label='z='+char_redshift[:-1], 
                    fontsize=14, fontweight='normal') 
            
        if ii==len(ThisRedshiftList)-1:
            plot_label_three_models (subplot, xlim, ylim, position='top_left')
                        
    #endfor
        
    plt.tight_layout()
    plt.savefig('./fig/plots_redfraction_color_cut.pdf')
    pdf.savefig()
    plt.close()
#endif redfraction_color_cut
