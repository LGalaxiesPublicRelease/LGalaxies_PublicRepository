# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

"""
read_snap
read_tree
redshift_to_time
select_current_redshift
def stellar_mass_with_err
grayify_cmap
"""
    
    
          
""" 
self.capacity *= 4
newdata = np.zeros((self.capacity,))
newdata[:self.size] = self.data
self.data = newdata
self.data[self.size] = x
     
data = self.data[:self.size]
return np.reshape(data, newshape=(len(data)/5, 5))"""


import numpy as np
from plots_input import *
import matplotlib.pyplot as plt

def read_snap(folder,FirstFile,LastFile,
              props,template,RedshiftsToRead,FullRedshiftList):    
    """ Reads L-Galaxy output files.
    Returns: (nTrees,nHalos,nTreeHalos,gals)
    Inputs: (folder,file_prefix,FirstFile,LastFile,props,template)
    props - list of properties to return
    template - structure dtype definition from database """     
    nTrees = 0
    nGals = 0    
    nTreeHalos = np.array([],dtype=np.int32)
    
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)  
            
    SnapshotList=np.array([],dtype=np.int32)
    
    #read only headers to figure out total nGals
    print ("\n\nReading Headers\n")
    for iredshift in range(0,len(FullRedshiftList)-1):
        if RedshiftsToRead[iredshift]:              
            for ifile in range(FirstFile,LastFile+1):
                char_redshift="%0.2f" % FullRedshiftList[iredshift]
                filename = folder+'/'+'SA_z'+char_redshift+"_"+"%d"%(ifile)               
                f = open(filename,"rb")
                
                this_nTrees =  np.fromfile(f,np.int32,1)
                nTrees += this_nTrees
                this_nGals = np.fromfile(f,np.int32,1)
                nGals += this_nGals
                
            print ("z=", char_redshift," nGals = ",nGals)  
               
    gals = np.zeros(nGals,dtype=filter_dtype)
  
    print ("\n")
    offset=0
    for iredshift in range(0,len(FullRedshiftList)-1):
        if RedshiftsToRead[iredshift]:  
            print ("\nReading redshift: ", FullRedshiftList[iredshift], "\n")
            for ifile in range(FirstFile,LastFile+1):
                char_redshift="%0.2f" % FullRedshiftList[iredshift]
                filename = folder+'/'+'SA_z'+char_redshift+"_"+"%d"%(ifile)
                #print(filename)
                f = open(filename,"rb")
                
                this_nTrees =  np.fromfile(f,np.int32,1)
                nTrees += this_nTrees
                this_nGals = np.fromfile(f,np.int32,1)
                nGals += this_nGals
                print ("File ", ifile," nGals = ",this_nGals)  
                
                addednTreeHalos = np.fromfile(f,np.int32,this_nTrees)
                nTreeHalos = np.append(nTreeHalos,addednTreeHalos)
                full_this_gals = np.fromfile(f,template,this_nGals) # all properties
                this_gals = np.zeros(this_nGals,dtype=filter_dtype) # selected props
                
                for prop in template.names:
                    if props[prop]:
                        this_gals[prop] = full_this_gals[prop]
                              
                gals[offset:offset+this_nGals] = this_gals[:]    
                offset+=this_nGals
                f.close()           
            #endfor
        #endif    
        #assign snapshot of current redshift given by the last galaxy on the last file 
        SnapshotList=np.append(SnapshotList,gals['SnapNum'][offset-1])       
    #endfor
    
    return (gals, SnapshotList)




def read_tree(folder,FirstFile,LastFile,
              props,template):    
    """ Reads L-Galaxy output files.
    Returns: (nTrees,nHalos,nTreeHalos,gals)
    Inputs: (folder,file_prefix,FirstFile,LastFile,props,template)
    props - list of properties to return
    template - structure dtype definition from database """   
    nGals = 0    
    
    filter_list = []
    for prop in props:
        if props[prop]:
            filter_list.append((prop,template[prop]))
    filter_dtype = np.dtype(filter_list)  
            
    SnapshotList=np.array([],dtype=np.int32)
    
    #read only headers to figure out total nGals
    print ("\n\nReading Headers\n")
    for ifile in range(FirstFile,LastFile+1):       
        filename = folder+'/'+'SA_galtree_'+"%d"%(ifile)               
        f = open(filename,"rb")
        one = np.fromfile(f,np.int32,1)
        nbytes = np.fromfile(f,np.int32,1)
        this_nGals = np.fromfile(f,np.int32,1)
        nGals += this_nGals               
    gals = np.zeros(nGals,dtype=filter_dtype)
    
    print("TotNgals=",nGals)
    print ("\n")
    
    offset=0
    for ifile in range(FirstFile,LastFile+1):         
        filename = folder+'/'+'SA_galtree_'+"%d"%(ifile)             
        f = open(filename,"rb")
        one = np.fromfile(f,np.int32,1)
        nbytes = np.fromfile(f,np.int32,1)  
        nskip=nbytes/4-3
        this_nGals = np.fromfile(f,np.int32,1)      
        nGals += this_nGals       
        print ("File ", ifile," nGals = ",this_nGals)      
        ib=np.fromfile(f,np.float32,int(nskip))          
       
        full_this_gals = np.fromfile(f,template,this_nGals) # all properties      
        this_gals = np.zeros(this_nGals,dtype=filter_dtype) # selected props
               
       
        for prop in template.names:
            if props[prop]:
                this_gals[prop] = full_this_gals[prop]
                              
        gals[offset:offset+this_nGals] = this_gals[:]    
        offset+=this_nGals
        f.close()           
    #endfor
   
    return (gals)

def redshift_to_time (z):
    Tyr = 977.8    ;# coefficent for converting 1/H into Gyr                               
    WM = 0.315
    WV = 0.685
    H0=67.3
    h = H0/100.
    WR = 4.165E-5/(h*h)   ;# includes 3 massless neutrino species, T0 = 2.72528            
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n=1000        ; # number of points in integrals                                        
    a=0
    for i in range(0, n-1):
        a = az*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot
    
    zage = az*age/n
    age_Gyr = (Tyr/H0)*zage
    
    return (age_Gyr)
#end redshift_to_time

def select_current_redshift(G_MR, ThisRedshiftList, ii):
    
    found_redshift=0
                           
    for jj in range(0, len(FullRedshiftList)):           
        if(ThisRedshiftList[ii]<1.):
            if round(FullRedshiftList[jj],1)==round(ThisRedshiftList[ii],1): 
                sel= (G_MR['SnapNum']==FullSnapshotList[jj])
                found_redshift=1                  
        else:    
            if round(FullRedshiftList[jj],0)==round(ThisRedshiftList[ii],0): 
                sel= (G_MR['SnapNum']==FullSnapshotList[jj])
                found_redshift=1 
                    
    if found_redshift==0:
        sys.exit("redshift:",ThisRedshiftList[ii],"needed for stellar mass function not read.") 
        
    return (sel)

#end  select_current_redshift


def stellar_mass_with_err(G0_MR, Hubble_h, redshift):
    
    mass= np.log10(G0_MR['StellarMass']*1.e10*Hubble_h) + np.random.randn(len(G0_MR['StellarMass']))*0.08*(1+redshift)

    return mass

#end get_stellar_mass


def median_and_percentiles (bin, xmin, xmax, x_variable, y_variable): 

    min_x=min(x_variable[x_variable > (-1e20)])
    if(min_x > xmin):
        xmin=min_x
    
    Nbins=(xmax-xmin)/bin+1

    median=np.zeros(Nbins, np.float32)
    pc16=np.zeros(Nbins, np.float32) 
    pc84=np.zeros(Nbins, np.float32)  
    x_min=xmin-bin/2.
    
    for ii in range(0,int(Nbins)): 

        x_variable_sel=x_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
        y_variable_sel=y_variable[(x_variable > (x_min+(ii*bin))) & (x_variable < (x_min+(ii+1.0)*bin))] 
       
        if(len(x_variable_sel) > 0):           
            median[ii]=np.median(y_variable_sel) 
            y_sorted = np.sort(y_variable_sel)
            pc16[ii] = y_sorted[16*len(y_variable_sel)/100]      
            pc84[ii] = y_sorted[84*len(y_variable_sel)/100]  
    #endfor

    x_binned=np.arange(Nbins)*((x_min+(Nbins*bin))-(x_min+(Nbins*0.0)))/(Nbins*1.)+x_min+bin/2.
 
    return (x_binned, median, pc16, pc84)

#end median_and_percentiles


def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]
    
    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)

#end grayify_cmap


def plot_label_three_models (subplot, xlim, ylim, position):
 
    if position=='top_left':
        x1=0.15        
        x21=0.04
        x22=0.13
        
        previous_model1_y1=0.9
        previous_model1_y2=0.92       
        previous_model2_y1=0.83
        previous_model2_y2=0.85        
        this_model_y1=0.76
        this_model_y2=0.78
        
    if position=='bottom_left':
        x1=0.15        
        x21=0.04
        x22=0.13
        
        previous_model1_y1=0.2
        previous_model1_y2=0.22       
        previous_model2_y1=0.13
        previous_model2_y2=0.15        
        this_model_y1=0.06
        this_model_y2=0.08
        
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=x1, y_percentage=previous_model1_y1, color='black', xlog=0, ylog=0, 
                label=prefix_previous_model1, fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=x21, y_percentage=previous_model1_y2, color='red', x2_percentage=x22, 
                xlog=0, ylog=0, linestyle=linestyle_previous_model1, linewidth=2)  
                    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=x1, y_percentage=previous_model2_y1, color='black', xlog=0, ylog=0, 
                label=prefix_previous_model2, fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=x21, y_percentage=previous_model2_y2, color='red', x2_percentage=x22, 
                xlog=0, ylog=0, linestyle=linestyle_previous_model2, linewidth=2)
                    
    plot_label (subplot, 'label', xlim, ylim, 
                x_percentage=x1, y_percentage=this_model_y1, color='black', xlog=0, ylog=0, 
                label=prefix_this_model, fontsize=10, fontweight='normal') 
    plot_label (subplot, 'line', xlim, ylim,
                x_percentage=x21, y_percentage=this_model_y2, color='red', x2_percentage=x22, 
                xlog=0, ylog=0, linestyle='-', linewidth=2)
                
#endf plot_label_three_models
            
                
def plot_label (subplot, label_type, xlim, ylim, x_percentage, y_percentage, color, 
                x2_percentage=0., xlog=0, ylog=0, label='', linestyle='-', linewidth=2, 
                fontsize=16, fontweight='normal', sym='o', sym_size=5, err_size=0.1):    
    
     if xlog==0 & ylog==0:
      
         if label_type=='label':
             x=xlim[0]+(xlim[1]-xlim[0])*x_percentage
             y=ylim[0]+(ylim[1]-ylim[0])*y_percentage             
             subplot.text(x,y,label, fontsize=fontsize, fontweight=fontweight)
         else:
             if label_type =='line':
                 x1=xlim[0]+(xlim[1]-xlim[0])*x_percentage
                 x2=xlim[0]+(xlim[1]-xlim[0])*x2_percentage
                 y=ylim[0]+(ylim[1]-ylim[0])*y_percentage
                 subplot.plot([x1,x2],[y,y],color=color,linestyle=linestyle, linewidth=2)
                
             else:
                 if label_type=='symbol':
                     x=xlim[0]+(xlim[1]-xlim[0])*x_percentage
                     y=ylim[0]+(ylim[1]-ylim[0])*y_percentage                     
                     subplot.errorbar(x,y, yerr=err_size, fmt=sym, markersize=sym_size, color=color)
                  
 
#end plot_label


# <codecell>


# <codecell>


