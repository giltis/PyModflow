# -*- coding: utf-8 -*-
"""
Tools for working with MODFLOW packages.

Created on Mon Dec 14 16:42:06 2015
@author: wzell
"""

import numpy as np
import pandas as pd
import os
import flopy.utils.binaryfile as bf

import matplotlib.pyplot as plt

def process_drains(budget_binary,baseflow_sim_file,ibound,discharge_ibound_map=None,drain_discharge_plot=None):
    '''Aggregates the drain discharge by zone and writes to file. IMPROVEMENTS TO MAKE:
    Save the plot of active drains to file for each time step.  Use pdf pages or animate.'''

    ibound = ibound[0,:,:]  # Only using the land surface to discretize to zones
    zones  = [int(izone) for izone in np.unique(ibound) if izone > 0]
    
    cbc = bf.CellBudgetFile(budget_binary)
    nper = cbc.kstpkper[-1][1]
    drns = cbc.get_data(text='DRAINS',full3D=True)

    with open(baseflow_sim_file,'w') as fout:
        
        for iper in range(nper):
            for izone in zones:
                
                iname = izone
                
                if (discharge_ibound_map is not None):
                    if (izone in discharge_ibound_map):
                        iname = discharge_ibound_map[izone]
                    else:
                        continue
                
                izone_drns = drns[iper][0,:,:][ibound == izone]
                izone_sum = np.sum(izone_drns) * -1           
                fout.write('%-15s%20.6e\n' %(iname,izone_sum))
    
    return

def get_bc_dict(bc_stages,nrow,ncol,delr,delc,nper=1,trans=None,decay_rate=None,hk=None,conductance_dy=None):
    '''Returns a dictionary with key=stress period and values=list of lists['Lay','Row','Col','Stage','Conductance'].
    Read by the flopy .drn and .ghb package constructors.
    bc_stages = (nrow,ncol) array with bc cells = land surface elevation and all other cells = np.nan'''
    
    # Generate the arrays of row,col indices    
    xx,yy = np.meshgrid(np.arange(ncol),np.arange(nrow))
    icols,irows = xx.ravel(),yy.ravel()
    cell_area = (delr * delc)    
    
    # The boundary stage is equal to the land surface elevation at each cell    
    stages = bc_stages.ravel()
    
    if (trans is not None):
        
        # Derive the conductance from the transmissivity field.
        # This is the tool used for the national model and its derivatives.
        bc_trans = trans.copy()

        if (conductance_dy is not None):
            multiplier = (decay_rate * cell_area)/(2 * np.log(2) * conductance_dy)
        else:
            multiplier = (decay_rate ** 2) * (cell_area/(2 * (np.log(2)) ** 2))
    
        # Mask the transmissivity array so that its nans correspond to the nans
        # in the bc_stages array.    
        bc_trans[np.isnan(bc_stages)] = np.nan
        conductances = (bc_trans * multiplier).ravel()
        
    elif (hk is not None):
        
        # Derive the conductance from the conductivity field.
        bc_hk = hk.copy()
        
        multiplier = cell_area * (1./(conductance_dy))
        
        # Mask the conductivity array so that its nans correspond to the nans
        # in the bc_stages array.    
        bc_hk[np.isnan(bc_stages)] = np.nan
        conductances = (bc_hk * multiplier).ravel()
        
    # Organize with a dataframe
    bc_df = pd.DataFrame(np.vstack((irows,icols,stages,conductances)).T,columns=['Row','Col','Stage','Conductance'])
    bc_df['Lay'] = 0
    bc_df = bc_df[['Lay','Row','Col','Stage','Conductance']]
    
    # Remove the NaN conductances (which should = inactive cells)
    bc_df = bc_df.dropna()
    
    bc_list = bc_df.values.tolist()    
    bc_dict = {}
    
    for i in range(nper):
        bc_dict[i] = bc_list

    return bc_dict

def write_openclose_layers(a,file_root=None,file_ext=None,fmt=None,write_dir=None):
    '''Writes a separate 2D array to file for each layer in a (potentially) 3D array. This is required
    because flopy does not currently support the EXTERNAL array read option.'''   
    
    if (a.ndim > 2):    # 3D array           
        for ilay in range(np.shape(a)[0]):            
            ifile = os.path.join(write_dir,'{}_{}.{}'.format(file_root,ilay+1,file_ext))
            np.savetxt(ifile,a[ilay,:,:],fmt=fmt)
            
    else:   # For one layer: return a single file name rather than a list of files       
        ifile = os.path.join(write_dir,'{}_1.{}'.format(file_root,file_ext))
        np.savetxt(ifile,a,fmt=fmt)
            
    return
    
def fetch_openclose_layers(fetch_dir,fetch_root,nlay,fetch_ext,start_dir=os.getcwd()):
    '''Gets the list of RELATIVE PATHS for files that contain the layer information for 3D array. This is required
    because flopy does not currently support the EXTERNAL array read option.'''   
    
    file_list = []
    for ilay in range(nlay):
        
        ifile = os.path.join(fetch_dir,fetch_root + '_' + str(ilay+1) + '.' + fetch_ext)
        if os.path.exists(ifile):
            file_list.append(os.path.relpath(ifile,start_dir))
        else:
            quit('%s does not exist.' %(ifile))
            
    return file_list