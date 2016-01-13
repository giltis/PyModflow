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
from pymodflow.pyobs import obs_tools

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
        bc_trans = trans[0,:,:].copy()

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
	
def update_param_fields(update_these_params=None):
    '''Ingests the current parameter estimates and generates updated
    parameter fields. 'update_which_params' is a dictionary with key=param name
    and values={'pval_file':str,'zone_array':real array}.'''

    model_param_field_dict = {}
    
    for iparam in update_these_params:
            
        print '    %s . . .' %(iparam)
        
        ipval_file = update_these_params[iparam]['pval_file']
        izone_array = update_these_params[iparam]['zone_array']
        
        _,ipval_zone_dict = obs_tools.pval_to_dict(ipval_file)         # Ingest the parameter information
        ifield = obs_tools.pval_to_field(ipval_zone_dict,izone_array)       # Update the parameter field and write the layer files

        model_param_field_dict[iparam] = ifield

    return model_param_field_dict
    
def update_bc_dicts(update_which_bc=None,ibound=None,landsurface=None,landuse=None):
    '''Compiles the updated boundary condition dictionaries and returns
    a dictionary of dictionaries. 'update_which_bc is a dict with
    key=bc_type ('DRN','GBH','RCH') and value=None or array (e.g, 'DRN'
    and 'GHB' require an array that serves as the basis of the
    conductance - e.g, the hk or trans arrays; 'RCH' requires the rch array).'''

    model_bc_dict = {}
        
    landsurface_to_bc = np.copy(landsurface)
    landsurface_to_bc[ibound[0,:,:] == 0] = np.nan    
        
    conductance_dy = lay_thick[0]
    print 'Conductancy dy = %.3f' %(conductance_dy)
    
    # Mask the land surface elevation array in order to reduce that array to only those
    # cells to which each boundary condition should be applied
    
    # The .ghb cells are specified using the ghb_landuse_dict
    ghb_stages = np.ma.MaskedArray(landsurface_to_bc,~np.in1d(landuse,ghb_landuse_dict.keys()))
    ghb_stages = ghb_stages.filled(np.nan)
    
    # The .drn cells are all active, non-constant head cells EXCEPT for those
    # cells specified as .ghb cells
    drn_stages = landsurface_to_bc
    drn_stages[np.isfinite(ghb_stages)] = np.nan # Exclude the ghb cells from the drn field      

    if ('DRN' in update_which_bc):
        drain_hk = update_which_bc['DRN']
        drain_hk = drain_hk[0,:,:]
        drn_dict = mfpt.get_bc_dict(drn_stages,nrow,ncol,delr,delc,nper=1,hk=drain_hk,conductance_dy=conductance_dy)
        model_bc_dict['DRN'] = drn_dict
    
    if ('GHB' in update_which_bc):
        ghb_hk = update_which_bc['GHB']
        ghb_hk = ghb_hk[0:,:,:]
        ghb_dict = mfpt.get_bc_dict(ghb_stages,nrow,ncol,delr,delc,nper=1,hk=ghb_hk,conductance_dy=conductance_dy)
        model_bc_dict['GHB'] = ghb_dict

    if ('RCH' in update_which_bc):
        rch_update = update_which_bc['RCH']
        rch_dict = {0:rch_update}
        model_bc_dict['RCH'] = rch_dict
    
    return model_bc_dict

def add_hob_to_namefile(nam_file=None,hob_file=None,heads_sim_file=None,hob_unit=60,heads_sim_unit=61):
    '''Adds the .hob package and the corresponding simulated heads to the namefile.'''
    
    with open(nam_file,'a') as fout:
        
        fout.write('%-14s%-3i%-s\n' %('HOB',hob_unit,hob_file))
        fout.write('%-14s%-3i%-s\n' %('DATA',heads_sim_unit,heads_sim_file))
        
    return

def reformat_sim_heads_file(heads_order,heads_fin,heads_fout):
    '''Reads the MODFLOW-generated simulated heads output and reformats.'''
    
    idf = pd.read_table(heads_fin,delim_whitespace=True)
    idf.columns = ['sim','obs','name']
    idf['name'] = pd.Categorical(idf['name'],heads_order)
    idf = idf.sort_values(by='name')
    
    with open(heads_fout,'w') as fout:
        for idx,irow in idf.iterrows():            
            iname,isim,iobs = irow['name'],irow['sim'],irow['obs']
            fout.write('%-15s%15.6e%15.6e\n' %(iname,isim,iobs))
        
    return

def get_dissolved_nitrate_dict(fertilizer_dir=None,rch_dict=None,new_loading_df=False,download_nass=False,\
                               nass_csv=None,loading_csv=None,nitrate_atm_csv=None,county_name=None,\
                               cdl=None,cdl_apply_list=None):
    '''Returns a dictionary with key=year (datetime object) and values=(nrow,ncol)
    array of dissolved nitrate concentrations.'''
    
    if (new_loading_df == True):
        county_df = county_nitrogen.get_county_nitrogen(fertilizer_dir=fertilizer_dir,download_nass=download_nass,\
                                                        nass_csv=nass_csv,loading_csv=loading_csv,nitrate_atm_csv=nitrate_atm_csv,\
                                                        county_name=county_name)
    else:
        county_df = pd.DataFrame.from_csv(county_nitrogen_csv)
        
    nitrate_array_dict,nitrate_avg_df = apply_nitrogen.apply_county_nitrogen(county_df,rch_dict=rch_dict,cdl=cdl,cdl_apply_list=cdl_apply_list,moving_average_window=1,\
                                                                             nrow=nrow,ncol=ncol,delr=delr,delc=delc,nitrate_rch_pickle=nitrate_rch_pickle)
     
    return nitrate_array_dict