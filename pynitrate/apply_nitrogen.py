# -*- coding: utf-8 -*-
"""
Processes the fertilizer time series dataframe from generate_nitrate_loading.py
and spatially distributes it per land use flags from the user specified NCDL dataset.
Saves a pickled dictionary with keys = year and values = (nrow,ncol) array of loading
MASS (i.e., not concentration) rates.  This pickled dictionary is read by
process_mpend_nitrate.py

Created on Tue May 19 17:32:12 2015
@author: Wesley Zell
"""

import numpy as np
import os,pickle
import pandas as pd
from scipy.stats import itemfreq
from datetime import datetime
from scipy.ndimage import gaussian_filter

import matplotlib.pyplot as plt
from matplotlib import gridspec

# ===================================
# --- START SCRIPT PARAMETER SET ----

# Put the unit conversions here in order
# to avoid import of constants with import of module
def get_unit_conversions(cell_area):
    
    acre_per_cell = 4046.86/cell_area    
    
    # Convert kg/(acre*year) to mg/(acre*year)
    convert_mass = 1000000. 
    
    # Convert m/day recharge to l/(acre*year) = days/year * m^2/cell * l/m^3 * cell/acre
    convert_volume = 365.25 * cell_area * 1000. * 1./acre_per_cell
    
    return convert_mass,convert_volume

# --- STOP SCRIPT PARAMETER SET -----
# ===================================

def summarize_cdl(cdl,cdl_xls='C:\\Data_National\\NASS\\CroplandDataLayer_Key.xls'):
    '''Returns a dataframe that summarizes cropland datalayer content
    for the input array (e.g., the model resolution cdl array). 'cdl_xls' is
    a locally stored legend for the cropland datalayer.'''

    cdl_key = pd.read_excel(cdl_xls)[['Code','Class Name']]

    cdl_sum = cdl_key.loc[np.unique(cdl)]
    cdl_sum['Count'] = cdl_sum['Code'].apply(lambda x: np.sum(cdl == int(x)))
    
    return cdl_sum

def plot_landuse():
    '''Plots one year of the NCDL and filtered application locations.'''
    
    plot_landuse = np.copy(loading_ts_dict[2008])
    gauss_landuse = gaussian_filter(loading_ts_dict[2008],pvl_dict['GaussianSigma'])
    
    plot_landuse[landuse_array == 0] = np.nan
    gauss_landuse[gauss_landuse == 0] = np.nan
    
    # Configure the plots
    fig = plt.figure(figsize=(10,10),dpi=300)
    gs = gridspec.GridSpec(1,2)
    
    ax_cdl = plt.subplot(gs[0:1])
    ax_cdl.imshow(plot_landuse,interpolation='nearest')
    ax_cdl.set_ylim(-10,150)
    ax_cdl.set_xlim(-10,ncol)
    ax_cdl.set_xticks([])
    ax_cdl.set_yticks([])
    ax_cdl.invert_yaxis()
    ax_cdl.set_xlabel('2008 NCDL\nLand Use',fontsize=fontsize)
    
    ax_gauss = plt.subplot(gs[1:2])
    ax_gauss.imshow(gauss_landuse,interpolation='nearest')
    ax_gauss.set_ylim(-10,150)
    ax_gauss.set_xlim(-10,ncol)
    ax_gauss.set_xticks([])
    ax_gauss.set_yticks([])
    ax_gauss.invert_yaxis()
    ax_gauss.set_xlabel('Filtered Loading Rates\n(kg/acre)',fontsize=fontsize)
    
    plt.suptitle(script_name,fontsize=fontsize * 0.8)
    plt.tight_layout()
    plt.savefig(plot_fout,format='png')
    
    return

def plot_applied_loading(county_df,nitrate_avg_df,fontsize=8,start_year=1930,stop_year=2005,plot_fout=None):
    '''Generates a composite plot of crop exports, agricultural inputs,
    and net recharging nitrogen input to the model domain.'''

    script_name = os.path.basename(__file__)

    fig,axes = plt.subplots(3,figsize=(12,6),dpi=300,sharex=True)    
    ax1,ax2,ax3 = axes
    
    for iax in axes:
        iax.tick_params(axis='both',labelsize=fontsize * 0.75)
    
    # Plot the crop data
    ax1.plot(county_df.index,county_df['CORN, GRAIN - ACRES HARVESTED'] * 1/1000,'y',lw=1,label='Corn')
    ax1.plot(county_df.index,county_df['SOYBEANS - ACRES HARVESTED'] * 1/1000,'g',lw=1,label='Soybeans',ls='--')
    ax1.plot(county_df.index,county_df['WHEAT - ACRES HARVESTED'] * 1/1000,'k',lw=1,label='Wheat',ls=':')
    ax1.set_ylabel('Harvested\nAcres\n(acres x 10$^3$)',fontsize=fontsize)
    
    # Plot the county-level aggregated fertilizer sales and poultry manure use
    ax2.plot(county_df.index,county_df['Gross_Ag_In (kg)'] * 1./1000000.,label='Fertilizer + Poultry Inputs',c='y',lw=1.5)
    ax2.plot(county_df.index,county_df['Total_N_Export (kg)'] * 1./1000000.,label='Corn + Soybean + Wheat Exports',c='g',lw=1,ls='--')
    ax2.set_ylabel('Agricultural\nNitrogen\n(kg x 10$^6$)',fontsize=fontsize)
    
    # Plot the recharge time-series
    ax3.plot(nitrate_avg_df.index,nitrate_avg_df['MaxAvgRchNitrate'],'r',lw=1.5,label='High Loading Scenario')
    ax3.plot(nitrate_avg_df.index,nitrate_avg_df['MinAvgRchNitrate'],'b',lw=1.5,label='Low Loading Scenario',ls='--')
    ax3.set_ylabel('Recharging\nAg Nitrate\n(mg/L)',fontsize=fontsize)
    ax3.set_ylim([0,75])
    
    plt.suptitle(script_name,fontsize=0.25*fontsize)    
    
    for iax,idx in zip(axes,['a','b','c']):
        iax.legend(loc='upper left',fontsize=0.75*fontsize,frameon=False)
        iax.set_title('(' + idx + ')',fontsize=fontsize*0.75)
    
    plt.xlabel('Year',fontsize=fontsize,labelpad=2)
    plt.xlim(datetime(start_year,01,01),datetime(stop_year,01,01))
    #plt.tight_layout()
    plt.savefig(plot_fout,format='png',dpi=300)
    
    return

def get_fraction_landuse(loading_df,landuse_array):
    '''Count the number of acres of each land use represented in the land use raster
    and determine what fraction that represents for the reported planted acres
    for the associated CDL year.  Note that itemfreq returns a tuple of (value,count).'''

    fraction_landuse_dict = dict()
    for i in itemfreq(landuse_array):
        if int(i[0]) in cdl_keep_dict:
            
            icode,icells = i
            iname = cdl_keep_dict[int(icode)]
            imodel_acres = icells * acres_per_cell
            
            itotal_county_acres = loading_df.loc[datetime(cdl_year,01,01)][iname + '_Acres']
            
            fraction_landuse_dict[iname] = imodel_acres/itotal_county_acres
    
    return fraction_landuse_dict

def get_fert_rate(irow):
    '''Not used by current scheme.  Returns the base rate of fertilization for a given year.  Invokes
    assumed distribution of fertilizer rates based on coefficients specified in the
    fertilizer_fraction_dict.  This rate must be muliplied by respective coefficent
    to generated the kg/acre fertilization rate for a particular crop in a particular year.'''
    
    fert_fract_dict = {'Corn':0.80,'Soybean':0.20}  # These must sum to 1    
    
    ifert_mass = irow['Fertilizer (kg)']
    icorn_acres,isoy_acres = irow['Corn_Acres'],irow['Soybean_Acres']
    icorn_coeff,isoy_coeff = fert_fract_dict['Corn'],fert_fract_dict['Soybean']
    
    base_rate = ifert_mass/(icorn_acres*icorn_coeff + isoy_acres*isoy_coeff)
    
    return base_rate

def mask_landuse(cdl,cdl_apply_list):
    '''Masks the landuses that DO NOT receive agricultural
    nitrogen inputs.'''

    cdl_apply = np.ma.MaskedArray(cdl,~np.in1d(cdl,cdl_apply_list)).filled(np.nan)
    
    return cdl_apply

def apply_loading_ts(loading_df,cdl_apply=None,nrow=None,ncol=None):
    '''Spatially distributeds the (possibly smoothed) time series
    of nitrogen inputs.'''
    
    idict = {}
    idict['max'] = {}
    idict['min'] = {}
    
    for iyear,irow in loading_df.iterrows():
        
        iwet = irow['WetDep (kg/acre)']
        imax = irow['Max_MovingAvg_Mass_Rate']
        imin = irow['Min_MovingAvg_Mass_Rate']
                
        iwet_dep = np.ones((nrow,ncol)) * iwet
                      
        imax_array = np.copy(iwet_dep) 
        imax_array[np.isfinite(cdl_apply)] += imax
            
        imin_array = np.copy(iwet_dep) 
        imin_array[np.isfinite(cdl_apply)] += imin
        
        idict['max'][iyear] = imax_array
        idict['min'][iyear] = imin_array

    return idict  

def dissolve_loading_ts(mass_rate_ts_dict,rch_dict,convert_mass,convert_volume):
    '''Dissolve the nitrate mass rate time series in the recharge time series.
    COULD BE ADJUSTED TO ACCOMODATE TRANSIENT RECHARGE. Currently assumes
    steady state.  Returns: (i) a dictionary with key=downscale method and
    value={year:array} and (ii) a dataframe of average recharging values.'''   
    
    idict = {}
    idict['max'] = {}
    idict['min'] = {}

    idf = pd.DataFrame(columns=['Year','MaxAvgRchNitrate','MinAvgRchNitrate'])
      
    ss_rch = rch_dict[0]
    
    icount = 0
    for iyear in mass_rate_ts_dict['max']:
        
        imax = mass_rate_ts_dict['max'][iyear]
        imin = mass_rate_ts_dict['min'][iyear]
        
        idict['max'][iyear] = (imax * convert_mass)/(ss_rch * convert_volume)
        idict['min'][iyear] = (imin * convert_mass)/(ss_rch * convert_volume)
                
        idf.loc[icount] = [iyear,np.mean(imax),np.mean(imin)]
        icount += 1
    
    idf = idf.set_index('Year')
    idf = idf.sort_index()
    
    return idict,idf
        
###############################################################################
def apply_county_nitrogen(county_df,rch_dict=None,cdl=None,cdl_apply_list=None,moving_average_window=1,\
                          nrow=None,ncol=None,delr=None,delc=None,nitrate_rch_pickle=None):
    '''Main function. Returns a time series of (nrow,ncol) arrays downscaled
    from the county-level loading dataframe (generated by county_nitrogen.py).

    Two downscaling methods (i.e., versions of the areal rate) are available.
    max_areal_rate = calculate using only the reported corn acreage.
    min_areal_rate = calculate using the total reported (corn + soybean + wheat) acreage.
    
    This function accepts a moving average window argument that will smooth
    the loading time series. Default is no smoothing (moving_average_window=1).'''
###############################################################################
    
    # Get the county level loads
    loading_df = county_df.copy()
    
    # Calculate the maximum and minimum loading rates and generate a time series of loading arrays
    loading_df['Max_Mass_Rate'] = loading_df['Net_Ag_In (kg)']/loading_df['CORN - ACRES PLANTED']
    loading_df['Min_Mass_Rate'] = loading_df['Net_Ag_In (kg)']/loading_df['Total_Acres']
    
    loading_df['Max_MovingAvg_Mass_Rate'] = pd.rolling_mean(loading_df['Max_Mass_Rate'],moving_average_window)
    loading_df.loc[loading_df['Max_MovingAvg_Mass_Rate'] < 0,'Max_MovingAvg_Mass_Rate'] = 0
    
    loading_df['Min_MovingAvg_Mass_Rate'] = pd.rolling_mean(loading_df['Min_Mass_Rate'],moving_average_window)
    loading_df.loc[loading_df['Min_MovingAvg_Mass_Rate'] < 0,'Min_MovingAvg_Mass_Rate'] = 0
    
    cdl_apply = mask_landuse(cdl,cdl_apply_list)
    mass_rate_ts_dict = apply_loading_ts(loading_df,cdl_apply,nrow=nrow,ncol=ncol)
    
    # Dissolve the nitrate load in recharge
    convert_mass,convert_volume = get_unit_conversions(delr*delc)
    nitrate_array_dict,nitrate_avg_df = dissolve_loading_ts(mass_rate_ts_dict,rch_dict,convert_mass,convert_volume)
    
    with open(nitrate_rch_pickle,'wb') as handle:
        pickle.dump(nitrate_array_dict,handle)
     
    return nitrate_array_dict,nitrate_avg_df