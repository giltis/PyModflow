# -*- coding: utf-8 -*-
"""
Tools for processing time series of atmospheric concentrations
into DISSOLVED MODEL UNITS.

Created on Tue Jun 10 21:13:08 2014
@author: Wesley Zell
"""

import os
from scipy import exp,log
from datetime import datetime, timedelta
import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import date2num

# ---------------------------
# --- START PARAMETER SET ---

# *** Specify the location of the atmospheric time series ***
default_atm_tracer_data_dir = 'C:\\Data_National\\AtmTracers'
default_atm_tracer_xls = os.path.join(default_atm_tracer_data_dir,'atm_tracers.xlsx')

mole_fraction_multiplier = 1e-12
to_pg_per_liter = 1e12  # Concentration conversion from g/liter to default model units = picograms/liter

# Plotting params
iformat,font_size,lw,ms,dpi = 'png',8,.5,3.5,300
matplotlib.rcParams.update({'font.size': font_size})
default_atm_tracer_plot = os.path.join(default_atm_tracer_data_dir,'atm_tracer.' + iformat)

# --- STOP PARAMETER SET ----
# ---------------------------

def model_conc(moles_per_liter,iFW,conc_multiplier):
    '''Converts moles per liter to model concentration units.  The base
    conversion: moles_per_liter * FW = g/liter.  The conc_multiplier argument
    subsequently converts g/liter to model units.'''
    
    return moles_per_liter * (iFW) * conc_multiplier
    
def henry_k(row,temp,salinity):
    '''Generates Henry's constant as function of chemical parameters and temperature and salinity.'''
       
    return exp(row.A1 + row.A2*(100/temp) + row.A3*log(temp/100) + salinity*(row.B1 + row.B2*(temp/100) + row.B3*((temp/100)**2)))

def molar_conc(mole_fraction,ihenry_k,temp,excess_air,altitude_meters,salinity):
    '''Converts mole fraction to moles per liter water.'''
    
    barometric_pressure = exp(-(altitude_meters/8300)) 
    vapor_pressure_h2o = exp(24.4543 - 67.4509*(100/temp) - 4.8489*log(temp/100) - 0.000544*salinity)
    
    return (ihenry_k*mole_fraction*(barometric_pressure - vapor_pressure_h2o) + (excess_air*mole_fraction)/22414)

def normalize_tracer_input_df(df):
    '''Returns a tracer input dataframe with each species time series normalized
    on a scale of 0 to 1.'''
    
    for itracer in ['SF6','CFC113','CFC11','CFC12','Trit']:
        
        iseries = df[itracer]
        iseries_max = iseries.max()
        df[itracer] = df[itracer]/iseries_max
        
    return df

def convert_decimal_year_to_datetime(decimal_year):
    '''Converts a float representation of a year to a datetime object.'''

    year = int(decimal_year)
    rem = decimal_year - year
    
    dt_base = datetime(year, 1, 1)
    dt_obj = dt_base + timedelta(seconds=(dt_base.replace(year=dt_base.year + 1) - dt_base).total_seconds() * rem)

    return dt_obj

def convert_years_to_plot_date(years):
    '''Converts a list of integer years to matplotlib date format.'''
    
    years = [str(x) for x in years]
    years = [datetime.strptime(x,'%Y') for x in years]
    years = [date2num(x) for x in years]
    
    return years

# ----------

def plot_tracer_signals(full_atm_df,annual_atm_df,plot_fig):
    '''Saves the ATMOSPHERIC (not dissolved) tracer input signals to a figure.
    Assumes that the input dataframe uses the following units: pptv SF6 and CFC; TU Tritium.'''
    
    full_atm_df.index.name = 'Date'
    full_atm_df = full_atm_df.reset_index()
    full_atm_df = full_atm_df.rename(columns={'index':'Date'})
    full_atm_df.Date = full_atm_df.Date.apply(convert_decimal_year_to_datetime)
    full_atm_df.Date = pd.to_datetime(full_atm_df.Date)
    full_atm_df.Date = full_atm_df.Date.apply(date2num)
    
    annual_atm_df.index.name = 'Date'
    annual_atm_df = annual_atm_df.reset_index()
    annual_atm_df.Date = annual_atm_df.Date.apply(date2num)
    
    xticks = convert_years_to_plot_date(np.arange(1940,2020,10))
    xlim = convert_years_to_plot_date([1940,2010]) 
    
    fig = plt.figure(figsize=(5,4))
    ax1 = fig.add_subplot(111)       
    ax2 = ax1.twinx()
    
    sols = ['CFC12','CFC11','CFC113','SF6','Trit']
    plot_sols = ['CFC-12','CFC-11','CFC-113','SF$_6$','$^3$H']
    mfcs = ['w','k','w','','']
    mss  = ['o','o','^','','']
    
    for isol,iplot_sol,mfc,ms in zip(sols,plot_sols,mfcs,mss):
        
        full_isol,full_date = full_atm_df[isol],full_atm_df.Date
        annual_isol,annual_date = annual_atm_df[isol],annual_atm_df.Date
        
        if (isol == 'Trit'):
            ax2.plot_date(full_date,full_isol,ls='-',lw=2,c='0.75',marker=None)
            ax2.plot_date(annual_date,annual_isol,'k-',lw=1.5,label=iplot_sol)

        elif (isol == 'SF6'):
            full_isol = full_isol * 100
            annual_isol = annual_isol * 100
            ax2.plot_date(annual_date,annual_isol,'k--',lw=1.5,label=iplot_sol)
        
        else:
            ax1.plot_date(annual_date,annual_isol,'k-',marker=ms,mfc=mfc,lw=1.5,label=iplot_sol,markevery=20)

    ax1.set_ylabel('Atmospheric CFC Concentration (pptv)')
    ax2.set_ylabel('Atmospheric SF$_6$ Concentration (pptv x 100)\nAtmospheric Tritium Concentration (TU)')
    plt.xticks(xticks) 
    plt.xlim(xlim)
    plt.xlabel('Year')
    
    # Ask matplotlib for the plotted objects and their labels in order to combine to single legend
    lines1,labels1 = ax1.get_legend_handles_labels()
    lines2,labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2,labels1 + labels2,loc='upper left',frameon=False)   
    
    plt.tight_layout()
    plt.savefig(plot_fig,format=iformat,dpi=dpi)
    
    return

def get_dissolved_tracer_df(tracer_xls=default_atm_tracer_xls, \
               excess_air=None,recharge_temp_celsius=None,altitude_meters=None,salinity = 0, \
               plot_flag=False,plot_fig=default_atm_tracer_plot,conc_multiplier=to_pg_per_liter):
    '''Returns a time series of dissolved tracer concentrations. The
    'conc_multiplier' argument converts g/liters to DISSOLVED MODEL UNITS (default
    dissolved model units for CFCs and SF6 are picograms/liter). NOTE: THIS FUNCTION
    DEPENDS ON PARAMETERS THAT ARE DEFINED IN THIS SCRIPT.'''

    temp = recharge_temp_celsius + 273.15  # Convert Celsius to Kelvin
    
    # Read in the Henry's K Parameters and the atmospheric concentration time series
    henrys_df = pd.io.excel.read_excel(tracer_xls,'HenrysKParams',header=0)
    CFC_SF6_pptv_df = pd.io.excel.read_excel(tracer_xls,'CFC_SF6_Atmosphere',header=0)
    Trit_TU_df = pd.io.excel.read_excel(tracer_xls,'Tritium_Atmosphere',header=0)
    
    # Calculate Henry's K for each solute given the specified temperature and salinity
    henrys_df['Henrys_K'] = henrys_df.apply(henry_k,axis=1,args=(temp,salinity))
    
    # Convert the pptv time series to mole_fraction time series then
    # convert the mole_fraction time series to moles per liter    
    CFC_SF6_dissolved_df = CFC_SF6_pptv_df.copy()    
    
    for itracer in ['SF6','CFC113','CFC11','CFC12']:
    
        iFW = henrys_df[henrys_df.TracerName == itracer]['FW']               
        ihenry_k = henrys_df[henrys_df.TracerName == itracer]['Henrys_K']                                   # extract the appropriate Henrys K
        CFC_SF6_dissolved_df[itracer] = CFC_SF6_dissolved_df[itracer] * mole_fraction_multiplier            # convert pptv to mole fraction
        CFC_SF6_dissolved_df[itracer] = CFC_SF6_dissolved_df[itracer].apply(molar_conc,args=(ihenry_k,temp,excess_air,altitude_meters,salinity))    # convert mole fraction to moles per liter
        CFC_SF6_dissolved_df[itracer] = CFC_SF6_dissolved_df[itracer].apply(model_conc,args=(iFW,conc_multiplier))

    # Create the annually-averaged input series that will be used by the calibration process.    
    # If multiple atmospheric observations in a single year, compute yearly average
    # and reduce to integer years. NOTE: the aggregation converts the join key (i.e., 'Date') to the index
    CFC_SF6_dissolved_annual_df = CFC_SF6_dissolved_df.copy()
    CFC_SF6_dissolved_annual_df.Date = CFC_SF6_dissolved_annual_df.Date.apply(np.round).astype(int).astype(str) # Convert each date to a year
    CFC_SF6_dissolved_annual_df = CFC_SF6_dissolved_annual_df.groupby('Date').agg({'CFC11':'mean','CFC12':'mean','CFC113':'mean','SF6':'mean'})
    CFC_SF6_dissolved_annual_df.index = pd.to_datetime(CFC_SF6_dissolved_annual_df.index,format='%Y')

    Trit_TU_annual_df = Trit_TU_df.copy()
    Trit_TU_annual_df.Date = Trit_TU_annual_df.Date.apply(np.round).astype(int).astype(str) # Convert each date to a year
    Trit_TU_annual_df = Trit_TU_annual_df.groupby('Date').agg({'Trit':'mean'})
    Trit_TU_annual_df.index = pd.to_datetime(Trit_TU_annual_df.index,format='%Y')   

    CFC_SF6_pptv_annual_df = CFC_SF6_pptv_df.copy()
    CFC_SF6_pptv_annual_df.Date = CFC_SF6_pptv_annual_df.Date.apply(np.round).astype(int).astype(str) # Convert each date to a year
    CFC_SF6_pptv_annual_df = CFC_SF6_pptv_annual_df.groupby('Date').agg({'CFC11':'mean','CFC12':'mean','CFC113':'mean','SF6':'mean'})
    CFC_SF6_pptv_annual_df.index = pd.to_datetime(CFC_SF6_pptv_annual_df.index,format='%Y')    

    Trit_TU_df = Trit_TU_df.set_index('Date')
    CFC_SF6_pptv_df = CFC_SF6_pptv_df.set_index('Date')
    
    # Generate three dataframes for the following time series: 
    # Annually-averaged series of dissolved concentrations (for model input),
    # full series of atmospheric concentrations, and annually-averaged series
    # of atmospheric concentrations
    full_atm_df = CFC_SF6_pptv_df.join(Trit_TU_df).fillna(0)
    annual_atm_df = CFC_SF6_pptv_annual_df.join(Trit_TU_annual_df).fillna(0)
    
    if (plot_flag == True):
        plot_tracer_signals(full_atm_df,annual_atm_df,plot_fig)
    
    annual_dissolved_df = CFC_SF6_dissolved_annual_df.join(Trit_TU_annual_df).fillna(0)
            
    return annual_dissolved_df