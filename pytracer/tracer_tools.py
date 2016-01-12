# -*- coding: utf-8 -*-
"""
Tools for working with atmospherically-derived tracers.
Note that the script for getting the atmospheric 

Created on Thu Feb 20 11:41:41 2014
@author: Wesley Zell
"""

import numpy as np
import pandas as pd
import os

# -----------------------------

# def normalize_tracer(tracer_df)

        # # Normalize the concentration and standard deviation measurements.
        # if (normalize_tracer_flag == True):
                        
            # imax_value = dissolved_input_df[jsol].max()     # Normalization multiplier = largest value in tracer atmospheric input series
           
            # if (isol == 'HelTrit'): # This is required to apply the multiplier to the HelTrit data
                # isol_df[iconc_name] = isol_df[iconc_name] * 1./imax_value
                # isol_df[istd_name] = isol_df[istd_name] * 1./imax_value                        
    
            # else:
                # isol_df[jconc_name] = isol_df[jconc_name] * 1./imax_value
                # isol_df[jstd_name] = isol_df[jstd_name] * 1./imax_value



#####################################
# START CHESTER TRACER DATA FUNCTIONS
#####################################

def get_conversion_factor(isol):
    '''Returns the multiplier necessary to convert the reported solute concentration
    units to model units.'''
    
    if isol == 'SF6':                                   # (SF6 reported as fmol/L)    
        FW = 146.0504
        conversion_factor = FW * (.001)                 # Converts fmol/L to pg/L

    if isol in ['CFC11','CFC12','CFC113']:              # (CFCs reported as pg/kg)    
        conversion_factor = 1                           # = pg/kg to pg/L

    if isol in ['Trit','HelTrit']:    
        conversion_factor = 1
        
    return conversion_factor
    
def get_chester_tracer_df(tracer_data_xls):
    '''Reads the Upper Chester atmospheric tracer information to a dataframe from the spreadsheet
    that synthesizes the Plummer Busenberg data.'''

    plum_bus_df = pd.io.excel.read_excel(tracer_data_xls,'RawData',header=0)

    # Identify the solutes represented in the dataframe;
    # Remove the columns containing interpreted age; these are indicated
    # by the string 'age'
    column_headers = list(plum_bus_df)
    
    # Extract the solute concentrations    
    tob_columns = [x for x in column_headers if 'age' not in x]
    tob_df = plum_bus_df[tob_columns]

    # Extract the AGES as reported by Plummer and Busenberg
    age_columns = [x for x in column_headers if 'Name' in x or 'Date' in x or 'age' in x]
    plummer_age_df = plum_bus_df[age_columns]
    
    return tob_df,plummer_age_df

def split_species(tob_df=None,tracer_data_dir=None,min_tracer_uncertainty=None):
    '''Writes a separate .csv for each tracer species.'''
        
    column_headers = list(tob_df)
    solutes = list(set([x.split('_')[0] for x in column_headers if ('Name' not in x and 'Date' not in x and 'age' not in x)]))
    
    for isol in solutes:

        iconc_name = isol + '_conc'    
        istd_name = iconc_name + '_std'
        
        # Collapse the dataframe to observations of the solute of interest
        isol_df = tob_df[['SiteName','Date',iconc_name,istd_name]].copy()
        isol_df = isol_df.dropna(subset=[iconc_name])
        isol_df['SiteName'] = isol_df['SiteName'].astype(str)       
              
        # Convert the concentration and standard deviation measurements.        
        iconvert = get_conversion_factor(isol)  # This converts to model units
        isol_df[iconc_name] = isol_df[iconc_name] * iconvert
        isol_df[istd_name] = isol_df[istd_name] * iconvert
    
        # These adjustments required since the HelTrit and Trit data both use
        # the same Tritium input signal
        if (isol == 'HelTrit'):
            jsol = 'Trit'
            jconc_name = iconc_name.replace('HelTrit','Trit')
            jstd_name = istd_name.replace('HelTrit','Trit')
        else:
            jsol = isol            
            jconc_name = iconc_name
            jstd_name = istd_name    
     
        # The standard deviation, if provided and non-zero, will be used to generate the
        # observation weight during parameter estimation. If no standard deviation
        # is provided, the average standard deviation is used. PEST and UCODE require non-zero weights;
        # therefore, if provided standard deviation is 0, it is replaced with the minimum standard deviation for all
        # measurements of that species.
        isol_df[istd_name] = isol_df[istd_name].fillna(isol_df[istd_name].mean())                       # Replace NaN with mean
        isol_df[istd_name] = isol_df[istd_name].fillna(isol_df[iconc_name] * min_tracer_uncertainty)    # If all NaNs (as with Tritium), replace with user-specified ratio of measurement
    
        # Get the minimum nonzero standard deviation  
        istd_list = isol_df[istd_name].values
        istd_min = np.min(istd_list[np.nonzero(istd_list)])
        isol_df[istd_name] = isol_df[istd_name].replace(to_replace=[0],value=istd_min)  # Replace 0 with minimum
    
        # Write the solute information to csv
        print 'Writing %s data to %s.' %(isol,os.path.join(tracer_data_dir,isol + '.csv'))
        isol_df.columns = ['station_nm','Date',isol,'Std']
        isol_df.to_csv(os.path.join(tracer_data_dir,isol + '.csv'))
    
    return

def write_chester_csvs(tracer_data_xls=None,tracer_data_dir=None,min_tracer_uncertainty=0.05):
    '''Writes the Plummer Busenber tracer data for the Upper Chester to species
    specific .csv files that are subsequently read by the observation utilities
    (e.g., write_pest.py). CSV VALUES ARE WRITTEN IN MODEL UNITS (i.e., the conversion
    to model units happens in this script).  No normalization in this function.'''
    
    # Write the Plummer/Busenberg information to dataframe
    tob_df,age_df = get_chester_tracer_df(tracer_data_xls)

    # Split the dataframe by species and write to .csv
    split_species(tob_df,tracer_data_dir,min_tracer_uncertainty)

    return

#####################################
# STOP CHESTER TRACER DATA FUNCTIONS
#####################################