# -*- coding: utf-8 -*-
"""
Main function and associated helper functions for generate the nitrogen
loading time series (agricultural inputs + atmospheric wet deposition)
for user-specified county.

Created on Fri Oct 03 15:22:29 2014
@author: Wesley Zell
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
from pymodflow.pydata import nass_tools as nass
from datetime import datetime

# ==================================
# --- START SCRIPT PARAMETER SET ---

fips_dict =      {'Kent':{'State':'MD','FIPS':'24029'}}
fips_dict.update({'QueenAnnes':{'State':'MD','FIPS':'24035'}})
fips_dict.update({'Talbot':{'State':'MD','FIPS':'24041'}})
fips_dict.update({'Caroline':{'State':'MD','FIPS':'24011'}})

uptake_dict = {'CORN, GRAIN':0.90,'CORN, SILAGE':9.7,'SOYBEANS':1.1,'WHEAT':1.5}    # lbs nitrogen per bushel grain or ton silage; from Murrell (2008)
for i in uptake_dict:
    uptake_dict[i] = uptake_dict[i] * (1./2.216) # Convert lbs to kg

# Constants
hectares_in_acre = 1./2.47105
lb_to_kg = 1./2.216

def get_time_bounds():
    '''Defines the dates that govern input/output time series.
    Wrapping the dates in this function in order to avoid automatic
    import with script import.'''

    start_year,fertilizer_start_year,stop_year = 1930,1935,2005
    
    return start_year,fertilizer_start_year,stop_year

# --- STOP SCRIPT PARAMETER SET ----
# ==================================

def usgs_45to85(idir=None,isubdir='USGS_1945_1985_Fertilizer',ifips=None):
    '''HELPER FUNCTION. Nitrogen Fertilizer Use: 1945-1985 (Alexander and Smith, 1990).
    Data stored as .dat files with leading space in records.'''

    alexander_smith_dir = os.path.join(idir,isubdir)
    alexander_smith_infiles = [x for x in os.listdir(alexander_smith_dir) if x.endswith('DAT')]
    
    # Load the first file to the dataframe    
    ifile = alexander_smith_infiles[0]
        
    iyears = ifile.split('.')[1]
    iyear_start,iyear_stop = int(iyears.split('-')[0]),int(iyears.split('-')[1]) + 1
    icolumns = ['FIPS'] + list(np.arange(iyear_start,iyear_stop,1))
    
    itemp_df = pd.read_table(os.path.join(alexander_smith_dir,ifile),delim_whitespace=True,header=None)
    itemp_df.columns = icolumns

    alexander_df = itemp_df
    
    # Merge the subsequent files to the dataframe        
    for ifile in alexander_smith_infiles[1:]:
    
        iyears = ifile.split('.')[1]
        iyear_start,iyear_stop = int(iyears.split('-')[0]),int(iyears.split('-')[1]) + 1
        icolumns = ['FIPS'] + list(np.arange(iyear_start,iyear_stop,1))
        
        itemp_df = pd.read_table(os.path.join(alexander_smith_dir,ifile),delim_whitespace=True,header=None)
        itemp_df.columns = icolumns
        
        alexander_df = pd.merge(alexander_df,itemp_df,on='FIPS')
    
    alexander_df = alexander_df.set_index('FIPS').transpose()
    alexander_df.index.names = ['Year']
    alexander_df = alexander_df[int(ifips)]     # Reduce to only the FIPS of interest
    alexander_df.name = 'Fertilizer (kg)'
    alexander_df = pd.DataFrame(alexander_df)

    return alexander_df

def usgs_82to01(idir=None,isubdir='USGS_1982_2001_Fertilizer',ifin='Nutrient_Inputs_1982-2001.NitrogenTables.xls'):
    '''HELPER FUNCTION. Nitrogen Fertilizer Use: 1982-2001 (Ruddy, Lorenz, and Mueller, 2006).
    Data stored as multiple tables (i.e., worksheets) in Excel file.
    Sheetnames = ['Farm','NonFarm','Livestock_Confined','Livestock_Unconfined','Atm_Deposition'].'''

    ruddy_dir = os.path.join(idir,isubdir)
    ruddy_fin = os.path.join(ruddy_dir,ifin)
    
    ruddy_farm_df = pd.io.excel.read_excel(ruddy_fin,sheetname='Farm',columns=0,index_col=0)
    
    return ruddy_farm_df

def construct_fips_string(irow):
    '''HELPER FUNCTION. Constructs the fips string from separate state
    and county fips identifiers. Problematic rows may be ignored.'''
    
    try:
        return str(int(irow['FIPS_ST'])) + str(0) + str(int(irow['FIPS_CO']))
    except:
        return '00000'

def usgs_87to06(idir=None,isubdir='USGS_1987_2006_Fertilizer',ifin='tblFarmNonfarmCountyNitrogen.xlsx',\
                isheet='tblFarmNonfarmCountyNitrogen',ifips=None):
    '''HELPER FUNCTION. Nitrogen Fertilizer Use: 1987-2006 (Gronberg and Spahr, 2012).
    Data originally stored in Access database; exported to Excel sheet due to
    problems accessing the MS Access database using 64-bit Python.
    Single table for nitrogen inputs: 'tblFarmNonfarmCountyNitrogen'.
    Attributes: 'FIPS_ST','FIPS_CO','STATE','CO','farmN1987','nonfN1987',etc. (1987-2006).'''    

    gronberg_dir = os.path.join(idir,isubdir)
    gronberg_fin = os.path.join(gronberg_dir,ifin)
    
    # Read the county level farm and nonfarm fertilizer sales 1987-2006 to a dataframe
    gronberg_df = pd.read_excel(gronberg_fin,isheet)
    
    # Derive the fips string from the individual state/co fips columns
    # and set as the index
    gronberg_df['FIPS'] = gronberg_df.apply(construct_fips_string,axis=1)
    gronberg_df = gronberg_df.set_index('FIPS')
    
    # Split the farm and nonfarm information into separate series   
    all_columns = gronberg_df.columns
    farm_columns = [x for x in all_columns if 'farm' in x]
    nonfarm_columns = [x for x in all_columns if 'nonf' in x]
    
    gronberg_farm_df = gronberg_df[farm_columns]
    gronberg_nonfarm_df = gronberg_df[nonfarm_columns]
    
    # Reduce the column names to years (this is required for plotting)    
    farm_columns = [x.split('N')[1] for x in farm_columns]
    nonfarm_columns = [x.split('N')[1] for x in nonfarm_columns]
    
    gronberg_farm_df.columns = farm_columns
    gronberg_nonfarm_df.columns = nonfarm_columns

    gronberg_farm_df = gronberg_farm_df.transpose()
    gronberg_farm_df = gronberg_farm_df[ifips]
    
    gronberg_farm_df.index.names = ['Year']
    gronberg_farm_df.name = 'Fertilizer (kg)'
    gronberg_farm_df = pd.DataFrame(gronberg_farm_df)
   
    return gronberg_farm_df
    
def sanford_45to04(idir=None,isubdir='SanfordPope',ifin='SanfordPopeNitrateInputs.xlsx',ifips=None):
    '''HELPER FUNCTION. Nitrogen Fertilizer Use: 1945-2004 (Sanford and Pope, 2012).
    Data stored in Excel Spreadsheet. Two loading tables (sheetnames = 'Fertilizer_Loading',
    'Poultry_Loading'. Attributes = 'Year','24029','24035' (i.e., FIPS for 
    Kent Co and Queen Anne Co, respectively).
    Sheetnames = ['Fertlizer_Loading','Poultry_Loading','Concentration']'''

    sanford_dir = os.path.join(idir,isubdir)
    sanford_fin = os.path.join(sanford_dir,ifin)
    
    # Rename the columns for the county of interest
    sanford_fertilizer_df = pd.read_excel(sanford_fin,sheetname='Fertilizer_Loading',columns=0).rename(columns={int(ifips):'Fertilizer (kg x 10^6)'})
    sanford_poultry_df = pd.read_excel(sanford_fin,sheetname='Poultry_Loading',columns=0).rename(columns={int(ifips):'Poultry (kg x 10^6)'})
    
    # Create a new dataframe that totals the fertilizer and poultry loads
    county_df = pd.merge(sanford_fertilizer_df,sanford_poultry_df,on='Year')
    county_df['Total_Load'] = county_df['Fertilizer (kg x 10^6)'] + county_df['Poultry (kg x 10^6)']

    return county_df

def get_fips(county_name):
    '''HELPER FUNCTION. Returns fips information.'''
    
    st_name = fips_dict[county_name]['State']
    fips = fips_dict[county_name]['FIPS']
    st_fips,co_fips = int(fips.split('0')[0]),int(fips.split('0')[1])
    
    return fips,st_name,st_fips,co_fips
    
def plot_county_nitrogen(idf=None,county_name=None,plot_fout=None,fontsize=8):
    '''Plots the nitrate loading time series for the user-specified county.'''
    
    plot_from_script = os.path.basename(__file__)
    fips,st_name,st_fips,co_fips = get_fips(county_name)
    start_year,fertilizer_start_year,stop_year = get_time_bounds()
    
    plt.figure()
    plt.plot(idf.index,idf['CORN, GRAIN - YIELD, MEASURED IN BU / ACRE'])
    plt.title('Corn Yield')
    plt.show()
    
    plt.figure()
    plt.plot(idf.index,idf['CORN, GRAIN - PRODUCTION, MEASURED IN BU'])
    plt.title('Corn Production')
    plt.show()
    
    fig,ax1 = plt.subplots(figsize=(8,6),dpi=300)
    ax2 = ax1.twinx()
    
    for iax in [ax1,ax2]:
        iax.tick_params(axis='both',labelsize=fontsize)
    
    # Plot the county-level aggregated fertilizer sales and poultry manure use
    ax1.plot(idf.index,idf['Gross_Ag_In (kg)'] * 1./1000000.,label='Fertilizer + Poultry Inputs',c='y',lw=1.5)
    #ax1.plot(idf.index,idf['Poultry (kg)'] * 1./1000000.,label='Kent Co. Fertilizer Sales',c='y',lw=1.5)
    ax1.plot(idf.index,idf['Total_N_Export (kg)'] * 1./1000000.,label='Estimated Corn + Soybean + Wheat Export',c='g',lw=1.5)
    ax1.set_ylabel('Nitrogen Inputs (kg x 10$^6$)',fontsize=fontsize)
    
    # Plot the number of acres planted for corn and soybeans
    
    ax2.plot(idf.index,idf['CORN, GRAIN - ACRES HARVESTED'] * 1/1000,'b--',lw=1,label='Acres HARVESTED in Corn')
    ax2.plot(idf.index,idf['SOYBEANS - ACRES HARVESTED'] * 1/1000,'b:',lw=1,label='Acres HARVESTED in Soybeans')
    #ax2.plot(idf.index,idf.Total_Acres * 1/1000,'r:',lw=1,label='Acres planted in Soybeans or Corn')
    #ax2.plot(idf.index,idf['CropFraction'] * 20,'k--')
    ax2.set_ylabel('Cultivated area (acres x 10$^3$)',fontsize=fontsize)
    
    # Ask matplotlib for the plotted objects and their labels in order to combine to single legend
    lines1,labels1 = ax1.get_legend_handles_labels()
    lines2,labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2,labels1 + labels2,loc='upper left',fontsize=fontsize,frameon=False)
    
    plt.xlabel('Year')
    plt.xlim(datetime(start_year,01,01),datetime(stop_year,01,01))
    plt.title('%s\nLand Surface Nitrogen Loading for %s County, %s' %(plot_from_script,county_name,st_name),fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(plot_fout,format='png',dpi=300)
    
    return

#########################################################
def get_county_nitrogen(fertilizer_dir=None,use_nass_cache=True,nass_csv=None,loading_csv=None,nitrate_atm_csv=None,county_name=None):
    '''Main function. Returns (as a dataframe) a time series (annual time steps) of COUNTY-LEVEL
    nitrogen inputs (from fertilizer and atmospheric wet deposition data) and
    nitrogen exports (from NASS records of crop harvest and export).'''
#########################################################

    fips,st_name,st_fips,co_fips = get_fips(county_name)
    start_year,fertilizer_start_year,stop_year = get_time_bounds()

    # -----------------------
    # --- NITROGEN INPUTS ---
    # -----------------------
    
    # Retrieve the county level fertilizer information
    # from various data sources. Note that functions (above)
    # allow for retrieval of multiple data sources, some of which
    # are redundant and not used here.  This script creates a composite
    # time series for the specified county using Alexander/Smith + Gronberg/Spahr
    
    alexander_df = usgs_45to85(idir=fertilizer_dir,ifips=fips)
    gronberg_farm_df = usgs_87to06(idir=fertilizer_dir,ifips=fips)
    county_df = pd.concat([alexander_df,gronberg_farm_df])
    
    # 1986 is missing from time series.  Fill with linear interpolation.
    county_df.index = county_df.index.astype(int)
    county_df.loc[1986] = np.nan
    county_df = county_df.sort_index()
    county_df = county_df.interpolate()
    
    # Retrieve the atmospheric deposition data.
    atm_df = pd.DataFrame.from_csv(nitrate_atm_csv)
    atm_df = atm_df[['Year','NO3']]
    atm_df['NO3'] = atm_df['NO3'] * hectares_in_acre
    atm_df = atm_df.rename(columns={'NO3':'WetDep (kg/acre)'})
    atm_df = atm_df.set_index('Year')
    
    # Retrieve the poultry manure information compiled by Sanford and Pope
    poultry_df = sanford_45to04(idir=fertilizer_dir,ifips=fips)
    poultry_df = poultry_df[['Year','Poultry (kg x 10^6)']]
    poultry_df['Poultry (kg x 10^6)'] = poultry_df['Poultry (kg x 10^6)'] * 1000000.
    poultry_df = poultry_df.rename(columns={'Poultry (kg x 10^6)':'Poultry (kg)'})
    poultry_df = poultry_df.set_index('Year')

    # Add the poultry and the atmospheric deposition to the dataframe
    county_df = pd.merge(county_df,atm_df,left_index=True,right_index=True,how='left')
    county_df = pd.merge(county_df,poultry_df,left_index=True,right_index=True,how='left')
    
    # Convert index to datetime    
    county_df.index = pd.to_datetime(county_df.index,format='%Y')
    
    # ------------------------
    # --- NITROGEN EXPORTS ---
    # ------------------------
    
    # Retrieve the county level NASS statistics on acres planted and harvested
    
    if (use_nass_cache == True):
        nass_df = pd.DataFrame.from_csv(nass_csv)
    else:
        nass_df = nass.nass_to_dataframe(start_year,co_fips)
        nass_df.to_csv(nass_csv)
    
    nass_df = nass_df[nass_df['commodity_desc'].isin(['CORN','SOYBEANS','WHEAT'])][['year','Value','short_desc']].rename(columns={'year':'Year'})
    
    col_strings_1 = ['WHEAT - ACRES PLANTED','CORN - ACRES PLANTED','SOYBEANS - ACRES PLANTED','CORN, SILAGE - PRODUCTION, MEASURED IN TONS']
    col_strings_temp = ['- ACRES HARVESTED','- PRODUCTION, MEASURED IN BU','- YIELD, MEASURED IN BU / ACRE']
    col_strings_2 = ['CORN, GRAIN ' + x for x in col_strings_temp]
    col_strings_3 = ['SOYBEANS ' + x for x in col_strings_temp]
    col_strings_4 = ['WHEAT ' + x for x in col_strings_temp]
    col_strings = col_strings_1 + col_strings_2 + col_strings_3 + col_strings_4
    
    for istring in col_strings:
            
        temp_df = nass_df[nass_df['short_desc'] == istring].rename(columns={'Value':istring})[['Year',istring]]
        
        if 'crop_df' not in locals():
            crop_df = temp_df
    
        else:
            crop_df = pd.merge(crop_df,temp_df,on='Year',how='outer')
    
    # Fix the comma delineation; this includes replacing the '(D)' flag (indicating data that was
    # not posted).  This seems to particularly (exclusively?) apply to Talbot county CENSUS data
    # for 2007 and 2012; for our commodities of interest (soybeans, wheat, and corn), these flags
    # get replaced by the SURVEY data for those same years
    
    # --- start helper function ---
    def reformat_floats(irow):
        
        if '(D)' in str(irow):
            return np.nan
        
        else:
            return float(str(irow).replace(',',''))
    # --- stop helper function ----
        
    for icol in crop_df.columns:
        if (icol == 'Year'):
            continue
        
        crop_df[icol] = crop_df[icol].apply(reformat_floats)
    
    # Some manipulations are needed in order to deal with the effects of duplicate years (presumably due
    # to census + other form of information; or as an artifact of the merge) and set the dates to first day of year
    crop_df['Year'] = crop_df['Year'].astype(int)
    crop_df = crop_df.groupby('Year').agg(np.mean)
    crop_df = crop_df.sort_index()
    crop_df.index = pd.to_datetime(crop_df.index,format='%Y')
    
    # Some data is missing; fill with the minimum reported data, except for soybean
    # acreage: start time series with minimum reported value and linearly interpolate
    for icol in col_strings:
        
        if 'SOYBEAN' in icol:
            crop_df[icol].iloc[0] = crop_df[icol].min() # Assigns the minimum value in the series to the first year
            crop_df[icol] = crop_df[icol].interpolate()
        
        else:
            crop_df[icol] = crop_df[icol].fillna(crop_df[icol].min())
            
    # Calculate the nitrogen uptake by corn, silage, soybeans, and wheat.
    # NOTE: this assumes that silage nitrogen leaves the catchment as livestock
    crop_df['Corn_Export'] = crop_df['CORN, GRAIN - PRODUCTION, MEASURED IN BU'] * uptake_dict['CORN, GRAIN']
    crop_df['Silage_Export'] = crop_df['CORN, SILAGE - PRODUCTION, MEASURED IN TONS'] * uptake_dict['CORN, SILAGE']
    crop_df['Soybean_Export'] = crop_df['SOYBEANS - PRODUCTION, MEASURED IN BU'] * uptake_dict['SOYBEANS']
    crop_df['Wheat_Export'] = crop_df['WHEAT - PRODUCTION, MEASURED IN BU'] * uptake_dict['WHEAT']
    crop_df['Total_N_Export (kg)'] = crop_df['Corn_Export'] + crop_df['Soybean_Export'] + crop_df['Wheat_Export'] + crop_df['Silage_Export']
    
    # Total acreage
    crop_df['Total_Acres'] = crop_df['CORN - ACRES PLANTED'] + crop_df['SOYBEANS - ACRES PLANTED'] + crop_df['WHEAT - ACRES PLANTED']
    
    # ------------------------------
    # --- NET = INPUTS - OUTPUTS ---
    # ------------------------------
    
    county_df = pd.merge(county_df,crop_df,left_index=True,right_index=True,how='outer')
    
    # There is more data reported from NASS than from fertilizer use (which is only reported 1945 to 2006).
    # Supplement fertilizer and wet deposition loading series with linear interpolation and save to csv.
    # Set the poultry input to 
    for icol in ['Fertilizer (kg)','WetDep (kg/acre)']:
        county_df[icol].loc[[datetime.strptime(str(x),'%Y') for x in range(start_year,fertilizer_start_year)]] = 0
    
    county_df['Poultry (kg)'] = county_df['Poultry (kg)'].fillna(county_df['Poultry (kg)'].min())
    county_df['Poultry (kg)'].iloc[-1] = county_df['Poultry (kg)'].max()
    county_df = county_df.interpolate()
    
    county_df['Gross_Ag_In (kg)'] = county_df['Fertilizer (kg)'] + county_df['Poultry (kg)']
    county_df['Net_Ag_In (kg)'] = county_df['Gross_Ag_In (kg)'] - county_df['Total_N_Export (kg)']
    county_df['CropFraction'] = county_df['Total_N_Export (kg)']/county_df['Gross_Ag_In (kg)']
    
    # Cache the loading to a csv
    county_df.to_csv(loading_csv)
    
    return county_df