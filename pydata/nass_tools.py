# -*- coding: utf-8 -*-
"""
Created on Tue Oct 07 10:45:51 2014

NASS PARAMETER NAMES:
--------------------
state_name
domain_desc
prodn_practice_desc
end_code
week_ending
location_desc
country_code
year
county_code
statisticcat_desc
agg_level_desc
county_ansi
CV (%)
short_desc
asd_code
commodity_desc
asd_desc
domaincat_desc
freq_desc
country_name
load_time
state_fips_code
class_desc
county_name
group_desc
Value
sector_desc
unit_desc
watershed_code
reference_period_desc
zip_5
state_ansi
util_practice_desc
begin_code
source_desc
congr_district_code
state_alpha
watershed_desc
region_desc

@author: Wesley Zell
"""

import pandas as pd
import numpy as np

# --- START PARAMETER SET ---

NASS_api = '1B149FA6-B96A-3B9F-98C5-C554D6FF64DD'

# --- STOP PARAMETER SET ----

def remove_commas(ival):
    '''Remove the comma separator and convert to a float.'''
    
    try:
        return float(ival.replace(',',''))
    except:
        return np.nan

def nass_to_dataframe(start_year,co_fips):
    '''Reads time series from NASS database to a pandas dataframe.'''
    
    url = 'http://quickstats.nass.usda.gov/api/api_GET/?key=' + NASS_api + '&year__GE=' + str(start_year) + '&state_alpha=MD&county_code=0' + str(co_fips) 
    df = pd.io.json.read_json(url,orient='split')
    
    df['Value'] = df['Value'].apply(remove_commas)

    return df
    
