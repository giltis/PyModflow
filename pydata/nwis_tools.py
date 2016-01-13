# -*- coding: utf-8 -*-
"""
Reads NWIS data and site information for the model area to separate
site and data CSVs.  Data is aggregated and merged with site information
in write_pest.py.

Data dataframe column names:
----------------------------
agency_cd                  - Agency Code
site_no                    - USGS site number
sample_dt                  - Begin date
sample_tm                  - Begin time
sample_end_dt              - End date
sample_end_tm              - End time
sample_start_time_datum_cd - Time datum
tm_datum_rlbty_cd          - Time datum reliability code
coll_ent_cd                - Agency Collecting Sample Code
medium_cd                  - Sample Medium Code
tu_id                      - Taxonomic unit code
body_part_id               - Body part code
parm_cd                    - Parameter code
remark_cd                  - Remark code
result_va                  - Parameter value
val_qual_tx                - Result value qualifier code
meth_cd                    - Method code
dqi_cd                     - Data-quality indicator code
rpt_lev_va                 - Reporting level
rpt_lev_cd                 - Reporting level type
lab_std_va                 - Lab standard deviation
anl_ent_cd                 - Analyzing entity code

Site dataframe column names:
----------------------------
site_no         -- Site identification number
station_nm      -- Site name
site_tp_cd      -- Site type
c
coord_acy_cd    -- Latitude-longitude accuracy
dec_coord_datum_cd -- Decimal Latitude-longitude datum
alt_va          -- Altitude of Gage/land surface
alt_acy_va      -- Altitude accuracy
alt_datum_cd    -- Altitude datum
aqfr_cd         -- Local aquifer code
well_depth_va   -- Well depth
hole_depth_va   -- Hole depth
gw_count_nu     -- Field water-level measurements count

Parameters of potential interest:
---------------------------------
72000   Altitude of land surface, feet
72001   Depth of hole, feet below land surface datum
30210   Depth to water level, below land surface datum (LSD), meters
72002   Depth to top of water
99137   Nitrate, water, in situ, milligrams per liter as nitrogen
72008   Depth of well, feet below land surface datum
618     Nitrate, water, filtered, milligrams per liter as nitrogen
71851   Nitrate, water, filtered, milligrams per liter
620     Nitrate, water, unfiltered, milligrams per liter as nitrogen
63149   Sulfur hexafluoride, water, unfiltered, recoverable, femtograms per kilogram
72015   Depth to top of sample interval, feet below land surface datum
72016   Depth to bottom of sample interval, feet below land surface datum
99121   Nitrate, water, filtered, field, milligrams per liter as nitrogen
72003   Depth to bottom of water
630     Nitrate plus nitrite, water, unfiltered, milligrams per liter as nitrogen
631     Nitrate plus nitrite, water, filtered, milligrams per liter as nitrogen
72019   Depth to water level, feet below land surface
00300   Dissolved oxygen, mg/L
99114   Iron(II), water, filtered, field, milligrams per liter
99113   Sulfate, water, filtered, field, milligrams per liter

Created on Tue Oct 20 18:00:46 2015
@author: wzell
"""
import fiona,urllib2,json
import pandas as pd
from pymodflow.pygrid import grid_util as gu

# --- START SCRIPT PARAM SET ---

global nwis_root
nwis_root = 'http://waterservices.usgs.gov/nwis/'

# --- STOP SCRIPT PARAM SET ----

def get_bbox(domain_shp):
    '''Returns the bounding box as a list of coordinates in same projection
    as model domain shape.'''
    
    with fiona.open(domain_shp,'r') as extent:
        bbox = extent.bounds
        
    return bbox
    
def get_bbox_string(bbox):
    '''Returns a string of bbox coordinates for insertion in NWIS call.'''
    
    return ('%10.7f,%10.7f,%10.7f,%10.7f' %(bbox))

def df_from_url(iurl,sep=None,comment='#',header=0,is_wq=False):
    '''Reads a USGS waterservices rdb response to a dataframe.'''
    
    irdb = urllib2.urlopen(iurl)   
    idf = pd.read_table(irdb,sep=sep,comment=comment,header=header)
    if (is_wq == False):
        idf = idf.drop(idf.index[0])    # Drop the description of the field widths    
    
    return idf

# =============================================================================

def get_station_df(site_url,site_type,keep_station_cols=None,sep=None):
    ''' Returns a dataframe of site information.  Current methodology only
    gets site information for sites with head or discharge information.
    Consequently, this function is not used with calls to the water quality web services,
    and the results from water quality web services are subsequently reduced to only those
    sites with head or discharge information.'''
    
    print '\nDownloading NWIS station information for sites of type %s.' %(site_type)
    station_df = df_from_url(site_url,sep)

    # Remove any dataframe rows that were read from blank rows in the file    
    station_df['site_no'] = station_df['site_no'].convert_objects(convert_numeric=True)    
    station_df = station_df.dropna(subset=['site_no'])   
    
    if (keep_station_cols is not None):    
        station_df = station_df[keep_station_cols]
    
    # Groundwater stations
    if (site_type == 'GW'):
        # Further clean the data:
        # Remove rows with alt_va = NaN or well_depth_va = NaN and
        # convert strings to floats
        station_df = station_df.dropna(subset = ['alt_va','well_depth_va'])
    
        for icol in [['dec_long_va','dec_lat_va','alt_va','well_depth_va']]:
            station_df[icol] = station_df[icol].convert_objects(convert_numeric=True)

        # Add column: elevation of the well = land surface elevation - well_depth
        station_df['well_elev'] = station_df['alt_va'] - station_df['well_depth_va']
   
    # Discharge stations
    if (site_type == 'ST'):    
        for icol in [['dec_long_va','dec_lat_va']]:
            station_df[icol] = station_df[icol].convert_objects(convert_numeric=True)
    
    # Remove the spaces in the station names (required for use of station names by UCODE)
    # and clean up the station information dataframe
    station_df['station_nm'] = station_df['station_nm'].apply(lambda x: x.replace(" ",""))
    
    print 'Found %i stations.' %(len(station_df))
    print 'NOTE: All stations cached to .csv and must be filtered by subsequent operations.\n'
        
    return station_df

# =============================================================================
# Enough variety between head, discharge, and water quality formats to require
# type-specific helper functions.
# =============================================================================

def get_head_results(data_url,sep=None):
    '''HELPER FUNCTION. Returns a dataframe of groundwater level observations.'''
    
    data_df = df_from_url(data_url,sep)
    
    # Remove any dataframe rows that were read from blank rows in the file
    # (or, in some cases, are comment lines and headers for additional sites)
    data_df['site_no'] = data_df['site_no'].convert_objects(convert_numeric=True)    
    data_df = data_df.dropna(subset=['site_no'])   
    
    # Convert the date column to a datetime object
    data_df['lev_dt'] = data_df['lev_dt'].apply(pd.to_datetime)
    data_df = data_df.rename(columns={'station_nm':'SiteName'})
        
    print 'Found %i groundwater level observations.' %(len(data_df))
    print 'NOTE: All observations cached to .csv and must be filtered by subsequent operations.\n'
    
    return data_df

def get_discharge_results(data_url,sep=None):
    '''HELPER FUNCTION. Returns a dataframe of discharge observations.'''
    
    data_df = df_from_url(data_url,sep)
    
    print data_df.columns
        
    # Remove any dataframe rows that were read from blank rows in the file
    # (or, in some cases, are comment lines and headers for additional sites)
    data_df['site_no'] = data_df['site_no'].convert_objects(convert_numeric=True)    
    data_df = data_df.dropna(subset=['site_no'])
    
    # Convert the date column to a datetime object    
    data_df['datetime'] = data_df['datetime'].apply(pd.to_datetime)
    data_df = data_df.rename(columns={'02_00060_00003':'Discharge (ft3/s)'})
        
    print 'Found %i discharge observations.' %(len(data_df))
    print 'NOTE: All observations cached to .csv and must be filtered by subsequent operations.\n'
    
    return data_df
  
def get_wq_results(data_url,sep=','):
    '''HELPER FUNCTION. Returns a dataframe of water quality observations.'''
    
    data_df = df_from_url(data_url,sep,is_wq=True)
    
    data_df = data_df[['MonitoringLocationIdentifier','ActivityStartDate','ResultMeasureValue']]
    data_df.columns = ['site_no','Date','Value']
    data_df['site_no'] = data_df['site_no'].apply(lambda x: x.replace('USGS-',''))
    data_df['Date'] = data_df['Date'].apply(pd.to_datetime)
    
    # Remove any dataframe rows that were read from blank rows in the file
    # (or, in some cases, are comment lines and headers for additional sites)
    data_df['site_no'] = data_df['site_no'].convert_objects(convert_numeric=True)    
    data_df = data_df.dropna(subset=['site_no'])
    
    return data_df

# =============================================================================
# The following functions called by model preparation scripts.
# =============================================================================

def get_gw(domain_shp,head_data_csv,head_station_csv,nwis_gw_format='rdb',model_epsg=None,sep='\t'):
    '''Retrieves NWIS groundwater levels and the associated station information
    and writes them to .csv.'''

    bbox = get_bbox(domain_shp)   
    bbox = gu.bbox_to_latlon(bbox,model_epsg)
    bbox_string = get_bbox_string(bbox)
        
    # Configure the NWIS request
    gw_data_root = nwis_root + 'gwlevels/?'
    gw_site_root = nwis_root + 'site/?'
    param_codes = '72019' # Other possible GW Level param codes = 62610,62611,72019,72150,72020
    gw_site_type = 'GW'
    period = 'P30000D'  # NOTE: Default period for webservices API = 'most recent'
          
    # For convenience, reduce the number of columns in the station dataframe
    keep_station_cols = ['site_no', 'station_nm','dec_lat_va','dec_long_va','alt_va','well_depth_va','hole_depth_va','site_tp_cd']
        
    # Configure the urls
    site_url = gw_site_root + 'format=' + nwis_gw_format + '&bBox=' + bbox_string + '&siteType=' + gw_site_type + '&siteOutput=expanded'
    data_url = gw_data_root + 'format=' + nwis_gw_format + '&bBox=' + bbox_string + \
                    '&parameterCd=' + param_codes + '&siteType=' + gw_site_type + '&period=' + period

    station_df = get_station_df(site_url,gw_site_type,keep_station_cols,sep=sep)
    data_df = get_head_results(data_url,sep=sep)

    station_df.to_csv(head_station_csv)
    data_df.to_csv(head_data_csv)
        
    return

def get_discharge(domain_shp,discharge_data_csv,discharge_station_csv,nwis_discharge_format='rdb',model_epsg=None,sep='\t'):
    '''Retrieves NWIS discharge and the associated station information and writes them to .csv.'''
    
    bbox = get_bbox(domain_shp)   
    bbox = gu.bbox_to_latlon(bbox,model_epsg)
    bbox_string = get_bbox_string(bbox)
            
    # Configure the NWIS request
    dv_data_root = nwis_root + 'dv/?'
    dv_site_root = nwis_root + 'site/?'
    param_codes = '00060'   # Discharge in ft3/s   
    discharge_site_type = 'ST'
    period = 'P30000D'  # NOTE: Default period for webservices API = 'most recent'
    
    # For convenience, reduce the number of columns in the station dataframe
    keep_station_cols = ['site_no', 'station_nm','dec_lat_va','dec_long_va','alt_va','site_tp_cd']
        
    # Configure the urls
    site_url = dv_site_root + 'format=' + nwis_discharge_format + '&bBox=' + bbox_string + '&siteType=' + discharge_site_type + '&siteOutput=expanded'
    data_url = dv_data_root + 'format=' + nwis_discharge_format + '&bBox=' + bbox_string + \
                    '&parameterCd=' + param_codes + '&siteType=' + discharge_site_type + '&period=' + period
    
    station_df = get_station_df(site_url,discharge_site_type,keep_station_cols,sep=sep)
    data_df = get_discharge_results(data_url,sep=sep)

    station_df.to_csv(discharge_station_csv)
    data_df.to_csv(discharge_data_csv)                   
                           
    return

def get_wq(domain_shp,data_csv,station_df,nwis_format='csv',nwis_param=None,model_epsg=None):
    '''Downloads NWIS (NOT STORET - need separate call) water quality data to .csv.
    Note that current workflow only includes sites for which groundwater level or discharge
    information exists; as a result, this function does not retrieve station
    information but instead uses the station_df argument to filter observations.'''
    
    bbox = get_bbox(domain_shp)   
    bbox = gu.bbox_to_latlon(bbox,model_epsg)
    bbox_string = get_bbox_string(bbox)
   
    # Water quality services uses different format string for bbox call
    bbox_string = bbox_string.replace(',','%2C')
            
    # Configure the web services request
    dv_data_root = 'http://waterqualitydata.us/Result/search?'
    data_url = dv_data_root + 'mimeType=' + nwis_format + '&bBox=' + bbox_string
    
    if nwis_param is not None:
        data_url = data_url + '&pCode=' + nwis_param  
    
    wq_df = get_wq_results(data_url,sep=',')

    # Water quality webservices does not allow filtering by site type.
    # Consequently need to filter using user-specified station information.
    wq_df = pd.merge(station_df,wq_df,on='site_no',how='inner')
    wq_df = wq_df.dropna(subset=['Value'])
    
    wq_df.to_csv(data_csv)
            
    return wq_df

# -------------------------------------------
# --- DEPRECATED OR UNDEVELOPED FUNCTIONS ---
# -------------------------------------------

def get_json_df(iurl):
    '''Reads a USGS waterservices json response to a dataframe.
    NOTE: NOT FULLY DEVELOPED. Current version partially reads the response
    from gw services/levels. Modification required to read site information.'''

    ijson = urllib2.urlopen(iurl).read()
    ijson = json.loads(ijson)

    data_df = pd.DataFrame(columns=['SiteID','Date','Lat','Lon','Srs','Value','Param'])
    
    icount = 0
    for record in ijson['value']['timeSeries']:
        
        iloc = record['sourceInfo']
        iid  = iloc['siteCode'][0]['value']
        igeo = iloc['geoLocation']['geogLocation']
        isrs,ilat,ilon = igeo['srs'],igeo['latitude'],igeo['longitude']
    
        ival = record['values'][0]['value'][0]['value']
        idate = record['values'][0]['value'][0]['dateTime'].split('T')[0]
        iparam = record['variable']['variableCode'][0]['value']
        
        data_df.loc[icount,:] = [iid,idate,ilat,ilon,isrs,ival,iparam]
        icount += 1
    
    return data_df

def count_comment_lines(ifile):
    '''Counts the number of lines to be skipped when dataframe is created. Necessary
    due to apparent pandas malfunction of comment kwarg. NOTE: Not currently
    used due to improved functionality of the pandas 'comment' kwarg.'''

    skip_count = 0
    with open(ifile,'r') as fin:
        
        for iline in fin:
            
            if (iline.startswith('#') == True):
                skip_count += 1
                
            else:
                header = iline.split('\t')  # Capture this for use as header (problems with header kwarg)
                skip_count += 1
                break
        
    fin.close()        
    return skip_count,header

def separate_params(df,param_codes):
    '''Separates parameter codes into separate columns.  This is required in order
    to plot the supporting parameter information'''

    df_list = []    # This is a list of dataframes that will be merged
    
    # Generate a separate dataframe for each user-specified parameter
    for icode in param_codes:
        
        i_df = df[df[param_col_header] == icode]
        i_df[icode] = i_df['result_va']
        i_df = i_df[['site_no','TobDate',icode]]
        df_list.append(i_df)
    
    # Merge the separate dataframes    
    j_df = df_list[0]
    
    if (len(df_list) > 1):
        
        for i_df in df_list[1:]:
            
            j_df = pd.merge(j_df,i_df,on=['site_no','TobDate'],how='outer')
        
    return j_df

def create_param_dict(data_file,param_list_flag,param_keywords):
    '''Creates a dictionary with key=param_id and value=param_name. NOTE: this
    is a filtered dictionary that only includes parameters with descriptions
    that include user-specified keywords and phrases.  Useful for determining
    which parameters to include in the param_code_list that is used to filter
    the data dataframe. THIS NEEDS TO BE MODIFIED IN ORDER TO OPERATE OFF
    OF THE CSV/DATAFRAME RATHER THAN THE TEXT FILE.'''

    param_dict = dict()  
    with open(data_file,'r') as fin:
        
        # Read to the beginning of the list of parameter ids
        for iline in fin:
            if param_list_flag in iline:
                break
        
        # Continue reading until the next blank line in the file
        for iline in fin:
            
            if (len(iline.split()) == 1):   # blank line still includes '#'
                break
    
            for iword in param_keywords:
                
                if iword in iline:
                    
                    iparam = iline.split('-')[1].strip()
                    iparam_id = int((iline.split('-')[0]).strip('#').strip())
                    param_dict[iparam_id] = iparam
            
        fin.close()
    return param_dict