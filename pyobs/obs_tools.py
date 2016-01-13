# -*- coding: utf-8 -*-
"""
Reads a raster description of a parameter field and writes a PEST parameter
template file.

Capitalization indicates PEST parameter handle (e.g., 'T_11').
Lower case indicates MODFLOW parameter name (e.g., 't_11')

Created on Tue Nov 10 08:57:21 2015
@author: wzell
"""

import rasterio,fiona,pyproj
from shapely.geometry import mapping,shape,Point
from fiona.crs import from_epsg
from pymodflow.pygrid import grid_util as gu
import numpy as np
import pandas as pd
from itertools import product
import string


from pandas.stats.api import ols

# --- START PARAMETER SET ----

min_dtw = 10                # Keep all wells shallower than this   
max_dtw = 300               # max depth to water (m) allowed for water table identification
max_trend = 0.05            # modeled trend (m/year) in head change from linear regression      
max_r2 = 0.10               # fit of modeled trend in head change from linear regression
std_mult = 0.10             # fraction of depth to water allowed for standard deviation of head observations
sat_thick_mult = 0.50       # fraction of (depth of well - depth to water) compared to depth of well
land_surface_buffer = 20    # sat thick (dow-dtw) condition only applies if well is deeper than this (m)   

http_dict = {'GW':['http://nwis.waterdata.usgs.gov/nwis/gwlevels/?site_no=','&agency_cd=USGS&amp'],\
             'ST':'http://waterdata.usgs.gov/usa/nwis/uv?'}

# --- Global Functions ---

def hob_file_to_df(hob_file,skiprows=3,header_start='#'):
    '''Reads a MODFLOW .hob file to a dataframe. This is useful for extracting
    the .hob observation write order and using it as a reference list for 
    writing additional simulated output and PEST instruction files.  Note that this
    assumes that there is a single header line that contains column identifiers.'''
    
    # First extract the column identifiers
    with open(hob_file,'r') as fin:        
        # Get the first line
        header_cols = fin.readlines()[0]
        # Remove the comment marker
        if header_cols.startswith(header_start):
            header_cols = header_cols[1:]
        # Make it a list
        header_cols = header_cols.split()
    
    # Then read the dataframe
    hob_df = pd.read_table(hob_file,skiprows=skiprows,names=header_cols,delim_whitespace=True)
    hob_df.columns = header_cols
    
    # Make the observation name the index
    hob_df = hob_df.set_index(header_cols[0])
    
    return hob_df

def pval_to_field(pval_dict,zones_array,constant_key=0):
    '''Maps values to an array given an array of flags and a dictionary
    with key=flag and value=value. For zones_array = NaNs, generates
    a constant array.'''

    print '\t',pval_dict

    if (np.all(np.isnan(zones_array))):
        ivalues = np.ones(np.shape(zones_array)) * pval_dict[constant_key]
    
    else:
        ivalues = np.zeros((np.shape(zones_array)))    
        for k,v in pval_dict.iteritems():
            ivalues[zones_array==k] = v

    return ivalues

def pval_to_dict(pval_file,split_zone_flags=True):
    '''Read a parameter value file and returns UP TO two dictionaries.
    (i) For all cases, returns a dictionary with key=param_name and value=param_value.
    (ii) If split_zone_flags=True, returns a dictionary with key=zone flag and
    value=param_value (if split_zone_flags=False, returns None for second dict).'''
    
    pval_df = pd.read_table(pval_file,delim_whitespace=True,header=None)    
    pval_df.columns = ['Parameter','Value']
    pval_name_dict  = dict(zip(pval_df['Parameter'],pval_df['Value']))
    pval_zone_dict  = None
    
    if (split_zone_flags == True):    
        pval_df['Flag'] = pval_df['Parameter'].apply(lambda x: int(x.split('_')[1]))    
        pval_zone_dict = dict(zip(pval_df['Flag'],pval_df['Value']))
    
    return pval_name_dict,pval_zone_dict

def trim_hob_name(iname):
    '''Reduces the observation name to 10 characters.'''
    
    if ('Localnumber,' in iname):
        iname = iname.split(',')[1]
        return iname[0:(np.min([10,len(iname)]))]
        
    else:
        return iname[0:(np.min([10,len(iname)]))]
    
def write_mf_hob(hob_df,hob_file,hob_param_dict):
    '''Writes the head observations from a dataframe to a MODFLOW .hob file.'''

    print_df = hob_df.copy()

    # Unpack the hob parameter dictionary and supplement the dataframe as necessary       
    MOBS,MAXM,IUHOBSV,HOBDRY,TOMULTH,IREFSP,TOFFSET = [hob_param_dict[k] for k in ['MOBS','MAXM','IUHOBSV','HOBDRY','TOMULTH','IREFSP','TOFFSET']]
    
    print_df['IREFSP'],print_df['TOFFSET'] = IREFSP,TOFFSET
    
    # Reduce the dataframe to only those columns that will be written to the hob file
    print_df = print_df[['station_nm','Layer','Row','Column','IREFSP','TOFFSET','ROFF','COFF','HeadObs','lev_va','Model_DTW','Easting','Northing','GlobZTop','GlobZBot','alt_va','ModelTop']]
    
    # Trim the observation name as necessary and
    # transform (lay,row,col) to MODFLOW indexing    
    print_df['station_nm'] = print_df['station_nm'].apply(trim_hob_name)
    for icol in ['Layer','Row','Column']:
        print_df[icol] = print_df[icol].apply(lambda x: x + 1)
    
    # HOB formatting
    hob_fout_formats = ['%-10s','%6i','%6i','%6i','%8.1f','%8.1f','%7.2f','%7.2f','%7.1f','%10.1f','%10.1f','%10.1f','%10.1f','%10.1f','%10.1f','%10.1f','%10.1f\n']        
           
    with open(hob_file,'w') as fout:
    
        fout.write('#%-9s%6s%6s%6s%8s%8s%7s%7s%7s%10s%10s%10s%10s%10s%10s%10s%10s\n' 
            % ('ObsName','Lay','Row','Col','IREFSP','TOFFSET','ROFF','COFF','Mean','NWIS_DTW','Model_DTW','Easting','Northing','GlobZtop','GlobZbot','GageElev','ModelElev'))
            
        fout.write('%i %i %i %i %i\n' %(len(hob_df),MOBS,MAXM,IUHOBSV,HOBDRY))
        fout.write('%i\n' %(TOMULTH))
        
        # Iterate through the observations (i.e., the .hob rows)
        for idx,irow in print_df.iterrows():
            
            # Iterate through the descriptors (i.e., the .hob columns)
            for jvalue,jformat in zip(irow,hob_fout_formats):   # NOTE: the hob_fout_formats are established at the beginning of this script                    
                fout.write(jformat %(jvalue))
                
    return

def write_ucode_hob(hob_df,ucode_obs_file,ucode_hob_instruction_file,ucode_param_dict):
    '''Writes the head observations from a dataframe to a UCODE .obs file and .instruction file.'''

    # Unpack the ucode parameter dictionary
    STATFLAG,GROUPNAME = [ucode_param_dict[k] for k in ['STATFLAG','GROUPNAME']]

    # Reduce the dataframe to only those columns that will be written to the UCODE file
    ucode_df = hob_df[['station_nm','HeadObs','StdError']]
    ucode_df['STATFLAG'] = STATFLAG
    ucode_df['GROUPNAME'] = GROUPNAME
    
    # Write the observation data file
    with file(ucode_obs_file,'w') as fout:
        
        fout.write('BEGIN OBSERVATION_DATA TABLE\n')
        fout.write('NROW=%-5iNCOL=%-5iCOLUMNLABELS\n' %(len(hob_df),5))
        fout.write('%-10s%10s%10s%10s%10s\n' \
            %('OBSNAME','OBSVALUE','STATISTIC','STATFLAG','GROUPNAME'))

        # Iterate through the observations (i.e., the .hob rows)
        for idx,irow in ucode_df.iterrows():
            
            # Iterate through the descriptors (i.e., the .hob columns)
            for jvalue,jformat in zip(irow,['%-10s','%10.2f','%10.2f','%10s','%10s']):                   
                fout.write(jformat %(jvalue))
            fout.write('\n')

        fout.write('END OBSERVATION_DATA')
        
    # Write the instruction file           
    with file(ucode_hob_instruction_file,'w') as fout:
        
        fout.write('jif @\n')
        fout.write('StandardFile 1 1 %5i\n' %(len(hob_df)))
        
        for idx,irow in ucode_df.iterrows():
            
            fout.write('%-10s\n' %(irow['station_nm']))    

    return

def write_ucode_tob(tob_df,ucode_tob_file,ucode_tob_instruct_file,isol=None):
    '''Writes the transport observations from a dateframe to a UCODE .obs file and .instruction file.'''

    nconcobs = tob_df.count()

    # Open the output files and write the headers
    with open(ucode_tob_file,'w') as ucode_tob_fout:  # This is the UCODE observation file
        
        ucode_tob_fout.write('BEGIN OBSERVATION_DATA TABLE\n')
        ucode_tob_fout.write('NROW=%i NCOL=5 COLUMNLABELS\n' %(nconcobs))
        ucode_tob_fout.write('%-20s%-15s%-15s%-15s%-15s%-15s\n' %('OBSNAME','OBSVALUE','STATISTIC','STATFLAG','GROUPNAME','DATE'))            
                    
        with open(ucode_tob_instruct_file,'w') as ucode_instruct_fout:    # This is the UCODE instruction file
    
            ucode_instruct_fout.write('jif @\n')
            ucode_instruct_fout.write('StandardFile 2 2 %10i\n' %(nconcobs))
    
            # ITERATE THROUGH OBSERVATIONS FOR THIS SOLUTE.
            # icount used to distinguish multiple observations of same solute at single site
            icount,oldsite = 0,0
            for idx,irow in tob_df.iterrows():
                
                isite = irow['SiteName']
                
                # Name and number the observation        
                if (isite != oldsite):
                    icount = 0
                itob = isite + '_' + isol + '_' + str(icount)
                oldsite = isite        
                icount += 1
                
                ival = irow['Value']
                istd = irow['Std']
    
                ucode_tob_fout.write('%-20s%-15.4e%-15.4f%-15s%-15s%-15s\n' %(itob,ival,istd,'SD','age',irow.Date.strftime('%Y-%m-%d')))    
                ucode_instruct_fout.write('%s\n' %(itob))
    
            ucode_tob_fout.write('END OBSERVATION_DATA')    
        
    return

def write_mt3_tob(model_name,tob_df,tob_file,isol=None,sim_start_date=None):
    '''Writes the transport observations from a dateframe to an MT3DMS .obs file and .instruction file.'''

    nconcobs = tob_df.count()

    # Assign the package parameters
    cscale,ioutcobs,iconclog,iconcintp = 1,1,0,0
    icomp,roff,coff,iweight = 1,0,0,1
    inConcObs,inFluxObs,inSaveObs = 51,0,0
    
    with open(tob_file,'w') as tob_fout:  # This is the MT3DMS .tob file
        
        tob_fout.write('%-5i%-5i%-80i' %(nconcobs,0,0))
        tob_fout.write('#Item1: MaxConcObs,MaxFluxObs,MaxFluxCells\n')
        
        tob_fout.write('%-15s%-5i%-5i%-65i' %(model_name,inConcObs,inFluxObs,inSaveObs))
        tob_fout.write('#Item2: OUTNAM,inConcObs,inFluxObs,inSaveObs\n')
        
        tob_fout.write('%-5i%-10.3f%-5i%-5i%-65i' %(nconcobs,cscale,ioutcobs,iconclog,iconcintp))
        tob_fout.write('#Item3: nConcObs,CScale,iOutCobs,iConcLOG,iConcINTP\n')

        # ITERATE THROUGH OBSERVATIONS FOR THIS SOLUTE.
        # icount used to distinguish multiple observations of same solute at single site
        icount,oldsite = 0,0
        for idx,irow in tob_df.iterrows():
            
            isite = irow['SiteName']
            idate_days = (irow['Date'] - sim_start_date).days # Gets the difference in days between the two datetime objects
            
            # Name and number the observation        
            if (isite != oldsite):
                icount = 0
            itob = isite + '_' + isol + '_' + str(icount)
            oldsite = isite        
            icount += 1
            
            ival = irow['Value']
            istd = irow['Std']
            iweight = 1./istd
                                        
            # Write the information from this row of the dataframe to the appropriate files
            tob_fout.write('%-20s%-5i%-5i%-5i%-5i%-15.5e%-5.1f%-5.1f%-5.2f%-20.5e' %(itob,irow.Lay,irow.Row,irow.Col,icomp,idate_days,roff,coff,iweight,ival))
            tob_fout.write('#Item4: COBSNAM,Lay,Row,Col,iComp,TimeObs,Roff,Coff,weight,COBS\n') 
    
    return

def process_filtered_df(idf,station_df,exclusion_dict):
    '''Processes the filtered (i.e., the keep and exclude) dataframes
    and re-associates them with the station information.'''

    # Reduce the observation to mean values plus associated statistics
    grouped = idf.groupby('site_no')
    mean_df = grouped.mean()
    sem_df  = grouped.sem().rename(columns={'lev_va':'StdError'})
    nobs_df = grouped.count().rename(columns={'lev_va':'NObs'})
    
    filtered_df = mean_df.merge(sem_df,left_index=True,right_index=True).merge(nobs_df,left_index=True,right_index=True)
    filtered_df = pd.merge(filtered_df,station_df.set_index('site_no'),how='left',left_index=True,right_index=True)
    
    exclusion_df = pd.DataFrame.from_dict(exclusion_dict,orient='index')
    exclusion_df.columns= ['Violation']
    filtered_df = pd.merge(filtered_df,exclusion_df,left_index=True,right_index=True,how='left')
    filtered_df['Violation'].fillna('No Violation',inplace=True)
    
    filtered_df['HeadObs'] = filtered_df['alt_va'] - filtered_df['lev_va']
    filtered_df['Model_DTW'] = filtered_df['ModelTop'] - filtered_df['HeadObs']
    
    return filtered_df

def filter_obs(hob_df,filter_water_table,filter_stresses):
    '''Applies depth and measurement variability criteria to the head observations.'''
   
    exclusion_dict = {}
    for iname,idf in hob_df.groupby('site_no'):
        
        idow  = idf['well_depth_va'].mean()
        idtw  = idf['lev_va'].mean()
        istd  = idf['lev_va'].std()
    
        if (filter_water_table == True): # Reduce the dataframe to only those sites that are likely to measure the water table
                
                # Keep all wells shallower than minimum depth
                if (idtw <= min_dtw):
                    continue
                
                # Criteria: Exclude depths to water that are deeper than the likely water table elevation
                if (idtw > max_dtw):
                    # print 'Too deep: ',iname,idtw
                    exclusion_dict[iname] = 'dtw = %i > max dtw' %(idtw)
                    continue
                
                # Criteria: Exclude wells that are likely measuring a confined aquifer
                if ((idow - idtw) > (sat_thick_mult * idow) and (idow > land_surface_buffer)):
                    # print 'Likely in a confined aquifer: ',iname,idtw
                    exclusion_dict[iname] = 'sat thick %i > %0.2f * dow' %((idow-idtw),sat_thick_mult)
                    continue
            
        if (filter_stresses == True):
            
                # Keep all wells shallower than minimum depth
                if (idtw <= min_dtw):
                    continue
                
                # Criteria: Exclude wells with a potential trend over time.
                # Perform ordinary least squares and test both the slope
                # and R^2 of the best fit
                idf['date_delta'] = (idf['lev_dt'] - idf['lev_dt'].min()) / np.timedelta64(1,'D')
                imodel = ols(y=idf['lev_va'],x=idf['date_delta'])                 
                itrend = imodel.beta['x'] * 365.  # The slope of the ols fit converted to length/year
                ir2 = imodel.r2                   # The R^2 value for the ols fit
                
                if ((itrend > max_trend) and (ir2 > max_r2)):
                    # print 'Apparent temporal trend: ',iname,itrend
                    exclusion_dict[iname] = 'Apparent temporal trend = %0.2f m/year > %0.2f & R^2 = %0.2f > %0.2f' %(itrend,max_trend,ir2,max_r2)
                    continue
                
                if  (istd > (std_mult * idtw)):
                    # print 'Excessive measurement variability: ',iname,istd
                    exclusion_dict[iname] = 'Measurement std = %4.2f > %0.2f * dtw' %(istd,std_mult)
                    continue
                
    return exclusion_dict

def hob_from_csv(station_df,data_csv,min_nobs=10,agg='mean',screen_length=1,convert2metric=True,filter_water_table=True,filter_stresses=True):
    '''Returns dataframe of merged station information and groundwater measurements
    that were previously cached from NWIS.'''
    
    data_df = pd.DataFrame.from_csv(data_csv)

    if (convert2metric == True):
        for icol in ['lev_va','sl_lev_va']:
            data_df[icol] *= 0.3048  # Convert feet to meters

    data_df['lev_dt'] = pd.to_datetime(data_df['lev_dt'])   # Convert the date string to a datetime object
    
    # Reduce the dataframe to those sites that meet the min_nobs criteria
    enough_observations = data_df['site_no'].value_counts()     # This is a series ...
    enough_observations = enough_observations[enough_observations >= min_nobs].index.values  # ... converted to a list of sites
    data_df = data_df[data_df['site_no'].isin(enough_observations)]
    
    # Merge the station information and the data to a single dataframe. The head observation
    # is computed by subtracting the depth to water (recorded with the measurement) 
    # from the gage altitude (recorded with the station information)
    hob_df = pd.merge(data_df[['site_no','lev_dt','lev_va']],station_df,on='site_no',how='left')    
    hob_df['HeadObs'] = hob_df['alt_va'] - hob_df['lev_va']
    hob_df['Model_DTW'] = hob_df['ModelTop'] - hob_df['HeadObs']

    # Drop sites with a negative groundwater elevation
    # and add a few columns for later convenience
    hob_df = hob_df[hob_df['HeadObs'] >= 0]
    hob_df = hob_df.sort_values(by='station_nm')    
    
    # Get a dictionary of observations that fail to meet depth 
    # and/or measurement variability criteria
    exclusion_dict = filter_obs(hob_df,filter_water_table,filter_stresses)
     
    keep_df = hob_df[~hob_df['site_no'].isin(exclusion_dict.keys())][['site_no','lev_va']]
    exclude_df = hob_df[hob_df['site_no'].isin(exclusion_dict.keys())][['site_no','lev_va']]
    
    keep_df = process_filtered_df(keep_df,station_df,exclusion_dict)
    exclude_df = process_filtered_df(exclude_df,station_df,exclusion_dict)

    return keep_df,exclude_df

def write_hob_shapefile(hob_df,shp_fout,model_epsg):
    '''Writes the groundwater observation locations to a 2D point shapefile.'''
    
    schema = {'geometry':'Point','properties':{'Site':'str','SiteName':'str',\
                                            'Head (m)':'float','DOW (m)':'float','NWIS_DTW (m)':'float',\
                                            'DOW-DTW (m)':'float','ModelLand (m)':'float',\
                                            'WellElev (m)':'float','NWISLand (m)':'float',\
                                            '(Lay,Row,Col)':'str','NWIS_Link':'str','Violation':'str'}}

    with fiona.collection(shp_fout, "w", "ESRI Shapefile",crs=from_epsg(model_epsg),schema=schema) as output:
        
        for index,irow in hob_df.iterrows():
            
            isite = str(index)
            ilink = http_dict['GW'][0] + str(isite) + http_dict['GW'][1]
            iname = irow['station_nm']
            ihob = irow['HeadObs']
            idow = irow['well_depth_va']
            idtw = irow['lev_va']
            iwellelev = irow['well_elev']
            inwis_land = irow['alt_va']
            imodel_land = irow['ModelTop']
            irowcol = '(%i,%i,%i)'%(irow['Layer'],irow['Row'],irow['Column'])
            iexclude = irow['Violation']
            
            point = Point(irow['Projected_X'],irow['Projected_Y'])
            output.write({'geometry': mapping(point),                            
                        'properties':{'Site':isite,'SiteName':iname,\
                                        'Head (m)':ihob,'DOW (m)':idow,'NWIS_DTW (m)':idtw,\
                                        'DOW-DTW (m)':(idow-idtw),'ModelLand (m)':imodel_land,\
                                        'WellElev (m)':iwellelev,'NWISLand (m)':inwis_land,\
                                        '(Lay,Row,Col)':irowcol,'NWIS_Link':ilink,'Violation':iexclude}})
            
    return

def write_discharge_shapefile(discharge_df,shp_fout,model_epsg,discharge_label):
    '''Writes the discharge measurement locations to a 2D point shapefile.'''
    
    schema = {'geometry':'Point','properties':{'Site':'int','SiteName':'str',\
                                            'GageElev':'float','NObs':'int',\
                                            discharge_label:'float','StdError':'float',\
                                            'NWIS_Link':'str'}}

    with fiona.collection(shp_fout, "w", "ESRI Shapefile",crs=from_epsg(model_epsg),schema=schema) as output:
        
        for index,irow in discharge_df.iterrows():
            
            isite = int(index)
            isite = str(isite).zfill(8) # NWIS discharge ids minimum 8 characters
            
            iname = irow['station_nm']
            igage = irow['alt_va']
            inobs = irow['NObs']
            imean = irow[discharge_label]
            istd  = irow['StdError']
            ilink = http_dict['ST'] + str(isite)
            
            point = Point(irow['Projected_X'],irow['Projected_Y'])
            output.write({'geometry': mapping(point),                            
                        'properties':{'Site':isite,'SiteName':iname,\
                                            'GageElev':igage,'NObs':inobs,\
                                            discharge_label:imean,'StdError':istd,\
                                            'NWIS_Link':ilink}})
            
    return

def get_pest_instruct_order(instruct_file):
    '''Reads the ordered list of observations from a PEST instruction file.'''
    
    obs_list = []
    with open(instruct_file,'r') as fin:
        
        all_lines = fin.readlines()
        
        # Extract the observation marker from the header line
        marker = all_lines[0].split()[1]        
        for iline in all_lines[1:]:            
            obs_list.append(iline.split('!')[1].strip())
        
    return obs_list

def write_pest_instruct(obs_list,obs_instruct_file,pif_flag='@'):
    '''Writes the PEST instruction file.'''

    with open(obs_instruct_file,'w') as fout:
        
        fout.write('pif ' + pif_flag + '\n')
        for iobs in obs_list:
            
            fout.write('l1 w !' + str(iobs) + '!\n')
    
    return

# -------------------------------
# Transport Observation Functions 
# -------------------------------
def tob_from_csv(station_df,data_csv):
    '''Returns dataframe of merged station information and groundwater measurements
    that were previously cached as NWIS-like (e.g., for the Upper Chester tracer
    information the individual tracer species were split into separate .csv files).'''

    data_df = pd.DataFrame.from_csv(data_csv)

    # Merge the station information and the data to a single dataframe.
    tob_df = pd.merge(data_df,station_df,on='station_nm',how='left')
    tob_df = tob_df.set_index('site_no')
    
    return tob_df

def enumerate_tob(isol_df,isol=None,obsname_col='ObsName'):
    '''Appends the observation name with additional information that identifies
    the solute and distinguishes between multiple observations of the same solute
    at a single site. Returns a list of updated observation labels.'''
    
    enumerated_obsnames = []
    icount,oldsite = 0,0
    for idx,irow in isol_df.iterrows():
                
        isite = irow[obsname_col]
        
        # Name and number the observation        
        if (isite != oldsite):
            icount = 0
        itob = isite + '_' + isol + '_' + str(icount)
        
        enumerated_obsnames.append(itob)        
        
        oldsite = isite        
        icount += 1
        
    return enumerated_obsnames

# -------------------------------
# Discharge Observation Functions 
# -------------------------------
def discharge_from_csv(station_df,data_csv,data_col='Discharge (ft3/s)',convert2metric=True,convertseconds2days=True,min_nobs=10,agg='mean',name_map=None):
    '''Returns a dataframe of merged station information and discharge measurements
    that were previously cached as NWIS-like .csv files.'''
    
    data_df = pd.DataFrame.from_csv(data_csv)

    if (convert2metric == True):
        data_df['Discharge (m3/s)'] = data_df[data_col] * 0.0283168 # Convert feet^3 to meters^3
        data_col = 'Discharge (m3/s)'
        
        if (convertseconds2days == True):
            data_df['Discharge (m3/day)'] = data_df[data_col] * 86400.
            data_col = 'Discharge (m3/day)'
    
    # Reduce the dataframe to those sites that meet the min_nobs criteria
    enough_observations = data_df['site_no'].value_counts()     # This is a series ...
    enough_observations = enough_observations[enough_observations >= min_nobs].index.values  # ... converted to a list of sites
    data_df = data_df[data_df['site_no'].isin(enough_observations)]
    
    # Reduce the observation to mean values plus associated statistics
    grouped = data_df[['site_no',data_col]].groupby('site_no')
    mean_df = grouped.mean()
    sem_df  = grouped.sem().rename(columns={data_col:'StdError'})
    nobs_df = grouped.count().rename(columns={data_col:'NObs'})
    
    aggregated_df = mean_df.merge(sem_df,left_index=True,right_index=True).merge(nobs_df,left_index=True,right_index=True)
    discharge_df = pd.merge(aggregated_df,station_df.set_index('site_no'),left_index=True,right_index=True)
    
    if (name_map is not None):
        discharge_df['site_no'] = discharge_df.index
        discharge_df['station_nm'] = discharge_df['site_no'].apply(lambda x: name_map[x])

    return discharge_df,data_col

class Obs(object):
    '''This class aggregates the grid geometry and the translates the observation
    locations to the grid.'''

    # ---------------------------------
    # The init generates the model grid
    # and reads the parameter zones from a text file (default) or raster.
    # ---------------------------------
    def __init__(self,modelname='',model_dir='',dx_dy=None,domain_shp='',model_epsg=0,ibound3D_file='',modeltop_file='',lay_thick=None):

        self.__name = modelname
        self.model_dir = model_dir
        self.delr,self.delc = dx_dy,dx_dy
        self.top = np.genfromtxt(modeltop_file)
        self.nrow,self.ncol = np.shape(self.top)
        self.ibound = gu.twoD_to_threeD(np.genfromtxt(ibound3D_file),(self.nrow,self.ncol))
        self.bottoms,self.midpoints = gu.create_grid_datums(self.top,lay_thick)
        self.lay_thick = lay_thick
        self.model_epsg = model_epsg
        
        # Read the bounding box extents from the model domain shapefile
        with fiona.open(domain_shp,'r') as extent:    
            self.bbox = extent.bounds
            
        return

    def station_to_modelgrid(self,station_csv,screen_length=1,bbox_is_latlon=False,convert2metric=True,is_discharge=False):
        '''Returns a dataframe of station information that includes projected model
        grid coordinates. If the model domain shape is already projected,
        need_to_project = False (default).'''
        
        print 'Reading groundwater station data and translating to the model grid.\n'        
        
        station_df = pd.DataFrame.from_csv(station_csv)
        
        if (convert2metric == True):
            convert_cols = ['alt_va']
            if (is_discharge == False): # For GW measurements . . .
                convert_cols = convert_cols + ['well_depth_va','hole_depth_va','well_elev']
            for icol in convert_cols:
                station_df[icol] *= 0.3048  # Convert feet to meters
        
        model_epsg = self.model_epsg        
        epsg_projected = pyproj.Proj('+init=EPSG:' + str(model_epsg))
        
        if (bbox_is_latlon == True):            
            w_lon,s_lat,e_lon,n_lat = self.bbox
            model_epsg = self.model_epsg

            # Configure the re-projection object and project the lat/lon coords            
            x_nw,y_nw = epsg_projected(w_lon,n_lat)
            x_sw,y_sw = epsg_projected(w_lon,s_lat)
            
        else:
            x_w,y_s,x_e,y_n = self.bbox
            x_nw,y_nw = x_w,y_n
            x_sw,y_sw = x_w,y_s
        
        # First return the projected distances from the NW corner of the active model (i.e., MODFLOW grid orientation.)
        station_df['Projected_X'],station_df['Projected_Y'],station_df['MODFLOW_X'],station_df['MODFLOW_Y'] = \
            zip(*station_df.apply(lambda row: gu.latlon_to_modelxy((row['dec_lat_va'],row['dec_long_va']),(x_nw,y_nw),epsg_projected),axis=1))
            
        # Now return the projected distances from the SW corner of the active model (i.e., easting and northing)        
        _,_,station_df['Easting'],station_df['Northing'] = \
            zip(*station_df.apply(lambda row: gu.latlon_to_modelxy((row['dec_lat_va'],row['dec_long_va']),(x_sw,y_sw),epsg_projected),axis=1))
            
        # The NWIS download occasionally captures stations that are outside the active model area.  grid_util.latlon_to_modelxy
        # returns these as model_xy = np.nan. Filter any such stations
        station_df = station_df.dropna(subset=['MODFLOW_X','MODFLOW_Y'])
        
        if (is_discharge == False): # For GW measurements . . . 

            # Translate the observation locations to the model grid
            station_df['Row'],station_df['Column'],station_df['ROFF'],station_df['COFF'] = \
                zip(*station_df.apply(lambda row: gu.modelxy_to_rowcol(row['MODFLOW_X'],row['MODFLOW_Y'],self.delr,self.delc),axis=1))
        
            station_df['ModelTop'],station_df['Layer'],station_df['LocalZ'] = \
                zip(*station_df.apply(gu.globalz_to_layer,args=(self.top,self.bottoms,self.midpoints,self.ibound,self.lay_thick),axis=1))
            
            # Remove the observations in inactive cells
            station_df = station_df.dropna(subset=['Layer'])
            
            # A few arithmetic transformations for later convenience . . .
            station_df['GlobZTop'],station_df['GlobZBot'] = station_df['well_elev'],station_df['well_elev'] - screen_length
        
        return station_df
        
# --------------------------
# PEST Specific Methods 
# --------------------------    

    def get_params(self,param_name_root,zones_array=None):
        '''Returns a dictionary of parameters with key=param name and value=zone.
        If zones_array=None, generates a spatially_constant parameter. If zones_array
        is provided, enumerates a separate parameter for each zone.'''
        
        if (zones_array is None):
            idict = {param_name_root + '0':'Constant'}
            
        else:       
            # Read the parameter zones        
            param_array = zones_array
    
            # Reduce the parameter array to only active model cells (i.e., IBOUND > 0)
            if (param_array.ndim == 2):  # For 2D arrays only check the top layer of the IBOUND
                param_array[self.ibound[0,:,:] == 0] = np.nan        
            else:        
                param_array[self.ibound == 0] = np.nan
            param_array[param_array == 0] = np.nan
            
            param_zones = list(np.unique(param_array[np.isfinite(param_array) == True]))
            param_zones = [int(i) for i in param_zones]
            
            idict = {}
            for i in param_zones:
                idict[param_name_root + str(i)] = i
        
        return idict   

    def write_pest_params(self,param_dict=None,tpl_file=None,pval_file=None,init_dict=None,ipar_init=1,fmt='%10.2f'):
        '''Writes a list of parameters to template (*.tpl) and parameter value
        (*.pval) files.'''
        
        print 'Writing template and parameter value files for the following parameters:'
        print param_dict.keys() 
    
        with open(tpl_file,'w') as tpl_out,open(pval_file,'w') as pval_out:
            tpl_out.write('ptf $\n')
            
            for ipar_name in param_dict:
                
                if (init_dict is not None):
                    ipar_init = init_dict['PARVAL1'][ipar_name]
                
                tpl_out.write('%-30s$%20s$\n' %(ipar_name,ipar_name))
                pval_out.write('%-30s %20f\n' %(ipar_name,ipar_init)) 
        
        return
                    
    def get_ppt(self,ppt_nx,ppt_ny,ppt_dz,zone_array,ppt_name_root='pt',param_name_root='k_',ipar_init=1):
        '''Generates a dataframe of pilot point locations. Note that grid_util.globalxy_to_localxy
        uses the northwest corner of the model grid as the origin.  However, the PEST 
        pilot points file assumes the SW corner of the model grid as the origin (i.e., 
        it prints coordinates as northings and eastings.)'''

        # Generate the pilot point grid in the row/col plane
        x_vector = np.linspace(0,self.delr * self.ncol,ppt_nx+1)[:-1]
        x_vector = np.add(x_vector,0.5*(x_vector[1]-x_vector[0]))
        
        y_vector = np.linspace(0,self.delc * self.nrow,ppt_ny+1)[:-1]
        y_vector = np.add(y_vector,0.5*(y_vector[1]-y_vector[0]))
        
        ppt_df = pd.DataFrame(columns=['PPT_Name','Param_Name','Easting','Northing','Elevation','Zone','Value'])
        icount = 0
        
        for ix,jy in product(x_vector,y_vector):
            
            # row and column returned from this function in PYTHON indexing
            irow,icol,_,_ = gu.modelxy_to_rowcol(ix,(self.nrow*self.delc - jy),self.delr,self.delc)
        
            # Generate the vertical dimension of the pilot point grid
            # This algorithm assumes that for any (row,col), no active cells are found
            # beneath any inactive cells.
            z_0 = self.midpoints[irow,icol][0]
            z_vector = np.arange(z_0,z_0 - np.sum(self.layer_depths),-(abs(ppt_dz)))
            
            # The PEST 10-character limit precludes
            # writing row,col, and layer identifiers into the pilot point name.
            # As an alternative, I've simply labeled them A,B,C, etc., starting
            # with the top pilot point in a given (row,col)
            for iz,iletter in zip(z_vector,string.uppercase[:(len(z_vector))]):
    
                ilay,_ = gu.globalz_to_localz(irow,icol,iz,self.midpoints,self.bottoms,self.layer_depths)
                if (gu.check_if_active(irow,icol,ilay,self.ibound) == True):
                    
                    izone = zone_array[irow-1,icol-1,ilay-1]
                    ppt_name = ppt_name_root + str(irow) + '_' + str(icol) + iletter
                    param_name = param_name_root + str(irow) + '_' + str(icol) + iletter
                    
                    ppt_df.loc[icount,:] = [ppt_name,param_name,ix,jy,iz,izone,ipar_init]
                    icount += 1
        
        return ppt_df
        
    def ppt_to_shapefile(self,ppt_df,shp_fout):
        '''Writes the (projected) pilot point locations to a 2D point shapefile.'''
        
        # Configure the re-projection object and project the lat/lon coords        
        w_lon,s_lat,e_lon,n_lat = self.bbox        
        epsg_projected = pyproj.Proj('+init=EPSG:' + str(self.model_epsg),preserve_units=True)
        
        # Determine the location of the projected coordinates origin (i.e.,
        # whether the model_x and model_y should be added to or subtracted from
        # the projected coordinates at the origin).        
        x_sw,y_sw = epsg_projected(w_lon,s_lat)
        x_se,y_se = epsg_projected(e_lon,s_lat)
        x_nw,y_nw = epsg_projected(w_lon,n_lat)
        
        if (x_sw < x_se and y_sw < y_nw):
            ppt_df['Global_X'] = ppt_df['Easting'] + x_sw
            ppt_df['Global_Y'] = ppt_df['Northing'] + y_sw
            
        else:
            print 'Need additional case for position of projection origin relative to model grid.'

        # Configure the shapefile
        schema = {'geometry':'Point','properties':{'PPT_Name':'str',\
                                                'Easting':'float','Northing':'float',\
                                                'Elevation':'float'}}
    
        with fiona.collection(shp_fout, "w", "ESRI Shapefile",crs=from_epsg(self.model_epsg),schema=schema) as output:
            
            for idx,row in ppt_df.iterrows():
                
                iname = row['PPT_Name']
                ieast,inorth = row['Easting'],row['Northing']
                ielev = row['Elevation']
                
                point = Point(row['Global_X'],row['Global_Y'])
                output.write({'properties':{'PPT_Name':iname,\
                                                'Easting':ieast,'Northing':inorth,\
                                                'Elevation':ielev},\
                                'geometry': mapping(point)})
                
        return        

    def write_ppt_params(self,ppt_df,tpl_file,ppt_file,param_name_root,ipar_init=1,add_params=None,ppt_replace_flag='$'):
        '''Write the active pilot points to file.'''
        
        with open(ppt_file,'w') as ppt_fout, open(tpl_file,'w') as tpl_fout:

            tpl_fout.write('ptf ' + ppt_replace_flag + '\n')
            
            for idx,irow in ppt_df.iterrows():
                
                ppt_name,param_name = irow['PPT_Name'],irow['Param_Name']
                ieast,inorth = irow['Easting'],irow['Northing']
                ielev = irow['Elevation']
                izone = irow['Zone']
                ival = ipar_init

                ppt_fout.write('%-12s%12.2f%12.2f%12.2f%5i%2s%20.5e%2s\n' %(ppt_name,ieast,inorth,ielev,izone,'',ival,''))
                tpl_fout.write('%-12s%12.2f%12.2f%12.2f%5i%2s%20s%2s\n' %(ppt_name,ieast,inorth,ielev,izone,ppt_replace_flag,param_name,ppt_replace_flag))
        
        return
		
# ----------------------
# UNDEVELOPED FUNCTIONS
# ----------------------

def separate_params(df):
    '''Separates parameter codes into separate columns.  This is required in order
    to plot the supporting parameter information. IMPORTED FROM CHESTER NITRATE 
	PROCESSING. Possibly redundant or defunct.'''

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
	
def write_tob_file(df):
    '''Writes the dataframe to a ucode observation file.'''

    nconcobs = len(df)

    # Open the output files and write the headers
    with open(tob_file + isuffix,'w') as tob_fout:  # This is the MT3DMS .tob file
        
        tob_fout.write('%-5i%-5i%-80i' %(nconcobs,0,0))
        tob_fout.write('#Item1: MaxConcObs,MaxFluxObs,MaxFluxCells\n')
        
        tob_fout.write('%-15s%-5i%-5i%-65i' %(model_root,inConcObs,inFluxObs,inSaveObs))
        tob_fout.write('#Item2: OUTNAM,inConcObs,inFluxObs,inSaveObs\n')
        
        tob_fout.write('%-5i%-10.3f%-5i%-5i%-65i' %(nconcobs,cscale,ioutcobs,iconclog,iconcintp))
        tob_fout.write('#Item3: nConcObs,CScale,iOutCobs,iConcLOG,iConcINTP\n')        
        
        with open(ucode_tob_file + isuffix,'w') as ucode_tob_fout:  # This is the UCODE observation file
            
            ucode_tob_fout.write('BEGIN OBSERVATION_DATA TABLE\n')
            ucode_tob_fout.write('NROW=%i NCOL=5 COLUMNLABELS\n' %(nconcobs))
            ucode_tob_fout.write('%-20s%-15s%-15s%-15s%-15s%-15s\n' %('OBSNAME','OBSVALUE','STATISTIC','STATFLAG','GROUPNAME','DATE'))            
                        
            with open(ucode_tob_instruct_file + isuffix,'w') as ucode_instruct_fout:    # This is the UCODE instruction file

                ucode_instruct_fout.write('jif @\n')
                ucode_instruct_fout.write('StandardFile 2 2 %10i\n' %(nconcobs))

                # ITERATE THROUGH OBSERVATIONS FOR THE USER SPECIFIED PARAM/SOLUTE CODE (e.g., nitrate).
                # icount used to distinguish multiple observations of same solute at single site
                
                for iname,igroup in df.groupby('SiteName'):                
                    icount,oldsite = 0,0
                    for idx,irow in igroup.iterrows():
                        
                        isite = irow.SiteName
                        idate_days = (irow.TobDate - sim_start_date).days # Gets the difference in days between the two datetime objects
                        
                        # Name and number the observation        
                        if (isite != oldsite):
                            icount = 0
                        itob = isite + '_' + isuffix + '_' + str(icount)
                        oldsite = isite        
                        icount += 1
                        
                        iconc_value = irow.result_va
                        istd_value = irow.lab_std_va
                        
                        if (np.isnan(istd_value) == True):
                            istd_value = 1
                                                    
                        # Write the information from this row of the dataframe to the appropriate files
                        tob_fout.write('%-20s%-5i%-5i%-5i%-5i%-15.5e%-5.1f%-5.1f%-5.2f%-20.5e' %(itob,irow.Lay,irow.Row,irow.Col,icomp,idate_days,roff,coff,iweight,iconc_value))
                        tob_fout.write('#Item4: COBSNAM,Lay,Row,Col,iComp,TimeObs,Roff,Coff,weight,COBS\n')                    
        
                        ucode_tob_fout.write('%-20s%-15.4e%-15.4f%-15s%-15s%-15s\n' %(itob,iconc_value,istd_value,'SD','nitrate',irow.TobDate.strftime('%Y-%m-%d')))
        
                        ucode_instruct_fout.write('%s\n' %(itob))
    
                ucode_tob_fout.write('END OBSERVATION_DATA')
    
    return