# -*- coding: utf-8 -*-
"""
Tools for generating MODFLOW framework arrays (i.e., IBOUND, modeltop elevations).

Created on Tue Jan 12 14:05:38 2016
@author: wzell
"""

import os,fiona,subprocess
from shapely.geometry import mapping,shape
import rasterio
from shapely.ops import unary_union
import numpy as np

# --- START MODULE PARAMETER SET ---

# Configure the cache for gdal subprocess calls
default_cachemax = 4000
default_cache_config = cache_config = ['--config','GDAL_CACHEMAX',str(default_cachemax)]

# --- STOP MODULE PARAMETER SET ----

def del_dir_contents(dir_to_delete,files_to_keep=None):
    '''Delete the files in the specified directory.'''
    
    for ifile in os.listdir(dir_to_delete):
        
        if files_to_keep is not None:
            if files_to_keep in ifile:
                continue
        
        os.remove(os.path.join(dir_to_delete,ifile))
    
    return

# === START IBOUND DEFINITION FUNCTIONS ====

def get_extents_from_huc(huc_data_shp=None,extents_output_shp=None,extents_huc_list=None):
    '''Extracts a user-specified HUC or list of HUCs from the national dataset and writes it
    to a shapefile. 'huc_data_shp'=shapefile that includes the huc polygons
    that will be extracted.'''
       
    extents_huc_scale = len(extents_huc_list[0])
    huc_field = 'HUC' + str(extents_huc_scale)
    
    with fiona.open(huc_data_shp) as vin:
        schema = vin.schema
        crs = vin.crs
        driver = vin.driver        
    
    # Reduce the extract schema to only the huc id field    
    schema['properties'] = {huc_field:'str'}
    
    # Now write the model domain shapefile     
    with fiona.open(huc_data_shp) as vect_in:
        polygon_list = []        
        for feature in vect_in:
            if (feature['properties'][huc_field] in extents_huc_list): 
                polygon_list.append(shape(feature['geometry']))
                merged = unary_union(polygon_list)

    with fiona.open(extents_output_shp,'w',driver=driver,crs=crs,schema=schema) as extract_out:
        extract_out.write({'geometry': mapping(merged),'properties':{huc_field:'Merged'}})
    
    return
    
def get_ibound_from_huc(huc_data_shp=None,ibound_output_shp=None,ibound_huc_scale=None,extents_huc_list=None):
    '''Reduces the shapefile used for the ibound zonation scheme to only those
    polygons (i.e., HUCs) that are in the active model area.  For higher resolution
    zonation (i.e., using smaller HUCs for the IBOUND zonation scheme) this is much
    faster than clipping the zonation shapefile before rasterizing.'''

    print '\nReducing the IBOUND zonation shapefile to the active model area.\n' 
    
    zone_id_field = 'HUC' + str(ibound_huc_scale)
    
    with fiona.open(huc_data_shp,'r') as vin:
        
        driver,crs,schema = vin.driver,vin.crs,vin.schema
        schema['properties'] = {zone_id_field:'str'}
        
        with fiona.open(ibound_output_shp,'w',driver=driver,crs=crs,schema=schema) as vout:
        
            for feature in vin:
                izone = int(feature['properties']['HUC' + str(ibound_huc_scale)])
                check_izone = str(izone).zfill(int(ibound_huc_scale))
                
                for ihuc in extents_huc_list:
                    if (check_izone.startswith(ihuc)):
                        
                        igeometry = shape(feature['geometry'])
                        vout.write({'geometry': mapping(igeometry),'properties':{zone_id_field:izone}})

    return

def check_existing_ibound(ibound_raster):
    '''Checks for an existing copy of the IBOUND raster.'''
    
    if (os.path.isfile(ibound_raster) == True):
        del_raster = None
        while (del_raster not in ['Y','N','y','n']):
            del_raster = raw_input('Do you want to overwrite the existing IBOUND raster? Y/N ')
        if (del_raster in ['Y','y']):
            print '\nAttempting to delete the raster.'
            print 'Note that you may need to release/delete the file in ArcMap/Catalog.\n'
            os.remove(ibound_raster)
            return
        else:
            print 'Stopping.'
            quit() 

# === STOP IBOUND DEFINITION FUNCTIONS ===

# === START VECTOR DATA FUNCTIONS ========

def reproject_shp(shp_fin=None,shp_fout=None,epsg_fout=None):
    '''Reprojects a shapefile.'''
    
    reproject_cmd = ['ogr2ogr','-f','ESRI Shapefile',shp_fout,shp_fin,'-t_srs','EPSG:%s' %(epsg_fout)]
    subprocess.call(reproject_cmd)
    
    return
    
def shp_to_model(shp_fin,raster_fout,grid_dx=None,grid_dy=None,rasterize_field=None,\
                  no_data=0,cache_config=default_cache_config,burn_constant=None):
    '''Rasterizes a PROJECTED shapefile. Default: rasterize based upon user-specified shapefile
    field. If burn_constant is provided, burn that constant. If bounds are provided,
    clip the raster to those bounds.'''

    print '\nRasterizing shapefile: %s' %(shp_fin)
    print 'Writing output raster to: %s\n' %(raster_fout)

    layer_name = os.path.basename(shp_fin).replace('.shp','')
    rasterize_cmd = ['gdal_rasterize'] + cache_config

    if (burn_constant is not None):
        rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-burn',str(burn_constant),'-of','GTiff','-ot','Float64','-l',layer_name,'-tr',str(grid_dx),str(grid_dy),shp_fin,raster_fout]
    else:            
        rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-a',rasterize_field,'-of','GTiff','-ot','Float64','-l',layer_name,'-tr',str(grid_dx),str(grid_dy),shp_fin,raster_fout]
    
    subprocess.call(rasterize_cmd)

    return

# === STOP VECTOR DATA FUNCTIONS =========

# === START RASTER DATA FUNCTIONS ========   

def get_raster_extents(raster_fin):
    
    with rasterio.open(raster_fin,'r') as src:
        
        # calculate extent of raster
        x_min = src.transform[0]
        x_max = src.transform[0] + src.transform[1]*src.width
        y_min = src.transform[3] + src.transform[5]*src.height
        y_max = src.transform[3]
        
    return x_min,y_min,x_max,y_max
    
def raster_to_model(data_fin,clipped_temp,raster_fout,bounds=None,delr=None,delc=None,
                    model_epsg=None,resample='average',\
                    cache_config=default_cache_config,cachemax=default_cachemax):
    '''Clips and resamples a raster to the model grid.'''
    
    x_min,y_min,x_max,y_max = bounds
    
    # Clip the raster to the model domain bound
    print '\nClipping the raster to the model domain.\n'
    clip_dem_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-te',str(x_min),str(y_min),str(x_max),str(y_max),data_fin,clipped_temp]
    subprocess.call(clip_dem_cmd)
    
    # Resample the raster to the model resolution
    print '\nResampling the raster to model grid resolution.\n'
    resample_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-r',resample,'-tr',str(delc),str(delr),\
        clipped_temp,raster_fout]
    subprocess.call(resample_cmd)

    return

# === STOP RASTER DATA FUNCTIONS ======== 

def write_raster_array(raster_in,file_out,fmt=None,multiplier=1):
    '''Writes an individual raster to file.'''
    
    with rasterio.open(raster_in,'r') as r:        
        iarray = r.read()[0,:,:] * multiplier
            
    np.savetxt(file_out,iarray,fmt=fmt)
        
    return 

def read_raster(raster_in,band=0):
    '''Reads a raster to array.'''
    
    with rasterio.open(raster_in,'r') as r:
        a = r.read()[band,:,:]
        no_data = r.nodata
        
    return a,no_data
    
def get_coastal_ibound(ibound=None,dem=None,coast=None,coast_nodata=None):
    '''Maps the coastline raster to the ibound raster and returns an (nrow,ncol)
    ibound array with the coastline cells = -1 (i.e., constant head).'''

    start_heads = dem
                         
    # Reduce the coastline cells to only those cells adjacent to the active model area.
    # This may be required because the coastline raster returns coastline within the 
    # active model bounding box rather than only at the perimeter of the active shapefile.
    ibound[ibound == 0] = np.nan
    rows,cols = np.where(coast != coast_nodata)
    adj_coast = np.empty(np.shape(ibound))

    for irow,icol in zip(rows,cols):
        
        n,nw,w,sw,s,se,e,ne = (irow-1,icol),(irow-1,icol-1),(irow,icol-1),(irow+1,icol-1),\
                              (irow+1,icol),(irow+1,icol+1),(irow,icol+1),(irow-1,icol+1)
                             
        for icheck in [n,nw,w,sw,s,se,e,ne]:                    
            try:
                # If the coastal cell is adjacent to an active cell, keep it . . .
                if (ibound[icheck] > 0):
                    adj_coast[irow,icol] = True
                    break                                
            except:                        
                if (irow == 0) and (icheck in [n,nw,ne]): continue
                if (irow == np.shape(ibound)[0]-1) and (icheck in [s,sw,se]): continue   
                if (icol == 0) and (icheck in [w,nw,sw]): continue
                if (icol == np.shape(ibound)[1]-1) and (icheck in [e,ne,se]): continue
                adj_coast[irow,icol] = False
   
    # Coordinate the IBOUND constant head flags and the starting heads constant head values             
    ibound[adj_coast == True] = -1 
    start_heads[ibound == -1] = 0

    inactive_idx = np.isnan(ibound)
    ibound[inactive_idx] = 0
            
    return ibound,start_heads