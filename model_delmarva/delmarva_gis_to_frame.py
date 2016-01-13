# -*- coding: utf-8 -*-
"""
Generates MODFLOW framework arrays (i.e., IBOUND, modeltop elevations) at
user-specified resolution for Delmarva models.

INPUT:
Required:
1. Active model boundary (domain shapefile)
2. Elevation data (georeferenced raster - e.g., DEM)
Optional:
3. Model geographic zone identifier (domain shapefile)
4. Model parameter zone identifier (shapefile - e.g., SSURGO shapefile) 

OUTPUT:
1. (Zoned) IBOUND array
2. Parameter zone array (e.g., recharge potential per soil drainage class)
3. (Unconfined) model top array (note that this may need further processing for
scenarios in which the simulated model top != land surface.)

Created on Wed Oct 14 09:38:20 2015
@author: wzell
"""

import numpy as np
import os,fiona,subprocess
from shapely.geometry import mapping,shape
import rasterio

from model_config import dx_dy,domain_shp,data_dir,frame_data_dir,model_epsg,\
                         landsurface_elev_file,landsurface_ibound_file,soils_file,landuse_file

# --- START PARAMETER SET ---

grid_dx,grid_dy = dx_dy,dx_dy
no_data = 0
dem_resample = 'bilinear'        # GDAL resampling options: 'average'|'mode'|'cubic'|etc. (http://www.gdal.org/gdalwarp.html)
landuse_resample = 'mode'

ibound_field = 'GRIDCODE'   # This is the field on the model domain shapefile that specifies the IBOUND zone

# Configure the cache for gdal subprocess calls
cachemax = 4000        
cache_config = ['--config','GDAL_CACHEMAX',str(cachemax)]

# -------------------
# Paths and filenames
# ------------------- 
dem_dir = os.path.join(data_dir,'LIDAR')
dem_data = os.path.join(dem_dir,'UC_LIDAR.tif')

soils_dir = os.path.join(data_dir,'SSURGO')
soils_data = os.path.join(soils_dir,'UC_SSURGO.shp')
project_soils = False
soils_param_field = 'drclasswet'
param_map = {'Somewhat excessively drained':1,'Well drained':1,'Moderately well drained':2,'Poorly drained':3,'Very poorly drained':4,None:9}

landuse_dir = 'C:\\Data_State\\MD\\MD_NLCD'
landuse_data = os.path.join(landuse_dir,'nlcd_md_utm18.tif')

# Temporary files (created and deleted)
scratch_dir = os.path.join(data_dir,'Scratch')

ibound_raster_temp = os.path.join(scratch_dir,'temp_ibound.tif')
dem_clipped = os.path.join(scratch_dir,'dem_clipped.tif')
soils_reproject = os.path.join(scratch_dir,'soils_reproject.shp')
soils_clipped = os.path.join(scratch_dir,'soils_clipped.shp')
soils_parameterized = os.path.join(scratch_dir,'soils_clipped_temp.shp')

# Output to model
ibound_raster = os.path.join(frame_data_dir,'IBOUND.tif')
dem_raster = os.path.join(frame_data_dir,'Elevation.tif')
soils_raster = os.path.join(frame_data_dir,'Soils.tif')
landuse_raster = os.path.join(frame_data_dir,'LandUse.tif')

# --- STOP PARAMETER SET ----

def del_scratch():
    '''Delete the temporary files in the scratch workspace.'''
    
    for ifile in os.listdir(scratch_dir):
        os.remove(os.path.join(scratch_dir,ifile))
    
    return



# ---- !!! SCRIPT STARTS HERE !!! ----

# -----------------------------------------------
# ---- SPECIFY MODEL GRID WITH IBOUND RASTER ----
# -----------------------------------------------

# Confirm that any existing IBOUND rasters are to be overwritten
check_existing_ibound(ibound_raster)

# Rasterize the PROJECTED IBOUND information.
# The bounds of this IBOUND becomes the clipping bounds for subsequent datasets
write_ibound_raster(domain_shp)
x_min,y_min,x_max,y_max = get_raster_extents(ibound_raster)

# -------------------------------------
# ---- RESAMPLE DATA TO MODEL GRID ----
# -------------------------------------

# ----------------------
# Land surface elevation
# ----------------------

# Clip the DEM to the model domain
print '\nClipping the DEM to the model domain.\n'
clip_dem_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
    ['-te',str(x_min),str(y_min),str(x_max),str(y_max),dem_data,dem_clipped]
subprocess.call(clip_dem_cmd)
dem_data = dem_clipped

# Resample the DEM to the model resolution
print '\nResampling the DEM to model grid resolution.\n'
resample_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
    ['-r',dem_resample,'-tr',str(grid_dx),str(grid_dy),dem_data,dem_raster]
subprocess.call(resample_cmd)

# ----------------------
# Soils
# ----------------------

if (project_soils == True):
    # Project the soils shapefile
    print '\nRe-projecting the surficial geology.\n'
    reproject_cmd = ['ogr2ogr','-f','ESRI Shapefile',soils_reproject,soils_data,'-t_srs','EPSG:%s' %(model_epsg)]
    subprocess.call(reproject_cmd)
    soils_data = soils_reproject

# Clip to the model domain
print '\nClipping the soils data to the model domain.\n'
clip_cmd = ['ogr2ogr','-f','ESRI Shapefile',soils_clipped,soils_data,'-clipsrc',str(x_min),str(y_min),str(x_max),str(y_max)]
subprocess.call(clip_cmd)

# Create a new shapefile with an attribute field that maps the soil parameter descriptor to an integer value.  This is 
# required in the event that the soil parameter field contains strings (which cannot be rasterized)
with fiona.open(soils_clipped,'r') as vin:
    src_crs = vin.crs
    src_driver = vin.driver
    src_schema = vin.schema
    
vout_schema = src_schema
vout_schema['properties'] = {soils_param_field:'str','ParamFlag':'int'}

with fiona.open(soils_clipped,'r') as vin, \
     fiona.open(soils_parameterized,'w',driver=src_driver,crs=src_crs,schema=vout_schema) as vout:    
    
    param_list = []
    for feature in vin:
        
        iparam = feature['properties'][soils_param_field]
        if iparam not in param_list:
            param_list.append(iparam)
            print iparam

        vout.write({'geometry':mapping(shape(feature['geometry'])),'properties':{soils_param_field:iparam,'ParamFlag':param_map[iparam]}})
        
# Rasterize the soils
print '\nRasterizing the soils map.\n'
layer_name = os.path.basename(soils_parameterized).replace('.shp','')
rasterize_cmd = ['gdal_rasterize'] + cache_config
rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-a','ParamFlag','-l',layer_name,'-tr',str(grid_dx),str(grid_dy),soils_parameterized,soils_raster]
subprocess.call(rasterize_cmd)

# ----------------------
# Land Use
# ----------------------

# NOTE: In some cases may need to reproject this raster before clipping
# Clip the DEM to the model domain

print '\nClipping the land use raster to the model domain.\n'
clip_landuse_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
    ['-te',str(x_min),str(y_min),str(x_max),str(y_max),'-r',landuse_resample,'-tr',str(grid_dx),str(grid_dy),landuse_data,landuse_raster]
subprocess.call(clip_landuse_cmd)

# --------------------

# Write the rasters to text arrays
for fin,fout,ifmt in zip([dem_raster,ibound_raster,soils_raster,landuse_raster],\
                         [landsurface_elev_file,landsurface_ibound_file,soils_file,landuse_file],\
                         ['%10.2f','%5i','%5i','%5i']):
    
    with rasterio.open(fin,'r') as src:        
        iarray = src.read()[0,:,:]
        np.savetxt(fout,iarray,fmt=ifmt)

# Clean up the temporary files
del_scratch()