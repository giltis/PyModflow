# -*- coding: utf-8 -*-
"""
Generates MODFLOW framework arrays (i.e., IBOUND, modeltop elevations) at
user-specified resolution.

INPUT:
Required:
1. Active model boundary (shapefile - e.g., national boundary or HUC4 boundary)
2. Elevation data (georeferenced raster - e.g., DEM)
Optional:
3. Model geographic zone identifier (shapefile - e.g., HUC12 boundaries)
4. Model parameter zone identifier (shapefile - e.g., surficial geology) 

OUTPUT:
1. (Zoned) IBOUND array
2. Parameter zone array (e.g., HK zoned per surficial geology)
3. (Unconfined) model top array (note that this may need further processing for
scenarios in which the simulated model top != land surface.)
Optional:
4. Arrays required for FloPy package construction.

NB: The current version has some flexibility for specifying the model domain
shapefile for purpose of generating the IBOUND raster. However, the function
call that generates the land surface elevation ASSUMES THAT THE LARGEST ACTIVE
MODEL AREA = HUC2; based upon this assumption it uses the elevation data
provided by HUC2 with the NHDPlus dataset.  This assumption is hard-coded into
the conditional 'if (model_huc_scale > 2)'.

Created on Wed Oct 14 09:38:20 2015
@author: wzell
"""

import os,fiona,subprocess
from shapely.geometry import mapping,shape
import rasterio
from shapely.ops import unary_union
import numpy as np
import matplotlib.pyplot as plt

import scipy.ndimage as ndimage          
from scipy import stats

# --- START PARAMETER SET ---

# Configure the cache for gdal subprocess calls
cachemax = 4000        
cache_config = ['--config','GDAL_CACHEMAX',str(cachemax)]

# --- STOP PARAMETER SET ----

def reduce_zone_shp(zone_base_shp,ibound_shp,zone_huc_scale,model_huc_list):
    '''Reduces the shapefile used for the ibound zonation scheme to only those
    polygons (i.e., HUCs) that are in the active model area.  For higher resolution
    zonation (i.e., using smaller HUCs for the IBOUND zonation scheme) this is much
    faster than clipping the zonation shapefile before rasterizing.'''

    print '\nReducing the IBOUND zonation shapefile to the active model area.\n' 
    
    zone_id_field = 'HUC' + str(zone_huc_scale).zfill(2)
    
    with fiona.open(zone_base_shp,'r') as vin:
        
        driver,crs,schema = vin.driver,vin.crs,vin.schema
        schema['properties'] = {zone_id_field:'str'}
        
        with fiona.open(ibound_shp,'w',driver=driver,crs=crs,schema=schema) as vout:
        
            for feature in vin:
                izone = int(feature['properties']['HUC' + str(zone_huc_scale)])
                check_izone = str(izone).zfill(int(zone_huc_scale))
                
                for ihuc in model_huc_list:
                    if (check_izone.startswith(ihuc)):
                        
                        igeometry = shape(feature['geometry'])
                        vout.write({'geometry': mapping(igeometry),'properties':{zone_id_field:izone}})

    return

def write_active_envelope(domain_base_shp=None,domain_extract_shp=None,huc_field=None,model_huc_list=None):
    '''Extracts a user-specified HUC or list of HUCs from the national dataset and writes it
    to a shapefile.'''

    # Check if the specified model_huc is a single huc or a list.
    # If a single huc, put into a list    
    try:
        assert isinstance(model_huc_list,list)
    except:
        model_huc_list = [model_huc_list]
    
    with fiona.open(domain_base_shp) as vin:
        schema = vin.schema
        crs = vin.crs
        driver = vin.driver        
    
    # Reduce the extract schema to only the huc id field    
    schema['properties'] = {huc_field:'str'}
    
    # Now write the model domain shapefile     
    with fiona.open(domain_base_shp) as vect_in:
        polygon_list = []        
        for feature in vect_in:
            if (feature['properties'][huc_field] in model_huc_list): 
                polygon_list.append(shape(feature['geometry']))
                merged = unary_union(polygon_list)

    with fiona.open(domain_extract_shp,'w',driver=driver,crs=crs,schema=schema) as extract_out:
        extract_out.write({'geometry': mapping(merged),'properties':{huc_field:'Merged'}})
    
    return

def del_scratch(scratch_dir):
    '''Delete the temporary files in the scratch workspace.'''
    
    for ifile in os.listdir(scratch_dir):
        os.remove(os.path.join(scratch_dir,ifile))
    
    return

def write_ibound_raster(unprojected_shp,projected_shp,projected_raster,delr,delc,src_epsg=4269,model_epsg=5070,discretize_ibound=True,no_data=0):
    '''Writes the IBOUND raster to disk. Default reprojection is from GCS NAD (epsg 4269) to Albers Equal Area (epsg 5070).
    For a single ibound flag set discretize_ibound = False.'''

    if (discretize_ibound == True):   # Discretize the IBOUND by HUC

        # Get the name of the field that identifies the HUC
        with fiona.open(unprojected_shp,'r') as vin:
            for feature in vin:
                huc_string = [x for x in feature['properties'].keys() if 'HUC' in x][0]
        
        # Project the HUC zones
        print '\nRe-projecting the IBOUND zonation scheme.\n'
        reproject_cmd = ['ogr2ogr',projected_shp,unprojected_shp,'-s_srs','EPSG:%s' %(src_epsg),'-t_srs','EPSG:%s' %(model_epsg)]
        subprocess.call(reproject_cmd)
        
        # Rasterize the projected HUC zones
        print '\nRasterizing the IBOUND zonation scheme.\n'
        layer_name = os.path.basename(projected_shp).replace('.shp','')
        rasterize_cmd = ['gdal_rasterize'] + cache_config
        rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-a',huc_string,'-l',layer_name,'-tr',str(delc),str(delr),projected_shp,projected_raster]
        subprocess.call(rasterize_cmd)

        
    else:   # Use a constant value = 1
    
        # Project the model domain
        print '\nRe-projecting the model domain.\n'
        reproject_cmd = ['ogr2ogr',projected_shp,unprojected_shp,'-s_srs','EPSG:%s' %(src_epsg),'-t_srs','EPSG:%s' %(model_epsg)]
        subprocess.call(reproject_cmd)    
    
        # Rasterize the projected model domain
        print '\nRasterizing the model domain.\n'    
        layer_name = os.path.basename(projected_shp).replace('.shp','')
        
        rasterize_cmd = ['gdal_rasterize'] + cache_config
        rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-burn',str(1),'-l',layer_name,'-tr',str(delc),str(delr),projected_shp,projected_raster]
        subprocess.call(rasterize_cmd)
            
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

def get_raster_extents(raster_fin):
    
    with rasterio.open(raster_fin,'r') as src:
        
        # calculate extent of raster
        x_min = src.transform[0]
        x_max = src.transform[0] + src.transform[1]*src.width
        y_min = src.transform[3] + src.transform[5]*src.height
        y_max = src.transform[3]
        
    return x_min,y_min,x_max,y_max
    
def raster_to_model(data,clipped,raster,bounds,delr,delc,model_epsg=5070,resample='average'):
    '''Clips and resamples a raster to the model grid.'''
    
    x_min,y_min,x_max,y_max = bounds
    
    # Clip the DEM to the model domain bound
    print '\nClipping the raster to the model domain.\n'
    clip_dem_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-te',str(x_min),str(y_min),str(x_max),str(y_max),data,clipped]
    subprocess.call(clip_dem_cmd)
    
    # Resample the DEM to the model resolution
    print '\nResampling the DEM to model grid resolution.\n'
    resample_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-t_srs','EPSG:%s' %(model_epsg),'-r',resample,'-tr',str(delc),str(delr),clipped,raster]
    subprocess.call(resample_cmd)

    return
    
def shp_to_model(shp_data,data_field,shp_clipped,shp_raster,bounds,delr,delc,model_epsg=5070,no_data=0,resample='average',burn_value=1):
    '''Reprojects and rasterizes a shapefile to the model domain.'''
    
    x_min,y_min,x_max,y_max = bounds

    # Project the shapefile
    print '\nRe-projecting shapefile and clipping to the model domain: %s\n' %(shp_data)
    reproject_cmd = ['ogr2ogr','-f','ESRI Shapefile',shp_clipped,shp_data,'-t_srs','EPSG:%s' %(model_epsg),'-clipdst',str(x_min),str(y_min),str(x_max),str(y_max)]
    subprocess.call(reproject_cmd)
    
    print '\nRasterizing the shapefile: %s.\n' %(shp_data)
    layer_name = os.path.basename(shp_clipped).replace('.shp','')
    rasterize_cmd = ['gdal_rasterize'] + cache_config   
    
    if (data_field == 'burn'):        
        rasterize_cmd = rasterize_cmd + ['-burn',str(burn_value)]

    else:        
        rasterize_cmd = rasterize_cmd + ['-a',data_field]
        
    rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-l',layer_name,'-tr',str(delc),str(delr),'-te',str(x_min),str(y_min),str(x_max),str(y_max),shp_clipped,shp_raster]
    subprocess.call(rasterize_cmd)
        
    return

def nan_to_coastal(window):
    '''Finds gaps between the active model (IBOUND > 0) and the coastal
    constant head boundary and converts them to an active cell or a constant head.'''

    # Use the mode of the IBOUND zones in the search window to assign the IBOUND zone
    imode = stats.mode(window)[0][0][0]

    if (imode > 0):
        return imode
    else:
        return np.nan
    
def write_framework_arrays(ibound_raster,ibound3D_file,dem_raster,landsurface_file,modeltop_file,cell_bottoms3D_file,lay_thick,starting_heads_file,coastline_raster=None,no_data=-9999):
    '''Writes the IBOUND array to file. If the DEM raster is provided, then
    all DEM elevations <= zero are set as constant head (sea level = 0) cells
    in the IBOUND array.'''
    
    with rasterio.open(ibound_raster,'r') as i,\
         rasterio.open(dem_raster,'r') as d:
               
        ibound = i.read()[0,:,:]
        nrow,ncol = np.shape(ibound) 
        
        dem = d.read()[0,:,:]        
        dem = dem * 0.01    # Convert cm to meters 
        start_heads = dem
    
        if coastline_raster is not None:
            
            # Map the coastal cells to the ibound array            
            with rasterio.open(coastline_raster,'r') as c:

                coast = c.read()[0,:,:]
                         
            # Reduce the coastline cells to only those cells adjacent to the active model area.
            # This may be required because the coastline raster returns coastline within the 
            # active model bounding box rather than only at the perimeter of the active shapefile.
            ibound[ibound == 0] = np.nan
            rows,cols = np.where(coast != c.nodata)
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
    
    botms = dem - lay_thick[0]
    if (len(lay_thick) > 1):        
        for ilay in len(lay_thick):
            ibotm = dem - lay_thick[ilay]
            ibotm[inactive_idx] = no_data
            botms = np.stack((botms,ibotm))           

    np.savetxt(ibound3D_file,ibound,fmt='%15i')
    np.savetxt(landsurface_file,dem,fmt='%15.3f')
    np.savetxt(modeltop_file,dem,fmt='%15.3f')
    np.savetxt(cell_bottoms3D_file,botms,fmt='%15.3f')
    np.savetxt(starting_heads_file,start_heads,fmt='%15.3f')
        
    return
    
def write_raster_array(raster_in,file_out,fmt=None,multiplier=1):
    '''Writes an individual raster to file.'''
    
    with rasterio.open(raster_in,'r') as r:        
        iarray = r.read()[0,:,:] * multiplier
            
    np.savetxt(file_out,iarray,fmt=fmt)
        
    return 