# -*- coding: utf-8 -*-
"""
This script generates and writes the model-resolution rasters/ASCII files
of hydrogeologic unit elevations from the Maryland Coastal
Plain Aquifer Information (MDCPAI) shapefiles.

PREREQUISITES: See list of variables imported from model_config.
Note that this includes the PROJECTED model domain shapefile, which identifies
the boundaries of the active model area.

Created on Tue Nov 24 17:06:39 2015
@author: wzell
"""
import numpy as np
import os,subprocess,fiona,rasterio
from shapely.geometry import mapping,shape

from scipy import interpolate

from model_config import data_dir,delr,delc,model_epsg,domain_shp,Kzone_dict,frame_data_dir,strat_fin_prefix

# ---------------------------
# --- START PARAMETER SET ---

# Consider making these parameters accessible/arguments to a main function
convert_strat = False
process_units = [iunit for iunit in Kzone_dict.keys() if iunit not in ['Surficial','AquiaCU','Hrnrstwn']]

# Local features added to the MDCPAI for the Upper Chester;
# thicknesses in model length units (m)
hrnstwn_thick = 15
aquiaCU_thick = 7

feet_to_meters = 0.3048

# Paths and files
scratch_dir = os.path.join(data_dir,'Scratch')
model_raster_dir = os.path.join(frame_data_dir,'Hydrogeo_Rasters')
strat_dir = os.path.join(data_dir,'MDCPAI')
fips_dir = os.path.join(strat_dir,'FIPS')
utm_dir = os.path.join(strat_dir,'UTM18N_Meters')
regional_dem = os.path.join(utm_dir,'Extent_DEM_utm18N_meters.tif')

extent_sfx,contour_sfx = '_ext.shp','_con.shp'
elev_field = 'CONTOUR'
no_data = -9999

# Configure the cache for gdal subprocess calls
interpolate_string = 'invdist:power=2.0:smoothing=1.0'
cachemax = 4000        
cache_config = ['--config','GDAL_CACHEMAX',str(cachemax)]

# --- STOP PARAMETER SET ----
# ---------------------------

def del_dir_contents(dir_to_delete,files_to_keep=None):
    '''Delete the files in the specified directory.'''
    
    for ifile in os.listdir(dir_to_delete):
        
        if files_to_keep is not None:
            if files_to_keep in ifile:
                continue
        
        os.remove(os.path.join(dir_to_delete,ifile))
    
    return

def reproject_shp(src_shp,reproj_shp,model_epsg):
    '''Reprojects the contour shapefile to model projection.'''
    
    print '\nRe-projecting the shapefile: %s \n' %(src_shp)
    reproject_cmd = ['ogr2ogr','-f','ESRI Shapefile',reproj_shp,src_shp,'-t_srs','EPSG:%s' %(model_epsg)]
    subprocess.call(reproject_cmd)
    
    return

def contour_to_metric(shp,save_dir):
    '''Converts contours FIPS >> meters.'''
    
    save_shp = os.path.join(save_dir,os.path.basename(shp))
    
    with fiona.open(shp,'r') as src:
        crs,driver,schema = src.crs,src.driver,src.schema
    
    save_schema = schema
    save_schema['properties'] = {elev_field:'float'}
    
    with fiona.open(shp,'r') as vin,\
         fiona.open(save_shp,'w',driver=driver,crs=crs,schema=save_schema) as vout:
    
        for feature in vin:
            icontour = feature['properties'][elev_field] * feet_to_meters
            vout.write({'geometry':mapping(shape(feature['geometry'])),'properties':{elev_field:icontour}})
    
    return

def strat_to_utm(fips_dir,utm_dir,scratch_dir,contour_sfx,extent_sfx,model_epsg) :
    '''Reprojects the MDCPAI extent and contour files to the model epsg and converts
    the elevation field to metric.'''

    for ifile in os.listdir(fips_dir):
        
        # If an extents file, write directly to the utm_dir
        if ifile.endswith(extent_sfx):
            
            src_shp = os.path.join(fips_dir,ifile)
            reproj_shp = os.path.join(utm_dir,os.path.basename(ifile))
            reproject_shp(src_shp,reproj_shp,model_epsg)
    
        # If a contours file, write to scratch before updating the elevation field
        if ifile.endswith(contour_sfx):
    
            src_shp = os.path.join(fips_dir,ifile)
            reproj_shp = os.path.join(scratch_dir,os.path.basename(ifile))
            reproject_shp(src_shp,reproj_shp,model_epsg)
            
            contour_to_metric(reproj_shp,utm_dir)
        
    return

def get_output_paths(scratch_dir,iunit):
    
    # Generated rasters    
    contour_raster = os.path.join(scratch_dir,iunit + '_con.tif')
    dem_raster = os.path.join(scratch_dir,iunit + '.tif')
    envelope_raster = os.path.join(scratch_dir,iunit + '_envelope.tif')
    
    return contour_raster,dem_raster,envelope_raster

def get_input_paths(idir,iunit):
    
    # Input files
    extent_shp = os.path.join(idir,iunit + extent_sfx)
    contour_shp = os.path.join(idir,iunit + contour_sfx)
    contour_layer_name = os.path.basename(contour_shp).replace('.shp','')
    
    return extent_shp,contour_shp,contour_layer_name

def interpolate_contours(contour_shp,layer_name,dem,elev_field):
    '''Interpolates the contours to a DEM.'''

    with fiona.open(contour_shp,'r') as shp:
        
        x0,y0,x1,y1 = shp.bounds
        icol = int(np.round(np.abs(x1-x0))/delr)
        irow = int(np.round(np.abs(y1-y0))/delc)

    print '\nInterpolating the unit DEM: %s \n' %(contour_shp)
    grid_cmd = ['gdal_grid'] + cache_config + ['-txe',str(x0),str(x1),'-tye',str(y0),str(y1),'-outsize',str(icol),str(irow)]
    grid_cmd = grid_cmd + ['-zfield',elev_field,'-a',interpolate_string,'-of','GTiff','-ot','Float64','-l',layer_name,contour_shp,dem]
    subprocess.call(grid_cmd)
    
    return

def burn_contours(contour_shp,layer_name,contour_raster,elev_field):
    '''Rasterizes the contour lines.'''
    
    print '\nRasterizing the contour lines: %s \n' %(contour_shp)
    layer_name = os.path.basename(contour_shp).replace('.shp','')
    rasterize_cmd = ['gdal_rasterize'] + cache_config
    rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-a',elev_field,'-of','GTiff','-ot','Float64','-l',layer_name,'-tr',str(delr),str(delc),contour_shp,contour_raster]
    subprocess.call(rasterize_cmd)
    
    return

def split_extents(extent_shp):
    '''Extracts the outcrop and subcrop polygons from the extents shapefile.
    Returns PATHS to the extracted features and NAMES of the associated layers.'''
    
    outcrop_shp = os.path.basename(extent_shp.replace(extent_sfx,'_outcrop.shp'))    
    subcrop_shp = os.path.basename(extent_shp.replace(extent_sfx,'_subcrop.shp'))
    extents_shp = os.path.basename(extent_shp.replace(extent_sfx,'_extents.shp'))
    
    outcrop_shp = os.path.join(scratch_dir,outcrop_shp)
    subcrop_shp = os.path.join(scratch_dir,subcrop_shp)
    extents_shp = os.path.join(scratch_dir,extents_shp)
    
    outcrop_layer = os.path.basename(outcrop_shp).replace('.shp','')
    subcrop_layer = os.path.basename(subcrop_shp).replace('.shp','')
    extents_layer = os.path.basename(extents_shp).replace('.shp','')
       
    with fiona.open(extent_shp,'r') as vin:
        src_driver,src_crs,src_schema = vin.driver,vin.crs,vin.schema
    
    with fiona.open(extent_shp,'r') as vin,\
         fiona.open(outcrop_shp,'w',driver=src_driver,crs=src_crs,schema=src_schema) as outcrop,\
         fiona.open(subcrop_shp,'w',driver=src_driver,crs=src_crs,schema=src_schema) as subcrop,\
         fiona.open(extents_shp,'w',driver=src_driver,crs=src_crs,schema=src_schema) as extents:
        
        for feature in vin:
            
            itype = feature['properties']['FType']
            iwrite = {'geometry':mapping(shape(feature['geometry'])),'properties':feature['properties']}
            
            if (itype == 'Outcrop Display'):
                outcrop.write(iwrite)
                
            if (itype == 'Subcrop Display'):
                subcrop.write(iwrite)
                
            if ('Extent' in itype):
                extents.write(iwrite)
                
    return outcrop_shp,subcrop_shp,extents_shp,outcrop_layer,subcrop_layer,extents_layer

def raster_to_array(iunit,iarray,Kzone_dict):
    '''Writes the model raster to an ASCII array.'''
    
    iarray[np.isnan(iarray)] = no_data
    ifout = strat_fin_prefix + str(Kzone_dict[iunit])
    np.savetxt(ifout,iarray,fmt='%10.2f')
    
    print '\n\n\n ******\n'
    print iunit,np.shape(iarray)
    print '\n******\n\n\n'
    
    return

def clip_to_model(full_unit_raster):
    '''Clips a DEM for a hydrogeologic unit to the model domain.'''
    
    model_unit_raster = full_unit_raster.replace('.tif','_model.tif')
    print 'Clipping the unit raster to the model domain: %s\n' %(model_unit_raster)
    clip_unit_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-te',str(model_x_min),str(model_y_min),str(model_x_max),str(model_y_max),'-of','GTiff','-ot','Float64','-tr',str(delr),str(delc),'-dstnodata',str(no_data),full_unit_raster,model_unit_raster]
    subprocess.call(clip_unit_cmd)
    
    return model_unit_raster

# ------------------
# SCRIPT STARTS HERE
# ------------------

# Confirm that previous versions of the output rasters and any contents of the 
# scratch directory should be overwritten
delete_rasters = raw_input('Do you want to overwrite the model resolution hydrogeological rasters? Y/N ')
if (delete_rasters in ['Y','y']):
    del_dir_contents(model_raster_dir)
else:
    print 'Stopping.'
    quit()

delete_scratch = raw_input('Do you want to overwrite the scratch directory contents? Y/N ')
if (delete_scratch in ['Y','y']):
    del_dir_contents(scratch_dir)
else:
    print 'Stopping.'
    quit()

if (convert_strat == True):
    print 'Converting MDCPAI shapefiles to model coordinates.\n'
    strat_to_utm(fips_dir,utm_dir,scratch_dir,contour_sfx,extent_sfx,model_epsg)
    del_dir_contents(scratch_dir)

# Get the extents of the model domain
with fiona.open(domain_shp,'r') as vin:
    model_x_min,model_y_min,model_x_max,model_y_max = vin.bounds

# Map the bottom of the surficial aquifer directly from the surficial contours
envelope_shp,contour_shp,contour_layer_name = get_input_paths(utm_dir,'Surficial')
surficial_contours,surficial_dem,surficial_envelope = get_output_paths(model_raster_dir,'Surficial')
interpolate_contours(contour_shp,contour_layer_name,surficial_dem,elev_field)

# Clip the surficial to the model domain and save as an ASCII array
model_surficial_raster = clip_to_model(surficial_dem)
with rasterio.open(model_surficial_raster,'r') as src:        
    model_surficial_array = src.read()[0,:,:]
    raster_to_array('Surficial',model_surficial_array,Kzone_dict)

# Now iterate through the hydrogeologic units.
for iunit in process_units:
    
    print iunit

    # Get the paths and extract the outcrop and subcrop
    envelope_shp,contour_shp,contour_layer_name = get_input_paths(utm_dir,iunit)
    contour_raster,dem_raster,envelope_raster = get_output_paths(scratch_dir,iunit)    
    outcrop_shp,subcrop_shp,extent_shp,outcrop_layer,subcrop_layer,extent_layer = split_extents(envelope_shp)
    outcrop_raster,subcrop_raster = outcrop_shp.replace('.shp','.tif'),subcrop_shp.replace('.shp','.tif')

    # Get a rasterized version of the unit extents in order to later mask the numpy interpolation operation
    envelope_layer = os.path.basename(envelope_shp.replace('.shp',''))    
    rasterize_cmd = ['gdal_rasterize'] + cache_config
    rasterize_cmd = rasterize_cmd + ['-a_nodata',str(no_data),'-burn',str(1),'-l',envelope_layer,'-tr',str(delr),str(delc),envelope_shp,envelope_raster]
    subprocess.call(rasterize_cmd)
    
    with rasterio.open(envelope_raster,'r') as src:
        unit_x_min,unit_y_min,unit_x_max,unit_y_max = src.bounds

    # Burn the contours for this unit to a raster
    contour_raster_temp1 = contour_raster.replace('.tif','_temp1.tif')
    contour_raster_temp2 = contour_raster.replace('.tif','_temp2.tif')
    burn_contours(contour_shp,contour_layer_name,contour_raster_temp1,elev_field)

    # Clip the contours to the polygon for only the non-subcrop subsurface extents (i.e., the total extents
    # of the hydrogeologic unit - (outcrop + subcrop)).
    clip_contours_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-cutline',extent_shp,'-of','GTiff','-ot','Float64','-cl',extent_layer,'-tr',str(delr),str(delc),'-dstnodata',str(no_data),contour_raster_temp1,contour_raster_temp2]
    subprocess.call(clip_contours_cmd)    
    
    # Reframe the contours within a bounding box that includes the total extents for this hydrogeologic unit
    # (i.e., including the subcrop and outcrop). This bounding box serves as the reference bounding box for the 
    # later numpy operation that merges and interpolates the outcrop, subcrop, and contour rasters for this unit.
    reframe_contours_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-te',str(unit_x_min),str(unit_y_min),str(unit_x_max),str(unit_y_max),'-of','GTiff','-ot','Float64','-tr',str(delr),str(delc),'-dstnodata',str(no_data),contour_raster_temp2,contour_raster]
    subprocess.call(reframe_contours_cmd)

    # Clip the regional DEM to the outcrop and then reframe to the reference bounding box.
    print '\nClipping the land surface to the outcrop.\n'    
    outcrop_temp = outcrop_raster.replace('.tif','_temp.tif')
    clip_dem_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-cutline',outcrop_shp,'-of','GTiff','-ot','Float64','-cl',outcrop_layer,'-tr',str(delr),str(delc),'-dstnodata',str(no_data),regional_dem,outcrop_temp]
    subprocess.call(clip_dem_cmd)
    
    reframe_dem_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-te',str(unit_x_min),str(unit_y_min),str(unit_x_max),str(unit_y_max),'-of','GTiff','-ot','Float64','-cl',outcrop_layer,'-tr',str(delr),str(delc),'-dstnodata',str(no_data),outcrop_temp,outcrop_raster]
    subprocess.call(reframe_dem_cmd)

    # Clip the bottom of the surficial aquifer to the subcrop and then reframe to the reference bounding box.
    print '\nClipping the bottom of the surficial aquifer to the subcrop.\n'
    subcrop_temp = subcrop_raster.replace('.tif','_temp.tif')
    clip_dem_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-cutline',subcrop_shp,'-of','GTiff','-ot','Float64','-cl',subcrop_layer,'-tr',str(delr),str(delc),'-dstnodata',str(no_data),surficial_dem,subcrop_temp]
    subprocess.call(clip_dem_cmd)
    
    reframe_dem_cmd = ['gdalwarp'] + cache_config + ['-wm',str(cachemax/2)] + \
        ['-te',str(unit_x_min),str(unit_y_min),str(unit_x_max),str(unit_y_max),'-of','GTiff','-ot','Float64','-cl',outcrop_layer,'-tr',str(delr),str(delc),'-dstnodata',str(no_data),subcrop_temp,subcrop_raster]
    subprocess.call(reframe_dem_cmd)

    # Create a virtual mosaic raster of the contours + outcrop + subcrop
    mosaic_raster = contour_raster.replace('.tif','_mosaic.vrt')
    vrt_raster_cmd = ['gdalbuildvrt','-tr',str(delr),str(delc),mosaic_raster,contour_raster,outcrop_raster,subcrop_raster]
    subprocess.call(vrt_raster_cmd)
    
    # WRITE THE RASTERS TO NUMPY ARRAYS AND INTERPOLATE    
    with rasterio.open(mosaic_raster,'r') as src_raster,\
         rasterio.open(envelope_raster,'r') as mask_raster:

        # Capture the raster specs in order to later write the numpy array to a raster
        kwargs = mask_raster.meta
        
        mosaic = src_raster.read()[0,:,:]
        mask = mask_raster.read()[0,:,:]
        irow,icol = np.shape(mosaic)

        # Mask everything outside of the extents domain in order to restrict
        # the interpolation to only pixels within the data extent of this
        # hydrogeologic unit (i.e., not within all pixels of the bounding box).        
        mosaic[mask == no_data] = np.nan
        
        # Generate boolean arrays specifying where we have and need data
        need_interp = (mosaic == no_data)
        
        mosaic[mosaic == no_data] = np.nan
        have_data = np.isfinite(mosaic)

        # Re-cast the 2D boolean arrays to arrays containing the pixel (x,y) coordinates (dim=(number of points,2)).
        # These are the coordinate arrays required by the interpolation function.       
        xx,yy = np.meshgrid(range(icol),range(irow))
        have_data_xy = np.vstack( (np.ravel(xx[have_data]),np.ravel(yy[have_data])) ).T
        need_interp_xy = np.vstack( (np.ravel(xx[need_interp]),np.ravel(yy[need_interp])) ).T
                
        # The locations with data as a 1D array with data listed in the same order as the coordinates in data_xy
        elev_data = np.ravel(mosaic[have_data])

        # Interpolate and map the arrays of (x,y) coordinates back to a 2D array        
        print '\nInterpolating the unit contours.\n'        
        result = interpolate.griddata(have_data_xy,elev_data,need_interp_xy,method='linear')
        for (j,i),z in zip(need_interp_xy,result):
            mosaic[i,j] = z
        
        # Now that the interpolation is complete, the only remaining NaNs should be
        # outside the data extent of the unit. Re-insert the specified no_data value
        # before returning these to rasters.
        mosaic[mosaic == np.nan] = no_data
                
        # WRITE THE NUMPY ARRAYS TO A RASTER
        full_unit_raster = os.path.join(model_raster_dir,iunit + '.tif')
        print 'Writing the interpolated unit to raster: %s\n' %(full_unit_raster)
        with rasterio.open(full_unit_raster,'w',**kwargs) as dst:
            
            dst.write_band(1,mosaic)
        
    model_unit_raster = clip_to_model(full_unit_raster)
    
    # SAVE THE UNIT RASTER AS AN ASCII ARRAY
    with rasterio.open(model_unit_raster,'r') as src:
        
        model_unit_array = src.read()[0,:,:]
        raster_to_array(iunit,model_unit_array,Kzone_dict)       
       
        # Add the local features, which in the case of the Upper Chester is
        # the Hornerstown Aquifer and Aquia Confining Unit       
        if (iunit == 'SevernCU'):
            
            hrnrstwn_array = model_unit_array + hrnstwn_thick
            aquiaCU_array = hrnrstwn_array + aquiaCU_thick
            
            raster_to_array('Hrnrstwn',hrnrstwn_array,Kzone_dict)
            raster_to_array('AquiaCU',aquiaCU_array,Kzone_dict)
       
    del_dir_contents(scratch_dir)