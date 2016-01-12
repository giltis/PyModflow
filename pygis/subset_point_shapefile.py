# -*- coding: utf-8 -*-
"""
Reads an UNPROJECTED (i.e., lat,lon) point shapefile and creates a PROJECTED subset.
APPLICATION: used to generate the pour point shapefile that is used
in watershed delineation.

vin,vout: 'vector_in','vector_out'

Created on Mon Nov 23 12:21:36 2015
@author: wzell
"""

import os,fiona,pyproj
from fiona.crs import from_epsg
from shapely.geometry import mapping,shape

from mfpytools import grid_util as gu

from model_config import model_epsg,data_dir

# --- START PARAMETER SET ---

shp_epsg = 4269

work_dir = os.path.join(data_dir,'Gages')
shp_fin  = os.path.join(work_dir,'Synoptic_Stations.shp')
shp_fout = os.path.join(work_dir,'Model_Gages.shp')

key_field = 'site_no'   # select on this field
select_dict = {'1493010':13,'1493190':17,'1493250':18,'1493195':19,'1493500':1,'1493112':20}

# --- STOP PARAMETER SET ----

# Define the coordinate system transform
epsg_projected = pyproj.Proj('+init=EPSG:' + str(model_epsg),preserve_units=True)

with fiona.open(shp_fin,'r') as vin:
    
    source_driver = vin.driver
    source_schema = vin.schema

# Add an additional schema field that maps the IBOUND zone number to the gage information
save_schema = source_schema
save_schema['properties'] = {'IBOUND':'int'}

with fiona.open(shp_fin) as vin, \
     fiona.open(shp_fout,'w',driver=source_driver,crs=from_epsg(model_epsg),schema=save_schema) as vout:
         
         for ifeature in vin:
             ikey = str(int(ifeature['properties']['site_no']))
             
             if (ikey in select_dict):
                 
                izone = select_dict[ikey]
                
                igeometry = shape(ifeature['geometry'])

                ilon,ilat = mapping(igeometry)['coordinates']
                ix,iy,_,_ = gu.latlon_to_modelxy((ilat,ilon),(0,0),epsg_projected)
                out_geometry = {'type':'Point','coordinates':(ix,iy)}
           
                vout.write({'geometry':out_geometry,'properties':{'IBOUND':izone}})