# -*- coding: utf-8 -*-
"""
Generates a PEST 3D pilot points file.
Columns = name, easting, northing, elevation, zone number, value

NB: In most (and perhaps all) implementations of the pilot point
file, the pilot point file is accompanied by a PEST template file that
directs PEST to update the 'value' field with the parameter value generated
by the current optimization iteration.  For cases in which pilot points are
used to calibrate the conductivity field, value is likely = log(conductivity).

WORKFLOW:
pilot points -> PPK2FAC3D -> FAC2REAL3D -> conductivity field -> MODFLOW LPF PACKAGE
pilot points -> PEST template file -> PEST control file (parameter data section)

NOTE: The current version of this script is built on the dissertation version
of the model_config workflow.  Adjustments may need to be made in order
to accomodate integration of flopy utilities.

OPTIONS TO ADD:
(i) restrict pilot points to user-specified zone
(ii) add additional pilot point grids with variable resolution

Created on Thu Oct 08 16:08:01 2015
@author: wzell
"""

import numpy as np
import os
from itertools import product
from mfpytools import grid_util as gu
import string

import matplotlib.pyplot as plt

from model_config import nrow,ncol,nlay,dx_dy,layer_depths,IBOUND3D_file,modeltop_file,zonelay_fout

# --- START PARAMETER SET ---

ppt_name_root = 'pt'
delr,delc = dx_dy,dx_dy
nx,ny = 8,8     # The x and y dimensions of the pilot point grid
ppt_dz = 75     # The vertical distance between pilot points. The first pilot point is placed midpoint in the top layer.

ppt_fout = os.path.join('PEST\PEST_Input','HK.pts')
ppt_tpl_fout = ppt_fout.replace('pts','tpl')
ppt_replace_flag = '$'
param_name_root = 'HK'
mu,sigma = 10,2

# --- STOP PARAMETER SET ----

# Generate the model grid
ibound = np.genfromtxt(IBOUND3D_file)
ibound = gu.twoD_to_threeD(ibound,(nrow,ncol))
modeltop = np.genfromtxt(modeltop_file)[:nrow,:]
bottoms,midpoints = gu.create_grid_datums(modeltop_file,nrow,ncol,layer_depths)
zones_3d = gu.create_layer_3D(nrow,ncol,nlay,zonelay_fout)

# Generate the pilot point grid in the row/col plane
x_vector = np.linspace(0,delr * ncol,nx+1)[:-1]
x_vector = np.add(x_vector,0.5*(x_vector[1]-x_vector[0]))

y_vector = np.linspace(0,delc * nrow,ny+1)[:-1]
y_vector = np.add(y_vector,0.5*(y_vector[1]-y_vector[0]))

# Evaluate the pilot point locations and write the active pilot points to file.
# Note that grid_util.globalxy_to_localxy
# uses the northwest corner of the model grid as the origin.  However, the PEST
# pilot points file uses the SW corner of the model grid as the origin (i.e.,
# it prints coordinates as northings and eastings.)

with open(ppt_fout,'w') as ppt_fout, open(ppt_tpl_fout,'w') as tpl_fout:
    
    tpl_fout.write('ptf ' + ppt_replace_flag + '\n')

    plot_x,plot_y = [],[]
    
    for ix,jy in product(x_vector,y_vector):
        
        # row and column returned from this function in MODFLOW indexing
        irow,icol,_,_ = gu.globalxy_to_localxy(ix,(nrow*dx_dy - jy),dx_dy,dx_dy,nrow,ncol)
    
        # Generate the vertical dimension of the pilot point grid
        # This algorithm assumes that for any (row,col), no active cells are found
        # beneath any inactive cells.
        z_0 = midpoints[irow-1,icol-1][0]
        z_vector = np.arange(z_0,z_0 - np.sum(layer_depths),-ppt_dz)
        
        # The PEST 10-character limit precludes
        # writing row,col, and layer identifiers into the pilot point name.
        # As an alternative, I've simply labeled them A,B,C, etc., starting
        # with the top pilot point in a given (row,col)
        for iz,iletter in zip(z_vector,string.uppercase[:(len(z_vector))]):
    
            ilay,_ = gu.globalz_to_localz(irow,icol,iz,midpoints,bottoms,layer_depths)
            if (gu.check_if_active(irow,icol,ilay,ibound) == True):
                
                izone = zones_3d[irow-1,icol-1,ilay-1]
                ppt_name = ppt_name_root + str(irow) + '_' + str(icol) + iletter
                param_name = param_name_root + str(irow) + '_' + str(icol) + iletter
                
                ival = np.log(np.random.normal(mu,sigma))                
                
                #name, easting, northing, elevation, zone number, value
                ppt_fout.write('%-12s%12.2f%12.2f%12.2f%5i%2s%20.5e%2s\n' %(ppt_name,ix,jy,iz,izone,'',ival,''))
                tpl_fout.write('%-12s%12.2f%12.2f%12.2f%5i%2s%20s%2s\n' %(ppt_name,ix,jy,iz,izone,ppt_replace_flag,param_name,ppt_replace_flag))
                
                plot_x.append(icol)
                plot_y.append(irow)
                
            else:
                break

print '\n%i pilot points written to %s\n' %(len(plot_x),ppt_fout)
        
plt.figure()
plt.imshow(ibound[:,:,0])
plt.plot(plot_x,plot_y,'yo')
plt.show()