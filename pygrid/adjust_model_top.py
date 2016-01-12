# -*- coding: utf-8 -*-
"""
This script generates an array to be used as the model top elevation
for a confined simulation using the simulated water table elevation and
the land surface elevation. The elevation
of each cell (i,j) is equal to the water table elevation for that cell
unless the water table is above the land surface (e.g., with an actively
discharging drain), in which case the elevation is set equal to the land surface elevation.

Created on Tue Sep 10 16:15:17 2013
@author: Wesley Zell
"""

import numpy as np
import pylab as plt
from scipy import ndimage
from mfpytools import binaryfile as bf
from mfpytools import grid_util as gu
from mfpytools import read_binaries as rb
from model_config import *

def max_gradients(array,string,landsurface_IBOUND):
    
    [nrow,ncol] = np.shape(array)
    
    gradients_array = np.zeros((nrow,ncol))    
    
    for irow in range(nrow):
        for jcol in range(ncol):
            
            if ( landsurface_IBOUND[irow,jcol] != 0 ):
        
                max_gradient = 0   
            
                adj_rows = [irow,irow,irow-1,irow+1] #,irow-1,irow-1,irow+1,irow+1]
                adj_cols = [jcol-1,jcol+1,jcol,jcol] #,jcol-1,jcol+1,jcol-1,jcol+1]
                   
                for check in range(len(adj_rows)):
                            
                    check_row = adj_rows[check]
                    check_col = adj_cols[check]
                    
                    if ( check_row >= 0 and check_row < 300
                    and check_col >= 0 and check_col < ncol ):
                        
                        if ( landsurface_IBOUND[check_row,check_col] != 0 ):                
                                                
                            igradient = np.abs(model_top[irow,jcol] - model_top[check_row,check_col])
                            max_gradient = np.max((igradient,max_gradient))
                
                gradients_array[irow,jcol] = max_gradient
    
    grad_row = np.where(gradients_array == np.nanmax(gradients_array))[0][0]
    grad_col = np.where(gradients_array == np.nanmax(gradients_array))[1][0]
    location = (grad_row,grad_col)    
    
    print_string = '%10s: Max active gradients = %3.1f at %s' %(string,np.nanmax(gradients_array),location)
    print print_string    
    
    return gradients_array,print_string

# ------------

# Read in the land surface and IBOUND arrays
landsurface_elev = np.genfromtxt(landsurface_elev_fin)
landsurface_IBOUND = np.genfromtxt(landsurface_ibound_fin)

if (raw_input('Use land surface as model top basis? (Y/N)\n (N = Use simulated water table): ') == 'Y'):
    use_land = True
    model_top = landsurface_elev
    
else:
    use_land = False
    print 'Using the water table as the model top basis.'
    print 'No flow areas will be filled with land surface elevations before filtering.\n'
 
    # Read in the cbc file and extract the boundary conditions to arrays
    cbc_dict = rb.read_budget(cbc_file,nrow,ncol,nlay)
    drains = cbc_dict['DRAINS'][:,:,0]
    #ghb_dict = cbc.get_data(kstp=1, kper=1,text='HEAD DEP BOUNDS')
    #ghb = gu.find_active_bc(ghb_dict,(nrow,ncol))
    
    # Read in the heads file and prepare it as the base for filtering
    heads = rb.read_heads(simulated_heads,nrow,ncol,nlay)[:,:,0]

    model_top = np.copy(heads)
    
    # Replace the no-flow flag in inactive cells with land surface elevation
    for i in range(nrow):
        for j in range(ncol):
            
            if (model_top[i,j] < 0):            
                model_top[i,j] = landsurface_elev[i,j]

gradients,print_string = max_gradients(model_top,'Unfiltered  ',landsurface_IBOUND)

plt.figure()
plt.imshow(model_top)
plt.title(print_string)
plt.savefig(os.path.join(flowmodel_output_dir,'Unfiltered'))

for step in range(nsteps):
    
    # Use the gaussian filter to smooth, but reset those cells
    # that are actively discharging drains to the landsurface elevations
    # (this is in order to retain valleys)
    gaussian = ndimage.filters.gaussian_filter(model_top,sigma)
    
    for i in range(nrow):
        for j in range(ncol):
            
            model_top[i,j] = np.min((gaussian[i,j],landsurface_elev[i,j]))
    
    gradients,print_string = max_gradients(model_top,('Filter # %3i' %(step+1)),landsurface_IBOUND)

    plt.figure()
    plt.imshow(model_top)
    plt.title(print_string)
    plt.savefig(os.path.join(flowmodel_output_dir,'Step_' + str(step)))
    
print 'Writing new model top to %s.' %(adj_model_top_file)
np.savetxt(adj_model_top_file,model_top,fmt='%9.2f')