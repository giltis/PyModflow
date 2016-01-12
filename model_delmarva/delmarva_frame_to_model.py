""" This program reads elevation data for the unit_elevsigraphic units (in the form of files unit_elevs01, unit_elevs02, etc.)
and sorts them into model grid cells """

import numpy as np
from pygrid import grid_util as gu
from pygrid import framework_tools as frame
from model_config import * # Global variables here

import matplotlib.pyplot as plt

# -----------------

def del_frame_dir(model_frame_dir):
    '''Deletes the contents of the model framework directory.'''
    
    for ifile in os.listdir(model_frame_dir):
        file_to_delete = os.path.join(model_frame_dir,ifile)
        if os.path.isfile(file_to_delete):
            os.unlink(file_to_delete)
    
    return

def clean_ibound(ibound):
    '''Removes isolated active cells from the IBOUND array.'''
    
    nlay,nrow,ncol = np.shape(ibound)
    clean_ibound = np.copy(ibound)
    for ilay in range(nlay):

        this_layer = clean_ibound[ilay,:,:]
        has_neighbor = np.zeros(this_layer.shape,bool)
        
        has_neighbor[:,1:]  = np.logical_or(has_neighbor[:,1:], this_layer[:,:-1] > 0)  # left
        has_neighbor[:,:-1] = np.logical_or(has_neighbor[:,:-1], this_layer[:,1:] > 0)  # right
        has_neighbor[1:,:]  = np.logical_or(has_neighbor[1:,:], this_layer[:-1,:] > 0)  # above
        has_neighbor[:-1,:] = np.logical_or(has_neighbor[:-1,:], this_layer[1:,:] > 0)  # below

        this_layer[np.logical_not(has_neighbor)] = 0
        ibound[ilay,:,:] = this_layer

    return ibound

# -----------------

if (raw_input('Would you like to use the land surface as the model top? (Y/N): ') in ['Y','y','yes','Yes','YES']):    
    model_top_fin = landsurface_elev_file    
else:   
    model_top_fin = adj_model_top_file
print 'Using %s for model top.\n' %(model_top_fin)

print '\nReading basic framework arrays.\n'
ibound = np.genfromtxt(landsurface_ibound_file)
ibound = np.tile(ibound,(nlay,1))
ibound = gu.twoD_to_threeD(ibound,(nrow,ncol))
ibound[np.isnan(ibound) == True] = 0

top = np.genfromtxt(model_top_fin)
landsurface = np.genfromtxt(landsurface_elev_file)

# Create 3D array of unit elevations
unit_elevs = frame.create_unit_elevs(nrow,ncol,nunit,min_elev,max_elev,\
                                   strat_fin_prefix,nodata)

# Create arrays of cell bottom elevations and cell midpoint elevations for current model top
print 'Calculating cell bottoms and midpoints.\n'
bottoms,midpoints = gu.create_grid_datums(model_top_fin,lay_thick)

# --- START MAIN ASSIGNMENT LOOP --------------------

# Using the datum elevation for each model grid cell
# assign a unit_elevsigraphic unit to that cell

cell_units = np.zeros((nlay,nrow,ncol))

for ilay in range(nlay):
    for irow in range(nrow):
        for icol in range(ncol):
            
            itop = top[irow,icol]
            iland = landsurface[irow,icol]
            imid = midpoints[ilay,irow,icol]
            
            stack = np.copy(unit_elevs[:,irow,icol])
                         
            # unit_elevs unit tops that are higher than the bottom of the surficial aquifer
            # are interpolation artifacts and may not actually exist;
            # a unit is assumed not to exist if its top is higher than the bottom
            # of the Pensauken AND if its top is higher than at least one other unit
            # whose top is also higher than the bottom of the Pensauken.
            # The lowest unit top elevation that is higher than the bottom of
            # the Pensauken is set equal to the bottom of the Pensauken.
            
            # find all the elements that are higher than the Pensauken
            # and send them to an index array
            TempA = np.ones((nunit)) * 999
            for strat_layer in range(1,nunit):
                if ( stack[strat_layer]  >= stack[0] ):
                    TempA[strat_layer] = stack[strat_layer]
            
            # set the units that are non-minimum nonzero elements in the Temp array
            # to -8888 in the stack array
            minvalA = np.min(TempA)
            for strat_layer in range(1,nunit):
                if ( TempA[strat_layer] != minvalA ) and ( stack[strat_layer] >= stack[0]):
                    stack[strat_layer] = -8888
                    
            # set the remaining value in the stack array that is greater than 
            # the bottom of the Pensauken equal to the bottom of the Pensauken
            for strat_layer in range(1,nunit):
                if ( stack[strat_layer] >= stack[0] ):
                    stack[strat_layer] = stack[0]
                    index = strat_layer					# keep of track of any values that change           
            
            # if all unit tops are less than the Pensauken (i.e., if there is a gap
            # between the the top of the topmost unit and the bottom of the Pensauken) then
            # activate the unit above the highest active unit by setting its top
            # elevation equal to the bottom of the Pensauken

            if (nunit > 1):
                maxval = np.max(stack[1:nunit])
                if ( maxval < stack[0] ):
                    if ( stack[1] == maxval ): # automatically close gap with CalvertCU if CalvertCU is active
                        stack[1] = stack[0]
                    else:
                        for strat_layer in range(2,nunit):
                            if ( stack[strat_layer] == maxval ):
                                maxunit = strat_layer
                        for strat_layer in range(2,nunit):
                            if ( strat_layer == maxunit ):
                                stack[strat_layer - 1] = stack[0]
            
            # find all the elements that are higher than the land surface
            # and send them to a temp array
            TempB = np.ones((nunit)) * 999
            for strat_layer in range(nunit):
                if ( stack[strat_layer] ) >= iland:
                    TempB[strat_layer] = stack[strat_layer]
                    
            # set the units that are non-minimum nonzero elements in the Temp array
            # equal to -7777 in the stack array (i.e., eliminate all units that are
            # above the land surface except for the lowest unit)
            minvalB = np.min(TempB)
            for strat_layer in range(nunit):
                if ( TempB[strat_layer] != minvalB ) and ( stack[strat_layer] >= iland):
                    stack[strat_layer] = -7777
            
            # set the remaining value in the stack array that is greater than 
            # the land surface elevation equal to the land surface
            for strat_layer in range(0,nunit):
                if ( stack[strat_layer] > itop ):
                    stack[strat_layer] = itop
                    if ( strat_layer >= 1 ) and ( np.max(stack[:strat_layer]) == iland ):
                        stack[np.where(stack[:strat_layer] == itop)] = -7777      
                    
            # sort the stack in order to compare the unit elevations to the
            # grid cell elevations
            
            Sorted_Stack = sorted(stack,reverse=True)
            
            unit_elevs_val = Sorted_Stack[0]
                    
            # find the lowest unit_elevs unit that the cell elevation is below
            for strat_layer in range(1,nunit):
                if ( imid < Sorted_Stack[strat_layer] ):
                    unit_elevs_val = Sorted_Stack[strat_layer]
            
            # determine the unsorted unit_elevs layer number associated with that sorted value
            unit_elevs_index = int(np.max(np.where(stack == unit_elevs_val))) + 1
            
            if ( imid >= unit_elevs[0,irow,icol] ):
                if ( imid <= itop ):
                    unit_elevs_index = 1
            
            cell_units[ilay,irow,icol] = unit_elevs_index
                
            del unit_elevs_index
   
    print 'Finished assigning unit codes to Layer %2i' % (ilay + 1)

# --- END MAIN ASSIGNMENT LOOP --------------------

if (change_bottom_flag == True):   
    bottom_unit = Kzone_dict[bottom_unit_name]    
    print '\nApplying the %s as the bottom of the model.' %(bottom_unit_name)
    ibound = frame.change_model_bottom(ibound,cell_units,bottom_unit)
    unit_labels = unit_labels[0:bottom_unit]

# Convert isolated active cells to inactive cells. This is required to
# avoid numerical issues (i.e., unrealistically isolated heads) in the head solution.
ibound = clean_ibound(ibound)
              
# Delete the framework directory and write the various MODFLOW framework input files
del_frame_dir(model_frame_dir)
np.savetxt(ibound3D_file,gu.threeD_to_twoD(ibound),fmt='%10i')
np.savetxt(modeltop_file,top,fmt='%5.2f')
np.savetxt(bottoms3D_file,gu.threeD_to_twoD(bottoms),fmt='%10.3f')
np.savetxt(geo_zones3D_file,gu.threeD_to_twoD(cell_units),fmt='%6i')
np.savetxt(starting_heads3D_file,np.tile(top,(nlay,1)),fmt='%10.3f')

if (plot_geozones == True):
    
    izones = np.copy(cell_units[:,:,:])
    izones[ibound == 0] = np.nan
    print 'Saving geozone plots to %s.\n' %(geozones_plot_file)
    frame.plot_layers(geozones_plot_file,izones,cbar_labels=unit_labels)    
    
if (plot_xsect_flag == True):
    frame.plot_cross_sections(top,cell_units,layer_depths,xsect_row_pdf,xsect_col_pdf,unit_labels)