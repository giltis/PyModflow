

import numpy as np
import os
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.ioff()

# -------------------------------------

def get_dtw(landsurface_file,heads_file):
    '''Returns a 2D array of the depth to water at each (row,col).'''
    
    land_surface = np.genfromtxt(landsurface_file)    
    
    heads = bf.HeadFile(heads_file)
    heads = heads.get_alldata()
    
    return land_surface - heads[-1,0,:,:]

def create_unit_elevs(nrow,ncol,nunit,min_elev,max_elev,fin_prefix,nodata):
    """Reads the 2D elevation arrays for each hydrogeologic unit and
    constructs a 3D array.  Note that unit 1 is defined by BOTTOM elevation,
    while all other units are defined by TOP elevation."""

    print 'Creating unit elevation arrays.\n'
    unit_elev = np.zeros((nunit,nrow,ncol))

    for i in range(nunit):
        
        filename_in = fin_prefix + str(i + 1)
    
        unit_elev[i,:,:] = np.genfromtxt(filename_in)
        check = np.ma.masked_less(unit_elev[i,:,:],nodata)   # Ignore no data
        check_min = np.min(check)
        check_max = np.max(check)
    
        if ( check_min < min_elev ):
            min_elev = check_min
    
        if ( check_max > max_elev ):
            max_elev = check_max

    return unit_elev

def plot_unit_elevs(unit_elev,unit_elev_fout,unit_labels,unit_contours,
                    unit_ticks,nodata,IBOUND,open_water,min_elev,max_elev):
    """Plots the elevations of the hydrogeologic units and saves to a pdf."""

    pp = PdfPages(unit_elev_fout)
    
    [nrow,ncol,nlay] = np.shape(unit_elev)
    nunit = len(unit_labels)
    
    for unit in range(nunit):
        plot_array = np.copy(np.ma.masked_less(unit_elev[:,:,unit],nodata))
        for i in range(nrow):
            for j in range(ncol):
                if ( IBOUND[i,j] == 0 ):
                    plot_array[i,j] = np.nan
        
        plt.figure()
        cmap = plt.cm.jet_r
        cmap.set_under('w',min_elev) 
        plt.contourf(plot_array,unit_contours,cmap=cmap,alpha=.90)
        plt.gca().invert_yaxis()
        plt.colorbar(ticks=unit_ticks,cmap=cmap)
        plt.imshow(open_water,cmap='binary')
        plt.imshow(IBOUND,cmap='gray',alpha=.75)

        if ( unit == 0):
            plt.title('%2s\nElevation of unit BOTTOM (ft above MSL)' 
                       % (unit_labels[unit]))
        else:
            plt.title('%2s\nElevation of unit TOP (ft above MSL)' 
                       % (unit_labels[unit]))
        pp.savefig()
    
    pp.close()
    
    return

def plot_layers(fout,array,cbar_labels=None,cbar_bins=20):
    """Plots layer information and saves to file."""

    nlay,nrow,ncol = np.shape(array)
    
    # If the cbar_labels are not provided, generate them from the data
    if (cbar_labels == None):       
        cbar_min,cbar_max = np.floor(np.nanmin(array)),np.ceil(np.nanmax(array))
        cbar_labels = list(np.linspace(cbar_min,cbar_max,cbar_bins))
        cbar_labels = [int(x) for x in cbar_labels]
        
    cbar_ticks = np.arange(1,len(cbar_labels) + 1)
    
    pp = PdfPages(fout)
    
    for ilay in range(nlay):
        plot_array = np.ma.masked_array(array[ilay,:,:],mask=[0])
        
        plt.figure()
        cmap = plt.cm.jet
        cmap.set_under('w')        
        ax = plt.imshow(plot_array,vmin=np.min(cbar_ticks),
                        vmax=np.max(cbar_ticks),cmap=cmap,interpolation='nearest')     
        cbar = plt.colorbar(ax,ticks=cbar_ticks,boundaries=cbar_ticks)
        cbar.ax.set_yticklabels(cbar_labels)
        cbar.ax.invert_yaxis()
        plt.title('Layer %2i' % (ilay+1))
        pp.savefig()
    
    pp.close()
    
    return

def create_grid_datums(model_top,layer_depths):
    """Constructs 3D arrays of cell bottoms and cell midpoints"""
    
    #[nrow,ncol] = np.shape(model_top)
    nlayer = len(layer_depths)

    bottoms = np.zeros((nrow,ncol,nlayer))
    imodel_top = model_top
    
    # The model-top file may have two instances of the model top elevations
    if (np.shape(model_top)[0] == 2 * nrow):
        imodel_top = model_top[0:nrow,:]
    
    bottoms[:,:,0] = imodel_top - layer_depths[0]

    for ilayer in range(1,nlayer):
        bottoms[:,:,ilayer] = bottoms[:,:,(ilayer-1)] - layer_depths[ilayer]

    midpoints = np.zeros((nrow,ncol,nlayer))
    midpoints[:,:,0] = imodel_top - (.5 * layer_depths[0])

    for ilayer in range(1,nlayer):
        for irow in range(nrow):
            for icol in range(ncol):
    
                midpoints[irow,icol,ilayer] = np.mean((bottoms[irow,icol,(ilayer-1)],
                                                       bottoms[irow,icol,ilayer]))

    return bottoms,midpoints

def write_zones_total(cell_units,strat_fout):
    """Write the stratigraphic codes for each cell to disk; all layers
    in ONE FILE."""
    
    nlayer = np.shape(cell_units)[2]    
    
    with file(strat_fout,'w') as outfile:
        for ilayer in range(nlayer):
            np.savetxt(outfile,cell_units[:,:,ilayer],fmt = '%3i')
    
    return

def write_zones_layers(cell_units,zonelay_fout):
    """Writes the stratigraphic codes for each layer to SEPARATE FILE."""

    extension = os.path.splitext(zonelay_fout)[1]

    nlayer = np.shape(cell_units)[2]

    for ilayer in range(nlayer):
       A = cell_units[:,:,ilayer]
       outfile = zonelay_fout.replace(extension,'Lay' + str(ilayer+1) + extension)
       np.savetxt(outfile,A,fmt = '%3i')
     
    return

def write_zones_namefile(zone_unit_start,zone_4name_fout,zonelay_fout,nlayer):
    """Writes the zone array portion of the name file."""

    with file(zone_4name_fout,'w') as outfile:
        iunit = zone_unit_start    
        for ilayer in range(nlayer):
            outfile.write('data   %-3i  %s\n' %(iunit,zonelay_fout + str(ilayer + 1)))
            iunit = iunit + 1
            
    return

def write_bottoms_total(cell_bottoms,cell_bottoms_fout):
    """Writes the cell bottom elevation for each cell to file."""

    nlayer = np.shape(cell_bottoms)[2]    
    
    with file(cell_bottoms_fout,'w') as outfile:
        for ilayer in range(nlayer):
            np.savetxt(outfile,cell_bottoms[:,:,ilayer],fmt = '%8.2f')
    
    return

def write_bottoms_layers(cell_bottoms,bottoms_fout):
    """Writes the cell bottom elevations for each layer to SEPARATE FILE."""

    extension = os.path.splitext(bottoms_fout)[1]

    nlayer = np.shape(cell_bottoms)[2]

    for ilayer in range(nlayer):
       A = cell_bottoms[:,:,ilayer]
       outfile = bottoms_fout.replace(extension,'Lay' + str(ilayer+1) + extension)
       np.savetxt(outfile,A,fmt = '%3i')
     
    return

def write_starting_heads(model_top,starting_heads_fout,nlayer):
    """Write the starting heads (using the land surface elevation) 
    for each cell to disk"""  
    
    with file(starting_heads_fout,'w') as outfile:
        for ilayer in range(nlayer):
            np.savetxt(outfile,model_top,fmt = '%8.2f')
            
    return

def write_IBOUND_3D(IBOUND_2D_array,IBOUND_3D_fout):
    """Write the 2D x nlayer IBOUND array to disk."""

    np.savetxt(IBOUND_3D_fout,IBOUND_2D_array,fmt = '%4i')
    
    return
    
def write_zone_file(zone_unit_start,ncol,nlayer,zone_fout):
    """Write the chester.zon file to disk."""
    
    with file(zone_fout, 'w') as outfile:
        outfile.write('%-4i\n' % (nlayer+1))
        iunit = zone_unit_start 
        for ilayer in range(nlayer):
            lay_name = 'ZoneLay' + str(ilayer + 1)
            outfile.write('%s\n' % (lay_name))
            outfile.write('EXTERNAL   %-3i  1   (%sI4)  -1\n' % (iunit,ncol))
            iunit = iunit + 1
        outfile.write('RechargeMap\n')
        outfile.write('EXTERNAL   %-3i  1   (%s)  -1\n' % (191,'FREE'))
    return

def change_model_bottom(IBOUND_array,unit_flags_array,bottom_layer):
    """Adjust the bottom of the active model area to a specified layer.
    All strat units INCLUDING THE UNIT DESIGNATED BY bottom_layer
    will be deactivated in the IBOUND array."""   
    
    active_cells_prior = np.sum(IBOUND_array > 0)
    
    IBOUND_array[unit_flags_array >= bottom_layer] = 0
    
    active_cells_post = np.sum(IBOUND_array > 0)
    
    active_cells_top_layer = np.sum(IBOUND_array[0,:,:] > 0)
    
    print 'Total number of cells = %i' %(unit_flags_array.size)
    print 'Number of active cells before reduction = %i' %(active_cells_prior)
    print 'Number of active cells after  reduction = %i' %(active_cells_post)
    print 'Number of active cells in top layer     = %i\n' %(active_cells_top_layer)
    
    return IBOUND_array

def plot_cross_sections(model_top,cell_units,layer_depths,row_fout,col_fout,unit_labels):

    [nrow,ncol,nlay] = np.shape(cell_units)
    model_depth = int(np.sum(layer_depths))
    max_elev = np.max(model_top)
    min_elev = np.min(model_top) - model_depth
    
    cbar_ticks = np.arange(1,len(unit_labels) + 2)
    
    plot_array = np.zeros((nrow,ncol,np.ceil(max_elev - min_elev))) # one layer per foot
    
    # Load the arrays that contain the zone codes and further discretize into one
    # foot layers
    
    layer_array = np.zeros((nrow,ncol,model_depth))
    
    for k in range(nlay):
        counter = k * int(layer_depths[k])
        for kk in range(int(layer_depths[k])):
            layer_array[:,:,counter] = cell_units[:,:,k]
            counter = counter + 1
    
    # Define the top layer of the ZoneArray as the max land surface
    # elevation; assign the top of each stack in the stratigraphic layer array
    # to the appropriate layer of the ZoneArray
    
    for row in range(nrow):
       for col in range(ncol):
          stack_top = model_top[row,col]
          start_layer = np.floor(max_elev - stack_top)
          for layer in range(model_depth):
              if ( layer == start_layer ):
                  stop_layer = start_layer + model_depth
                  plot_array[row,col,start_layer:stop_layer] = layer_array[row,col,:]
    
    multiplier = 10
    
    prnt_rows = PdfPages(row_fout)
    
    for i in range(nrow):
        
            rowslice = i        
            
            iplot = np.rot90(plot_array[rowslice,:,:],3)
            iplot = np.repeat(iplot,multiplier,axis = 1)
            iplot = np.fliplr(iplot)
            
            fig = plt.figure()
            cmap = plt.cm.jet
            cmap.set_under('w')        
            ax = plt.imshow(iplot,vmin=np.min(cbar_ticks),vmax=np.max(cbar_ticks),cmap=cmap)     
            plt.gca().set_yticks(np.arange(0,400,25))
            plt.gca().set_yticklabels([])
            x_ticks = (plt.gca().get_xticks())/multiplier
            x_ticks = x_ticks.astype(int)
            plt.gca().set_xticklabels(x_ticks)
            cbar = plt.colorbar(ax,ticks=cbar_ticks,boundaries=cbar_ticks)
            cbar.ax.set_yticklabels(unit_labels)
            cbar.ax.invert_yaxis()
            plt.xlabel('Column Number')
            plt.ylabel('Ticks on Vertical Scale = 25 feet')
            plt.title('Cross Section along Row %i\n<- West  East ->' %(rowslice + 1))
            prnt_rows.savefig()
            plt.close()
            print 'Printing Cross Section for Row %i' %(i+1)
    
    prnt_rows.close()
    
    prnt_cols = PdfPages(col_fout)   
    
    for j in range(ncol):
    	
            colslice = j        
            
            jplot = np.rot90(plot_array[:,colslice,:],3)
            jplot = np.repeat(jplot,multiplier,axis = 1)
            jplot = np.fliplr(jplot)
            
            fig = plt.figure()
            cmap = plt.cm.jet
            cmap.set_under('w')        
            ax = plt.imshow(jplot,vmin=np.min(cbar_ticks),vmax=np.max(cbar_ticks),cmap=cmap)
            plt.gca().set_yticks(np.arange(0,400,25))
            plt.gca().set_yticklabels([])
            x_ticks = (plt.gca().get_xticks())/multiplier
            x_ticks = x_ticks.astype(int)
            plt.gca().set_xticklabels(x_ticks)
            cbar = plt.colorbar(ax,ticks=cbar_ticks,boundaries=cbar_ticks)
            cbar.ax.set_yticklabels(unit_labels)
            cbar.ax.invert_yaxis()    
            plt.xlabel('Row Number')
            plt.ylabel('Ticks on Vertical Scale = 25 feet')
            plt.title('Cross Section along Column %i\n<- North  South ->' %(colslice + 1))
            prnt_cols.savefig()
            plt.close()
            print 'Printing Cross Section for Col %i' %(j+1)
    
    prnt_cols.close()
    
    return
    
def write_multiplier_array(multiplier_array,fout_root,ilayer):
    
    fout = fout_root + str(ilayer + 1)  # Adjust python indexing
    
    np.savetxt(fout,multiplier_array,fmt='%9.2f')
    
    return
    
    
    