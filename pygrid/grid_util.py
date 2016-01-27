# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 10:56:55 2013

@author: Wesley Zell
"""
import numpy as np
import pyproj
import xlsxwriter

# ----------------------------
# --- CALLED BY HEADOBS.PY ---
# ----------------------------
def latlon_to_modelxy((ilat,ilon),(x0,y0),epsg_projected,x_increases_east=True):
    '''Projects NWIS (lat,lon) to grid global coordinates.  Current version
    assumes no grid rotation.'''
    
    proj_x,proj_y = epsg_projected(ilon,ilat)
    
    if (proj_x < x0):
        return np.nan,np.nan,np.nan,np.nan
        
    else:
        model_x,model_y = (proj_x - x0),(proj_y-y0)
        return proj_x,proj_y,model_x,model_y
        
def bbox_to_latlon(bbox,model_epsg):
    '''Returns the lat/lon equivalent of a projected bounding box.'''

    epsg_unproject = pyproj.Proj(init='epsg:'+str(model_epsg))

    xmin,ymin,xmax,ymax = bbox
    west,south = epsg_unproject(xmin,ymin,inverse=True)
    east,north = epsg_unproject(xmax,ymax,inverse=True)
        
    return west,south,east,north

def create_grid_datums(top,lay_thick):
    """Constructs 3D arrays of cell bottoms and cell midpoints"""  
       
    nrow,ncol = np.shape(top)
    nlay = len(lay_thick)

    bottoms = np.zeros((nlay,nrow,ncol))
    bottoms[0,:,:] = top - lay_thick[0]

    for ilayer in range(1,nlay):
        bottoms[ilayer,:,:] = bottoms[(ilayer-1),:,:] - lay_thick[ilayer]

    midpoints = np.zeros((nlay,nrow,ncol))
    midpoints[0,:,:] = top - (.5 * lay_thick[0])

    for ilay in range(1,nlay):
        for irow in range(nrow):
            for icol in range(ncol):
    
                midpoints[ilay,irow,icol] = np.mean((bottoms[(ilay-1),irow,icol],
                                                       bottoms[ilay,irow,icol]))

    return bottoms,midpoints

def globalz_to_localz(irow,icol,globalz,midpoints,bottoms,lay_thick):
    '''Converts a global z value to a layer number with a local z ratio. MODPATH
    local z = 0 at bottom of cell and 1 at top of cell. NOTE: Expects
    PYTHON INDEXING for irow,icol and RETURNS PYTHON INDEXING.'''    
    
    stacked_midpoints = midpoints[:,irow,icol]
    ilay = (np.abs(stacked_midpoints - globalz)).argmin()
     
    cell_bottom = bottoms[ilay,irow,icol]
    localz = (globalz - cell_bottom)/lay_thick[ilay]
     
    return ilay,localz

def modelxy_to_rowcol(x,y,delr,delc):
    '''Converts the model x,y coordinates to (row,col) with local offsets.
    RETURNS PYTHON BASED INDEXING. Assumes the origin=NW corner of an unrotated
    grid.'''
    
    x,y = abs(x),abs(y)    
    
    icol = int(np.floor(x/delr))
    part_col = x - icol*delr
    COFF = (part_col - (delr/2))/delr
    
    irow = int(np.floor(y/delc))
    part_row = y - irow*delc
    ROFF = (part_row - (delc/2))/delc
    
    return irow,icol,ROFF,COFF    

def check_if_active(irow,icol,ilay,IBOUND_3D):
    '''Checks row,col index against IBOUND in order to filter inactive cells.
    EXPECTS PYTHON-BASED INDEXING.'''

    check = False
      
    if (IBOUND_3D[ilay,irow,icol] > 0):        
        check = True    
    
    return check

def globalz_to_layer(row,modeltop,bottoms,midpoints,ibound3D,layer_thick):
    '''Flag those observations that are outside the active model area. If
    active, return the layer and local z. Else return NaNs. RETURNS
    PYTHON BASED INDEXING.'''
    
    nrow,ncol = np.shape(modeltop)
    iname,irow,icol,iwellelev = row['station_nm'],row['Row'],row['Column'],row['well_elev']
    
    if (irow in range(nrow) and icol in range(ncol)):          
        ilay,localz = globalz_to_localz(irow,icol,iwellelev,midpoints,bottoms,layer_thick)               
    else:
        print 'Removed as inactive:%10s%5i%5i (MODFLOW Indexing)' %(iname,irow+1,icol+1)
        return np.nan,np.nan,np.nan
 
    if (check_if_active(irow,icol,ilay,ibound3D) == True):
        return modeltop[irow,icol],ilay,localz       
    else:
        print 'Removed as inactive:%10s%5i%5i%5i (MODFLOW Indexing)' %(iname,irow+1,icol+1,ilay+1)
        return np.nan,np.nan,np.nan

def modelxy_to_localxy(model_x,model_y,delr,delc,nrow,ncol):
    '''Converts model x and y coordinates to a row,column index and a local x
    and y. Local x and local y reported in MODPATH terms (local x INCREASES
    along increasing col index, local y DECREASES along increasing row index).
    Global X,Y origin at NW corner of finite difference grid.
    Returns row and column number in PYTHON indexing.'''
    
    jcol = int(np.floor(model_x/delr))
    local_x = np.mod(model_x,delr)/delr
           
    irow = int(np.floor(model_y/delc))
    local_y = 1 - np.mod(model_y,delc)/delc    # global y and row indices increase in opposite directions

    if (irow < 0):
        irow = 0
        local_y = 1
    
    if (irow > nrow-1):
        irow = nrow-1
        local_y = 0
        
    if (jcol < 0):
        jcol = 0
        local_x = 0
        
    if (jcol > ncol-1):
        jcol = ncol-1
        local_x = 1
    
    return irow,jcol,local_x,local_y

# -----------------------------------------------------------------------------

def arrays_to_xls(array_list,sheet_name_list,xls_packet_fout):
    '''Writes arrays to Excel spreadsheet.'''
    
    ibook = xlsxwriter.Workbook(xls_packet_fout,{'nan_inf_to_errors': True}) 
    
    # Write the arrays to a bundle of Excel sheets to facilitate array inspection
    for iarray,isheet_name in zip(array_list,sheet_name_list):
                                  
        if (iarray.ndim > 2):    # 3D array           
            for ilay in range(np.shape(iarray)[0]):
                jarray = iarray[ilay,:,:]
                isheet_name = isheet_name + '_' + str(ilay+1)
                isheet = ibook.add_worksheet(isheet_name)
                print 'Writing sheet: ',isheet_name
                for (irow,icol),ival in np.ndenumerate(jarray):
                    isheet.write(irow,icol,ival)
                
        else:   # For one layer        
            isheet = ibook.add_worksheet(isheet_name)
            print 'Writing sheet: ',isheet_name               
            for (irow,icol),ival in np.ndenumerate(iarray):
                isheet.write(irow,icol,ival)                                
    
    ibook.close()
    
    return

def clean_ibound(ibound):
    '''Removes isolated active cells from the IBOUND array.'''
    
    if (ibound.ndim > 2):
        nlay,nrow,ncol = np.shape(ibound)
    else:
        nrow,ncol = np.shape(ibound)
        nlay = 1
		
    clean_ibound = np.copy(ibound)

    for ilay in range(nlay):

        if (ibound.ndim > 2):
            this_layer = clean_ibound[ilay,:,:]
        else:
            this_layer = clean_ibound
            
        has_neighbor = np.zeros(this_layer.shape,bool)
        
        has_neighbor[:,1:]  = np.logical_or(has_neighbor[:,1:], this_layer[:,:-1] > 0)  # left
        has_neighbor[:,:-1] = np.logical_or(has_neighbor[:,:-1], this_layer[:,1:] > 0)  # right
        has_neighbor[1:,:]  = np.logical_or(has_neighbor[1:,:], this_layer[:-1,:] > 0)  # above
        has_neighbor[:-1,:] = np.logical_or(has_neighbor[:-1,:], this_layer[1:,:] > 0)  # below

        this_layer[np.logical_not(has_neighbor)] = 0
        
        if (ibound.ndim > 2):
            clean_ibound[ilay,:,:] = this_layer
        else:
            clean_ibound = this_layer

    return clean_ibound

def get_exp_decay(Y0,lam,x_array,asymptote=0):
    '''Returns an ndim array of exponentially-decayed values.'''
    
    return (Y0 - asymptote) * np.exp((-lam) * x_array) + asymptote

def localxy_to_globalxy(irow,icol,localx,localy,dx,dy):
    '''Converts a row,col tuple with localx,localy to a global xy coordinate.
    Origin is in NW corner of model domain.'''
    
    globalx = ((icol-1) * dx) + (localx * dx)
    globaly = ((irow-1) * dy) + ((1-localy) * dy)
    globaly = - globaly
    
    return globalx,globaly

def localz_to_globalz(irow,icol,ilay,localz,model_top,cell_bottoms):
    '''Converts a (row,col,lay) tuple with localz to a global z coordinate.'''

    if (ilay - 1 == 0): 
        lay_top = model_top[irow-1,icol-1]
    else:
        lay_top = cell_bottoms[irow-1,icol-1,ilay-2]

    lay_bot = cell_bottoms[irow-1,icol-1,ilay-1]        
        
    globalz = lay_bot + localz*(lay_top - lay_bot)
       
    return globalz

def get_water_table(i_array,cell_bottoms,layer_depths):
    '''Returns an (nrow,ncol) 2D array with array values = water table elevation
    and a 2D array with array values = layer in PYTHON INDEXING in which
    the water table occurs.'''
    
    nrow,ncol,nlay = np.shape(i_array)
    wt = np.zeros((nrow,ncol))
    wt_layers = np.zeros((nrow,ncol))
    local_zs = np.zeros((nrow,ncol,nlay))
    
    for i in range(nrow):
        for j in range(ncol):
    
            i_heads = np.array(i_array[i,j,:])            
            i_bottoms = np.array(cell_bottoms[i,j,:])            
            i_check = i_heads > i_bottoms
            iwt_layer = np.max(np.where(i_check == False)) + 1
            wt[i,j] = i_heads[iwt_layer]
            wt_layers[i,j] = iwt_layer

            local_zs[i,j] = (wt[i,j] - i_bottoms[iwt_layer]) * (1./layer_depths[iwt_layer])

    return wt,wt_layers,local_zs

#def write_watertable(wt,bottoms):
#    '''Write the water table characteristics to file.'''
#
#    with open('WaterTable','w') as fout:
#        
#        for iper in range(nperiods_per_year):
#            
#            # Calculate the saturated thickness as the difference between the water table
#            # and the bottom of the model. Note that this is total saturated thickness, not
#            # the sat thickness for a particular hydrogeologic layer, and it may be interrupted by confining units            
#            
#            isat_thck = wt[:,:,iper] - bottoms[:,:,-1]            
#            
#            fout.write('Stress Period %i\n' %(iper))
#            np.savetxt(fout,wt[:,:,iper],fmt='%5.2f')
#            fout.write('Mean saturated thickness = %5.2f\n' %(np.mean(isat_thck)))
#      
#    return


def get_water_table_info(heads_array,bottoms,nrow,ncol,nstp,layer_depths):
    '''Returns an (nrow,ncol,nstp) array with each value = to the layer (PYTHON INDEXING)
    in which the water table occurs.'''

    water_table = np.zeros((nrow,ncol,nstp))
    water_table_layers = np.zeros((nrow,ncol,nstp))
    local_zs = np.zeros((nrow,ncol,nlay,nstp))
    global_zs = np.zeros((nrow,ncol,nstp))
    
    for istp in range(nstp):
        
        iwt,iwt_layers = get_water_table(i_array,bottoms)
        
        i_localzs = water_table[:,:,:,istp] - bottoms   # This is the residual - note that it will compute for each layer but will only be true for the layer in which the water table occurs
        local_zs[:,:,:,istp] = i_localzs * (1./5.)
        local_zs[local_zs < 0] = 0
        local_zs[local_zs > 1] = 0
        global_zs[:,:,istp] = iarray

    return wt_layers,local_zs,global_zs

def depth_to_water(landsurface_fin,heads_fin,ibound_fin):
    '''Reads the landsurface elevations and the heads array and returns
    a 2D array of the depth to the water table.  Current version only
    works with a confined model such that heads in layer 1 = water table
    elevation.'''
    
    water_table = read_binaries.read_heads(heads_fin,nrow,ncol,nlay)[:,:,0]
    ibound_flags = np.genfromtxt(ibound_fin)[:,:]
    landsurface = np.genfromtxt(landsurface_fin)
    water_table_depth = np.zeros((nrow,ncol))
    water_table_depth[:,:] = np.nan
    
    for i in range(nrow):
        for j in range(ncol):
            
            if (ibound_flags[i,j] > 0):
                
                water_table_depth[i,j] = np.max([0,(landsurface[i,j] - water_table[i,j])])

    return water_table_depth

def create_retardation_factors(retardation_multiplier,landsurface_fin,heads_fin,ibound_fin):
    '''Returns a 2D array of retardation factors for the water table
    cells as a linear function of the depth to the water table.'''

    #water_table_depth = depth_to_water(landsurface_fin,heads_fin,ibound_fin)
    
    retardation_factors = np.ones((nrow,ncol)) * retardation_multiplier #water_table_depth * retardation_multiplier
    
    return retardation_factors
    
def read_pvl(pvl_file):
    '''Reads the .pvl file to a dictionary with key=paramname,value=paramvalue.'''
    
    pvl_dict = dict()
    
    with open(pvl_file,'r') as fin:

        next(fin)   # Skip the first line        
        for iline in fin:
            iparamname,iparamval = iline.split()
            pvl_dict[iparamname] = float(iparamval)
            
        return pvl_dict

def build_recharge_map(recharge_flags_fin,recharge_rates_fout,pvl_file,rch_file):
    '''Reads the appropriate files to build a 2D array of the recharge multipliers
    applied to the recharge rate. This function should be cleaned up to properly
    read and store .rch package parameters. (Add a function that reads the .rch
    file). NOTE: Requires that recharge_map already exist according to user-specs.'''
    
    rch_flags = np.genfromtxt(recharge_flags_fin)
    rch_rates = np.ones((nrow,ncol))
    
    pvl_dict = read_pvl(pvl_file)
    
    # Read the recharge parameter names from the .rch file
    # and create a dictionary where key=paramflag from recharge map,value=paramname
    rch_flag_dict = dict()
    
    with open(rch_file,'r') as fin:
    
        for iline in fin:
            
            if (iline.split()[0] in pvl_dict) and (len(iline.split()) > 1):
                iparamname = iline.split()[0]
                iparamflag = int(next(fin).split()[2])
                rch_flag_dict[iparamflag] = iparamname
                
        fin.close()
                
    # Translate the recharge map array (with integer flags) to a recharge
    # rate array
    for i in range(nrow):
        for j in range(ncol):
            
            iflag = int(rch_flags[i,j])
            iparamname = rch_flag_dict[iflag]            
            iparamval = pvl_dict[iparamname]
            rch_rates[i,j] = iparamval
    
    np.savetxt(recharge_rates_fout,rch_rates,fmt='%6.3e')

    return                

def get_hydraulics(heads_file,hob_dict,sim_obs_heads_file,landsurface_elev_fin,modeltop_file):
    '''Returns a dictionary with key=hob_name and value=(unweighted heads residual,
    water table depth below the land surface.'''

    hydraulics_dict = dict()
    
    heads = read_binaries.read_heads(heads_file,nrow,ncol,nlay,nper=1)
    
    for ihob in hob_dict:
        
        irow,icol = hob_dict[ihob][1],hob_dict[ihob][2] # Note that this is MODFLOW, not Python, indexing
        
        with open(sim_obs_heads_file,'r') as fin:
            
            for iline in fin.readlines():
                
                if (iline.split()[2] == ihob):
                    
                    isim,iobs = float(iline.split()[0]),float(iline.split()[1])
                    try:                    
                        iresid = (iobs - isim)/iobs
                    except:
                        iresid = 999
                    
        landsurface = np.genfromtxt(landsurface_elev_fin)
        model_top = np.genfromtxt(modeltop_file)
        
        iwatertable = heads[irow-1,icol-1,0]
        iland = landsurface[irow-1,icol-1]
        itop = model_top[irow-1,icol-1]
               
        idelta = iland - itop
        ivadose = iland - iwatertable
        isurplus = iwatertable - itop
        
        hydraulics_dict[ihob] = (iresid,idelta,ivadose,isurplus)
        
    return hydraulics_dict

def get_vadose_thickness(heads_file,landsurface_elev_fin):
    '''Returns the thickness of the simulated unsaturated zone. NOTE: this is a
    2D array of the distance between the landsurface and the simulated water table,
    NOT the distance between the landsurface and the model top (which, for 
    a confined model, should approximate the water table).'''
    
    heads = read_binaries.read_heads(heads_file,nrow,ncol,nlay,nper=1)
    watertable = np.squeeze(heads[:,:,0])
    landsurface = np.genfromtxt(landsurface_elev_fin)   
    
    return landsurface - watertable
        
def get_geozone(zones_3d,iloc_tuple,unit_labels):
    
    irow,icol,ilay = iloc_tuple
    izone = zones_3d[irow,icol,ilay]

    return unit_labels[int(izone-1)] 

def create_layer_3D(nrow,ncol,nlay,zonelay_fout):
    
    zones_3d = np.zeros((nrow,ncol,nlay))

    for k in range(nlay):
        
        izone = np.genfromtxt(zonelay_fout + str(k+1))
        zones_3d[:,:,k] = izone
        
    return zones_3d

def threeD_to_twoD(threeD_array):
    
    [nlay,nrow,ncol] = np.shape(threeD_array)
    
    twoD_array = np.zeros((nrow * nlay,ncol))
    
    for k in range(nlay):
        
        start_row = nrow * k
        stop_row  = nrow * (k + 1)
        twoD_array[start_row:stop_row,:] = threeD_array[k,:,:]
    
    return twoD_array

def twoD_to_threeD(twoD_array,(nrow,ncol)):
    
    twoD_rows = np.shape(twoD_array)[0]
    nlay = int(twoD_rows/nrow)
    
    threeD_array = np.zeros((nlay,nrow,ncol))
    
    for k in range(nlay):
        
        start_row = nrow * k
        stop_row  = nrow * (k + 1)
            
        threeD_array[k,:,:] = twoD_array[start_row:stop_row,:]
    
    return threeD_array
    
def thksat(model_top_fin,layer_bottoms_fin,heads,min_head,max_head):

    '''Creates a boolean 3D array that locates the cells in
    which the water table occurs and calculates the saturated thickness for
    each water table cell.
    
    model_top   = 2D ASCII array of model top elevations for all surface cells
    layer_bottoms  = 3D ASCII array of bottom elevations for all cells
    heads_object = 3D NUMPY array of heads
    
    Create the heads 3D array using the binaryfile module. Note that the indexing
    convention used by the binaryfile model is [nlay,nrow,ncol]'''    
    
    [nlay,nrow,ncol] = np.shape(heads)        
    
    model_top = np.genfromtxt(model_top_fin)    
    
    layer_bottoms = np.genfromtxt(layer_bottoms_fin)
    layer_bottoms = twoD_to_threeD(layer_bottoms,(nrow,ncol))
    
    watertable = np.zeros((nrow,ncol,nlay),dtype=bool)
    thksat = np.zeros((nrow,ncol))
    
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                
                head = heads[k,i,j]
                bottom = layer_bottoms[i,j,k]                
                
                if (k == 0):
                    top = model_top[i,j]
                    
                else:
                    top = layer_bottoms[i,j,k-1]
                
                if (head < top) and (head > bottom):
                    watertable[i,j,k] = True
                    thksat[i,j] = head - bottom
    
    return watertable,thksat

def lateral_connections(watertable):
    
    '''Determines the number of cells to which each watertable
    cell is connected (i.e., the number of adjacent cells that are NOT dry).
    
    watertable = 3D boolean numpy array'''
    
    [nrow,ncol,nlay] = np.shape(watertable)
    connections = np.zeros((nrow,ncol,nlay))
    
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                
                if (watertable[i,j,k] == True):
                    
                    adj_rows = [i,i,i-1,i+1]
                    adj_cols = [j-1,j+1,j,j]
                    count = 0                    
                    
                    for check in range(len(adj_rows)):
                        
                        check_row = adj_rows[check]
                        check_col = adj_cols[check]
                        
                        if (watertable[check_row,check_col,k] == True):
                            count += 1
                    
                    connections[i,j,k] = count
    
    return connections

def convert_flowpy_3d(flowpy_3D_array):
    '''Converts flowpy 3d array (nlay,nrow,ncol) to 3d array suitable for
    my other functions (nrow,ncol,nlay)'''
    
    [nlay,nrow,ncol] = np.shape(flowpy_3D_array)
    
    my_3D_array = np.zeros((nrow,ncol,nlay))
    
    for k in range(nlay):
        
        my_3D_array[:,:,k] = flowpy_3D_array[k,:,:]
    
    return my_3D_array

def convert_binary_heads(head_fin,head_fout,kstp,kper):
    '''Reads a binary head file and writes to ASCII.'''

    heads = bf.HeadFile(head_fin)
    heads = heads.get_data(kstp=kstp,kper=kper)
    heads = convert_flowpy_3d(heads)
    heads = threeD_to_twoD(heads)
    
    np.savetxt(head_fout,heads,fmt='%10.2f')
    
    return
  
def find_active_bc(bc_dict,(nrow,ncol)):
    '''Converts the ordered dictionary of a boundary condition
    (that is created by the call to binaryfile.py) to a 2D array of cells
    for which the boundary condition is active.'''
    
    active_bc = np.empty((nrow,ncol))
    active_bc[:,:] = np.nan

    for bc_cell in bc_dict:
        
        index = np.mod(bc_cell,(nrow*ncol))

        if (bc_dict[bc_cell] != 0):
            bc_row = np.floor(index/ncol) + 1
            bc_col = index - ((bc_row - 1) * ncol)
            active_bc[bc_row-1,bc_col-1] = bc_dict[bc_cell]
            
    return active_bc

def find_active_edge(IBOUND_2D):
    '''Reads the IBOUND array and returns an (nrow,ncol) array with
    all perimeter cells flagged 1, all other cells 0.'''
    
    [nrow,ncol] = np.shape(IBOUND_2D)
    edge_array = np.zeros((nrow,ncol))
    
    for i in range(nrow):
        for j in range(ncol):
            
            adj_rows = [i,i,i-1,i+1,i-1,i-1,i+1,i+1]
            adj_cols = [j-1,j+1,j,j,j-1,j+1,j-1,j+1]
            
            for check in range(len(adj_rows)):
                
                        check_row = adj_rows[check]
                        check_col = adj_cols[check]
                        
                        if ( check_row >= 0 and check_row < nrow
                        and check_col >= 0 and check_col < ncol):
                            
                            if (IBOUND_2D[check_row,check_col] == 0):
                                edge_array[i,j] = 1
    
    return edge_array
    
def read_dis(dis_fin):
    
    with open(dis_fin,'r') as f:

        dim_line = f.readlines()[1]        
        
        (nlay,nrow,ncol) = int(dim_line.split()[0]),int(dim_line.split()[1]),int(dim_line.split()[2])
        
    return nlay,nrow,ncol
       
def bc_to_2D_array(bc_fin):
    '''Returns a 2D array of model cells subject to a boundary condition.
    Current version works with the DRN and GHB packages.'''
    
    bc_array = np.zeros((nrow,ncol))
    bc_array[:,:] = np.nan
    
    with open(bc_fin,'r') as fin:
        
        for iline in fin.readlines():            
            if (len(iline.split()) == 6):
                irow,icol = int(iline.split()[1]),int(iline.split()[2])
                bc_array[irow,icol] = 1

    return bc_array
    
def get_rch_flags(rch_file,ss_rch_dict):
    '''Returns the recharge parameter names and their associated flags from the .rch file
    and returns a dictionary where key=paramname value=rch_flag.'''
    
    rch_flag_dict = dict()
    
    with open(rch_file,'r') as fin:
    
        for iline in fin:
            
            if (iline.split()[0] in ss_rch_dict) and (len(iline.split()) > 1):
                iparamname = iline.split()[0]
                iparamflag = int(next(fin).split()[2])
                rch_flag_dict[iparamname] = iparamflag
                
        fin.close()
        
    return rch_flag_dict

def get_rch_params(pvl_dict):
    '''Reads the recharge parameters and steady state values from the .pvl file
    and returns a dictionary.'''

    pvl_keys = pvl_dict.keys()
    rch_keys = [x for x in pvl_keys if 'Rch' in x]
    rch_dict = {rch_key:pvl_dict[rch_key] for rch_key in rch_keys}  
    
    return rch_dict
    

