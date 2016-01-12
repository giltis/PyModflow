# -*- coding: utf-8 -*-
"""
Created on Mon Mar 09 11:55:53 2015

@author: Wesley Zell
"""

import numpy as np
from pygrid import grid_util as gu
import itertools
from scipy.spatial.distance import cdist

# ============================================    
# === START REGISTRATION VOLUME FUNCTIONS ====
# ============================================

def get_disc_particles(iradius,iparticles_per_cell,idelr,idelc):
    '''Returns the calculated number of particles that should be distributed on
    a vogel disc of a user-specified radius.'''
        
    area_disc = 3.14159 * ((iradius) ** 2)
    area_cell = idelr * idelc
    
    particles_per_disc = np.max([(int(np.ceil(iparticles_per_cell * (area_disc/area_cell)))),iparticles_per_cell])
    
    return particles_per_disc

def vogel_disc((x0,y0),npoints,radius):
    '''Generates an evenly distributed set of points on a disc. Returns a list of
    x,y coordinates of points relative to user-specified center (x0,y0).'''

    radii = radius * (np.sqrt(np.arange(npoints) / float(npoints)))
     
    golden_angle = np.pi * (3 - np.sqrt(5))
    thetas = golden_angle * np.arange(npoints)
     
    points = np.zeros((npoints, 2))
    points[:,0] = np.cos(thetas)
    points[:,1] = np.sin(thetas)
    points *= radii.reshape((npoints, 1))
    
    # Translate the points per the specified center
    for point in points:
        point[0] = point[0] + x0
        point[1] = point[1] + y0
    
    return points

def dist_bt_points(points):
    '''Computes an approximate average closest distance between points
    distributed on a vogel disc.'''
    
    # Generate a symmetric matrix that includes all distances between all points
    all_dist = cdist(points,points)
    all_dist = np.ma.masked_where(all_dist == 0,all_dist)   # Don't include distance to self
    
    # Average over the minimum distance to another point for each point
    nearest_neighbors = np.min(all_dist,axis=1)  
    average_min_distance = np.mean(nearest_neighbors)
    
    return average_min_distance

def points_on_plane((x0,y0),iradius,ztop,zbottom,nparticles_per_cell,delc,lay_thick,print_specs=False):
    '''Returns a list of tuples (x,y of points) in MODEL GLOBAL coordinate system.
    Note that the nparticles_per_cell is used to approximate point density.'''
    
    # Calculate the point density and scale to size of registration plane
    cell_cross_sectional_area = delc * np.mean(lay_thick)
    point_density = nparticles_per_cell * 1./(cell_cross_sectional_area)

    start_x,stop_x = x0-iradius,x0+iradius
    reg_width = stop_x - start_x
    reg_height = ztop - zbottom
    reg_area = reg_width * reg_height
    
    points_in_plane = np.max((point_density * reg_area,nparticles_per_cell))
    
    print "Points in plane: %0.3e" %(points_in_plane)
    
    # Distribute the points as a grid on the plane. Find grid dimensions as
    # function of registration plane dimensions
    if (reg_width < 1):
        n_grid_cols = 2
    else:
        n_grid_cols = np.ceil(points_in_plane / reg_width)
    
    n_grid_rows = np.ceil(points_in_plane / n_grid_cols)  

    if (print_specs == True):
        print "point_density * reg_area: %0.3e" %(point_density * reg_area)
        print "N Particles per cell: %0.3e" %(nparticles_per_cell)
        print "Point density: %0.3e" %(point_density)
        print "Cell cross sectional area: %0.3e" %(cell_cross_sectional_area)
        print "Radius, top, bottom: %0.3e %0.3e %0.3e" %(iradius,ztop,zbottom)
        print "Area of plane: %0.3e" %(reg_area)
        print "(Width,Height) of plane: (%0.3e,%0.3e)" %(reg_width,reg_height)
        print "Number of (rows,columns) in plane: (%i,%i)" %(n_grid_rows,n_grid_cols)
    
    vert_particle_line = np.linspace(ztop,zbottom,n_grid_rows)
    horiz_particle_line = np.linspace(start_x,stop_x,n_grid_cols)
    
    points = []
    for ix in horiz_particle_line:
        points.append((ix,y0))
        
    return points, vert_particle_line
    
# ===========================================    
# === STOP REGISTRATION VOLUME FUNCTIONS ====
# ===========================================   

def get_well_startloc(start_df,top=None,bottoms=None,midpoints=None,min_particle_separation=1,\
                      reg_cylinder=None,particles_per_cell=None,delr=None,delc=None,lay_thick=None):
    '''Returns a dictionary of particle starting locations for input to MODPATH
    package writing functions. Removes duplicate station LOCATIONS from
    the dataframe that may include multiple observations at a particular station.
    If the specs for a registration volume are provided,
    distributes particles throughout the registration volume.  If no registration
    volume is provided, distributes particles in well screen.  (NOTE: Function
    currently only allows generation of a registration CYLINDER.) Returns PYTHON
    indexing.'''
    
    nlay,nrow,ncol = np.shape(bottoms)
    
    start_df = start_df.reset_index()
    start_df = start_df.drop_duplicates('index').set_index('index')
    start_df.index.names = ['station_nm']
    start_df = start_df.drop(['ObsName','TobDate'],axis=1)
    
    # If the dimensions of a registration cylinder are provided, generate a
    # vogel disc in order to determine the distance between points.  This will
    # be used to assign the number of vertical increments to the registration
    # cylinder.  Otherwise (i.e., for cases in which particles are only generated
    # from the well screen), use the min_particle_separation parameter to vertically
    # discretize the well screen.
    if (reg_cylinder is not None):        
        iradius,dz_top,dz_bottom = reg_cylinder['Radius'],reg_cylinder['Top'],reg_cylinder['Bottom']
        particles_per_disc = get_disc_particles(iradius,particles_per_cell,delr,delc)
        points_xy = vogel_disc((0,0),particles_per_disc,iradius)
        distance_between_points = np.max([dist_bt_points(points_xy),min_particle_separation])
    
    else:
        distance_between_points = min_particle_separation
        dz_top,dz_bottom = 0,0

    # Loop through all the observation locations
    particle_start_dict = dict()
    for iname,irow in start_df.iterrows():
        
        print 'Generating particle starting locations for %s.' %(iname)
        
        # Read the well location from the dataframe 
        x0,y0,ztop,zbottom = irow['ModelX'],irow['ModelY'],irow['GlobZtop'],irow['GlobZbot']
        
        # Update the upper and lower bounds on the starting particle locations
        ztop = ztop + dz_top
        zbottom = zbottom - dz_bottom
        nz_increments = int(np.ceil((ztop - zbottom)/distance_between_points))
                            
        # Distribute the particles throughout the well screen or the associated
        # registration volume. 'vert_particle_line' is a vector of elevations
        # at which the particles are set
        if (reg_cylinder is not None):
            if (nrow == 1): # For 2D case generate a plane rather than a cylinder
                points_xy,vert_particle_line = points_on_plane((x0,y0),iradius,ztop,zbottom,particles_per_cell,delc,lay_thick)
            else:
                points_xy = vogel_disc((x0,y0),particles_per_disc,iradius)
                vert_particle_line = np.linspace(zbottom,ztop,num=nz_increments)
                
        else:   # Only generate particles within the well screen
            points_xy = [(x0,y0)]
            vert_particle_line = np.linspace(zbottom,ztop,num=nz_increments)       
        
        # Convert the global elevations to layer number with localz
        for ivert in vert_particle_line:     
            
            for ipoint in points_xy:
                
                xx,yy = ipoint[0],ipoint[1]
                irow,jcol,local_x,local_y = gu.modelxy_to_localxy(xx,yy,delr,delc,nrow,ncol)
                
                # Do not place particles above the model top or below the model bottom
                if (ivert > top[irow,jcol]) or (ivert < bottoms[-1,irow,jcol]):
                    break
                
                ilay,local_z = gu.globalz_to_localz(irow,jcol,ivert,midpoints,bottoms,lay_thick)
            
                if iname not in particle_start_dict:
                    particle_start_dict[iname] = [(ilay,irow,jcol,local_x,local_y,local_z)]
                else:
                    particle_start_dict[iname].append((ilay,irow,jcol,local_x,local_y,local_z))

    return particle_start_dict
    
def get_landsurface_startloc(zoned_startloc_array,nx,ny,nz=1,part_group_root='IBOUND'):
    '''Returns a dictionary of particle starting locations for input to MODPATH
    package writing functions. CURRENT VERSION ASSUMES THAT PARTICLES START
    IN THE TOP LAYER. Returns PYTHON indexing.'''    
        
    nrow,ncol = np.shape(zoned_startloc_array)
    startloc_list = np.unique(zoned_startloc_array)
    
    # Assumes that cells tagged with '0' are inactive and does not
    # write starting locations for those cells.
    startloc_list = [int(x) for x in startloc_list if x != 0]
    
    particle_start_dict = dict()
    for izone in startloc_list:
        
        igroup = part_group_root + '_' + str(izone)
        particle_start_dict[igroup] = []
        
        for i,j,k in itertools.product(range(nrow),range(ncol),[0]):
            if (zoned_startloc_array[i,j] == izone):
                
                xs,ys,zs = np.linspace(0,1,nx+2)[1:-1],np.linspace(0,1,ny+2)[1:-1],np.linspace(0,1,nz+2)[1:-1]
                
                for local_x,local_y,local_z in itertools.product(xs,ys,zs):
                    
                    particle_start_dict[igroup].append((k,i,j,local_x,local_y,local_z))
                    
    return particle_start_dict
    
def read_mpstartloc(startloc_file,model_top,cell_bottoms,dx,dy):
    '''Reads the MODPATH starting locations file to a dictionary with key=groupname
    and values=(row,col,lay), where row=irow+(1-localy), col=icol+localx, and 
    lay=ilay+localz.  MAY NEED UPDATING.'''
    
    particle_dict = dict()    
    
    with open(startloc_file,'r') as fin:
        
        lines = fin.readlines()
        
        igroupcount = int(lines[1].split()[0])
        startline = 2        
        
        for igroup in range(igroupcount):
            
            igroupname = lines[startline].split()[0]
            ilocationcount = int(lines[startline+1].split()[0])
            
            for i in range(startline+2,startline+2+ilocationcount):
                
                iline = lines[i]           
                ilay,irow,icol = int(iline.split()[1]),int(iline.split()[2]),int(iline.split()[3])
                localx,localy,localz = float(iline.split()[4]),float(iline.split()[5]),float(iline.split()[6])
                
                globalx,globaly = gu.localxy_to_globalxy(irow,icol,localx,localy,dx,dy)
                globalz = gu.localz_to_globalz(irow,icol,ilay,localz,model_top,cell_bottoms)
                
                if igroupname not in particle_dict:
                    particle_dict[igroupname] = [(globalx,globaly,globalz)]
                else:
                    particle_dict[igroupname].append((globalx,globaly,globalz))
                
            startline += ilocationcount + 2
                  
    return particle_dict

def plot_mpstartloc(startloc_file,model_top,cell_bottoms,dx,dy):
    '''Reads the MODPATH starting locations file and generates a 3D plot
    of particle starting locations.  MAY NEED UPDATING.'''
    
    particle_dict = read_mpstartloc(startloc_file,model_top,cell_bottoms,dx,dy)
    
    for igroup in particle_dict:
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        xs,ys,zs = [],[],[]        
        
        for iparticle in particle_dict[igroup]:
            ys.append(iparticle[0])
            xs.append(iparticle[1])
            zs.append(iparticle[2])
                    
        ax.scatter(xs,ys,zs)

        plt.show()

    return