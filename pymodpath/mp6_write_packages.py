# -*- coding: utf-8 -*-
"""
Created on Mon Mar 03 13:52:07 2014

@author: Wesley Zell
"""

from pymodflow.pygrid import grid_util as gu

def write_startloc(particle_dict,startloc_file,mp_startloc_params):
    '''Writes a dictionary of particle starting locations (key=obsname,
    values=(ilay,irow,icol,localx,localy,localz)) to a MODPATH
    starting locations file. One group of particles written for each observation.'''
    
    print '\nWriting start location file for %i particle starting groups' %(len(particle_dict))
    print 'to %s' %(startloc_file)
    print 'Note that some of these observations may be located in the same MODFLOW cell.\n'    

    grid = mp_startloc_params['grid']
    inputstyle = mp_startloc_params['inputstyle']
    releasestarttime = mp_startloc_params['releasestarttime']
    releaseoption = mp_startloc_params['releaseoption']
       
    with open(startloc_file,'w') as fout:
        
        fout.write('%-55i%-20s\n' %(inputstyle,'#Item1: InputStyle'))
        fout.write('%-55i%-20s\n' %((int(len(particle_dict))),'#Item6: GroupCount'))
            
        for igroup in particle_dict:
            
            count = 0
            
            fout.write('%-55s%-20s\n' %(igroup,'#Item7: GroupName'))
            fout.write('%-10i%-6i%-6i%33s%-20s\n' %(len(particle_dict[igroup]),releasestarttime,releaseoption,'','#Item8: LocationCount,ReleaseStartTime,ReleaseOption'))
            
            for iparticle in particle_dict[igroup]:
                
                label = str(igroup) + '-' + str(count)  
                count += 1
                
                ilay,irow,icol = iparticle[0],iparticle[1],iparticle[2]                              
                localx,localy,localz = iparticle[3],iparticle[4],iparticle[5]
                fout.write('%-3i%-5i%-5i%-5i%-7.4f%-7.4f%-7.4f%15s\n' %(grid,ilay+1,irow+1,icol+1,localx,localy,localz,label))

    return

def write_mp6(mp6_file,ibound_files=None,porosity_files=None,hnoflo=None,hdry=None,mp6_budgetlabels=None,laytyp=None):
    '''Writes MODPATH .mp6 file. Assumes that the ibound and porosity files are passed as a list of relative paths
    and writes using 'OPEN/CLOSE' array control.'''
    
    with open(mp6_file,'w') as fout:
        
        fout.write('%-10.1f%-10.1f\n' %(hnoflo,hdry))
        fout.write('%-6i\n' %(len(mp6_budgetlabels)))
        
        for ilabel in mp6_budgetlabels:
            fout.write('%-10s\n%-5i\n' %(ilabel,mp6_budgetlabels[ilabel]['IFACE']))
        
        for ilaytyp in laytyp:
            fout.write('%-3i' %(ilaytyp))
        fout.write('\n')
        
        ilay = 0
        for ifile in ibound_files:
            fout.write('%-15s%-50s%5i%10s%5i Layer %3i IBOUND\n' %('OPEN/CLOSE',ifile,1,'(FREE)',-1,ilay+1))
            ilay += 1
        
        ilay = 0
        for ifile in porosity_files:
            fout.write('%-15s%-50s%5i%10s%5i Layer %3i Porosity\n' %('OPEN/CLOSE',ifile,1,'(FREE)',-1,ilay+1))
            ilay += 1
            
    return

def write_mpnam(mpnam_file,mpnam_dict):
    '''Writes the MODPATH .mpnam file.''' 
        
    with open(mpnam_file,'w') as fout:

        for itype in mpnam_dict:
            
            ifile = mpnam_dict[itype]['File']
            iunit = mpnam_dict[itype]['Unit']
            
            fout.write('%-20s%-10i%-40s\n' %(itype,iunit,ifile))
    
    return

def write_mpsim(mpsim_file=None,mpnam_relpath=None,mplist_relpath=None,mpendpoint_relpath=None,mpstartloc_relpath=None,mppathline_relpath=None,\
                mpsim_options=None,mp_referencetime=1,mp_stopzones=0,zone_files=None,retard_files=None):
    '''Writes the MODPATH .mpsim file. 'mpsim_options' is an ORDERED dictionary. 'zone_files' and 'retard_files' are lists of relative paths.'''

    # If mp_stopzones is a not a list, convert it to a list
    if ((type(mp_stopzones)==list) == False):
        ilist = []
        ilist.append(mp_stopzones)
        mp_stopzones = ilist
    
    with open(mpsim_file,'w') as fout:
        
        # Items 1 and 2: Modpath name and listing files
        fout.write('%s\n' %(mpnam_relpath))
        fout.write('%s\n' %(mplist_relpath))
        # Item 3: Option flags
        for i in mpsim_options:
            fout.write('%3i' %(mpsim_options[i]))
        fout.write('\n')
        # Item 4: Endpoint file
        fout.write('%s\n' %(mpendpoint_relpath))
        # Item 5: Pathline file        
        if (mpsim_options['simtype'] == 2):
            fout.write('%s\n' %(mppathline_relpath))
        # Item 8: Reference time
        fout.write('%3i\n' %(mp_referencetime))
        # Item 22: Particle starting locations file        
        fout.write('%s\n' %(mpstartloc_relpath))
        # Item 30: Stop zone
        for izone in mp_stopzones:
            fout.write('%-5i' %(izone))
        fout.write('\n')
        
        # Item 31: Zone flags
        ilay = 0
        for ifile in zone_files:
            fout.write('%-15s%-50s%5i%10s%5i Layer %3i Zone\n' %('OPEN/CLOSE',ifile,1,'(FREE)',-1,ilay+1))
            ilay += 1

        # Item 32:
        if (mpsim_options['retard'] == 2):
            ilay = 0
            for ifile in retard_files:
                fout.write('%-15s%-50s%5i%10s%5i Layer %3i RetardationFactors\n' %('OPEN/CLOSE',ifile,1,'(FREE)',-1,ilay+1))
                ilay += 1
                
    return