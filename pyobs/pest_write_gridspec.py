# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 10:09:23 2015

@author: wzell
"""
import numpy as np

nrow,ncol,nlay = 107,137,27
dx,dy = 250,250

delr = (np.ones((1,ncol)) * dy)
delc = (np.ones((nrow,1)) * dx).T

rotation = 0
easting = 0
spectype = 1

gridspec_file = 'C:\Delmarva\UC_250\NWT_SS_Unconfined\PEST\PEST_Input\chester.gridspec'
bottoms_root = 'C:\Delmarva\UC_250\NWT_SS_Unconfined\Framework\Lay'
modeltop_file = 'C:\Delmarva\UC_250\NWT_SS_Unconfined\Framework\chester_modeltop.twice'

with open(gridspec_file,'w') as fout:
    
    fout.write('%10i%10i%10i\n' %(nrow,ncol,nlay))
    fout.write('%15.7e%15.7e%10.2f\n' %(easting,nrow*dy,rotation))
    np.savetxt(fout,delr,fmt='%10.2f')
    np.savetxt(fout,delc,fmt='%10.2f')
    fout.write(str(spectype) + '\n')
    fout.write(modeltop_file + '\n')
    
    for i in range(nlay):
        fout.write(bottoms_root + str(i+1) + '.bottoms\n')