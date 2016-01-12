# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 17:04:20 2015

@author: wzell
"""

import numpy as np
import matplotlib.pyplot as plt


fin = 'PEST\PEST_Input\HK.out'
threshold = 1e20

with open(fin,'r') as fin:
    
    ncol,nrow,nlay = [int(i) for i in fin.readline().split()]
    conduct = np.empty((nrow,ncol,nlay))

    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
            
                conduct[i,j,k] = float(fin.readline())

conduct[conduct > threshold] = np.nan            
plt.imshow(conduct[:,:,0],interpolation='nearest')
        