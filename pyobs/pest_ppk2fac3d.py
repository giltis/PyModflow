# -*- coding: utf-8 -*-
"""
Executes the PEST PPK2FAC3D utility.

PPK2FAC3D computes kriging factors that are subsequently read by the FAC2REAL3D
utility that interpolates kriged values to MODFLOW grid.
Reads responses to user prompts from Excel notebook.

It is my current understanding that these factors depend only on spatial
information associated with the model framework (and the specified or derived 
spatial variation of conductivity or other parameter) and therefore need only
be computed a single time.  

FILE DEPENDENCIES:
-- settings.fig
-- Grid specification file that includes 3D information (i.e., nlay and
layer elevation filenames).
-- 3D pilot points file
-- structure file (which describes the geostatistical structure)

Utility user prompts:
1. Enter name of grid specification file:
2. Enter the name of the 3D pilot points file:
3. Enter filename base of layer zonal integer array files:

(Need to confirm that, if provided a filename base, it does not subsequently
prompt for more zone information - e.g., 'Press Enter if entire grid belongs
to single zone.')

4. Enter name of structure file:

5. FOR EACH ZONE WITHIN THE MODEL DOMAIN:
5a. Enter structure name for zone with integer value of 1:
5b. Perform simple or ordinary kriging [s/o]:
5c. Enter search radius in maximum horizontal elongation dirn:
5d. Enter search radius in minimum horizontal elongation dirn:
5e. Enter search radius in vertical dirn:
5f. Enter minimum number of points to use for interpolation:
5g. Enter maximum number of points to use for interpolation:

6. Enter name for interpolation factor file:
7. Is this a formatted or unformatted file? [f/u]:

Created on Fri Oct 09 10:20:25 2015
@author: wzell
"""
