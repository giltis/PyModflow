# -*- coding: utf-8 -*-
"""
Executes the PEST FAC2REAL3D utility.

FAC2REAL3D reads the kriging factors calculated by PPK2FAC3D and interpolates
parameter values listed in a 3D pilot points file to MODFLOW cell centers.

RETURNS: nlay files where each file contains a (nrow x ncol) array of parameter values.

FILE DEPENDENCIES:
-- settings.fig
-- Grid specification file that includes 3D information (i.e., nlay and
layer elevation filenames).
-- 3D pilot points file
-- structure file (which describes the geostatistical structure)

Utility user prompts:
1. Enter the name of the interpolation factor file:
2. Is this a formatted or unformatted file? [f/u]:
3. Enter name of 3D pilot points file [or press ENTER to confirm
that the pilot points file used by PPK2FAC3D is the correct pilot points file]:
4. Enter lower interpolation limit:
5. Enter upper interpolation limit:
6. Write outputs to single 3D table or multiple 2D real array files? [s/m]:

The following prompt assumes that 'm' was entered for 6. (See the documentation
for the 's' option.)
7. Enter filename base for output real array files:

8. Enter value for elements to which no interpolation takes place (e.g., cells
in zones with no spatial interpolation or cells that are further removed from
the nearest pilot point than the maximum search radius specified by PPK2FAC3D): 

Created on Fri Oct 09 10:51:45 2015
@author: wzell
"""

