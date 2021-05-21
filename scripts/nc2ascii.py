# -*- coding: utf-8 -*-
"""
Created on Fri May 21 12:48:24 2021

@author: smrak@bu.edu
"""

import numpy as np
import xarray, glob, os
import matplotlib.pyplot as plt

idir ='G:\\My Drive\\eclipse\\mask\\20170821\\aia\\'
wl  = 193
altkm = 150
EOFF = np.array(glob.glob(idir + "*{}*.nc".format(altkm)))

for f in EOFF:
    EOF = xarray.open_dataset(f)
    mask = EOF.of.values.ravel()
    ff = open(os.path.splitext(f)[0]+'.txt', 'w')
    np.savetxt(ff, mask, fmt='%.4f')
    ff.close()
    EOF.close()
#    break

#a = np.genfromtxt(os.path.splitext(EOFF[40])[0]+'.txt', dtype=float).reshape(180,360)
#plt.imshow(a, origin='lower')