# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 19:41:58 2020

@author: smrak@bu.edu
"""

import xarray, os, glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import cartomap.geogmap as gm
import subprocess
import cartopy.crs as ccrs

save = 0
date = '20170821'
wl = 193
galt = 150
root = root = "\\".join(os.getcwd().split('\\')[:-1])
EOFF = glob.glob(root + "\\data\\{}\\SDO193ephem_test\\{}A_{}*km*.nc".format(date,wl,galt))
odir = root + "\\data\\{}\\EOF\\{}_{}\\".format(date, wl, galt)
for f in EOFF:
    EOF = xarray.open_dataset(f)
    t = np.datetime64(EOF.time.values, 's').astype(datetime)
    fig, ax = gm.plotCartoMap(projection='plate', lon0=0,
                              figsize=[8,5],
                              lonlim=[-180,180], latlim=[-90,90],
                              title = '{}, Alt = {} km'.format(t, EOF.alt_km.values),
                              meridians=np.arange(-180,180.1,40), parallels=np.arange(-80,81,20),
                              background_color='grey')
#    lap = abs(ndimage.laplace(EOF.of.values))
#    lap[lap<0.001] = np.nan
#    LAP = ax.contourf(EOF.glon.values, EOF.glat.values, lap.T, cmap='terrain', 
#                      levels = np.linspace(0.008, 0.1, 40),# vmin=0.001, vmax=0.08,
#                      transform=ccrs.PlateCarree())
    OF = ax.pcolormesh(EOF.glon.values, EOF.glat.values, EOF.of.values, cmap='terrain',
                    vmin=0, vmax=1,
                    transform=ccrs.PlateCarree())
    OFC = ax.contour(EOF.glon.values, EOF.glat.values, EOF.of.values, colors='w', 
                    levels=np.arange(0.39, 1.1, 0.1),
                    transform=ccrs.PlateCarree())
#    plt.colorbar()
    ax.clabel(OFC,OFC.levels, inline=True)
    
    posn0 = ax.get_position()
    cax = fig.add_axes([posn0.x0+posn0.width+0.01, posn0.y0, 0.02, posn0.height])
    fig.colorbar(OF, cax=cax, label='EOF',format='%.2f')
#    plt.tight_layout()
    if save:
        if not os.path.exists(odir):
            subprocess.call('mkdir "{}"'.format(odir), timeout=2, shell=True)
        fig.savefig(odir+ "{}.png".format(t.strftime("%Y%m%d_%H%M")))
        plt.close(fig)
    else:
        break