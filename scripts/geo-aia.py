# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 13:26:19 2021

@author: smrak@bu.edu
"""

import xarray, glob, os
import subprocess
from cartomap import geogmap as gm
import cartopy.crs as ccrs
import numpy as np
from dateutil import parser
import matplotlib.pyplot as plt

save = 1
alt_km = 150
scale = 1.1
wl = 193
eta = 1

date = '20210610'
root = 'G:\\My Drive\\eclipse\\mask\\{}\\'.format(date)
root = 'C:\\Users\\smrak\\OneDrive - UCB-O365\\Projects\\eclipse_LASP\\{}\\'.format(date)
odir = root + 'diff\\'

Fgeo = np.array(sorted(glob.glob(root + "geo\\*_{}km_{}_{}.nc".format(alt_km, 'geo', scale))))
geo_times = np.array([parser.parse(os.path.split(f)[1][:14]) for f in Fgeo])

Faia = np.array(sorted(glob.glob(root + "aia\\*_{}km_{}*_{}.nc".format(alt_km, wl, eta))))
aia_times = np.array([parser.parse(os.path.split(f)[1][:14]) for f in Faia])

for i,f in enumerate(Faia):
    EOFaia= xarray.open_dataset(f).of.values
    glon = xarray.open_dataset(f).glon.values
    glat = xarray.open_dataset(f).glat.values
    idt = abs(geo_times-aia_times[i]).argmin()
    if abs(geo_times-aia_times[i]).min().seconds > 100:
        continue
    EOFgeo = xarray.open_dataset(Fgeo[idt]).of.values
    
    fig, ax = gm.plotCartoMap(projection='ortographic', lon0=-50, lat0=40,
                              figsize=[8,5],
                              lonlim=[-180,180], latlim=[-90,90],
                              title = aia_times[i],
                              meridians=np.arange(-180,180.1,40), parallels=np.arange(-80,81,20),
)#                              background_color='grey',)
#                              apex=0, mlat_levels=None, mlat_colors='b',
#                              mlat_labels=0)
    
    sza = xarray.open_dataset(f).sza.values
    night = np.copy(sza)
    
    aia = EOFaia
    aia[night>=90] = 0
    
    geo = EOFgeo
    geo[night>=90] = 0
    
    diff = geo-aia
    
    imd = ax.pcolormesh(glon, glat, diff, cmap='nipy_spectral',
                  vmin=-0.2, vmax=0.2,
                  transform=ccrs.PlateCarree())
#    ax.contour(glon, glat, aia, levels=np.linspace(0.2,1,15), colors='w',linewidths=0.75,
#               transform=ccrs.PlateCarree())
    
    posn0 = ax.get_position()
    cax = fig.add_axes([posn0.x0+posn0.width+0.01, posn0.y0, 0.02, posn0.height])
    fig.colorbar(imd, cax=cax, label='GEO-AIA',format='%.2f')
    
    if save:
        if not os.path.exists(odir):
            subprocess.call('mkdir "{}"'.format(odir), shell=True)
        fig.savefig(odir + 'geo{}-aia{}_{}.png'.format(scale,wl,aia_times[i].strftime("%Y%m%d%H%M")))
        plt.close(fig)
    else:
        break