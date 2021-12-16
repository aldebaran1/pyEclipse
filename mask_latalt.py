# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 13:01:16 2021

@author: smrak@bu.edu
"""
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from eclipse import utils, eio
from mpl_toolkits.axes_grid1 import make_axes_locatable

t0 = datetime(2021,6,10,9,0,0)
t1 = datetime(2021,6,10,10,0,0)
T = datetime(2021,12,4,8,15,0)
#tsdo = datetime(2021,12,4,7,30,0)
tsdo = None
T = datetime(2017,8,21,18,0)
dm = 50
ds = 0
srad_fact = 1.15

glon_center = [160, -20]
glon_center = [-100]

wl = 94

animation = 0
save = 0
folder = 'G:\\My Drive\\eclipse\\mask\\{}\\latalt\\lon{}\\'.format(T.strftime("%Y%m%d"), glon_center)

glat = [np.arange(-40, -90, -1), np.arange(-90, -50.1, 1)]
glat = [np.arange(0, 90.1, 1)]
if len(glat) > 1:
    for i in range(len(glat)):
        if i == 0:
            xglat = glat[i]
        else:
            xglat = np.hstack((xglat, glat[i]))
else:
    xglat = glat[0]
ghgt = np.arange(0, 400.1e3, 50e3)

if animation:
    times = []
    while t0 <= t1:
        times.append(t0)
        t0 += timedelta(minutes = dm, seconds = ds)
    times = np.array(times)
else:
    times = np.array([T])
    
if wl != 'geo':
    aiafolder = 'G:\\My Drive\\eclipse\\sdoaia\\'
    if tsdo is None:
        tsdo = times[0]
    SDO = eio.sunaia(folder=aiafolder, wl=wl, time=tsdo)

latsize = 0
for i in range(len(glat)):
    latsize += glat[i].size
OF = np.nan*np.ones((latsize, ghgt.size, times.size))
for tt, T in enumerate(times):
    print ("Processing {} {}/{}".format(T, tt+1, times.size))
    for i in range(len(glat)):
        if wl == 'geo':
            if i == 0:
                OF[0:glat[i].size,:,tt] = a = utils.mask_latalt_geo(T, glon0=glon_center[i], glat=glat[i], ghgt=ghgt, srad_fact=srad_fact)
            else:
                OF[glat[0].size:,:,tt] = a = utils.mask_latalt_geo(T, glon0=glon_center[i], glat=glat[i], ghgt=ghgt, srad_fact=srad_fact)
        else:
            if i == 0:
                OF[0:glat[i].size,:,tt] = a = utils.mask_latalt_sdo(SDO, T, glon=glon_center[i], glat=glat[i], ghgt=ghgt, wl=wl)
            else:
                OF[glat[0].size:,:,tt] = a = utils.mask_latalt_sdo(SDO, T, glon=glon_center[i], glat=glat[i], ghgt=ghgt, wl=wl)

plt.figure(figsize=[8,5])
plt.title(f'Time:{T}, glon:{glon_center}')
im = plt.pcolormesh(np.arange(latsize), ghgt/1e3, np.squeeze(OF).T, cmap='nipy_spectral', 
                    shading='gouraud')
plt.contour(np.arange(latsize), ghgt/1e3, np.squeeze(OF).T, levels=np.linspace(0.1,1,40), 
            colors='w', linewidths=0.5)
plt.ylabel('Altitutde [km]')
plt.xlabel('Geographic latitude')
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(mappable=im, cax=cax)
cbar.set_label('EOF')
plt.draw()

ticks = ax.get_xticklabels()
xticklabels = []
for tick in ticks:
    try:
        xticklabels.append(int(xglat[int(tick.get_text())]))
    except:
        pass
ax.set_xticklabels(xticklabels)
plt.draw()
