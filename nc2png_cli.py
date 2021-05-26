# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 19:41:58 2020

@author: smrak@bu.edu
"""

import xarray, os, glob
import numpy as np
import matplotlib.pyplot as plt
import cartomap.geogmap as gm
import subprocess
import cartopy.crs as ccrs
from scipy import ndimage
from dateutil import parser
from argparse import ArgumentParser

projection = 'ortographic'
clw = 0.5
mlat_levels = [-65, -75, 60, 75]
tlim = None

def main(idir=None, wl=None, alt_km=None, odir=None, tlim=None,
         laplace=0,auroral_oval=0,clabel=0,
         lat0 = None, lon0 = None, save=1):
    if odir is None:
        if laplace:
            odir = os.path.join(idir, "{}_{}_lap\\".format(wl, alt_km))
        else:
            odir = os.path.join(idir, "{}_{}\\".format(wl, alt_km))
    if not os.path.exists(odir):
        subprocess.call('mkdir "{}"'.format(odir), timeout=2, shell=True)
    assert os.path.exists(odir)
    
    EOFF = np.array(glob.glob(idir + "*_*{}km_*{}*.nc".format(alt_km,wl)))
    f_times = np.array([parser.parse(os.path.split(f)[1][:14]) for f in EOFF])
    if tlim is None:
        idt = np.ones(f_times.size,dtype=bool)
    else:
        tlim = [parser.parse(tlim[0]), parser.parse(tlim[1])]
        assert len(tlim) == 2
        idt = (f_times >= tlim[0]) & (f_times <= tlim[1])
    #
    F = EOFF[idt]
    #
    for i,f in enumerate(F):
        t = f_times[idt][i]
        save_fn = odir + "{}_{}_{}.png".format(t.strftime("%Y%m%d_%H%M%S"), alt_km, wl)
        EOF = xarray.open_dataset(f)
        fig, ax = gm.plotCartoMap(projection=projection, lon0=lon0, lat0=lat0,
                                  figsize=[8,5],
                                  lonlim=[-180,180], latlim=[-90,90],
                                  title = '{}, Alt = {} km'.format(t, EOF.alt_km.values),
                                  meridians=np.arange(-180,180.1,40), parallels=np.arange(-80,81,20),
                                  background_color='grey',
                                  apex=True, mlat_levels=mlat_levels, mlat_colors='b',
                                  mlat_labels=0)
            
        OF = ax.pcolormesh(EOF.glon.values, EOF.glat.values, EOF.of.values, cmap='gray',
                               vmin=0, vmax=1,
                               transform=ccrs.PlateCarree())
    #        
        if laplace:
            lap = abs(ndimage.laplace(EOF.of.values))
            lap[lap<0.001] = np.nan
            ax.contour(EOF.glon.values, EOF.glat.values, lap, cmap='terrain', 
                              transform=ccrs.PlateCarree())        
    
        try:
            OFC = ax.contour(EOF.glon.values, EOF.glat.values, EOF.of.values, colors='r', linewidths=clw, 
                             levels=np.linspace(0.0, 1.0, 25),
                                     transform=ccrs.PlateCarree())
            if clabel:
                ax.clabel(OFC, OFC.levels, inline=True)
        except:
            pass
        posn0 = ax.get_position()
        cax = fig.add_axes([posn0.x0+posn0.width+0.01, posn0.y0, 0.02, posn0.height])
        fig.colorbar(OF, cax=cax, label='EOF',format='%.2f')
        EOF.close()
        
        if save:
            fig.savefig(save_fn)
            plt.close(fig)
        else:
            break
    
if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('folder', help='input folder', type=str)
    p.add_argument('--wl', help='eclipse wavelength. If uniform type geo', default='geo')
    p.add_argument('--odir', help='Save images to:', default=None)
    p.add_argument('--altkm', help='altitude in km', default = 150, type=int)
    p.add_argument('--tlim', help='time limit, 2 args "start" "stop"', nargs=2, default=None, type=str)
    
    p.add_argument('--lon0', help='Map center: longitude default = 0', type=int, default=0)
    p.add_argument('--lat0', help='Map center: Latitude default = 0', type=int, default=0)
    
    p.add_argument('--laplace', help='Plot Laplace?', action='store_true')
    p.add_argument('--oval', help='Plot auroral oval?', action='store_true')
    p.add_argument('--clabel', help='Render countor labels?', action='store_true')
    
    P = p.parse_args()
    main(idir=P.folder, wl=P.wl, alt_km=P.altkm, odir=P.odir, tlim=P.tlim,
         laplace=P.laplace, auroral_oval=P.oval, clabel=P.clabel,
         lat0=P.lat0, lon0=P.lon0, save=1)
    