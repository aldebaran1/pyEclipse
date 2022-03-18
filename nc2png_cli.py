# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 19:41:58 2020

@author: smrak@bu.edu
"""

import xarray, os, glob, datetime
import numpy as np
import matplotlib.pyplot as plt
import cartomap.geogmap as gm
import subprocess
import cartopy.crs as ccrs
from scipy import ndimage
from dateutil import parser
from argparse import ArgumentParser

projection = 'ortographic'
mlat_levels = [-65, -75, 60, 75]

def main(idir=None, wl=None, alt_km=None, odir=None, tlim=None,
         laplace=0,auroral_oval=0,clabel=0, srad=1.0,
         lat0 = None, lon0 = None, save=1):
    cmax = 10 if clabel else 25
    clw = 1 if clabel else 0.5
    
    if auroral_oval:
        apx=1
    else:
        apx = 0
        
    if odir is None:
        if wl != 'geo':
            if laplace:
                odir = os.path.join(idir, "{}_{}_lap\\".format(wl, alt_km))
            else:
                odir = os.path.join(idir, "{}_{}\\".format(wl, alt_km))
        else:
            if laplace:
                odir = os.path.join(idir, "{}_{}_{}_lap\\".format(wl, srad, alt_km))
            else:
                odir = os.path.join(idir, "{}_{}_{}\\".format(wl, srad, alt_km))
    if not os.path.exists(odir):
        subprocess.call('mkdir "{}"'.format(odir), timeout=2, shell=True)
    assert os.path.exists(odir)
    
    if wl != 'geo':
        EOFF = np.array(sorted(glob.glob(idir + "*_{}km_{}*.nc".format(alt_km, wl))))
    else:
        EOFF = np.array(sorted(glob.glob(idir + "*_{}km_{}*_{}.nc".format(alt_km, wl, srad))))
    f_times = []
    for f in EOFF:
        try:
            f_times.append(parser.parse(os.path.split(f)[1][:14]))
        except:
            pass
    f_times = np.array(f_times)
    assert (f_times.size > 0)
    
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
        if wl != 'geo':
            try:
                t_sdo = EOF.time_sdo.values.astype('datetime64[s]').astype(datetime.datetime).strftime("%Y-%m-%d %H:%M")
            except: 
                t_sdo = EOF.time_image.values.astype('datetime64[s]').astype(datetime.datetime).strftime("%Y-%m-%d %H:%M")
            title = '{}, Alt = {} km\nT-Image = {}'.format(t, EOF.alt_km.values, t_sdo)
        else:
            title = '{}, Alt = {} km'.format(t, EOF.alt_km.values)
        fig, ax = gm.plotCartoMap(projection=projection, lon0=lon0, lat0=lat0,
                                  figsize=[8,5],
                                  lonlim=[-180,180], latlim=[-90,90],
                                  title = title,
                                  meridians=np.arange(-180,180.1,40), parallels=np.arange(-80,81,20),
                                  background_color='grey',
                                  apex=apx, mlat_levels=mlat_levels, mlat_colors='b',
                                  mlat_labels=0)

        sza = EOF.sza.values
        night = np.copy(sza)
        inight = (sza > 90)
        if laplace:
            lap = abs(ndimage.laplace(EOF.of.values))
            lap[lap<0.001] = np.nan
            lapm = 0.1
            lap[lap>lapm] = lapm
            night[inight] = 0
            night[~inight] = 1
            OF = ax.pcolormesh(EOF.glon.values, EOF.glat.values, EOF.of.values, cmap='gray',
                               vmin=0, vmax=1,
                               transform=ccrs.PlateCarree())
            
            try:
                ax.contour(EOF.glon.values, EOF.glat.values, lap, cmap='nipy_spectral', 
                            levels=np.linspace(0.005, lapm, 20),
                            transform=ccrs.PlateCarree())
            except:
                pass
            if clabel:
                try:
                    OFC = ax.contour(EOF.glon.values, EOF.glat.values, EOF.of.values, 
                                     cmap='jet', linewidths=clw, 
                                     levels=np.linspace(0.0, 1.0, cmax),
                                     transform=ccrs.PlateCarree())
                    if clabel:
                        ax.clabel(OFC, OFC.levels, inline=True)
                except:
                    pass
        else:
            OF = ax.pcolormesh(EOF.glon.values, EOF.glat.values, EOF.of.values, cmap='gray',
                               vmin=0, vmax=1,
                               transform=ccrs.PlateCarree())
            try:
                Zcont = EOF.of.values
                Zcont[inight] = np.nan
                OFC = ax.contour(EOF.glon.values, EOF.glat.values, Zcont, 
                                 colors='r', linewidths=clw, 
                                 levels=np.linspace(np.nanmin(Zcont), 1.0, cmax),
                                 transform=ccrs.PlateCarree())
                if clabel:
                    ax.clabel(OFC, OFC.levels, inline=True, fmt='%.2f')
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
    p.add_argument('--srad', help='altitude in km', default = 1.0, type=float)
    
    p.add_argument('--lon0', help='Map center: longitude default = 0', type=int, default=0)
    p.add_argument('--lat0', help='Map center: Latitude default = 0', type=int, default=0)
    
    p.add_argument('--laplace', help='Plot Laplace?', action='store_true')
    p.add_argument('--oval', help='Plot auroral oval?', action='store_true')
    p.add_argument('--clabel', help='Render countor labels?', action='store_true')
    
    P = p.parse_args()
    main(idir=P.folder, wl=P.wl, alt_km=P.altkm, odir=P.odir, tlim=P.tlim,
         laplace=P.laplace, auroral_oval=P.oval, clabel=P.clabel, srad=P.srad,
         lat0=P.lat0, lon0=P.lon0, save=1)
    
