# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:21:50 2021

@author: smrak@bu.edu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 11:47:20 2019

@author: smrak
"""
from datetime import datetime, timedelta
from eclipse import eio, utils
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os, platform, subprocess
from argparse import ArgumentParser
from dateutil import parser

if platform.system() == 'Windows':
    sep = '\\'
else:
    sep = '/'

RE = 6371 #km
#aiares = 4096
#srad_fact = 1.0
#

    
#startend = ['2017-8-21T19:00', '2017-8-21T20:30']
#glat = 10
#glon = -100
#alt_km = 150
#wl = 193

#odir = 'G:\\My Drive\\eclipse\\mask\\20170821\\fov\\{}\\'.format(wl)
cmap = 'gist_stern'

def main(startend =None, wl = 193, alt_km=None, glon=None,glat=None,
         tsdo = None, dt = 10,odir = None, sdodir = None, plot=None,
         srad_fact=1, clim = [0, 1000], cmap = 'gist_stern', animation = False):
    
    global sep, RE
    assert odir is not None
    
    fwddir = '{}E{}N{}H_{}{}'.format(glon, glat, int(alt_km),wl,sep) if glon > 0 else '{}E{}N{}H_{}{}'.format(360+glon, glat, int(alt_km),wl,sep)
    odir = os.path.join(odir, fwddir)
    
    times = utils.get_times(parser.parse(startend[0]), parser.parse(startend[1]), dm=dt)
    ghgt = alt_km*1e3
    
    if tsdo is None:
        tsdo = times[0]
    else:
        tsdo = parser.parse(tsdo)
    assert isinstance(tsdo, datetime)
    
    #SDO Image
    if wl != 'geo':
        if platform.system() == 'Windows':
            sdodir = 'G:\\My Drive\\eclipse\\sdoaia\\'
        if sdodir is None:
            sdodir = input("Type path to the sdoaia directory:\n")
        assert os.path.exists(sdodir)
        try: 
            D = eio.sunaia(folder=sdodir, wl=int(wl), time=tsdo)
            image = D['AIA'+wl].values
#            moon = plt.imshow(image, origin='lower',aspect='equal',
#                           vmin=clim[0],vmax=clim[1],
#                           extent=[-2047,2048,-2047,2048],
#                           cmap=cmap)
#            plt.show()
        except BaseException as e:
            raise (e)
#        exit()
        
    for i,T in enumerate(times):
        sun, moon = utils.objects(T, glon, glat, ghgt)
        horizon = (-np.arccos(RE / (RE + alt_km)) - sun.alt - sun.radius)
        sep = utils.separation(sun.az, sun.alt, moon.az, moon.alt) 
        azm = utils.azimuth(sun.az, sun.alt, moon.az, moon.alt)
        mx0, my0 = utils.rotate(sep, azm, 0.0, 0.0)
        
        if wl == 'geo':
            sr0 = sun.radius
            sr = sr0 * srad_fact
            of = utils.get_EOF(sr, moon.radius, mx0, my0)
            
            fig, ax = plt.subplots(figsize=[5,5])
            ax.set_title('Uniform eclipse, \n LON {}, LAT {}, at {} km\n{}'.format(glon,glat,ghgt/1e3,T))
            ax.set(xlim=(-sr0*2, sr0*2), ylim = (-sr0*2, sr0*2))
            s_im = plt.Circle((0, 0), sr0, color='y')
            m_im = plt.Circle((mx0, my0), moon.radius, color='k')
            ax.text(-sr0*1.9, 1.8*sr0, "Solar radii = {}sr\nEOF={:.2f}".format(srad_fact,of),color='k')
            ax.add_artist(s_im)
            ax.add_artist(m_im)
            ax.set_axis_off()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        
        else:
            wl = int(wl)
            if horizon*D.pixscale >= (image.shape[0]/2-image.shape[0]):
                hmask = utils.horizon_mask(horizon=horizon, selv=sun.alt, 
                                           imsdo=image, 
                                           pixscale=D.pixscale, y0=D.y0)
            else:
                hmask = np.ones_like(image)
            eta = utils.parallactic_angle(sun.az, sun.dec, glat)
                
            
            # Rotation of the moon if ~100x faster than rotation of the Sun for the parallactinc angle
            imrot=  ndimage.rotate(image, np.rad2deg(eta), reshape=False) # --- takes about 3seconds to compute
            
            mmask = utils.moon_mask(imrot.shape[0], mx0*D.pixscale + D.x0, 
                                    my0*D.pixscale + D.y0, 
                                    np.round(moon.radius, 8)*D.pixscale)
            mask = np.multiply(hmask, mmask)
            of =  np.nansum(np.multiply(imrot,mask)) / np.nansum(imrot)
            
            sunclip = np.multiply(imrot, mask)
            
            fig = plt.figure(figsize=(8,8))
            ax = plt.gca()
            ax.set_title('SDO AIA {}A, \n LON {}, LAT {}, at {} km\n{}'.format(wl,glon,glat,ghgt/1e3,T))
            moon = ax.imshow(sunclip, origin='lower',aspect='equal',
                       vmin=clim[0],vmax=clim[1],
                       extent=[-2047,2048,-2047,2048],
                       cmap=cmap)
            ax.set_axis_off()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.text(-2000, 1900, "EOF={:.3f}".format(of), color='w')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            
            cbar = plt.colorbar(mappable=moon, cax=cax)
            cbar.set_label('counts')
        plt.tight_layout()
        
        print ('Image {}/{}, t={}'.format(i+1, times.size, T))
        if not os.path.exists(odir):
            if platform.system() == 'Windows':
                subprocess.call('mkdir "{}"'.format(odir), shell=True, timeout=1)
            else:
                subprocess.call('mkdir -p {}'.format(odir), shell=True, timeout=1)
        plt.savefig(odir + "{}.png".format(T.strftime('%Y%m%d_%H%M%S')), DPI=200)
        if plot:
            plt.show(fig)
        plt.close(fig)
        
    if animation:
        if platform.system() == 'Windows':
            subprocess.call('python animation.py "{}" -f 5'.format(odir), shell=True, timeout=60)
        else:
            subprocess.call('python animation.py {} -f 5'.format(odir), shell=True, timeout=60)
        
if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('startend', help="strat end times", nargs=2)
    p.add_argument('glon', help="Geo. longitude", type=float)
    p.add_argument('glat', help="Geo. Latitude", type=float)
    p.add_argument('altkm', help="Geo. Altitude", type=float)
    p.add_argument('--tsdo', help = 'Timestamp for the SDO image. If None, it is <tstart>', type = str, default = None)
    p.add_argument('-o', '--odir', help = 'Output directory.', default=None, type = str)
    p.add_argument('--sdodir', help = 'Directory with SDOAIA images.', default=None, type = str)
    p.add_argument('--dt', help = 'time step - minutes', default = 5, type = float)
    p.add_argument('--wl', help = 'SDO wavelength in Angstrom. Default = 193. If uniform type "geo".', type = str, default='193')
    p.add_argument('--srad', help = 'Solar raii inflation factor. Defult=1.0.', type = float, default=1.0)
    p.add_argument('--clim', help = 'Colorbar limits / units caunts. Default 0 -- 1000', type = int, nargs=2, default = [0, 1000])
    p.add_argument('--cmap', help = 'Colormap. Defualt = gist_stern', type = str, default='gist_stern')
    p.add_argument('--animation', help = 'Create a .mp4 movie', action='store_true')
    p.add_argument('--plot', help = 'Plot indivitual images?', action='store_true')
    P = p.parse_args()
#    
    main(startend=P.startend, glon=P.glon, glat=P.glat, alt_km=P.altkm, 
         tsdo = P.tsdo, dt = P.dt, wl = P.wl, odir = P.odir, sdodir = P.sdodir,
         clim = P.clim, cmap = P.cmap, animation=P.animation, plot=P.plot, srad_fact=P.srad)