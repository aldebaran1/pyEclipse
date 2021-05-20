# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:00:33 2021

@author: smrak@bu.edu
"""

import os, subprocess, platform
import xarray
from eclipse import utils
from dateutil import parser
import numpy as np
from argparse import ArgumentParser


def main(startend, odir, glonlim=[-180,180], glatlim=[-90,90], alt_km=0, 
         dlon=None, dlat=None, dt=10, srad_fact=1,
         ):
    """
    dt = delta_t between t0 and t1 [minutes]
    glonlim, glatlim is a listwith [glon_min, glon_max]
    dlon, dlat is a spacing in lon and lat
    alt_km is alitutde in kilometer
    """
    if not os.path.exists(odir):
        if platform.system() == 'Windows':
            subprocess.call('mkdir "{}"'.format(odir), shell=True)
        elif platform.system() == 'Linux':
            subprocess.call('mkdir -p {}'.format(odir), shell=True)
        else:
            print ("Platform {} is currently not supported.".format(platform.system()))
            exit()
    
    assert os.path.exists(odir)

    times = utils.get_times(parser.parse(startend[0]), parser.parse(startend[1]), dm=dt)
    glon = np.arange(glonlim[0], glonlim[1], dlon)
    glat = np.arange(glatlim[0], glatlim[1], dlat)
    ghgt = alt_km * 1000
    
    for it,T in enumerate(times):
        save_fn = odir + "{}_{}km_{}_{}.nc".format(T.strftime("%Y%m%d%H%M%S"), int(alt_km), 'geo', srad_fact)
        if not os.path.exists(save_fn):
            print ("Processing {} // {}/{}".format(T, it+1, times.size))
            of = utils.mask_lonlat_geo(T=T, glon=glon, glat=glat, ghgt=ghgt, srad_fact=srad_fact)
            # TO XARRAY
            X = xarray.Dataset(
                {
                    "of": (("glat", "glon"), of.T),
                },
                {"glon": glon, "glat": glat}
            )
            X['time'] = T
            X['alt_km'] = alt_km
            X['srad_fact'] = srad_fact
            X.to_netcdf(save_fn)
            X.close()

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('startend', help='start and end times UTC e.g. "yyyy-mm-ddThh:mm"', nargs=2)
    p.add_argument('odir', help='directory to write downloaded FITS to', default='')
    p.add_argument('--tres', help='Time resolution in minutes. Default=15min', default=15, type=int)
    p.add_argument('--glon', help='start and end GLON, Default=-180 180', nargs=2, default=[-180, 180], type=int)
    p.add_argument('--glat', help='start and end GLAT, Default=-90 90', nargs=2, default=[-90, 90], type=int)
    p.add_argument('--dlon', help='GLON spacing, Default=1', type=int, default=1)
    p.add_argument('--dlat', help='GLAT spacing, Default=1', type=int, default=1)
    p.add_argument('--altkm', help='altitude in km', default = 350, type=int)
    p.add_argument('-s', '--sunradii', help='Solar radii inflation factor', default=1.0, type=float)
    
    P = p.parse_args()
    
    main(P.startend, odir=P.odir, dt = P.tres,
         glonlim=P.glon, glatlim=P.glat, dlon=P.dlon, dlat=P.dlat,
         alt_km= P.altkm, srad_fact=P.sunradii)
