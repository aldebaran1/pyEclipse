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

root = os.getcwd()
odir = root + 'plots//'

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
glon_grid, glat_grid = np.meshgrid(glon,glat)
for it,T in enumerate(times):
    save_fn = odir + "{}_{}km_{}_{}.nc".format(T.strftime("%Y%m%d%H%M%S"), int(alt_km), 'geo', srad_fact)
    if not os.path.exists(save_fn):
        print ("Processing {} // {}/{}".format(T, it+1, times.size))
        sza = utils.get_sza(T, glon=glon_grid, glat=glat_grid, alt_km=alt_km)
        of = utils.mask_lonlat_geo(T=T, glon=glon, glat=glat, ghgt=ghgt, srad_fact=srad_fact)
        # TO XARRAY
        X = xarray.Dataset(
            {
                "of": (("glat", "glon"), of.T),
                "sza": (("glat", "glon"), sza)
            },
            {"glon": glon, "glat": glat}
        )
        X['time'] = T
        X['alt_km'] = alt_km
        X['srad_fact'] = srad_fact
        X.to_netcdf(save_fn)
        X.close()

