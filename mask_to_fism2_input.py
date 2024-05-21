# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 21:30:55 2023

@author: sebastijan.mrak@colorado.edu
"""

import xarray as xr
import os, glob
import numpy as np
from dateutil import parser
import subprocess, platform

def main(folder, altkm:int, odir=None):
    files = np.unique(np.array(sorted(glob.glob(folder + f"*{altkm}km*_1.nc"))))
    file_dates = np.array([parser.parse(os.path.split(f)[1][:14]) for f in files])
    if odir is None:
        odir = folder + 'fism2_input' + os.sep
    if not os.path.exists(odir):
        if platform.system() in ('Linux', 'darwin'):
            subprocess.call(f'mkdir -p "{odir}"', shell=True)
        else:
            subprocess.call(f'mkdir "{odir}"', shell=True)
    for idf, d in enumerate(file_dates):
        print (f'{idf+1}/{files.size},  {d}')
        idt = np.isin(file_dates, d)
        tmp = files[idt]
        X = xr.Dataset()
        ofn = odir + f"{os.path.split(files[idf])[1][:20]}_alleof.nc"

        for iff, f in enumerate(tmp):
            if iff == 0:
                X['glon'] = xr.open_dataset(f).glon.values
                X['glat'] = xr.open_dataset(f).glat.values
                X['sza'] = (['glat', 'glon'], xr.open_dataset(f).sza.values)
            X[str(xr.open_dataset(f).wl.values)] = (['glat', 'glon'], xr.open_dataset(f).of.values)
        X.to_netcdf(ofn)
        
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('folder', help='Input folder', type=str)
    p.add_argument('--altkm', type=int, help='Height of the masks, Defaul=150km', default=150)
    p.add_argument('--odir', type=str, help='Path to output directory', default=None)
    
    P = p.parse_args()
    
    main(P.folder, P.altkm, P.odir)
