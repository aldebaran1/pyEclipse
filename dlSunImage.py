#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 14:21:05 2019

@author: smrak
"""

from sunpy.net import Fido, attrs as a
#from sunpy.net.vso import attrs as vso_attrs
import astropy.units as u
from numpy import ndarray, array
from typing import Union
#from pathlib import Path
import numpy as np
import os
from argparse import ArgumentParser
from subprocess import call
#odir = os.getcwd() + '/aia/'
from dateutil import parser
from datetime import timedelta, datetime
import wget
import urllib.request
from bs4 import BeautifulSoup

def dlSun(tlim: Union[str, list, ndarray] = [],
          instrument: str = 'aia',
          wl: int = 193,
          odir: str = '',
          first: bool = False):
    
    if 'suvi' in instrument:
        sn = int(instrument[-2:])
        instrument = 'suvi'
    
    if len(tlim) == 1:
        if instrument in ('aia', 'suvi'):
            tlim_dt = [parser.parse(tlim[0]), parser.parse(tlim[0])+timedelta(hours=1)]
        else:
            tlim_dt = [parser.parse(tlim[0])-timedelta(hours=parser.parse(tlim[0]).hour), parser.parse(tlim[0])+timedelta(hours=24-parser.parse(tlim[0]).hour)]
        
        first = True
    else:
        tlim_dt = [parser.parse(tlim[0]), parser.parse(tlim[1])]
        
    tlim = array([tlim_dt[0].strftime('%Y/%m/%d %H:%M'), tlim_dt[1].strftime('%Y/%m/%d %H:%M')])
    tlim = array(tlim) if isinstance(tlim, list) else tlim
    assert tlim.shape[0] == 2, 'tlim must have length of 2 (start, stop) (yyyy/mm/dd hh:mm)'
    
    if odir == "":
        odir = os.path.join('C:\\Users\\smrak\\Google Drive\\BU\\Projects\\Eclipse-Nsf\\data\\', instrument, "")
    if not os.path.exists(odir):
        call('mkdir "{}"'.format(odir), shell=True)
    
    attrs_time = a.Time(tlim[0], tlim[1])
    
    
    if instrument == 'suvi':
        suvimap = {'94': 'fe094', '131': 'fe131', '171': 'fe171', '195': 'fe195', '284': 'fe284', '304': 'he304'}
        suvimap = {'94': 'ci094', '131': 'ci131', '171': 'ci171', '195': 'ci195', '284': 'ci284', '304': 'ci304'}
        
#        import tarfile
        line = suvimap[str(wl)]
        # url = f'https://www.ncei.noaa.gov/data/goes-r-series-l2-operational-space-weather-products/access/goes{sn}/suvi/{tlim_dt[0].strftime("%Y")}/{tlim_dt[0].strftime("%m")}/{tlim_dt[0].strftime("%d")}/'
        url = f'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{sn}/l2/data/suvi-l2-{line}/{tlim_dt[0].strftime("%Y/%m/%d/")}'
#        url = f'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{sn}/suvi/{tlim_dt[0].strftime("%Y")}/{tlim_dt[0].strftime("%m")}/{tlim_dt[0].strftime("%d")}/'
        # url = f'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes{sn}/l1b/suvi-l1b-{line}/{tlim_dt[0].strftime("%Y")}/{tlim_dt[0].strftime("%m")}/{tlim_dt[0].strftime("%d")}/'
#        fnmask = f'OR_SUVI-L1b-{line}_G{sn}_{tlim_dt[0].strftime("%Y%j%H%M")}*fits.gz'
        filenames = []
        filenamedates = []
        with urllib.request.urlopen(url) as response:
            html = response.read().decode('ascii')
            soup = BeautifulSoup(html, 'html.parser')
            for link in soup.find_all('a'):
                if link.get('href') is not None and link.get('href')[:2] == 'dr':
                    filenames.append(link.get('href'))
                    filenamedates.append(datetime.strptime(link.get('href')[22:38], '%Y%m%dT%H%M%SZ'))
        fnd = np.array(filenamedates)
        idt = abs(fnd - tlim_dt[0]).argmin()
        f = np.array(filenames)[idt]
        try:
#            subprocess.call(f'wget {url} -P "{odir}" -A {fnmask}', shell=True)
            wget.download(url + f, out=odir)
            # import gzip
            # import shutil
            # import subprocess
            # with gzip.open(f'{odir}{f}', 'rb') as f_in:
                # with open(f'{odir}{f[:-3]}', 'wb') as f_out:
                    # shutil.copyfileobj(f_in, f_out)
            # subprocess.call(f'del /f "{odir}{f}"', shell=True)
        except Exception as e:
            print (e)
            # print (f"File doesn't exists:\n{url + np.array(filenames)[idt]}")
            # print ('Check if data exists for this date first on NOAA webpage')
            # print (url)
        
        # SunPy interface doesn't work
#        result = Fido.search(attrs_time,
#                         a.Instrument(instrument),
#                         a.Wavelength(wl * u.angstrom),
#                         a.goes.SatelliteNumber(sn),
#                         a.Level.two)
    else:
        result = Fido.search(attrs_time,
                         a.Instrument(instrument),
                         a.Wavelength(wl * u.angstrom))
    
    
        if first:
            print ('Downloading to {}:'.format(odir))
            print (result[0,0])
            Fido.fetch(result[0,0], path=odir)
        else:
            print ('Downloading {} elements.... to {}:'.format(result.file_num, odir))
            print (result)
            Fido.fetch(result, path=odir)
        print ('Downloading complete.')


def main():
    p = ArgumentParser()
    p.add_argument('startend', help='start/end times UTC e.g. "yyyy-mm-ddThh:mm"', nargs='+')
    p.add_argument('-o', '--odir', help='directory to write downloaded FITS to', default='')
    p.add_argument('-i', '--instrument', help='aia or eit (SDO or SOHO), default is SDO-AIA', default='aia')
    p.add_argument('-w', '--wl', help='Choose the wavelength. Default is 193A for the SDO', default = 193)
    p.add_argument('--first', help="download only the first image from the result list!", action='store_true')
    
    P = p.parse_args()
    
    dlSun(tlim = P.startend, instrument = P.instrument, wl=int(P.wl), odir = P.odir, first = P.first)


if __name__ == '__main__':
    main()

#
#tlim = ['2017/08/21 17:00', '2017/08/21 17:01']
#instrument = 'aia'
#wl = 193
#odir = os.path.join(os.getcwd(),'aia', "")
#dlSun(tlim = tlim, instrument=instrument, wl=wl, odir=odir)