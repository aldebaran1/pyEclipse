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
import os
from argparse import ArgumentParser
from subprocess import call
#odir = os.getcwd() + '/aia/'
from dateutil import parser
from datetime import timedelta

def dlSun(tlim: Union[str, list, ndarray] = [],
          instrument: str = 'aia',
          wl: int = 193,
          odir: str = '',
          first: bool = False):
    if len(tlim) == 1:
        tlim_dt = [parser.parse(tlim[0]), parser.parse(tlim[0])+timedelta(minutes=1)]
    
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