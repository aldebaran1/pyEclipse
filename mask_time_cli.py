# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:00:33 2021

@author: smrak@bu.edu
"""

import os, subprocess, platform
import xarray, datetime
import numpy as np
from eclipse import utils, eio
from dateutil import parser
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from argparse import ArgumentParser

def main(startend=None, glon=None, glat=None, alt_km=100, odir=None,
         wll=None, dt=10, aiafolder=None, tsdo=None, srad=1.0,save=1, plot=0,
         instrument=None,
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
    
    instrument = instrument.upper()
    assert instrument in ('AIA', 'EIT', 'SUVI')
    
    if aiafolder is None:
        import yaml
        try:
            stream = yaml.load(open(os.path.join(os.getcwd(), 'cfg', 'cfg.yaml'), 'r'), Loader=yaml.SafeLoader)
            aiafolder= stream.get('sdodir')
        except BaseException as e:
            raise(e)
    if len(wll) > 1 or np.squeeze(wll) != 'geo':
        if aiafolder is None or (not os.path.exists(aiafolder)):
            if platform.system() == 'Windows':
                aiafolder = 'G:\\My Drive\\eclipse\\sdoaia\\'
            else:
                aiafolder = input("type path to the aiafolder: \n")
        assert os.path.exists(aiafolder), "This folder doesn't exists"
    
    if tsdo is None:
        tsdo = startend[0]
    assert isinstance(parser.parse(tsdo), datetime.datetime), "Time format wrong"
    
    tlim = [parser.parse(startend[0]), parser.parse(startend[1])]
    ghgt = alt_km * 1000
    save_fn = os.path.join(odir, 
                           '{}lat_{}lon_{}_{}.nc'.format(glat,glon,
                                tlim[0].strftime("%Y%m%d%H%M%S"), 
                                tlim[1].strftime("%Y%m%d%H%M%S")), 
                           )
    
    if plot:
        plt.figure()
    
    D = {}
    for i,wl in enumerate(wll):
        if wl == 'geo':
            times, of = utils.eof_time(tlim[0], tlim[1], glon, glat, ghgt, dm=0, ds=dt, srad_fact=srad)
        else:
            if isinstance(wl, str):
                try:
                    wl = int(wl)
                    SDO = eio.load(folder=aiafolder, wl=wl, time=parser.parse(tsdo), instrument=instrument)
                    times, of = utils.eof_time_sdo(SDO=SDO, t0=tlim[0], t1=tlim[1], glon=glon, glat=glat, ghgt=ghgt, wl=wl, dm=0, ds=dt)
                except:
                    continue
    
        if plot:
            if wl == 'geo':
                plt.plot(times, of, label = f'{wl} - {srad}')
            else:
                plt.plot(times, of, label = wl)
            
        if i == 0:
            D['time'] = times
        D[str(wl)] = of
        
    
    if plot:
        plt.grid(axis='y')
        plt.legend()
        ax = plt.gca()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.show()
    if save:
        X = xarray.Dataset(D)
        X['glon'] = glon
        X['glat'] = glat
        X['alt_km'] = alt_km
        X.to_netcdf(save_fn)
        
        X.close()


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('startend', help='start and end times UTC e.g. "yyyy-mm-ddThh:mm"', nargs=2)
    p.add_argument('glon', help='Geodetic longitude', type=float)
    p.add_argument('glat', help='Geodetic latitude', type=float)
    p.add_argument('odir', help='directory to write netcdf file', type=str)
    p.add_argument('--altkm', help='Altitude in km. Defult=0', default = 0, type=int)
    p.add_argument('--tres', help='Time resolution in seconds. Default=60sec', default=60, type=int)
    p.add_argument('--tsdo', help='Time for the SDO image', default=None)
    p.add_argument('--sdodir', help='Folder to SDO images', default=None)
    p.add_argument('--instrument', help='Instrument, AIA,EIT, or SUVI. Default=AIA', default='AIA')
    p.add_argument('-w', '--wl', help='Choose the wavelengths. Comma separated, "geo" for geometric eclipse. Geo is default', type=str, default='geo', nargs='+')
    p.add_argument('--plot', help='Plot the lines?', action='store_true')
    p.add_argument('--save', help='Do not save the result', action='store_false')
    p.add_argument('--srad', help='Altitude in km. Defult=0', default = 1.0, type=float)
    P = p.parse_args()
    
    main(startend=P.startend, glon=P.glon, glat=P.glat, alt_km=P.altkm, odir=P.odir,
         wll=P.wl, dt=P.tres, aiafolder=P.sdodir, tsdo=P.tsdo, srad=P.srad,
         save=P.save, plot=P.plot,instrument=P.instrument
         )