# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:00:33 2021

@author: smrak@bu.edu
"""

import os, subprocess, platform
import xarray, datetime
from eclipse import utils, eio
from dateutil import parser
import matplotlib.pyplot as plt
from argparse import ArgumentParser

def main(startend=None, glon=None, glat=None, alt_km=100, odir=None,
         wll=None, dt=10, aiafolder=None, tsdo=None, save=1, plot=0
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
    
    if aiafolder is None:
        import yaml
        try:
            stream = yaml.load(open(os.path.join(os.getcwd(), 'cfg', 'cfg.yaml'), 'r'), Loader=yaml.SafeLoader)
            aiafolder= stream.get('sdodir')
        except BaseException as e:
            raise(e)
    if aiafolder is None or (not os.path.exists(aiafolder)):
        aiafolder = input("type path to the aiafolder: \n")
    assert os.path.exists(aiafolder)
    if tsdo is None:
        tsdo = startend[0]
    assert isinstance(parser.parse(tsdo), datetime.datetime)
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
            times, of = utils.eof_time(tlim[0], tlim[1], glon, glat, ghgt, dm=0, ds=dt)
        else:
            if isinstance(wl, str):
                try:
                    wl = int(wl)
                    SDO = eio.sunaia(folder=aiafolder, wl=wl, time=parser.parse(tsdo))
                    times, of = utils.eof_time_sdo(SDO, tlim[0], tlim[1], parser.parse(tsdo), glon, glat, ghgt, wl=wl, dm=0,ds=dt)
                except:
                    continue
    
        if plot:
            plt.plot(times, of, label = wl)
            plt.legend()
        if i == 0:
            D['time'] = times
        D[str(wl)] = of
        
        
    if plot:
        plt.grid(axis='y')
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
    p.add_argument('-w', '--wl', help='Choose the wavelengths. Comma separated, "geo" for geometric eclipse. Geo is default', type=str, default='geo', nargs='+')
    p.add_argument('--plot', help='Plot the lines?', action='store_true')
    p.add_argument('--save', help='Do not save the result', action='store_false')
    P = p.parse_args()
    
    main(startend=P.startend, glon=P.glon, glat=P.glat, alt_km=P.altkm, odir=P.odir,
         wll=P.wl, dt=P.tres, aiafolder=P.sdodir, tsdo=None, save=P.save, plot=P.plot
         )