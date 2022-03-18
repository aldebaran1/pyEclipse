# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:00:33 2021

@author: smrak@bu.edu
"""

import os, subprocess, platform
import xarray, yaml
from eclipse import utils, eio
from dateutil import parser
from datetime import datetime
import numpy as np
import concurrent.futures
from argparse import ArgumentParser

def _eof(sep, azm, eta, mrad, wl):
#    global SDO, imres, rad2pix, wl
    return utils.get_eof_mask_from_angles(SDO['AIA{}'.format(wl)].values, sep,azm,eta,mrad, x0, y0, imres=imres,pixscale=rad2pix)

def main(startend, odir, tsdo=None, glonlim=[-180,180], glatlim=[-90,90], alt_km=0, 
         wl=193, dlon=None, dlat=None, dt=10, j=1,
         aiafolder=None, instrument=None, parallactic_angle=False, sn=None
         ):
    global SDO, x0, y0, imres, rad2pix
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
        try:
            print (os.path.join(os.getcwd(), 'cfg', 'cfg.yaml'))
            stream = yaml.load(open(os.path.join(os.getcwd(), 'cfg', 'cfg.yaml'), 'r'), Loader=yaml.SafeLoader)
            aiafolder= stream.get('sdofolder')
        except BaseException as e:
            raise(e)
        print (aiafolder)
        if aiafolder is None or (not os.path.exists(aiafolder)):
                aiafolder = input("Type valid path to the sdoaia folder:\n")
    print (aiafolder)
    assert os.path.exists(aiafolder), "AIAFOLDER doesn't exists"

    if tsdo is None:
        tsdo = startend[0]
    times = utils.get_times(parser.parse(startend[0]), parser.parse(startend[1]), dm=dt)
    glon = np.arange(glonlim[0], glonlim[1], dlon)
    glat = np.arange(glatlim[0], glatlim[1], dlat)
    ghgt = alt_km * 1000
    SDO = eio.load(folder=aiafolder, wl=wl, time=parser.parse(tsdo), 
                   instrument=instrument, sn=sn)
    x0, y0 = SDO.x0, SDO.y0
    
    glon_grid, glat_grid = np.meshgrid(glon,glat)
    rad2pix = SDO.pixscale
    imres = int(SDO['AIA{}'.format(wl)].shape[0])
    critical_distance = imres/2 + imres*1.414
    
    for it,T in enumerate(times):
        save_fn = odir + "{}_{}km_{}_{}.nc".format(T.strftime("%Y%m%d%H%M%S"), int(alt_km), wl, int(parallactic_angle))
        if not os.path.exists(save_fn):
            print ("Processing {} // {}/{}".format(T, it+1, times.size))
            tt = datetime.now()
            of = np.nan * np.ones((glat.size, glon.size))
            sza = utils.get_sza(T,glon=glon_grid,glat=glat_grid,alt_km=alt_km)
            sep, azm, mrad = utils.get_angles(T, glon=glon_grid, glat=glat_grid, ghgt=ghgt)
            if parallactic_angle:
                eta = utils.get_parallactic_angle(T, glon=glon_grid, glat=glat_grid, ghgt=ghgt)
            else:
                eta = np.zeros(glon_grid.shape)
            mask = (sep*rad2pix > critical_distance)
            of[mask] = 1.0
            mask = (sza > 90)
            of[mask] = 0.0

            ix = np.where(np.isnan(of))
            print ("Loading for one iteration {}".format(datetime.now()-tt))
            with concurrent.futures.ThreadPoolExecutor(max_workers=j) as ex:
                of_worker = np.asarray([ex.submit(_eof, sep[ix[0][i],ix[1][i]], azm[ix[0][i],ix[1][i]], eta[ix[0][i],ix[1][i]], mrad[ix[0][i],ix[1][i]],wl) for i in range(ix[0].size)])
            print ("Total {}".format(datetime.now()-tt))
            for i in range(of_worker.size):
                of[ix[0][i],ix[1][i]] = of_worker[i].result()
            # TO XARRAY
            X = xarray.Dataset(
                {
                    "of": (("glat", "glon"), of),
                    "sza": (("glat", "glon"), sza),
                },
                {"glon": glon, "glat": glat}
            )
            X['time'] = T
            X['alt_km'] = alt_km
            X['wl'] = wl
            X['time_image'] = SDO.time
            X['instrument'] = instrument
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
    p.add_argument('--tsdo', help='Time for the SDO image', default=None)
    p.add_argument('--sdodir', help='Folder to SDO images', default=None)
    p.add_argument('-w', '--wl', help='Choose the wavelength. Default is 193A for the SDO', default = 193, type=int)
    p.add_argument('--altkm', help='altitude in km', default = 350, type=int)
    p.add_argument('--sn', help='GOES satellite number', default = 16, type=int)
    p.add_argument('-i', '--instrument', help='AIA, SUVI or EIT', default = 'aia', type=str)
    p.add_argument('-j', '--proc', help='Numer of parallel processes', default=10, type=int)
    p.add_argument('--eta', help='Utilize parallactic angle?', action='store_true')
    
    P = p.parse_args()
    
    main(P.startend, odir=P.odir, dt = P.tres, tsdo = P.tsdo,
         glonlim=P.glon, glatlim=P.glat, dlon=P.dlon, dlat=P.dlat,
         wl = P.wl, alt_km= P.altkm, aiafolder = P.sdodir, j=P.proc,
         instrument=P.instrument, sn=P.sn, parallactic_angle=P.eta)
