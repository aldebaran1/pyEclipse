# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:34:29 2019

@author: smrak
"""

from glob import glob
import numpy as np
from dateutil import parser
import xarray
import os
import urllib
from datetime import datetime
from typing import Union
from numpy import float64, pi, array, ones, arange
from astropy.io import fits
from scipy import ndimage

# https://ftp.iana.org/tz/tzdb-2019a/leap-seconds.list
LEAPS = arange(10, 38)
LEAP_dates = ['1 Jan 1972', '1 Jul 1972', '1 Jan 1973',
         '1 Jan 1974', '1 Jan 1975', '1 Jan 1976',
         '1 Jan 1977', '1 Jan 1978', '1 Jan 1979',
         '1 Jan 1980', '1 Jul 1981', '1 Jul 1982',
         '1 Jul 1983', '1 Jul 1985', '1 Jan 1988',
         '1 Jan 1990', '1 Jan 1991', '1 Jul 1992',
         '1 Jul 1993', '1 Jul 1994', '1 Jan 1996',
         '1 Jul 1997', '1 Jan 1999', '1 Jan 2006',
         '1 Jan 2009', '1 Jul 2012', '1 Jul 2015',
         '1 Jan 2017']

oneradian_arcsec = (180 * 3600) / pi
ddd2D = {'Jan': 1, 'Jul': 7}

def get_filename(folder, wl, time, instrument, sn=None):
    
    assert (os.path.exists(folder)), f'{folder} doesnt exists'
        
    if time is None:
        if instrument == 'aia':
            seekfor = '*{}a*.fits'.format(wl)
        elif instrument == 'eit':
            seekfor = 'efz*'
        elif instrument.lower() == 'xrt':
            seekfor = 'synop_XRT*'
        else:
            print ("Entered instrument is not on th elist of supported telescopes")
        try:
            fn = glob(folder + seekfor)[0]
            print ("File choosen: " + os.path.split(fn)[1])
        except Exception as e:
            raise (e)
    else:
        if isinstance(time, str):
            time = parser.parse(time)
        assert isinstance(time, datetime), 'Time must be in appropraite format (str or datetime)'
        
        if instrument in ('aia', 'AIA'):
            wlfilt = '*{}a*.fits'.format(wl)
            fnlist = np.array(glob(folder + wlfilt))
            names = [os.path.split(ff)[1] for ff in fnlist]
            if len(str(wl)) == 3:
                filedates = [n[14:27].replace('_', '-').upper() + n[27:33].replace('_', ':') for n in names]
            elif len(str(wl)) == 2:
                filedates = [n[13:26].replace('_', '-').upper() + n[26:32].replace('_', ':') for n in names]
            elif len(str(wl)) == 4:
                filedates = [n[15:28].replace('_', '-').upper() + n[28:34].replace('_', ':') for n in names]
            else:
                raise ('Wrong wavelength argument')
            fdate_dt = array([parser.parse(d) for d in filedates])
        elif instrument in ('eit', 'EIT'):
            wlfilt = 'efz{}*'.format(time.strftime('%Y%m%d'))
            fnlist = np.array(glob(folder + wlfilt))
            names = [os.path.split(ff)[1] for ff in fnlist]
            filedates = [f'{n[3:11]}T{n[12:]}' for n in names]
            fdate_dt = np.array([parser.parse(d) for d in filedates])
            
        elif instrument in ('suvi', 'SUVI'):
            fnlist = np.array(sorted(glob(folder + '*suvi*.fits')))
            suvi_wl = np.array([int(os.path.split(f)[1][13:16]) for f in fnlist])
            suvi_sn = np.array([int(os.path.split(f)[1][18:20]) for f in fnlist])
            fdate_dt = np.array([datetime.strptime(os.path.split(f)[1][22:37], '%Y%m%dT%H%M%S') for f in fnlist])
            
            idr = np.logical_and(suvi_sn==sn, suvi_wl==wl)
            fdate_dt = fdate_dt[idr]
            fnlist = fnlist[idr]
        elif instrument.lower() == 'xrt':
            fnlist = np.array(glob(folder + 'synop_XRT*.fits'))
            fdate_dt = np.array([datetime.strptime(os.path.split(f)[1][9:-7], '%Y%m%d_%H%M%S') for f in fnlist])
        
        if fdate_dt.size > 1:    
            idX = abs(fdate_dt - time).argmin()
            if (abs(fdate_dt - time)[idX].total_seconds()) > 24*60*60:
                    print ('No images available within 24 hours from the given time, ', time)
                    return 0
            fn = fnlist[idX]
        elif fdate_dt.size == 1:
            if (abs(fdate_dt[0] - time).total_seconds()) > 24*60*60:
                print ('No images available within 24 hours from the given time')
                return 0
            else:
                fn = fnlist[0]
        else:
            raise("No matches were found")
        print ("File choosen: " + os.path.split(fn)[1])
    return fn

def load(folder: str = None, 
           wl: int = 193, 
           time: Union[datetime, str] = None,
           sn: int = None,
           instrument: str = 'aia',
           xrt_treshold: int = 3) -> xarray.Dataset:
    
    assert (instrument in ('aia' , 'eit', 'suvi', 'AIA', 'EIT', 'SUVI', 'xrt', 'XRT')), "Currently we support SDO AIA, SOHO EIT, Hinode XRT, and GOES-R SUVI telescopes"
    
    if not os.path.isfile(folder):
        fn = get_filename(folder, wl, time, instrument, sn)
    else:
        fn = folder
    
    ix = 1 if instrument in ('aia', 'AIA', 'suvi', 'SUVI') else 0
    
    try:
        FITS = fits.open(fn)
    except Exception as e:
        raise (e)
    FITS[ix].verify('fix')
    # bitpix = FITS[ix].header['BITPIX']
    if instrument.lower() in ('aia', 'suvi') :
        imtime =  parser.parse(FITS[ix].header['DATE-OBS']) 
    # elif instrument.lower() in ('suvi'):
        # print(FITS[ix].header)
        # imtime = parser.parse(FITS[ix].header['DATE_OBS'])
    elif instrument.lower() in ('xrt', 'eit'):
        # imtime =  parser.parse(FITS[ix].header['DATE-OBS']) 
        # print (FITS[1].header)
        imtime = parser.parse(FITS[ix].header['DATE_OBS'])
    else:
        print ("Instrument not on the list of available data")
        return 0
            
    
    imx0 = FITS[ix].header['CRPIX1']
    imy0 = FITS[ix].header['CRPIX2']
    pixel2arcsec = FITS[ix].header['CDELT1']
    if instrument.lower() in ('aia', 'suvi', 'eit'): 
        wl = FITS[ix].header['WAVELNTH']
    elif instrument.lower() in ('xrt'):
        wl = 'X'
    else:
        print ("Instrument not on the list of available data")
        return 0
        
    #assert (FITS[ix].header['BITPIX'] == bitpix), 'Assure the data encoded with uint-16bit'
    
    im = float64(FITS[ix].data.squeeze())
    if instrument in ('eit', 'EIT'):
        im = ndimage.rotate(im, 180)
    # Set data outside the detector values to 0
    im[im < 0] = 0
    im[im > 2**16] = 0
    if instrument.lower() == 'xrt':
        im[im < xrt_treshold] = 0
        im = ndimage.median_filter(ndimage.median_filter(im, 3), 3)
    # Construct  xarray for easier data manipulation and access
    D = xarray.Dataset({f'{instrument.upper()}{wl}': (('x', 'y'), im)})
    D.attrs['time'] = imtime
    D.attrs['x0'] = imx0
    D.attrs['y0'] = imy0
    D.attrs['pxarcsec'] = pixel2arcsec
    D.attrs['pixscale'] = oneradian_arcsec / pixel2arcsec
    if instrument.lower() == 'xrt':
        D.attrs['mode'] = FITS[ix].header['EC_FW1_']
    
    return D

#def sunaia(folder: str = None, 
#           wl: int = 193, 
#           time: Union[datetime, str] = None):
#    if folder is None:
#        folder = os.path.join(os.getcwd(), 'aia', "")
#        
#    if time is None:
#        seekfor = '*{}a*.fits'.format(wl)
#        try:
#            fn = glob(folder + seekfor)[0]
#            print ("File choosen: " + os.path.split(fn)[1])
#        except Exception as e:
#            raise (e)
#    else:
#        if isinstance(time, str):
#            time = parser.parse(time)
#        assert isinstance(time, datetime), 'Time must be in appropraite format (str or datetime)'
#        
#        wlfilt = '*{}a*.fits'.format(wl)
#        fnlist = glob(folder + wlfilt)
#        names = [os.path.split(ff)[1] for ff in fnlist]
#        if len(str(wl)) == 3:
#            filedates = [n[14:27].replace('_', '-').upper() + n[27:33].replace('_', ':') for n in names]
#        elif len(str(wl)) == 2:
#            filedates = [n[13:26].replace('_', '-').upper() + n[26:32].replace('_', ':') for n in names]
#        elif len(str(wl)) == 4:
#            filedates = [n[15:28].replace('_', '-').upper() + n[28:34].replace('_', ':') for n in names]
#        else:
#            raise ('Wrong wavelength argument')
#        fdate_dt = array([parser.parse(d) for d in filedates])
#        idX = abs(fdate_dt - time).argmin()
#        if (abs(fdate_dt - time)[idX].total_seconds()) > 12*60*60:
#            raise ValueError ("Closest SDO image is more than 12 hours away.")
#        fn = fnlist[idX]
#        print ("File choosen: " + os.path.split(fn)[1])
#            
#    try:
#        imdata = fits.open(fn)
#    except Exception as e:
#        raise (e)
#    imtime = parser.parse(imdata[1].header['T_REC'])
#    imx0 = imdata[1].header['CRPIX1']
#    imy0 = imdata[1].header['CRPIX2']
#    pixel2arcsec = imdata[1].header['CDELT1']
#    wl = imdata[1].header['WAVELNTH']
#    
#    assert (imdata[1].header['BITPIX'] == 16), 'Assure the data encoded with uint-16bit'
#    
#    imdata[1].verify('fix')
#    im = float64(imdata[1].data.squeeze())
#    # Set data outside the detector values to 0
#    im[im < 0] = 0
#    im[im > 2**16] = 0
#    # Construct  xarray for easier data manipulation and access
#    D = xarray.Dataset({'AIA'+str(wl): (('x', 'y'), im)})
#    D.attrs['time'] = imtime
#    D.attrs['x0'] = imx0
#    D.attrs['y0'] = imy0
#    D.attrs['pxarcsec'] = pixel2arcsec
#    D.attrs['pixscale'] = oneradian_arcsec / pixel2arcsec
#    return D

def spaceTimeCorrections(year: int = None, month: int = None, day: int = None):
    """
    """
    def _get_last(doc):
        s = []
        for line in doc:
            s.append(sum(c.isdigit() for c in line[135:]))
        s = np.array(s)
        ixall = np.where(s > 0)[0]
        ix = ixall[-1]
        xp = float(doc[ix][135:145])
        yp = float(doc[ix][145:155])
        ut1utc = float(doc[ix][155:165])
        return xp, yp, ut1utc
    if year is None and month is None and day is None: 
        year = 2017
        month = 8
        day = 21
        print ('Default date is set for the 21 August 2017 eclipse')
    assert year is not None, 'set year argument "YYYY"'
    assert month is not None, 'set month argument "MM"'
    assert day is not None, 'set day argument "DD"'
    date = datetime(year, month, day)
    link = "ftp://ftp.iers.org/products/eop/rapid/standard/finals2000A.all"
    urldoc = urllib.request.urlopen(link)
    listdoc = urldoc.readlines()
    txtdoc = array([a.decode('ascii').replace('\n',"") for a in listdoc])
    if date.month >= 10 and date.day >= 10:
        match = date.strftime('%y%m%d')
        isin = np.char.find(txtdoc, match)
    elif date.month >= 10 and date.day < 10:
        match = "{} {}".format(date.strftime('%y%m'), date.strftime('%y%m').replace("0", ""))
        isin = np.char.find(txtdoc, match)
    elif date.month < 10 and date.day >= 10:
        match = "{} {}{}".format(date.strftime('%y'), date.strftime('%m').replace("0", ""), date.strftime('%d'))
        isin = np.char.find(txtdoc, match)
    else:
        match = "{} {} {}".format(date.strftime('%y'), date.strftime('%m').replace("0", ""), date.strftime('%d').replace("0", ""))
        isin = np.char.find(txtdoc, match)
    ix = np.where(isin > -1)[0]
    if ix.size == 1:
        ix = ix.item()
    elif ix.size > 1:
        for l in ix:
            if txtdoc[l][:6] == match:
                ix = l
    else:
        raise ("Somethng went wrong")
    line = txtdoc[ix]
    if sum(c.isdigit() for c in line[135:]) > 1:
        xp = float(line[135:145])
        yp = float(line[145:155])
        ut1utc = float(line[155:165])
    else:
        xp, yp, ut1utc = _get_last(txtdoc)

    YY = ones(LEAPS.shape[0], dtype=int)
    MM = ones(LEAPS.shape[0], dtype=int)
    for j, line in enumerate(LEAP_dates):
        YY[j] = int(line[5:])
        MM[j] = int(ddd2D[line[2:5]])
    
    idY = year >= YY
    leap = LEAPS[idY][-1]
    params = {'xp': xp, 'yp': yp, 'ut1utc': ut1utc, 'leap': leap}
    
    return params
#
