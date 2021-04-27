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

def sunaia(folder: str = None, 
           wl: int = 193, 
           time: Union[datetime, str] = None):
    if folder is None:
        folder = os.path.join(os.getcwd(), 'aia', "")
        
    if time is None:
        seekfor = '*{}a*.fits'.format(wl)
        try:
            fn = glob(folder + seekfor)[0]
            print ("File choosen: " + os.path.split(fn)[1])
        except Exception as e:
            raise (e)
    else:
        if isinstance(time, str):
            time = parser.parse(time)
        assert isinstance(time, datetime), 'Time must be in appropraite format (str or datetime)'
        
        wlfilt = '*{}a*.fits'.format(wl)
        fnlist = glob(folder + wlfilt)
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
        idX = abs(fdate_dt - time).argmin()
        
        fn = fnlist[idX]
        print ("File choosen: " + os.path.split(fn)[1])
            
    try:
        imdata = fits.open(fn)
    except Exception as e:
        raise (e)
    imtime = parser.parse(imdata[1].header['T_REC'])
    imx0 = imdata[1].header['CRPIX1']
    imy0 = imdata[1].header['CRPIX2']
    pixel2arcsec = imdata[1].header['CDELT1']
    wl = imdata[1].header['WAVELNTH']
    
    assert (imdata[1].header['BITPIX'] == 16), 'Assure the data encoded with uint-16bit'
    
    imdata[1].verify('fix')
    im = float64(imdata[1].data.squeeze())
    # Set data outside the detector values to 0
    im[im < 0] = 0
    im[im > 2**16] = 0
    # Construct  xarray for easier data manipulation and access
    D = xarray.Dataset({'AIA'+str(wl): (('x', 'y'), im)})
    D.attrs['time'] = imtime
    D.attrs['x0'] = imx0
    D.attrs['y0'] = imy0
    D.attrs['pxarcsec'] = pixel2arcsec
    D.attrs['pixscale'] = oneradian_arcsec / pixel2arcsec
    return D

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
#def parallacticAngle(sazm, sdec, glat):
#    sineta = np.sin(sazm) * np.cos(np.radians(glat)) / np.cos(sdec)
#    return -np.arcsin(sineta)



#def parallacticAngle(sazm, selv, glat):
#    DOUG ????
#    """
#    Compute solar parallactic angle. This is an angle between loacl (observer's)
#    zenith and sun's axis of rotation.
#    
#    Parameters
#    ----------
#    sazm: float
#          Solar azimuth angle
#    selev: float
#          Solar elevation angle
#    glat: float
#         geographic latitude in WSG84 [degree]
#    Returns
#    ---------
#    float: solar parallactic angle
#    """
#    
#    olat = np.deg2rad(glat)
#    a = -np.cos(olat) * np.sin(sazm)
#    b = np.sin(olat) * np.cos(selv) - (np.cos(olat) * np.sin(selv) * np.cos(sazm))
#    eta = np.arctan2(a,b)
#    
#    return eta