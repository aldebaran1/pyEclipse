# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 13:17:51 2021

@author: smrak@bu.edu
"""

import numpy as np
from datetime import datetime, timedelta
import ephem
import concurrent.futures
from scipy import ndimage
RE = 6371 #km

def get_sza(time, glon, glat, horizon=None, alt_km=None):
    
    if horizon is None:
        if alt_km is None:
            alt_km = 0
        elif isinstance(alt_km, (list, np.ndarray)):
            alt_km = alt_km[0]
        horizon = -np.degrees(np.arccos(RE/(RE + alt_km)))
    
    def _sza(x, y):
        
        obs = ephem.Observer()
        obs.lat = np.deg2rad(y)
        obs.lon = np.deg2rad(x)
        obs.date = ephem.Date(time)
        
        sun = ephem.Sun()
        sun.compute(obs)
        sza = 90 - np.degrees(sun.alt) + horizon
        return sza
    def _sza_time(t):
        
        obs = ephem.Observer()
        obs.lat = np.deg2rad(glat)
        obs.lon = np.deg2rad(glon)
        obs.date = ephem.Date(t)
        
        sun = ephem.Sun()
        sun.compute(obs)
        sza = 90 - np.degrees(sun.alt) + horizon
        return sza
    
    def _sza_all(t, x, y):
        obs = ephem.Observer()
        obs.lat = np.deg2rad(y)
        obs.lon = np.deg2rad(x)
        obs.date = ephem.Date(t)
        
        sun = ephem.Sun()
        sun.compute(obs)
        sza = 90 - np.degrees(sun.alt) + horizon
        return sza
    
    if isinstance(glon, np.ndarray) and isinstance(time, datetime):
        with concurrent.futures.ThreadPoolExecutor(max_workers=50) as ex:
            sza_worker = np.asarray([ex.submit(_sza, glon.ravel()[i], glat.ravel()[i]) for i in range(glon.size)])

        sza = np.nan*np.ones(glon.ravel().size)
        for i in range(sza_worker.size):
            sza[i] = sza_worker[i].result()
        sza = sza.reshape(glon.shape)
    
    elif isinstance(time, np.ndarray) and not isinstance(glon, np.ndarray):
        with concurrent.futures.ThreadPoolExecutor(max_workers=50) as ex:
            sza_worker = np.asarray([ex.submit(_sza_time, time[i]) for i in range(time.size)])

        sza = np.nan*np.ones(time.size)
        for i in range(sza_worker.size):
            sza[i] = sza_worker[i].result()
            
    elif isinstance(time, np.ndarray) and isinstance(glon, np.ndarray) and isinstance(glat, np.ndarray):
        sza = np.nan * np.ones(time.size)
        for i in range(time.size):
            sza[i] = _sza_all(time[i], glon[i], glat[i])
            
    else:
        sza = _sza_time(time)
    return sza

def get_angles(time, glon, glat, ghgt=0):
    
    def _angles(glon, glat):
        sun, moon = objects(time,glon,glat,ghgt)
        sun_moon_sep = separation(sun.az, sun.alt, moon.az, moon.alt)
        sun_moon_azimuth = azimuth(sun.az, sun.alt, moon.az, moon.alt)
        return sun_moon_sep,sun_moon_azimuth,moon.radius
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=50) as ex:
        angles_worker = np.asarray([ex.submit(_angles, glon.ravel()[i], glat.ravel()[i]) for i in range(glon.size)])

    sun_moon_separation_radian = np.nan*np.ones(glon.ravel().size)
    sun_moon_azimuth_radian = np.nan*np.ones(glon.ravel().size)
    moon_radius_radian = np.nan*np.ones(glon.ravel().size)
    for i in range(angles_worker.size):
        sun_moon_separation_radian[i] = angles_worker[i].result()[0]
        sun_moon_azimuth_radian[i] = angles_worker[i].result()[1]
        moon_radius_radian[i] = angles_worker[i].result()[2]
    
    return sun_moon_separation_radian.reshape(glon.shape), sun_moon_azimuth_radian.reshape(glon.shape), moon_radius_radian.reshape(glon.shape)

def get_parallactic_angle(time, glon, glat, ghgt):
    
    def _eta(glon, glat):
        sun, moon = objects(time,glon,glat,ghgt)
        return parallactic_angle(sun.az, sun.dec, glat)
    
    def _eta_times(t):
        sun, moon = objects(t, glon, glat, ghgt)
        return parallactic_angle(sun.az, sun.dec, glat)
    
    if isinstance(glon, np.ndarray) and isinstance(time, datetime):
        with concurrent.futures.ThreadPoolExecutor(max_workers=50) as ex:
            eta_worker = np.asarray([ex.submit(_eta, glon.ravel()[i], glat.ravel()[i]) for i in range(glon.size)])
        
        eta = np.nan*np.ones(glon.ravel().size)
        for i in range(eta_worker.size):
            eta[i] = eta_worker[i].result()
        
        return eta.reshape(glon.shape)
    
    elif isinstance(time, np.ndarray):
        eta = np.nan*np.ones(time.size)
        with concurrent.futures.ThreadPoolExecutor(max_workers=50) as ex:
            eta_worker = np.asarray([ex.submit(_eta_times, time[i]) for i in range(time.size)])
        for i in range(eta_worker.size):
            eta[i] = eta_worker[i].result()
        return eta
    else:
        print ('a)')
        eta = _eta(glon, glat)
        return eta
def get_eof_mask_from_angles(image, sep, azm, eta, mrad, x0, y0, imres, pixscale):
    mx0, my0 = rotate(sep, azm-eta, 0.0, 0.0)
    mask = moon_mask(imres, mx0*pixscale + x0, my0*pixscale + y0, np.round(mrad, 8)*pixscale)
    return np.nansum(np.multiply(image,mask)) / np.nansum(image)

def get_times(t0,t1,dm=0, ds=0):
    assert (dm>0) or (ds>0)
    times =[]
    while t0 <= t1:
        times.append(t0)
        t0 += timedelta(minutes = dm, seconds=ds)
    return np.array(times)

def moon_mask(N, mx1, my1, r):
    xaxis = np.arange(N) - mx1
    yaxis = np.arange(N) - my1
    mask = ((np.zeros((N, N)) + yaxis*yaxis).T + xaxis*xaxis) < r*r
    return ~mask

def horizon_mask(horizon, pixscale, selv, y0, imsdo):
    if np.degrees(horizon) > -3:
        hmask = np.ones(imsdo.shape[0])
        for i in range(hmask.size):
            tmp = selv + ((i - int(y0)) / pixscale)
            if tmp < horizon:
                hmask[i] = 0
        hmask = np.array([hmask,]*imsdo.shape[0]).T
    else:
        hmask = np.ones(imsdo.shape)
    return hmask

def parallactic_angle(sazm, sdec, glat):
    sineta = np.sin(sazm) * np.cos(np.radians(glat)) / np.cos(sdec)
    return -np.arcsin(sineta)

def rotate(sep, azm, x0, y0):
#    (elv,azm,olon,olat,xlon,xlat)
    sint = np.sin(sep)
    cost = np.cos(sep)
    sinp = np.sin(azm)
    cosp = np.cos(azm)
    sinl = np.sin(y0)
    cosl = np.cos(y0)
    
    y1 = np.arcsin( sinl*cost + cosl*sint*cosp )
    x1 = x0 + np.arctan2( sint*sinp, cosl*cost - sinl*sint*cosp )
    
    return x1, y1

def objects(T, glon, glat, ghgt):
    assert isinstance(T, datetime)
    observer = ephem.Observer()
    observer.lon = str(glon)
    observer.lat = str(glat)
    observer.elevation = ghgt
    observer.date = T
    observer.pressure = 0
    observer.horizon = -np.arccos(RE / (RE + ghgt/1e3))
    
    sun = ephem.Sun(observer)
    moon = ephem.Moon(observer)
    
    return sun, moon

def azimuth(sazm, selv, mazm, melv):
    coslt1 = np.cos(melv)
    sinlt1 = np.sin(melv)
    
    coslt0 = np.cos(selv)
    sinlt0 = np.sin(selv)
    
    cosl0l1 = np.cos(mazm-sazm)
    sinl0l1 = np.sin(mazm-sazm)
    
    cosc = (sinlt0 * sinlt1) + (coslt0 * coslt1 * cosl0l1)  # Cos(a to b)
    sinc = np.sqrt(1 - (cosc*cosc))
    
    if (abs(sinc) > 1e-7):# Small angle?
        cosaz = (coslt0 * sinlt1 - sinlt0 * coslt1 * cosl0l1) / sinc  # Azimuth
        sinaz = (sinl0l1 * coslt1)/sinc
    else:   # It is antipodal
        cosaz = 1
        sinaz = 0
    
    azm = np.arctan2(sinaz, cosaz)
    return azm

def separation(sazm, selv, mazm, melv):
    return np.round(ephem.separation((sazm, selv), (mazm, melv)), 8)

def get_EOF(sr, mr, mx0, my0):
    
    sx0 = 0
    sy0 = 0
    
    if mr > sr:
        r1 = mr
        r2 = sr
    else:
        r1 = sr
        r2 = mr
    d = np.sqrt(abs(sx0-mx0)**2 + abs(sy0-my0)**2)
    if d > (r1+r2)*1.05:
        of = 1
    elif d <= (r1-r2):
        of = 1 - ((np.pi * mr**2) / (np.pi * sr**2))
    else:
        d1 = (r1**2 - r2**2 + d**2) / (2*d)
        d2 = d - d1
        
        A = ( r1**2 * np.arccos(d1/r1) - (d1 * np.sqrt(r1**2 - d1**2)) ) + \
            ( r2**2 * np.arccos(d2/r2) - (d2 * np.sqrt(r2**2 - d2**2)) )

        of = 1 - (A / (np.pi * sr**2)) if A > 0 else 1
    
    return of

#def eclipse_geo_novas(T, glon, glat, ghgt=0, srad_fact=1, plot=0):
#    assert isinstance(T, datetime)
#    hour = T.hour + T.minute / 60 + T.second / 3600
#    sazm,selv,srad,sdec,mazm,melv,mrad,azm,sm_dist = eomff.eom(T.year,T.month,T.day,hour,glon,glat,ghgt,srad_fact)
#    
#    horizon = (-np.arccos(RE / (RE + ghgt/1e3)) - selv - srad) #/ srad
#    if horizon > 0:
#        return 0
#    else:
#        mx0, my0 = eomff.rotate(sm_dist, azm, 0.0, 0.0)
#        of = get_EOF(srad, mrad, mx0, my0)
#        return of
    
def mask_geo_ephem(T, glon, glat, ghgt=0, srad_fact=1, plot=0):
    sun, moon = objects(T, glon, glat, ghgt)
    srad = np.round(sun.radius * srad_fact, 8)
    mrad = np.round(moon.radius, 8)
    horizon = (-np.arccos(RE / (RE + ghgt/1e3)) - sun.alt - sun.radius) #/ srad
    if horizon > 0:
        return 0
    else:
        sm_dist = np.round(ephem.separation((sun.az, sun.alt), (moon.az, moon.alt)), 8)
        azm = azimuth(sun.az, sun.alt, moon.az, moon.alt)
        mx0, my0 = rotate(sm_dist, azm, 0.0, 0.0)
        of = get_EOF(srad, mrad, mx0, my0)
        return of

def mask_sdo_ephem(T, glon, glat, ghgt, x0, y0, imsdo, pixscale):
    sun, moon = objects(T,glon,glat,ghgt)
    horizon = (-np.arccos(RE / (RE + ghgt/1e3)) - sun.alt - sun.radius)
    sep = separation(sun.az, sun.alt, moon.az, moon.alt) 
    
    if horizon*pixscale >= (imsdo.shape[0]/2-imsdo.shape[0]):
        hmask = horizon_mask(horizon=horizon, selv=sun.alt, imsdo=imsdo, pixscale=pixscale, y0=y0)
    else:
        hmask = np.ones_like(imsdo)
    eta = parallactic_angle(sun.az, sun.dec, glat)
#    if (sep*pixscale) < (imsdo.shape[0]*1.4142+np.round(moon.radius, 8)*pixscale):
    if (sep*pixscale) < (imsdo.shape[0]+np.round(moon.radius, 8)*pixscale):
        
        azm = azimuth(sun.az, sun.alt, moon.az, moon.alt)
        # Rotation of the moon if ~100x faster than rotation of the Sun for the parallactinc angle
        # ndimage.rotate(imsdo, np.rad2deg(eta), reshape=False) --- takes about 3seconds to compute
        mx0, my0 = rotate(sep, azm-eta, 0.0, 0.0)
        mmask = moon_mask(imsdo.shape[0], mx0*pixscale + x0, my0*pixscale + y0, np.round(moon.radius, 8)*pixscale)
        mask = np.multiply(hmask, mmask)
        of =  np.nansum(np.multiply(imsdo,mask)) / np.nansum(imsdo)
    else:
        of =  np.nansum(np.multiply(imsdo,hmask)) / np.nansum(imsdo)
    return of

def eof_time(t0, t1, glon, glat, ghgt, srad_fact=1, dm=10, ds=0, mode='ephem'):
    times = get_times(t0,t1,dm=dm, ds=ds)
    if not isinstance(srad_fact, list):
        srad_fact = [srad_fact]
    
    OF = np.nan * np.zeros((times.size, len(srad_fact)))
    
    for sf, srad_factor in enumerate(srad_fact):
        for i,T in enumerate(times):
#            if mode == 'novas':
#                OF[i, sf] = eclipse_geo_novas(T, glon=glon, glat=glat, ghgt=ghgt, srad_fact=srad_factor)
#            else:
            OF[i, sf] = mask_geo_ephem(T, glon=glon, glat=glat, ghgt=ghgt, srad_fact=srad_factor)
    OF[OF<0] = 0
    return times, np.squeeze(OF)

def mask_latalt_geo(T, glat, ghgt, glon0=0, srad_fact=1, mode='novas'):
    assert isinstance(ghgt, np.ndarray)
    assert isinstance(glat, np.ndarray)
    assert isinstance(T, datetime)
    if not isinstance(srad_fact, list):
        srad_fact_list = [srad_fact]
    else:
        srad_fact_list = srad_fact
        
    OF = np.zeros((glat.size, ghgt.size, len(srad_fact_list)))
    for s, srad_fact in enumerate(srad_fact_list):
        for i, lat in enumerate(glat):
            for j, alt in enumerate(ghgt):
#                if mode == 'novas':
#                    OF[i,j,s] = eclipse_geo_novas(T=T, glon=glon0, glat=lat, ghgt=alt, srad_fact=srad_fact)
#                else:
                OF[i,j,s] = mask_geo_ephem(T=T, glon=glon0, glat=lat, ghgt=alt, srad_fact=srad_fact)
    return np.squeeze(OF)

def mask_lonlat_geo(T, glon, glat, ghgt=0, srad_fact=1, mode='novas'):
    
    assert isinstance(glon, np.ndarray)
    assert isinstance(glat, np.ndarray)
    assert isinstance(T, datetime)
    
    OF = np.zeros((glon.size, glat.size))
    for i, lon in enumerate(glon):
        for j, lat in enumerate(glat):
#            if mode == 'novas':
#                OF[i,j] = eclipse_geo_novas(T=T, glon=lon, glat=lat, ghgt=ghgt, srad_fact=srad_fact)
#            else:
            OF[i,j] = mask_geo_ephem(T=T, glon=lon, glat=lat, ghgt=ghgt, srad_fact=srad_fact)
    return np.squeeze(OF)

#%% SDO AIA
def eof_time_sdo(SDO, t0, t1, glon, glat, ghgt, wl=193, dm=10,ds=0):

    times = get_times(t0,t1,dm=dm,ds=ds)
    OF = np.ones(times.size)
    imsdo = SDO['AIA{}'.format(wl)].values
    if glat < 0:
        imsdo = ndimage.rotate(imsdo, 180)
    for i,T in enumerate(times):
        if (i+1)%10 == 0:
            print ("Processing {}/{}".format(i+1, times.size))
        OF[i] = mask_sdo_ephem(T, glon, glat, ghgt, SDO.x0, SDO.y0, imsdo, SDO.pixscale)
    return times, OF

def mask_lonlat_sdo(SDO, T, glon, glat, ghgt, wl=193, verbose=False):
    assert isinstance(glon, np.ndarray)
    assert isinstance(glat, np.ndarray)
    OF = np.ones((glon.size, glat.size))
    if verbose:
        c = 1
        C = glon.size*glat.size
    for i in range(glon.size):
        if verbose:
            t0 = datetime.now()
        for j in range(glat.size):
            if verbose:
                if c%1000 == 0:
                    print ("{}/{}".format(c,C))
            OF[i,j] =  mask_sdo_ephem(T, glon[i], glat[j], ghgt=ghgt, x0=SDO.x0, y0=SDO.y0, imsdo=SDO['AIA{}'.format(wl)].values, pixscale=SDO.pixscale)
            if verbose:
                c+=1
        if verbose:
            print(datetime.now()-t0)
    return OF

def mask_latalt_sdo(SDO, T, glon, glat, ghgt, wl=193, verbose=False):
    assert isinstance(ghgt, np.ndarray)
    assert isinstance(glat, np.ndarray)
    assert isinstance(glon, (int, float))
    
    OF = np.ones((glat.size, ghgt.size))
    if verbose:
        c = 1
        C = glon.size*glat.size
    for i in range(glat.size):
        if verbose:
            t0 = datetime.now()
        for j in range(ghgt.size):
            if verbose:
                if c%1000 == 0:
                    print ("{}/{}".format(c,C))
            OF[i,j] =  mask_sdo_ephem(T, glon, glat[i], ghgt=ghgt[j], x0=SDO.x0, y0=SDO.y0, imsdo=SDO['AIA{}'.format(wl)].values, pixscale=SDO.pixscale)
            if verbose:
                c+=1
        if verbose:
            print(datetime.now()-t0)
    return OF

# Spacecraft
def eof_satellite(times, glon, glat, ghgt, SDO=None, srad=1.0, wl='geo', verbose=False):
    OF = np.zeros(times.size)
    for i,T in enumerate(times):
        if verbose:
            if (i+1)%10 == 0:
                print ("Processing {}/{}".format(i+1, times.size))
        sza = get_sza(T, glon[i], glat[i], alt_km=ghgt[i]*1e3)
        if sza < 95:
            if SDO is not None:
                OF[i] = mask_sdo_ephem(T, glon[i], glat[i], ghgt[i]*1e3, SDO.x0, SDO.y0, SDO['AIA{}'.format(wl)].values, SDO.pixscale)
            else:
                OF[i] = mask_geo_ephem(T, glon[i], glat[i], ghgt[i]*1e3, srad_fact=srad)
    return times, OF