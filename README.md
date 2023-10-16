# pyEclipse

[![DOI](https://zenodo.org/badge/362197547.svg)](https://zenodo.org/badge/latestdoi/362197547)

This model was suported by National Science Fundation grant AGS-13019141 and by NASA grant 13018916. Additon of Hinode XRT is supported by AFOSR under grant 23RT1097 to CU Boulder and JHU APL.


PyEclipse is a pythonic software framework for computing solar eclipses. The unique feautes of this framework is a computation of eclipse penumbra at wavelengths captured by solar telescopes at various wavelengths. Momenteraly, we interface with Virtual Solar Observatory, using SunPy library, to retrieve Solar Dynamics Observatory (SDO) Atmospheric Imaging Assembley (AIA) images at Extreme UltraViolet (EUV) Wavelengths. SDO AIA images are high-resolution (4096x4096) and high cadence (1 minute per wavelength). We can download SOHO EIT as well using SunPy interface, however, the data quality is very poor. PyEclipse has a custon interface to noaa.ncei to download GOES SUVI L1b and L2 data. Only L2 data is good enough quality for precise eclipse modeling.*[Other telescopes can be added to the package, for instance \, GOES SXI, as loong as they are privided the metadata described below]*   

SDO AIA EUV telescolpe images coronal (and chromosphere) emissions at: 
- 94A (Angstrom, 10 A = 1 nm)
- 131A
- 171A
- 193A
- 211A
- 304A
- 335A    
- 1600 A (continuum line from the photosphere/transition region)
For more information about these coronal lines refer to Lemen et al (2012) [a]. The original works used NOVAS Fortran routines developed by Douglas Drob (NRL). The results wee published in Geophysical Research Letters [b,c,d]. This software framework is purely pythonic, all NOVAS Fortran routines are re-written in Python. The 2-D maps take a long time to compute, so the tedious *for*-loops are parallelized with `concurrent.futures.ThreadPoolExecutor()`. Running more scirpts in parallel will jam your CPUs [Just saying].

GOES-R SUVI telescope images coronal emissions at:
- 94A
- 131A
- 171A
- 195A
- 285A
- 304A

We are accressing the SUVI data via NOAA NCEI. This is problematic becuase NCEI has numerous paths to parts of the data at different data levels (from L1 to L2) and differnt data availability. The latest version of pyEclipse connects to SUVI L2 data whihc appears to be available for 2022 and 2023.

TO DO: add capability to download Honode XRT synoptic maps from https://solar.physics.montana.edu/HINODE/XRT/SCIA/synop_official/
Presetnly, the XRT images need to be downlaoded manually, but pyEclipse reads the data with eio.load('hinode_image.fits')

To download desired SDO AIA images in .fits format, 

`python dlSunImage.py  2017-8-21T18:00  2017-8-21T18:15 --wl 335 --odir /path/to/your/directory/`

This command will download AIA images at 335A taken between 18 and 18:10 UT to you output directory. If you want to download just the first image from the batch add `--first` flag in the command.

These level 1 image `.fits` files are annotated. We utilze *Astropy* to read `fits` files as `astropy.io.fits.open(filname.fits)` and read the important metadata from the header:

`timestamp = Data.header[T_REC] `   
`center_x = Data.header[CRPIX1] `   
`center_y = Data.header[CRPIX2] `   
`pixel2arcsec = Data.header[CDELT1]`   
`wavelength = Data.header[WAVELNTH]`

The model takes into acount the altitude of an observer and a horizon, which also changes with altitude. This altitude dependance assumes *the Earth as a sphere*, which limits the model's accuracy at higher altitudes >~300 km and high solar zenith angles. Atmospehric refraction is parameterazed in pyEphem as a funcition of temperature, pressure and altitude. All these parameters are set default.

### Setup
Install the software framework  
`python setup.py install`

Import the functionalities in your python codes as   
`from eclipse import utils,eio`

### Eclipse -- Observer's Point of View:

Compute the visual eclipse transion for an observer's geolocation: *glon, glat, galt*. The Eclipse POV can be computed using SDO AIA images or assuming a uniform Sun with adjustable solar radii. POV Eclipse images are computed in a horizontal coordinate systems, also known as Altitude/azimute system. In this orientation, obserers pespectice is fixed in time and geolocation causing apparent roration of the Sun due to an (changing) angle between observer's zenith and Sun's north. This sngle is known as the parallactic angle. The parallactic angle was verificated against Peater Meadows routine and software *Helio*[1].  

`python eclipse_animation.py 2017-8-21T16:00 2017-8-21T20:00 -100 40.7 0 --dt 5 --wl 171 --odir /path/to/your/directory/pov/ --animation`     
`--animation` creates a .mp4 file. This functionality requres *py-lapse* package [2]    
`--wl` specifies the desired wavelength. Choose `geo` to create POV for uniform eclipse.     
`--srad` spacifies the solar radii for uniform eclipse. Defualt is 1.0. Set it to 1.15, that is inflating the solar radii by 15% to get maximum elcipse of ~0.9. Depending on solar activity.

Example of a POV EUV eclipse at 19.3 nm doe 2017 Eclipse:
![a](https://github.com/aldebaran1/pyEclipse/blob/master/misc/Aug2017_pov.gif)

### Eclipse time series

Compute the eclipse occultation as a function of time for a desired geographic location (glon, glat, galt) and eclipse mode. Choose between SDO AIA wavelengths you downloaded previously or GEOmetricaly uniform sun with adjustable solar radii. The timeseries profile is saved into a netCDF file by default.

`
pyhton mask_time_cli.py 2017-8-21T16:00 2017-8-21T19:45 -100 40.8 0 /path/to/output/directory/mask.nc --altkm 0 --tres 60 --sdodir /path/to/sdoaia/images/ --wl 94 193 211 geo --plot --srad 1.15
`

`--plot` option plots the output in the console  
`--save` option does not save the output to data  
`--wl` takes >1 argument. `--wl 94 193 211 geo` computes 4 profiles in the given order

Example of a time series for August 2017 eclipse for parameters in the above command.
![b](https://github.com/aldebaran1/pyEclipse/blob/master/misc/time_series_2017.png)

### Eclipse mask 2-D / longitude-latitude

Compute a lat/lon eclipse mask at a given altitude and time. Routines for computeing SDO-AIA and geometric masks have their own scripts. An eclipse mask for a uniform eclipse is computed pretty fast, it takes ~30 minutes to compute mask for a duration of an eclipse with a 5 minute time resolution. Computation od SDO AIA mask is considerably slower, becuase the Eclipse occultation funcion is computed as a ratio between the sum of all (unocculted) SDO AIA pixels devided by occulted SDO AIA pixels. This requres computation on 2D matrices of a size of 4096x4096 = 16.8M pixels for each observer's location. We make this computation only if a horizon or the moon apperars on the image. Othewise the finction returns 1 for SZA<90 and 0 for SZA>90. Horizon is continous.

Compute a 2D eclipse mask at 18:00UT on the 21th August 2017 using 211A image in /path/to/sdoaia/images/aia_filename.fits. Output mask is saved as a netCDF file.     
`python sdomask_lonlat_cli.py 2017-8-21T18:00 2017-8-21T18:01 /path/to/output/aiamask/ --sdoaia /path/to/sdoaia/ --dlon 1 --dlat 1 --wl 211 --altkm 125`    
`--dlon, --dlat` argumets set lon and lat resolution in degrees    
`--tres` set at time resolution to compute masks between *time1 time2* arguments.
`--sdodir` specify the directory where you have you sdo aia images  

For a geometrically symmetric eclipse masks type:     
`python geomask_lonlat_cli.py 2017-8-21T18:00 2017-8-21T18:01 /path/to/output/geomask/ --sdoaia /path/to/sdoaia/ --dlon 1 --dlat 1 --sunradii 1.15 --altkm 125`    

Here are two examples:   
SDO AIA 94A Penumbra at 150 km   
![c](https://github.com/aldebaran1/pyEclipse/blob/master/misc/aug2017_aia94A.gif)      
Uniform eclipse with inflated solar radii by 10%    
![d](https://github.com/aldebaran1/pyEclipse/blob/master/misc/Aug2017_geo1.1.gif)      

### Lateral cuts thru a penumbra: Lat-alt projection

Compute a lateral cut tru a penumbra to analyze height dependence: tilt in the shadow as a function of latitude, spatial irregularities, etc. This functionality is not a command line interface, you need to adjust the input arguments in the code `geomask_latalt.py`

Here is an example profile:    
![e](https://github.com/aldebaran1/pyEclipse/blob/master/misc/lat_alt_2017.png)

### References
All details available in the following associated paper. Please cite this software and methodology as:
Mrak, S., Zhu, Q., Deng, Y., Dammasch, I. E., Dominique, M., Hairston, M. R., Nishimura, Y., & Semeter, J. (2022). Modeling Solar Eclipses at Extreme Ultra Violet Wavelengths and the Effects of Nonuniform Eclipse Shadow on the Ionosphere-Thermosphere system. Journal of Geophysical Research: Space Physics, e2022JA031058. https://doi.org/10.1029/2022JA031058


[a] Lemen et al. (2012), Doi:10.1007/s11207-011-9776-8     
[b] Huba nd Drob (2017), Doi:10.1002/2017GL073549    
[c] Mrak et al. (2018), Doi:10.1029/2017GL076771     
[d] Hairston et al. (2018), Doi:10.1029/2018GL077381    
[1] https://www.petermeadows.com/html/parallactic.html    
[2] https://github.com/aldebaran1/pyEclipse/py-lapse.git    
