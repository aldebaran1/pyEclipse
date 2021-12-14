# pyEclipse

PyEclipse is a pythonic software framework for computing solar eclipses. The unique feautes of this framework is a computation of eclipse penumbra at wavelengths captured by solar telescopes at various wavelengths. Momenteraly, we interface with Virtual Solar Observatory, using SunPy library, to retrieve Solar Dynamics Observatory (SDO) Atmospheric Imaging Assembley (AIA) images at Extreme UltraViolet (EUV) Wavelengths. SDO AIA images are high-resolution (4096x4096) and high cadence (1 minute per wavelength). *[Other telescopes can be added to the package, for instance SOHO EIT, GOES SXI, as loong as they are privided the metadata described below]*   
SDO AIA EUV telescolpe images coronal emissions at: 
- 94A (Angstrom, 10 A = 1 nm)
- 131A
- 171A
- 193A
- 211A
- 304A
- 335A
For more information about these coronal lines refer to Lemen et al (2012) [a].

To download desired SDO AIA images in .fits format, 

`python dlSunImage.py  2017-8-21T18:00  2017-8-21T18:15 --wl 335 --odir /path/to/your/directory/`

This command will download AIA images at 335A taken between 18 and 18:10 UT to you output directory. If you want to download just the first image from the batch add `--first` flag in the command.

These level 1 image `.fits` files are annotated. We utilze *Astropy* to read `fits` files as `astropy.io.fits.open(filname.fits)` and read the important metadata from the header:

`timestamp = Data.header[T_REC] `   
`center_x = Data.header[CRPIX1] `   
`center_y = Data.header[CRPIX2] `   
`pixel2arcsec = Data.header[CDELT1]`   
`wavelength = Data.header[WAVELNTH]`

### Eclipse -- Observer's Point of View:

Compute the visual eclipse transion for an observer's geolocation: *glon, glat, galt*. The Eclipse POV can be computed using SDO AIA images or assuming a uniform Sun with adjustable solar radii. POV Eclipse images are computed in a horizontal coordinate systems, also known as Altitude/azimute system. In this orientation, obserers pespectice is fixed in time and geolocation causing apparent roration of the Sun due to an (changing) angle between observer's zenith and Sun's north. This sngle is known as the parallactic angle. The parallactic angle was verificated against Peater Meadows routine and software *Helio*[1].  

`python eclipse_animation.py 2017-8-21T16:00 2017-8-21T20:00 -100 40.7 0 --dt 5 --wl 171 --odir /path/to/your/directory/pov/ --animation`     
`--animation` creates a .mp4 file. This functionality requres *py-lapse* package [2]    
`--wl` specifies the desired wavelength. Choose `geo` to create POV for uniform eclipse.     
`--srad` spacifies the solar radii for uniform eclipse. Defualt is 1.0. Set it to 1.15, that is inflating the solar radii by 15% to get maximum elcipse of ~0.9. Depending on solar activity.

Example of a POV EUV eclipse at 19.3 nm doe 2017 Eclipse:
![a](https://github.com/aldebaran1/pyEclipse/blob/master/misc/Aug2017_pov.gif)

[a] Lemen et al. (2012), Doi:10.1007/s11207-011-9776-8     
[1] https://www.petermeadows.com/html/parallactic.html    
[2] https://github.com/aldebaran1/pyEclipse/py-lapse.git    
