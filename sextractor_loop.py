from astropy.table import QTable
import numpy as np
from astropy.io import ascii
from astropy.table import Table
import numpy as np           # to define our table
from astropy.io import fits  # for fits image
import bz2
import numpy.ma as ma
from astropy import wcs
from pylab import *
import astropy.units as u
#from photutils.segmentation import SourceFinder
from photutils.morphology import data_properties
import sys

from astropy.wcs import WCS

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time

import astroscrappy

# Disable some annoying warnings from astropy
import warnings
from astropy.wcs import FITSFixedWarning
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter(action='ignore', category=FITSFixedWarning)
warnings.simplefilter(action='ignore', category=VerifyWarning)
# Also ignore FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)


from stdpipe import astrometry, photometry, catalogs, cutouts, templates, subtraction, plots, psf, pipeline, utils


d= "pairs_sextractor.dat"
data = np.loadtxt(d,  dtype={'names':('Index','ID','RA','Dec'),'formats': ('i','O','f', 'f')}, skiprows=1)

# then define what each column has
Index= data['Index']
ID=data['ID']
RA= data['RA']
Dec=data['Dec']

### galaxies that couldn´t be run ## 84,85,86,91,92,95,96,97,98

for i in range(len(ID)):
        try:
            target = ID[i]
            print ("Iter %s..." % Index[i])


        #target_url ='http://das.sdss.org/raw/4670/40/corr/5/fpC-004670-r5-0102.fit.gz'
        #dat = fits.open(target_url)

            hdu=fits.open(ID[i])
            dat = hdu[0].data

            filename = ID[i]
            image = fits.getdata(filename).astype(np.double)
            header = fits.getheader(filename)

            time = utils.get_obs_time(header, verbose=False)
            fname = header.get('FILTER')
            gain = header.get('GAIN')

            #print('Processing %s: filter %s gain %.2f at %s' % (filename, fname, gain, time))

            mask = image > 50000

            cmask, cimage = astroscrappy.detect_cosmics(image, mask, verbose=True)
            print('Done masking cosmics: %d pixels masked' % np.sum(cmask))
            mask |= cmask

            #plots.imshow(image)
            #plt.title('Pre-processed image');
            #plt.show()

            #plots.imshow(mask)
            #plt.title('Mask');
            #plt.show()





            # We will detect objects using SExtractor and get their measurements in apertures with 3 pixels radius
            obj = photometry.get_objects_sextractor(image, mask=mask, aper=3.0, gain=gain, edge=10)
            print(len(obj), 'objects detected')

            # Let's see 3 brightest objects
            print(obj[:3])

            # Rough estimation of average FWHM of detected objects, taking into account only unflagged (e.g. not saturated) ones
            fwhm = np.median(obj['fwhm'][obj['flags'] == 0])
            print('Average FWHM is %.1f pixels' % fwhm)

            ### Average FWHM is 4.5 pixels #####

            # We will pass this FWHM to measurement function so that aperture and background radii will be relative to it.
            # We will also reject all objects with measured S/N < 5
            obj = photometry.measure_objects(obj, image, mask=mask, fwhm=fwhm, gain=gain, aper=1.0, bkgann=[5, 7], sn=5, verbose=True)
            print(len(obj), 'objects properly measured')

            #print('obj:', obj.columns)
            ### 424 objects properly measured ###

            #####obj: <TableColumns names=('mag','magerr','flux','fluxerr','x','y','xerr','yerr','a','b','theta',
            #####'FLUX_RADIUS','fwhm','flags','bg','ra','dec','bg_local')>


            #print(*obj)

            #sys.exit()



            ###In our case, FITS header already contains WCS solution, so we will just load it.

            # Load initial WCS
            wcs = WCS(header)

            # Get the center position, size and pixel scale for the image
            center_ra,center_dec,center_sr = astrometry.get_frame_center(wcs=wcs, width=image.shape[1], height=image.shape[0])
            pixscale = astrometry.get_pixscale(wcs=wcs)

            print('Frame center is %.2f %.2f radius %.2f deg, %.2f arcsec/pixel' % (center_ra, center_dec, center_sr, pixscale*3600))
########Catalogue name may be any Vizier identifier, or one of supported shortcuts for popular choices (ps1, gaiadr2, gaiaedr3, usnob1, gsc, skymapper, apass, sdss, atlas, vsx etc).##
            #Let's get PanSTARRS objects brighter than r=18 mag
            cat = catalogs.get_cat_vizier(center_ra, center_dec, center_sr, 'sdss', filters={'rmag':'<22.5'})
            print(len(cat), 'catalogue stars')

            print (cat[:3])

            col_name = cat.columns
            #print('tb:', col_name)

            # Let's use SCAMP for astrometric refinement.
            wcs = pipeline.refine_astrometry(obj, cat, 5*pixscale, wcs=wcs, method='astropy', cat_col_mag='rmag', verbose=True)

            if wcs is not None:
                # Update WCS info in the header
                astrometry.clear_wcs(header, remove_comments=True, remove_underscored=True, remove_history=True)
                header.update(wcs.to_header(relax=True))

            m = pipeline.calibrate_photometry(obj, cat, sr=2/3600, cat_col_mag='rmag', cat_col_mag1='gmag', cat_col_mag2='rmag', max_intrinsic_rms=0.02, order=2, verbose=True)


            print(obj[:5])



            RA2 = obj['ra']
            DEC2 = obj['dec']
            x0 = obj['x']
            y0 = obj['y']
            a = obj['a']
            b = obj['b']
            pa = obj['theta']
            mag = obj['mag_calib']
            mag_error = obj['mag_calib_err']

            c1 = SkyCoord(RA[i], Dec[i], frame='fk5', unit=(u.degree, u.degree))
            f_out = open('%s.dat'%Index[i], 'w')
            print ("# Index ID RA DEC x0  y0  a b pa mag  magerr", file= f_out)
            for j in range (len(RA2)):
                c2 = SkyCoord(RA2[j], DEC2[j], frame='fk5', unit=(u.degree, u.degree))
                angular_sep = c1.separation(c2)  # Differing frames handled correctly
                if angular_sep.arcsec < 2:
                    print('as:', angular_sep.arcsec)
                    print('c2:', c2)
                    print(Index[i],ID[i],RA2[j], DEC2[j], x0[j], y0[j], a[j], b[j], pa[j], mag[j], mag_error[j],file=f_out)

        except:
            pass




        #t= (RA, DEC, x0, y0, a, b, pa)
        #ascii.write(t, '%s.csv'%filename[:-7] , overwrite=False)


        #for i in range(len(RA)):

            #f_out = open('output%s.dat' %Index, 'w')
            #print ("# RA DEC x y a b theta", file=f_out)

            #print(RA[i], DEC[i], x0[i], y0[i], a[i], b[i], pa[i], file=f_out)

    #t = QTable(RA, DEC, x0, y0, a, b, pa)

    #print ('t:', t)
