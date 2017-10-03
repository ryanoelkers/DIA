#this program will do the photometry on the reference frame to get stars 
#for the subtraction and to get stars for later

#if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50

#import the relevant libraries for basic tools
import numpy
import scipy
from scipy import stats
import scipy.ndimage as ndimage
import astropy
from astropy.stats import sigma_clipped_stats
import math
import time
from photutils import DAOStarFinder
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils import CircularAnnulus
import pyfits

#for reading in fits files
from astropy.io import fits

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists

#import relevant spline libraries
from scipy.interpolate import Rbf

###BEGIN UPDATE#####
#what is the exposure time?
exp_time = 60.0

#optimal radius
rad = 3.0

#useful directories
cdedir = '../Python/' #code directory
caldir = '../calib/' # directory to output the flux files and location of master frame

##END UPDATE####

#read in the master frame
mlist = fits.open(caldir+'master.fits')
mheader = mlist[0].header #get the header info
mast = mlist[0].data #get the image info
exp_time = pyfits.getval(caldir+'master.fits','EXPTIME')

#get the positions in the master
mean, median, std = sigma_clipped_stats(mast, sigma = 3.0, iters = 5)
daofind = DAOStarFinder(fwhm = 1.5, threshold = 5.*std)
sources = daofind(mast-median)

#set up for aperture photometry on the master frame
mean, median, std = sigma_clipped_stats(mast, sigma = 3.0, iters = 5)
positions = (sources['xcentroid'], sources['ycentroid'])

print 'Getting the photometry on the master frame.'

#do the aperture photometry
apertures = CircularAperture(positions, r = rad)
phot_table = aperture_photometry(mast, apertures, method = 'exact')

#set up for the background annuli
annulus_apertures = CircularAnnulus(positions, r_in = rad+2, r_out = rad+4)
bkgf_table = aperture_photometry(mast, annulus_apertures, method = 'exact')
bkg_mean = bkgf_table['aperture_sum']/annulus_apertures.area()

#convert to magnitudes
flx = phot_table['aperture_sum']-(bkg_mean*apertures.area())
flx_er = numpy.sqrt(phot_table['aperture_sum'])
x_pix = sources['xcentroid']
y_pix = sources['ycentroid']

#create the magnitudes from the flux
mag = 25.0-2.5*numpy.log10(flx)+2.5*numpy.log10(exp_time)
err = (2.5/numpy.log(10.))*(flx_er/flx)

#write the magnitudes to a file
output = open(caldir+'master.ap', 'w')
for ii in range(0, len(phot_table['id'])):
	output.write(str(ii)+','+str(x_pix[ii])+','+str(y_pix[ii])+','+str(mag[ii])+','+str(err[ii])+'\n')
output.close()

#write the fluxes to a file
output = open(caldir+'master.flux', 'w')
for ii in range(0, len(phot_table['id'])):
	output.write(str(ii)+','+str(x_pix[ii])+','+str(y_pix[ii])+','+str(flx[ii])+','+str(flx_er[ii])+'\n')
output.close()

#write the star list to a file
output = open(caldir+'starlist.txt', 'w')
for ii in range(0, len(phot_table['id'])):
	output.write(str(ii)+','+str(long(x_pix[ii]))+','+str(long(y_pix[ii]))+'\n')
output.close()
