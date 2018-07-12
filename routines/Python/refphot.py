#this program will do the photometry on the reference frame to get stars 
#for the subtraction and to get stars for later

#if you use this code, please cite Oelkers & Stassun 2018

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

#for reading in fits files
from astropy.io import fits
from astropy.wcs import WCS

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join, exists

###BEGIN UPDATE#####
camera = '2'
ccd = '2'

#useful directories
cdedir = '../code/master/' #code directory
caldir = '../code/master/fin/' # directory to output the flux files and location of master frame
clndir = '../clean/'#directory with the cleaned images
##END UPDATE####

#read in the master frame
mast, mheader = fits.getdata(caldir+camera+'_'+ccd+'_master_py.fits', header = True)

###POSSIBLE TO UNCOMMENT HERE IF YOU HAVE NO STAR LIST###

#get the positions in the master, if no star list is provided 
mean, median, std = sigma_clipped_stats(mast, sigma = 3.0, iters = 5)
daofind = DAOStarFinder(fwhm = 1.5, threshold = 5.*std)
sources = daofind(mast-median)
#set up for aperture photometry on the master frame
mean, median, std = sigma_clipped_stats(mast, sigma = 3.0, iters = 5)
positions = (sources['xcentroid'], sources['ycentroid'])
x = sources['xcentroid']
y = sources['ycentroid']

###END POSSIBLE UNCOMMENT###

###POSSIBLE TO UNCOMMENT HERE IF YOU HAVE A STAR LIST###
#get the positions from the star list if one is provided ###COMMENT HERE###
#ticid, tmag, ra, dec = numpy.loadtxt('ffi_test.csv', unpack = 1, delimiter = ',', skiprows =1)

#convert the ra/dec into positions
#os.chdir(rawdir) #changes to the raw image direcotory
#files = [f for f in glob.glob("*.fits") if isfile(join(rawdir, f))] #gets the relevant files with the proper extension
#files.sort()
#os.chdir(cdedir) #changes back to the code directory

#get wcs from first calibrated file
#w = WCS(clndir+files[0])
#x, y = w.all_world2pix(ra, dec, 0)
#positions = (x, y)
###END POSSIBLE UNCOMMENT###

print 'Getting the photometry on the master frame.'
#apertures to test for optimal apeture size
rads = [1,2,3,4,5,6]

#do the aperture photometry and find the optimal aperture
apertures = [CircularAperture(positions, r=r) for r in rads]
phot_table = aperture_photometry(mast, apertures, method = 'exact')
idx = 0

offset = numpy.zeros((len(rads),len(x)))
for ii in range(0, len(x)):
	#exclude stars near the edge of the frame
	if (x[ii] > 0) and (x[ii] < 2048) and (y[ii] > 0) and (y[ii] < 2048):
		dist = numpy.sqrt((x[ii]-x)**2+(y[ii]-y)**2) #exclude 'crowded' stars
		chk = numpy.where(dist < 6.)
		if (len(chk[0]) == 1):
			#subtract the next highest radius mag from the current mag
			for jj in range(1, len(rads)):
				if (phot_table[ii][jj+3] > 0) and (phot_table[ii][jj+2] > 0):
					mg1 = 25.-2.5*numpy.log10(phot_table[ii][jj+3])
					snr1 = phot_table[ii][jj+3]/numpy.sqrt(phot_table[ii][jj+3])
					mg0 = 25.-2.5*numpy.log10(phot_table[ii][jj+2])
					snr2 = phot_table[ii][jj+2]/numpy.sqrt(phot_table[ii][jj+2])
					if (numpy.isfinite(mg1)) and (numpy.isfinite(mg0)) and (snr1 > 100) and (snr2 > 100) and (snr1 < 1000) and (snr2 < 1000):
						offset[jj,ii] = numpy.abs(mg1-mg0)

chk_off = numpy.zeros(len(rads)-1)
chk_rad = rads[1:]
for ii in range(1, len(rads)):
	nzro_idx = numpy.nonzero(offset[ii,:]) #get the index of the zeros
	chk_off[ii-1] = numpy.median(offset[ii,nzro_idx])
fin_rads = numpy.arange(chk_rad[0],chk_rad[-1],0.1)
fin_off = numpy.interp(fin_rads,chk_rad,chk_off)

opt_rad = fin_rads[0] #starting aperture to make sure if everything fails, there is still an aperture
#select the smallest radius where the next highest radius has only a slight change in magnitude
for ii in range(1, len(fin_off)):
	diff = numpy.abs(fin_off[ii]-fin_off[ii-1])
	if (diff < 0.01):
		opt_rad = fin_rads[ii-1]
		break
print 'The optimal aperture size is '+str(opt_rad)+'.'

#do the aperture photometry
apertures = CircularAperture(positions, r = opt_rad)
phot_table = aperture_photometry(mast, apertures, method = 'exact')

#get the background of the image
cimg, clow, chigh = scipy.stats.sigmaclip(mast, low=2.5, high = 2.5) #do a 2.5 sigma clipping
bkg_mean = numpy.median(cimg) #determine the sky value
sig = numpy.std(cimg) #determine the sigma(sky)

#convert to magnitudes
flx = phot_table['aperture_sum']-(bkg_mean*(numpy.pi*opt_rad**2))
flx_er = numpy.sqrt(phot_table['aperture_sum'])
x_pix = x
y_pix = y

#create the magnitudes from the flux
#create the magnitudes from the flux
gd_flx = numpy.where(flx > 0) #check against any negative flux values
if len(gd_flx[0] > 0):
	flx = flx[gd_flx[0]]
	flx_er = flx_er[gd_flx[0]]
	x_pix = x_pix[gd_flx[0]]
	y_pix = y_pix[gd_flx[0]]
	ticid = ticid[gd_flx[0]]
mag = 25.0-2.5*numpy.log10(flx)
err = (2.5/numpy.log(10.))*(flx_er/flx)

#open and write the output files
output = open(caldir+camera+'_'+ccd+'_master_py.ap', 'w')
for ii in range(0, len(x_pix)):
	if (x_pix[ii] > 0) and (x_pix[ii] < 2048) and (y_pix[ii] > 0) and (y_pix[ii] < 2048) and (numpy.isnan(mag[ii]) != 1):
		output.write(str(long(ticid[ii]))+','+str(x_pix[ii])+','+str(y_pix[ii])+','+str(mag[ii])+','+str(err[ii])+'\n')
output.close()

#write the fluxes to a file
output = open(caldir+camera+'_'+ccd+'_master_py.flux', 'w')
for ii in range(0, len(x_pix)):
	if (x_pix[ii] > 0) and (x_pix[ii] < 2048) and (y_pix[ii] > 0) and (y_pix[ii] < 2048) and (numpy.isnan(mag[ii]) != 1):
		output.write(str(long(ticid[ii]))+','+str(x_pix[ii])+','+str(y_pix[ii])+','+str(flx[ii])+','+str(flx_er[ii])+'\n')
output.close()

#write the star list to a file
output = open(caldir+camera+'_'+ccd+'_starlist_py.txt', 'w')
for ii in range(0, len(x_pix)):
	if (x_pix[ii] > 0) and (x_pix[ii] < 2048) and (y_pix[ii] > 0) and (y_pix[ii] < 2048) and (numpy.isnan(mag[ii]) != 1):
		output.write(str(long(ticid[ii]))+','+str(long(x_pix[ii]))+','+str(long(y_pix[ii]))+'\n')
output.close()

print 'All done! See ya later alligator!'
