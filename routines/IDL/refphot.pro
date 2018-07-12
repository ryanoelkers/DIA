PRO refphot

;this program will do the photometry on the reference frame to get stars for the subtraction
;and to get stars for later. it will also tell you the best aperture to use for the reduction
;you will either need to give it a list of RA/Dec information for stars or you can have the program
;detect sources on its own.

;if you use this code, please cite Oelkers & Stassun 2018

;;;;UPDATE INFORMATION HERE;;;;;

;useful directories
cdedir = '../code/master/' ;directory where the code 'lives'
caldir = '../code/master/fin/' ; directory to output the flux files, typically in the location of master frame
clndir = '../clean/';directory where the cleaned and aligned data lives

;get image list in the clean directory
cd, clndir
spawn, 'ls *.fits', files
cd, cdedir

;which field are you looking at?
camera = '2';camera name
ccd = '2';ccd name
gain = 5.2; gain from header
x = [] & y = []

;if a list does exist, read in the list and get the photometry of the stars;;;UPDATE HERE IF YOU HAVE RA/DEC COORDS
readcol, 'ffi_test.csv', ticid, tmag, ra, dec, format = '(a, f, f, f)', /silent
img = readfits(clndir+files[0], rhead, /silent);update here to get position information

;convert the ra and dec to a pixel position, if you are using ra/dec rather than x/y
adxy, rhead, ra, dec, x, y
;;;END UPDATE INFORMATION HERE;;;;

;update the aperture size based on the images
apr = dindgen(32)/4.+2.;;set apertures from 2 to 10 in size
tor = [2,4];;dummy holder for the torus

;read in the reference file for photometry
ref = readfits(caldir+camera'_'+ccd+'_master_idl.fits', mhead, /silent)

print, 'Getting flux information from master frame at '+systime()+'.'

;find the sky information
sky, ref, skym, skys, /silent

;if no list of stars based on ra/dec was supplied the program finds stars on its own
if (n_elements(x) eq 0) then find, ref, x,y, flux, sharp, round, 3*skys, 1.88, [-1.,1.], [.2,1.], /silent

;remove stars outside of the frame
bd = where(x lt 0 or x gt 2048 or y lt 0 or y gt 2048)
remove, bd, ticid, tmag, ra, dec, x, y
s = sort(tmag)
ticid = ticid[s] & tmag = tmag[s] & ra = ra[s] & dec = dec[s] & x = x[s] & y = y[s]

;do the aperture photometry on those stars in both flux and magntiude space
aper, ref, x, y, mags, errap, sky, skyer, gain, apr, tor,[-100, 55000], /silent, setskyval = skym

;now check for the optimal aperture 
ok = where(mags[0,*] lt 30 and mags[0,*] gt 0 and errap[0,*] gt 0.001 and errap[0,*] le 0.01)

offset = fltarr(n_elements(apr), n_elements(ok))
for ii = 0l, n_elements(ok)-1 do begin
	dist = sqrt((x[ok[ii]]-x)^2+(y[ok[ii]]-y)^2);;check for isolated stars to do the aperture check
	bd = where(dist lt 6 and dist gt 0)
	if bd[0] eq -1 then begin
		for jj = 1l, n_elements(apr)-1 do begin
			offset[jj,ii] = abs(mags[jj,ii]-mags[jj-1,ii])
		endfor
	endif
endfor

;now find when the change in mags is < 0.001
prv = 1.
opt_apr = apr[0] ; set the initial aperture to be the smallest, in case something fails below
for ii = 1l, n_elements(apr)-1 do begin
	meanclip, offset[ii,*], chk, schk, clipsig = 2.5, maxiter = 1000
	if abs(chk-prv) lt 0.001 then begin
		opt_apr = apr[ii]
		break
	endif
	if abs(chk-prv) ge 0.001 then prv = chk;if the change is larger than 1mmag continue
endfor
print, 'The optimal aperture is '+strcompress(apr[ii], /remove_all)+' pixels.'

;do the aperture photometry on those stars in both flux and magnitude space
aper, ref, x, y, mags, errap, sky, skyer, gain, opt_apr, tor,[-100, 55000], /silent, setskyval = skym
aper, ref, x, y, flux, fluxer, fsky, fskyer, gain, opt_apr, tor, [-100, 55000], /silent, /flux, setskyval = skym

;get rid of bad photomerty
good = where(mags lt 30. and mags gt 0. and errap gt 0 and errap lt 0.1)

;get rid of possible bad stars
ticid = ticid[good] & tmag = tmag[good] & x = x[good] & y = y[good] & flux = flux[good] & fluxer = fluxer[good] & sky = sky[good]
skyer = skyer[good]& fsky = fsky[good] & fskyer = fskyer[good] & mags = mags[good] & errap = errap[good]

;write magnitude information
openw, 1, caldir+camera+'_'+ccd+'_master_idl.ap'
for i = 0L, n_elements(x) -1 do printf, 1, ticid[i], x[i], y[i], tmag[i], mags[i], errap[i], sky[i], skyer[i], format = '(i0, 1x, d0,1x, d0, 1x,d0,1x,d0,1x,d0,1x,d0, 1x, d0)'
close, 1
;write the flux information
openw, 1, caldir+camera+'_'+ccd+'_master_idl.flux'
for i = 0L, n_elements(x) -1 do printf, 1, ticid[i], x[i], y[i], flux[i], fluxer[i], fsky[i],fskyer[i], format = '(i0, 1x, d0,1x, d0,1x,d0,1x,d0,1x,d0, 1x, d0)'
close, 1
;write the star list of differencing 
openw, 1, caldir+camera+'_'+ccd+'_starlist_idl.txt'
for i = 0L, n_elements(x) -1 do printf, 1, ticid[i], x[i], y[i], format = '(i0,1x,d0,1x,d0)'
close, 1

print, 'Finished at '+systime()+'. See ya later alligator!'

END
