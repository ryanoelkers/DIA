PRO clean

;this program will apply the bias subtraction, flat fielding, subtract the gradient background and align to the first image in the list

;;;;UPDATE INFORMATION HERE;;;;;
;DO YOU WANT TO FLAT FIELD AND BIAS SUBTRACT?
biassub = 0 ; yes = 1 no = 0 to bias subtract
flatdiv = 0 ; yes = 1 no = 0 to flat field
align = 1; yes = 1 no = 0 to align based on coordinates

;useful directories
rawdir = '.../cal/' ;directory with the raw images
cdedir = '.../code/clean/' ;directory where the code 'lives'
caldir = 'N/A' ; directory with the calibration images such as bias & flat
clndir = '../clean/'; directory for the cleaned images to be output

;sample every how many pixels? usually 32x32 is OK but it can be larger or smaller
pix = 32 ; UPDATE HERE FOR BACKGROUND SPACING
axs = 2048 ; UPDATE HERE FOR IMAGE AXIS SIZE
;;;END UPDATE INFORMATION;;;;

;setup background subtraciton
lop=2*pix
sze= (axs/pix)*(axs/pix)+2*(axs/pix)+1

;read in the bias frame if it exists
if biassub eq 1 then begin
	bias = readfits(caldir+'bias.fits', bheader, /silent)
endif

;read in the flat field if it exists
if flatdiv eq 1 then begin
	flat = readfits(caldir+'flat.fits', fheader, /silent)
endif

;output a file with the image statistics, for those interested
openw, 1, 'images.dat'

;get the file list
cd, rawdir
spawn, 'ls *.fits', files
nfiles = n_elements(files)
cd, cdedir

;read in the first image to get the 'master' header
img0 = mrdfits(rawdir+files[0], 1, head0, /silent)
hextract, img0, head0, mast, mhead, 44, 2091, 0, 2047, /silent

for ii = 0L, nfiles-1 do begin
	;split the name to make a final file
	nmehold = strsplit(files[ii], '.', /extract)
	if biassub eq 0 and flatdiv eq 0 and align eq 0 then nmefin = nmehold[0]+'_s.fits'
	if biassub eq 1 and flatdiv eq 1 and align eq 0 then nmefin = nmehold[0]+'_sfb.fits'
	if biassub eq 0 and flatdiv eq 0 and align eq 1 then nmefin = nmehold[0]+'_sa.fits'
	if biassub eq 1 and flatdiv eq 1 and align eq 1 then nmefin = nmehold[0]+'_sfba.fits'

	;read in the image
	if file_test(clndir+nmefin) eq 1 then print, 'Image '+files[ii]+' already processed.'
	if file_test(clndir+nmefin) eq 0 then begin
		print, 'Beginning to clean '+files[ii]+' at '+systime()+'.'
		
		;read in the fits file
		img = mrdfits(rawdir+files[ii], 1,head, /silent)
		hextract, img, head, a, header, 44, 2091, 0, 2047, /silent ; this may need to be updated for cutouts based on camera

		;if you want to bias subtract and flat divide do it here
		if biassub eq 1 then a = (a-bias)
		if flatdiv eq 1 then a = a/flat

		;Calculate sky statistics
		kp = where(a gt -100. and a lt 25000.)
		q=a(kp)
		mmm,q,sky,sig,/silent
		sig0=sig
		print,'The sky counts are '+strcompress(sky, /remove_all)+' ADU per pixel.'

		;Sample "local" sky value
		x=fltarr(sze)
		y=fltarr(sze)
		v=fltarr(sze)
		s=fltarr(sze)
		nd=0l
	
		for i=0,axs,pix do begin
			for j=0,axs,pix do begin
		
				il=max([i-lop,0])
				ih=min([i+lop,axs-1])
				jl=max([j-lop,0])
				jh=min([j+lop,axs-1])
				c=a(il:ih,jl:jh)
				kp=where(c gt -100. and c lt 25000.)
				c=c(kp)


				mmm,c,mean,sig,/silent
				x(nd)=min([i,axs-1])
				y(nd)=min([j,axs-2])
				v(nd)=mean
				s(nd)=sig
				nd++

			endfor
		endfor

		rj=where(v lt 0 or s gt 300,nr,complement=kp,ncomplement=nk)
		    
		if (nr ge 1) then begin
			xgood=x(kp)
			ygood=y(kp)
			vgood=v(kp)
			sgood=s(kp)
			for j=0l,nr-1 do begin
				xbad=x(rj(j))
				ybad=y(rj(j))
				rd=sqrt((xgood-xbad)^2+(ygood-ybad)^2)
				sr=sort(rd)
				vnear=vgood(sr(0:9))
				ave=median(vnear)
				ave=ave(0)
				v(rj(j))=ave
			endfor
		endif

		;interpolate for rest of image and subtract
		res=grid_tps(x,y,v,ngrid=[axs,axs],start=[0,0],delta=[1,1])

		;remove the gradient and return the sky
		a = double(a) - res
		res = uint(res)
		medv = median(v)
		a = a+medv
		a = float(a)
		    
		;compute new sky statistics
		kp=where(a gt -1000.)
		q=a(kp)
		mmm,q,sky,sig,/silent
		sigf=sig
		skys=stddev(res)
		skyk=kurtosis(res)
	
		;now we want to align the images based on wcs

		;align to the first image	
		hastrom, a, header, newa, newheader, mhead, interp = 2, cubic = -0.5
		
		;add the sky statistics to the header file
		sxaddpar,newheader,'SKYM',sky
		sxaddpar,newheader,'SKYS',skys
		sxaddpar,newheader,'SKYK',skyk
		sxaddpar,newheader,'SIG0',sig0
		sxaddpar,newheader,'SIGF',sigf

		;add what we did to the header file
		if biassub eq 1 then sxaddpar, newheader, 'BIAS', 'yes'
		if flatdiv eq 1 then sxaddpar, newheader, 'FLAT', 'yes'
		if align eq 1 then sxaddpar, newheader, 'ALIGN', 'yes'
		sxaddpar, newheader, 'BKSUB', 'yes'
	
		;get the julian dates of all the images	
		obs1 = sxpar(newheader, 'TSTART', /silent)
		obs2 = sxpar(newheader, 'TSTOP', /silent)
		bjd = mean([obs1, obs2])

		;get the rough number of stars on the frame
		find, newa, x, y, flux, sharp, round, 3*skys, 1.88, [-1.,1.],[.2,1.], /silent
		
		;write the file to the drive
		writefits, clndir+nmefin,newa, newheader
		
		;print the information to a data file
		printf, 1, nmefin, bjd[0], sky, skys, n_elements(x), format = '(a, 1x, f11.6, 1x, f8.2, 1x, f7.2, 1x, i7)'
		print, 'Image '+files[ii]+' cleaned at '+systime()+'.'
	endif
endfor


print, 'All done at ', systime(),'. See ya later alligator.'
close, 1
END
