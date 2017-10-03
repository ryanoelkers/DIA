PRO clean

;this program will apply the bias subtraction,
;flat fielding and subtract the gradient background,
;it will also align the images to the first image in the list

;if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50

;;;;UPDATE INFORMATION HERE;;;;;
;DO YOU WANT TO FLAT FIELD AND BIAS SUBTRACT?
biassub = 1 ; yes = 1 no = 0 to bias subtract
flatdiv = 1 ; yes = 1 no = 0 to flat field

;useful directories
rawdir = '../raw/' ;directory with the raw images
cdedir = '../IDL/' ;code directory
caldir = '../calib/' ; directory with the calibration images such as bias & flat
clndir = '../clean/'; directory for the cleaned images

;sample every how many pixels? usually 32x32 is OK but it can be larger or smaller
pix = 128 ; UPDATE HERE FOR BACKGROUND SPACING
axs = 4096 ; UPDATE HERE FOR IMAGE AXIS SIZE
;;;END UPDATE INFORMATION;;;;

;setup background subtraciton
lop=2*pix
sze= (axs/pix)*(axs/pix)+2*(axs/pix)+1

;change into the image directory and get the directory listings
cd, rawdir
spawn, 'ls *.fits', files
nfiles = n_elements(files)
cd, cdedir

;read in the bias frame if it exists
if biassub eq 1 then begin
	bias = readfits(caldir+'bias.fits', bheader, /silent)
endif

;read in the flat field if it exists
if flatdiv eq 1 then begin
	flat = readfits(caldir+'flat.fits', fheader, /silent)
endif

;if no images are found in the directory then stop
if nfiles eq 0 then begin
	print, 'No images found.'
	stop
endif

if nfiles gt 0 then begin
	for ii = 0L, nfiles-1 do begin
		;split the name to make a final file
		nmehold = strsplit(files[ii], '.', /extract)
		if biassub eq 0 and flatdiv eq 0 then nmefin = nmehold[0]+'_s.fits'
    		if biassub eq 1 and flatdiv eq 1 then nmefin = nmehold[0]+'_sfb.fits'

		;read in the image
		if file_test(clndir+nmefin) eq 1 then print, 'Image '+files[ii]+' already processed.'
		if file_test(clndir+nmefin) eq 0 then begin
			print, 'Beginning to clean '+files[ii]+' at '+systime()+'.'
			
			;read in the fit file
			a = readfits(rawdir+files[ii], header,/silent)
			a = float(a)

			;if you want to bias subtract and flat divide do it here
			if biassub eq 1 then a = (a-bias)
			if flatdiv eq 1 then a = a/flatdiv

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
		
			;add the sky statistics to the header file
			sxaddpar,header,'SKYM',sky
			sxaddpar,header,'SKYS',skys
			sxaddpar,header,'SKYK',skyk
			sxaddpar,header,'SIG0',sig0
			sxaddpar,header,'SIGF',sigf
			if biassub eq 1 then sxaddpar, header, 'BIAS', 'yes'
			if flatdiv eq 1 then sxaddpar, header, 'FLAT', 'yes'
			sxaddpar, header, 'BKSUB', 'yes'

			;now do the alignment, as a default it will be to the first images
			if ii eq 0 then begin
				algn = a
			endif
			if ii gt 0 then begin
				dim = size(algn)
				f1 = fft(algn,/doub)
				f2 = fft(a,/doub)
				rat = f1*conj(f2)/abs(f1*f2)
				irat = fft(rat,/inv,/dou)
				xx = max(real_part(irat), xy0)
				xc = xy0 mod dim[1]
				yc = xy0/dim[1]
				if yc gt dim[2]/2 then yc = yc-dim[2]
				if xc gt dim[1]/2 then xc = xc-dim[1]
				a = shift(a, xc, yc)
			endif
			sxaddpar, header, 'ALIGN', 'yes'
			;write out the fit file
			writefits, clndir+nmefin, a, header
			
			print, 'Image '+files[ii]+' cleaned at '+systime()+'.'
		endif
	endfor
endif

print, 'All done at ', systime(),'. See ya later alligator.'

END
