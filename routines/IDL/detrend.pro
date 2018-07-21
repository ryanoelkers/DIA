PRO detrend

;this program will detrend stars based on similar magnitude

;if you use this code, please cite Oelkers et al. 2015, AJ, 149, 50 -- Oelkers et al. 2018 submitted

;;;;UPDATE INFORMATION HERE;;;;;

;useful directories
lcdir = '../lc/' ;the directory where the raw light curves 'live'
otdir = '../lc/detrend/' ;the directory where the detrended light curves will be output to
cdedir = '../code/phot/' ;directory where this code 'lives'

;read in the star information
readcol, lcdir+'rms_hr.dat', nme, clp, mg, rms, format = '(a, f, f, f)', /silent, count = nstars

;read in all light curves to decrease the detrending time
;make the light curve holder
big = fltarr(nstars, 1348)
s = sort(mg)

for ii =0l, nstars-1 do begin
	;read in the light curves based on magnitude
	readcol, lcdir+nme[s[ii]], j, m, e, format = '(f,f,f)', /silent, count=nobs
	big[ii,*] = m
	if (ii mod 1000 eq 0 and ii gt 0) then print, '1000 light curves read in at '+systime()
endfor

;now loop through the stars and build the trend
for ii = 0l, nstars-1 do begin
	;read in the star to detrend
	readcol, lcdir+nme[s[ii]], j, m, e, format = '(f,f,f)', /silent

	;find the +/-500 stars to use in magnitude space
	udx = 1000-ii
	if (udx lt 500) then udx = ii+500
	if (udx gt nstars-1) then udx = nstars-1
	ldx = ii-500
	if (ldx lt 0) then ldx = 0

	;initial check of rms values
	chk1 = stddev(m)
	;holder for the trend stars
	trd = fltarr(1000,1348)

	;loop through and check for >10% improvement in rms
	for jj = ldx, udx-1 do begin

		if jj ne ii then begin
			;second check of rms values
			chk2 = stddev(m-(big[jj,*]-median(big[jj,*])))
			if chk1/chk2 gt 1.1 then begin ;if it is better than 10%, include the star
				trd[jj-ldx,*] = big[jj,*]-median(big[jj,*])
			endif
		endif
	endfor

	trend = fltarr(1348)
	brk = [0]
	offs = fltarr(1348)
	;now loop through and median combine the trends, if they exist
	for jj = 0, 1347 do begin
		vv = where(trd[*,jj] ne 0)
		if vv[0] ne -1 then begin
			trend[jj] = median(trd[vv,jj])
			if jj gt 0 then offs[jj] = abs(trend[jj]-trend[jj-1])
		endif
	endfor
	
	;check to see if a trend was found to remove
	chk3 = where(trend ne 0)
	otrend = trend ;keep the old trend just in case

	;if a trend was found, begin attempting to scale
	if chk3[0] ne -1 then begin
		;check for large deviations between subsequent dat apoints
		meanclip, offs, mm, ss, clipsig = 3, maxiter = 1000
		brk = mm+5*ss
		tme = where(offs gt brk)
		tme = [0,tme]
		;if large devaitions exist, move through them and scale
		for jj = 0, n_elements(tme)-1 do begin
	
			lw = tme[jj]
			if jj lt n_elements(tme)-1 then upp = tme[jj+1]-1
			if jj eq n_elements(tme)-1 then upp = 1347

			mm = median(m[lw:upp])-median(m)
			mt = median(trend[lw:upp])

			off = mm-mt
			chk5 = stddev(m[lw:upp])
			chk4 = stddev(m[lw:upp] - (trend[lw:upp]+off))

			trend[lw:upp] = trend[lw:upp]+off

		endfor
		df1 = fltarr(1348)
		df2 = fltarr(1348)
		;check to make sure scaling didn't increase the nnoise
		for jj =1l, 1347 do begin
			df1[jj] = m[jj]-m[jj-1]
			df2[jj] = (m[jj]-trend[jj])-(m[jj-1]-trend[jj-1])
		endfor
		meanclip, df1, mm1, ss1, clipsig = 2.5, maxiter = 1000
		meanclip, df2, mm2, ss2, clipsig = 2.5, maxiter = 1000

		num1 = n_elements(where(abs(df1) gt mm1+5*ss1))
		num2 = n_elements(where(abs(df2) gt mm2+5*ss2))

		;if more breaks are found in the light curve, keep the old trend
		if num2 ge num1 then trend = otrend

		;output the cleaned light curve
		openw, 1, otdir+nme[s[ii]]
		for jj =0l, n_elements(j)-1 do printf, 1, j[jj],m[jj]-trend[jj],e[jj], format ='(f11.6, 1x, f9.6, 1x, f8.6)'
		close, 1
	endif
	
	;if no trend stars exist, the output the raw light curve
	if chk3[0] eq -1 then begin
		openw, 1, otdir+nme[s[ii]]
		for jj =0l, n_elements(j)-1 do printf, 1, j[jj],m[jj],e[jj], format ='(f11.6, 1x, f9.6, 1x, f8.6)'
		close, 1
	endif
	if (ii mod 1000 eq 0 and ii gt 0) then print, 'Working on the next 1000 curves at '+systime()

endfor
print, 'All done at '+strcompress(systime(),/remove_all)+' See ya later alligator!'
END
