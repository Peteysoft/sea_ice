@ice_model4
@plane_parallel
@eps_parm

@cver/ice_model_c

;goto, skip

pro ridged_ice_model, rfile, tfile, lfile, outfile, $
		ntest=ntest, niter=niter, doplot=doplot, $
		fov=fov, tice=temp, tsea=tsea, sea_sal=sea_sal, $
		rindex=rindex, eps_ice0=eps, eps_parm0=eps_parm, $
		eps_parm1=eps_parm1, eps_ice1=eps_t, $
		eps_parm2=eps_parm2, eps_ice2=eps_b, $
		dvary=dvary, smooth_iw=smooth_iw, $
		$
		tbv=tbv, tbh=tbh

  t1=systime(1)

  topofile="icetopo.txt"
  initfile="iceinit.txt"

  openw, lun, "ridged_model_error_log.txt", /get_lun

  ;step 1: read in the relevant data:
  t0='2007/3/13'
  read_police, rfile, t0, tbs, lonr, latr, tr, /cor, $
		zangle=zangle, angle0=-40., alt=zr
  read_thickness, tfile, dt, lont, latt, tt, flag
  read_laser, lfile, zl, lonl, latl, tl

  ;use our clever flight-track rotation program to align them to the x-axis:
  rotate_flight_track, lonr, latr, xr, yr, /new
  rotate_flight_track, lont, latt, xt, yt
  rotate_flight_track, lonl, latl, xl, yl

  ;right off the bat, we calculate the ice-water interface:
  if keyword_set(smooth_iw) then begin
    zave=smooth(zl, 11, /edge_truncate)
    ht=interpol(zave, xl, xt)-dt
  endif else begin
    ht=interpol(zl, xl, xt)-dt
  endelse

  ;stop

  ;***the fact that the data is laid out on three different grids 
  ;(depending on what it is) makes
  ;this program very confusing and prone to errors...

  ;initialize the simulation:

  if keyword_set(dvary) then begin
    eps_parm_t="eps_parm2a"
    eps_parm_b="eps_parm2b"
  endif

  ;complex permittivity is related to ice thickness:
  if n_elements(eps_parm) ne 1 then eps_parm="eps_parm2"
  if n_elements(eps) eq 0 then begin
    ;really need a more elegant way to go about this...
    call_procedure, eps_parm, dt, eps_re, eps_im
    eps=complex(eps_re, eps_im)
  endif else begin
    if n_elements(eps) eq 1 then eps=replicate(eps, n_elements(dt))
    eps_im=imaginary(eps)
    eps_re=float(eps)
  endelse

  if n_elements(eps_t) eq 0 then begin
    if n_elements(eps_parm1) eq 1 then begin
      call_procedure, eps_parm1, dt, eps_re_t, eps_im_t
      eps_t=complex(eps_re_t, eps_im_t)
    endif
  endif else begin
    eps_t=replicate(eps_t, n_elements(dt))
    eps_im_t=imaginary(eps_t)
    eps_re_t=float(eps_t)
  endelse

  if n_elements(eps_b) eq 0 then begin
    if n_elements(eps_parm2) eq 1 then begin
      call_procedure, eps_parm2, dt, eps_re_b, eps_im_b
      eps_b=complex(eps_re_b, eps_im_b)
    endif
  endif else begin
    eps_b=replicate(eps_b, n_elements(dt))
    eps_im_b=imaginary(eps_b)
    eps_re_b=float(eps_b)
  endelse

  freq=1.4

  if n_elements(fov) ne 1 then fov=13.16	;angular field of view

  ;calculate an effective field of view:
  ;fov_eff=atan(sqrt((data_r.z*tan(fov*!dtor))^2-data_r.dy^2)/data_r.z)

  fov1=fov*!dtor

  ;physical temperature of ice:
  if n_elements(temp) ne 1 then temp=273.

  ;complex refractive index of (brackish) sea:
  if n_elements(tsea) ne 1 then tsea=272.
  if n_elements(sea_sal) ne 1 then sea_sal=5.
  eps_sea=eps_oc2(freq, tsea, sea_sal)
  nsea=sqrt(eps_sea)

  ;pick a random starting point:
  bmargin=3000.		;depends on fov
  fmargin=2000.
  min_xt=min(xt, j1)
  max_xt=max(xt, j2)
  if j1 gt j2 then begin
    xt=reverse(xt)
    yt=reverse(yt)
    dt=reverse(dt)
    reverse_flag=1
    nt=n_elements(xt)
  endif
  min_xl=min(xl, j1)
  max_xl=max(xl, j2)
  if j1 gt j2 then begin
    xl=reverse(xl)
    yl=reverse(yl)
    zl=reverse(zl)
    ;stop
  endif
  ind=where(xr gt min_xt > min_xl+bmargin and xr lt max_xt < max_xl-fmargin, nr2)

  ;if we have the subscripts we want to use, use them:
  if n_elements(ntest) ne 1 then begin
    if n_elements(rindex) eq 0 then begin
      ntest=2000
      rindex=ind[long(randomu(seed, ntest)*nr2)]
    endif else begin
      ntest=n_elements(rindex)
      no_retry_on_err=1
      ;set error value for out-of-range subscripts:
      check=where(rindex gt max(ind) or rindex lt min(ind), cnt)
      if cnt gt 0 then rindex[check]=-1
    endelse
  endif else begin
    rindex=ind[long(randomu(seed, ntest)*nr2)]
  endelse
  if n_elements(niter) ne 1 then niter=1000

  tbv_meas=fltarr(ntest)
  tbh_meas=fltarr(ntest)

  tbv=fltarr(ntest)
  tbh=fltarr(ntest)

  tbvc=fltarr(ntest)
  tbhc=fltarr(ntest)

  verrc=fltarr(ntest)
  herrc=fltarr(ntest)

  tbvc2=fltarr(ntest)
  tbhc2=fltarr(ntest)

  verrc2=fltarr(ntest)
  herrc2=fltarr(ntest)

  tbv_pp=fltarr(ntest)
  tbv_pp=fltarr(ntest)
  tbh_pp=fltarr(ntest)

  tbv_ppe=fltarr(ntest)
  tbh_ppe=fltarr(ntest)

  verr=fltarr(ntest)
  herr=fltarr(ntest)

  thick=fltarr(ntest)
  thstd=fltarr(ntest)
  hstd=fltarr(ntest)

  zangarr=fltarr(ntest)

  tindex1=lonarr(ntest)
  tindex2=lonarr(ntest)

  ;auxiliary data:
  s=fltarr(ntest)
  lon=fltarr(ntest)
  lat=fltarr(ntest)
  wfactor=fltarr(ntest)
  swflag=fltarr(ntest)		;shallow water flag

  ;take a Gaussian-weighted mean of the measurement points:
  var=2.^2

  ;stop

  for i=0L, ntest-1 do begin
    if rindex[i] lt 0 then continue
    ;fov1=fov_eff[rindex[i]]
    ;fov1=fov*!dtor

    ;zangarr[i]=zangle[rindex[i]]
    s0=xr[rindex[i]]

    weight=exp(-(s0-xr)^2/var/2)		;bit of a waste summing all of the points...
    wt=total(weight)
    zangarr[i]=total(zangle*weight)/wt
    zangle1=zangarr[i]*!dtor

    z0=zr[rindex[i]]
    dx1=z0*tan(zangle1)
    x0=s0+dx1

    ;tbv_meas[i]=tbs[1, rindex[i]]
    ;tbh_meas[i]=tbs[0, rindex[i]]

    tbv_meas[i]=total(tbs[1, *]*weight)/wt
    tbh_meas[i]=total(tbs[0, *]*weight)/wt

    lon[i]=lonr[rindex[i]]
    lat[i]=latr[rindex[i]]
    s[i]=s0

    ;since IDL is too stupid to fit the entire surface, we have to use
    ;only a piece:
    indn=where(xl lt s0+fmargin and xl gt s0-bmargin, n)
    tindex1[i]=indn[0]
    if keyword_set(reverse_flag) then begin
      inter=tindex1[i]
      tindex1[i]=nt-tindex2[i]-1
      tindex2[i]=nt-inter-1
    endif

    indn2=where(xt lt s0+fmargin and xt gt s0-bmargin, n2)
    dx2=xt[indn2]-x0
    dy2=yt[indn2]-yr[rindex[i]]
    z2=z0^2
    ;offset angle from instrument dead-centre:
    r1=sqrt(dx2^2+dy2^2+z2)
    r2=sqrt(dx1^2+z2)
    theta=acos((dx1*dx2+z0^2)/r2/r1)
    ;calculated weights:
    w2=exp(-theta^2/fov1^2/2)

    ;plot, theta, w2
    ;stop

    tw2=total(w2)
    wfactor[i]=tw2
    swflag[i]=total(flag[indn2]*w2)/tw2

    ;average over instrument weighting function:
    thick[i]=total(w2*dt[indn2])/tw2
    thstd[i]=sqrt(total(w2*(dt[indn2]-thick[i])^2)/tw2)

    s2=xl[indn]
    ;because of the poor resolution of the data, s needs to be smoothed:
    m=regress(findgen(n), s2, const=b1)		;this is a really, really ugly hack
    s1=m[0]*findgen(n)+b1

    ;this is very, very confusing
    ;--here I switch around the variable names!!
    h=zl[indn]

    ;interpolate the thickness to the laser altimetry data:
    ;d=interpol(dt, xt, s1, /spline)
    ;z=h-d
    z=interpol(ht, xt, s1, /spline)

    ;sometimes the fitting fails (use linear interpolation):
    indf=where(finite(z) ne 1, cntf)
    if cntf gt 0 then z[indf]=interpol(ht, xt, s1[indf])

    ;stop

    ;need to interpolate the complex permittivity data as well:
    eps_im1=interpol(eps_im, xt, s1)
    eps_re1=interpol(eps_re, xt, s1)
    eps1=complex(eps_re1, eps_im1)

    if n_elements(eps_b) ne 0 then begin
      eps_im2=interpol(eps_im_b, xt, s1)
      eps_re2=interpol(eps_re_b, xt, s1)
      eps2=complex(eps_re2, eps_im2)
    endif else begin
      eps2=eps1
    endelse
    nice2=sqrt(eps2)

    if n_elements(eps_t) ne 0 then begin
      eps_im0=interpol(eps_im_t, xt, s1)
      eps_re0=interpol(eps_re_t, xt, s1)
      eps0=complex(eps_re0, eps_im0)
    endif else begin
      eps0=eps1
    endelse
    nice0=sqrt(eps0)

  ;goto, skip20
    ;dy0=min(abs(dy2), j)
    ;offset=dy0/r2
    ;different method of calculating this offset:
    offset=abs(sin(min(theta)))

    ;goto, skip1

    t1=systime(1)

    ;supplementally, we call the c routines:
    openw, lun1, topofile, /get_lun
    printf, lun1, n_elements(s1)
    printf, lun1, transpose([[s1], [h], [z], [eps_re1], [eps_im1]]), format="(5f)"
    free_lun, lun1

    openw, lun1, initfile, /get_lun
    printf, lun1, 1
    printf, lun1, x0, z0, -zangle1, offset
    free_lun, lun1

    spawn, "cver/ice_main "+topofile+" "+initfile, results
    print, results
    tbvc1=0.
    tbhc1=0.
    verrc1=0.
    herrc1=0.
    field=strsplit(results[0], "+/-", /extract)
    reads, field[0], tbvc1
    reads, field[1], verrc1 
    field=strsplit(results[1], "+/-", /extract)
    reads, field[0], tbhc1
    reads, field[1], herrc1

    tbvc[i]=tbvc1
    tbhc[i]=tbhc1
    verrc[i]=verrc1
    herrc[i]=herrc1

    skip1:

    t2=systime(1)
    ;print, t2-t1, "s"
    ;stop

    ;goto, skip2

    ice_model_c, s1, h, z, x0, z0, -zangle1, fov1, $
	    	tbv2, tbh2, verr2, herr2, fail=fail, $
		eps=eps0, freq=freq, tice=temp, offset=offset, $
		plot=doplot, niter=niter, nsea=nsea, tw=tsea, $
		sw=sea_sal

    tbvc2[i]=tbv2
    tbhc2[i]=tbh2
    verrc2[i]=verr2
    herrc2[i]=herr2

    skip2:

    if tbvc1 gt 300 or tbvc1 lt 150 or tbhc1 gt 300 or tbhc1 lt 100 then stop

    t2a=systime(1)
    print, t2a-t2, "s"

    print, x0, z0
    tbv[i]=ice_model(s1, h, z, x0, z0, -zangle1, fov1, pol="v", fail=fail, $
		n1=nice0, n2=nice2, freq=freq, tice=temp, offset=offset, $
		err=err, plot=doplot, niter=niter, nsea=nsea, tw=tsea)
    if fail ne 0 then begin
      printf, lun, !error_state	   
      printf, lun, rfile
      printf, lun, tfile
      printf, lun, lfile
      printf, lun, "EMIRAD measurement index=", rindex[i]
      if keyword_set(no_retry_on_error) then begin
        rindex[i]=-1
      endif else begin
        rindex[i]=ind[long(randomu(seed, 1)*nr2)]
        i=i-1
      endelse
      continue
    endif
    verr[i]=err

    tbh[i]=ice_model(s1, h, z, x0, z0, -zangle1, fov1, pol="h", fail=fail, $
		n1=nice0, n2=nice2, freq=freq, tice=temp, offset=offset, $
		err=err, plot=doplot, niter=niter, nsea=nsea, tw=tsea)
    if fail ne 0 then begin
      printf, lun, !error_state
      printf, lun, rfile
      printf, lun, tfile
      printf, lun, lfile
      printf, lun, "EMIRAD measurement index=", rindex[i]
      if keyword_set(no_retry_on_error) then begin
        rindex[i]=-1
      endif else begin
        rindex[i]=ind[long(randomu(seed, 1)*nr2)]
        i=i-1
      endelse
      continue
    endif
    herr[i]=err

    t3=systime(1)
    print, t3-t2a, "s"
  ;skip20:

    ;most basic plane-parallel model, taking just zenith angle
    ;and average ice thickness:
    eps1a=total(w2*eps[indn2])/tw2
    three_layer, freq, zangle1, thick[i], temp, eps1a, tbv1, tbh1, eps_oc=eps_sea, tw=tsea
    tbv_pp[i]=tbv1[0]
    tbh_pp[i]=tbh1[0]

    ;plane-parallel model taking an ensemble average
    ;for different thicknesses:
    zangle2=atan(sqrt(dx2^2+dy2^2), z0)
    three_layer, freq, zangle2, dt[indn2], temp, eps[indn2], tbv1, tbh1, eps_oc=eps_sea, tw=tsea
    tbv_ppe[i]=total(w2*tbv1)/tw2
    tbh_ppe[i]=total(w2*tbh1)/tw2

    ;stop
  
  endfor

  ;ps_end
  t2=systime(1)

  print, "Elapsed time:", t2-t1

  save, filename=outfile, tbv, tbh, verr, herr, tbv_meas, tbh_meas, $
		thick, thstd, zangarr, tbv_pp, tbh_pp, $
		tbv_ppe, tbh_ppe, s, lon, lat, wfactor, $
		eps, tsea, temp, freq, fov, niter, $
		rfile, tfile, lfile, rindex, tindex1, tindex2, $
		swflag

  free_lun, lun

  stop

end


