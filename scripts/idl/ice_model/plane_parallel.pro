@ice_model_util
@eps_calc

;dirt simple plane parallel radiative transfer model:

;most basic version, takes only the complex permittivities:
;2008-11-06 PM: cleaned things up by converting to matrix operations...
pro ppeps, nu, zangle, dz, t, eps, bound, tbv, tbh, theta=theta, $
		weight_h=weight_h, avet_h=avet_h, weight_v=weight_v, avet_v=avet_v, $
		raw_weight=raw_weight
;+*****************************************************************************
;		PPEPS
;******************************************************************************
;
; usage:	ppeps, nu, zangle, dz, t, eps, bound, tbv, tbh
;
; purpose:	Plane-parallel radiative-transfer model with discontinuous interfaces.
;		Given a profile comprised of layer-thicknesses, temperatures
;		and complex permittivities, solves for upwelling- and downwelling-
;		brightness temperatures.
;
; input:	nu	Frequency in GHz
;
;		zangle	Zenith angle in radians.
;
;		dz	Layer thickneses in m.  Should be from top to bottom.
;
;		t	Temperatures in K.
;
;		eps	Complex permittivities.
;
;		bound	Indices of boundaries between layers for which the
;			reflection coefficients should be calculated
;			(discontinuous interfaces).  0 is boundary between 
;			layer 0 and layer 1.
;
; output:	tbv	Vertically polarized brightness temperatures for all layers.  
;			First half are upwelling, second half, downwelling.
;
; 		tbh	Horizontally polarized brightness temperatures for all layers.  
;			First half are upwelling, second half, downwelling.
;
; keywords:	theta	Returns array of beam angles.
;
;		weight_h	Returns the temp. weighting function.
;		weight_v	Returns the temp. weighting function.
;
;		avet_h	Returns the temperature averaged over the weighting
;		avet_v	function.
;
; dependencies:	ice_model_util.pro
;
; author: 	Peter Mills (pmills@iup.physik.uni-bremen.de)
;
; history:	Formally documented 2009-6-02
;
;-*****************************************************************************

  nz=n_elements(eps)

  ;calculate refractive indices and absorption coefficients:
  n=sqrt(eps)
  ;it makes a pretty big difference when we convert to real values and whether
  ;we take the absolute value or simply discard the imaginary component...
  n1=sqrt(float(eps))
  c=0.3
  abs_coef=4*!pi*imaginary(n)*nu/c
  
  ;calculate refraction angles and Fresnel coefficients:
  nbound=n_elements(bound)
  ;theta1=fltarr(nz)

  theta=asin(n1[0]*sin(zangle)/n1)

  ;try it the recurive way and compare:
  ;theta1=theta
  ;theta1[0]=zangle
  ;for i=1, nz-1 do theta1[i]=asin(n[i-1]*sin(theta[i-1])/n[i])
  ;print, theta
  ;print, theta1

  ;stop

  costhet=cos(theta)
  ;print, theta, costhet

  ;at what point do we convert to real values??
  ;do we do so by taking the absolute value or by discarding the imaginary component??

  ;rh=(n1[bound]*costhet[bound]-n1[bound+1]*costhet[bound+1])/ $
	;	(n1[bound]*costhet[bound]+n1[bound+1]*costhet[bound+1])
  ;rv=(n1[bound]*costhet[bound+1]-n1[bound+1]*costhet[bound])/ $
	;	(n1[bound]*costhet[bound+1]+n1[bound+1]*costhet[bound])

  rh=(n[bound]*costhet[bound]-n[bound+1]*costhet[bound+1])/ $
		(n[bound]*costhet[bound]+n[bound+1]*costhet[bound+1])
  rv=(n[bound]*costhet[bound+1]-n[bound+1]*costhet[bound])/ $
		(n[bound]*costhet[bound+1]+n[bound+1]*costhet[bound])

  ;real angles and reflection coefficients for each level
  ;to make computation easier:
  rh1=fltarr(nz)
  rv1=fltarr(nz)

  rh1[bound]=abs(rh*rh)
  rv1[bound]=abs(rv*rv)
  ;print, rv1, rh1

  ;transmission coefficient at each layer
  trans=exp(-abs_coef*dz/cos(abs(theta)))
  ;emission temperature of each layer:
  tem=t*(1-trans)

  ;construct the matrix:

  a1=fltarr(nz*2, nz*2)		;vertical polarisation
  a2=fltarr(nz*2, nz*2)		;horizontal polarisation

  ;more efficient way to populate the matrix:
  ind2=lindgen(nz*2)
  a1[ind2, ind2]=1		;diagonal is 1
  a2[ind2, ind2]=1		;diagonal is 1

  ind1=lindgen(nz-1)

  ;coefficients for upwelling brightness temperatures:
  a1[ind1+1, ind1]=-(1-rv1[0:nz-2])*trans[0:nz-2]
  a1[ind1+nz, ind1]=-rv1[0:nz-2]*trans[0:nz-2]
  a2[ind1+1, ind1]=-(1-rh1[0:nz-2])*trans[0:nz-2]
  a2[ind1+nz, ind1]=-rh1[0:nz-2]*trans[0:nz-2]

  ;coefficients for downwelling brightness temperatures:
  a1[ind1+nz, ind1+nz+1]=-(1-rv1[0:nz-2])*trans[1:nz-1]
  a1[ind1+1, ind1+nz+1]=-rv1[0:nz-2]*trans[1:nz-1]
  a2[ind1+nz, ind1+nz+1]=-(1-rh1[0:nz-2])*trans[1:nz-1]
  a2[ind1+1, ind1+nz+1]=-rh1[0:nz-2]*trans[1:nz-1]

  b=fltarr(nz*2)
  b[nz-1]=t[nz-1]		;sea temperature
  b[nz]=t[0]			;sky temperature

  ;emission temperature of each layer:
  b[ind1+1]=tem[1:nz-1]		;upwelling
  b[ind1+nz+1]=tem[1:nz-1]	;downwelling

  ;solve the sparse matrix:
  ;(using non-sparse methods--note to self...)
  a=a1
  ludc, a, index
  tbv=lusol(a, index, b)
  a=a2
  ludc, a, index
  tbh=lusol(a, index, b)

  if arg_present(weight_h) or arg_present(avet_h) then begin
    trans[0]=0
    ;weight_h=invert(a2)*rebin([1-trans, 1-trans], 2*nz, 2*nz)
    ;w=weight_h[1, [lindgen(nz-2)+1, lindgen(nz-2)+nz+1]]
    ;if keyword_set(raw_weight) ne 1 then begin
      ;weight_h=reform(weight_h[1, lindgen(nz-1)+1]+weight_h[1, lindgen(nz-1)+nz+1])
    ;endif
    ;remove top and bottom layers from calculation:
    ;t1=t[1:nz-2]
    ;avet_h=total(w*[t1, t1])/total(w)
    ;stop

    ;since matrix inversion doesn't work, we simply set each temperature in turn to 1
    ;all others to 0 and solve:
    weight_h=fltarr(nz*2)
    for i=0, nz-1 do begin
      bt=fltarr(nz*2)
      bt[i]=1-trans[i]
      bt[i+nz]=1-trans[i]
      a=a2
      ludc, a, index
      weight_h[i]=(lusol(a, index, bt))[0]
      print, lusol(a, index, bt)
    endfor
   
  endif

  if arg_present(weight_v) or arg_present(avet_v) then begin
    trans[0]=0
    weight_v=invert(a1)*rebin([1-trans, 1-trans], 2*nz, 2*nz)
    w=weight_v[1, [lindgen(nz-2)+1, lindgen(nz-2)+nz+1]]
    if keyword_set(raw_weight) ne 1 then begin
      weight_v=reform(weight_v[1, lindgen(nz-1)+1]+weight_v[1, lindgen(nz-1)+nz+1])
    endif
    ;remove top and bottom layers from calculation:
    t1=t[1:nz-2]
    avet_v=total(w*[t1, t1])/total(w)
    ;stop
  endif

  ;stop
end

;dirt simple plane parallel radiative transfer model:
;nu	= frequency
;zangle = zenith angle
;dz	= layer thickness
;t	= temperature (only used for ice)
;type	= 0=air, 1=snow, 2=compacted snow, 3=fy, 4=my, 5=water
;		(currently only snow vs. ice distinquished
;rho	= density (currently only used for snow)
;s	= salinity (currently only used for ice)
;bound	= boundaries between distinct layers (indices)

pro plane_parallel1, nu, zangle, dz, t, type, rho, s, bound, $
		tbv, tbh, theta=theta, eps=eps, mix=model, ice_type=code, $
		no_rat=no_rat, eps_calc=eps_func
  ;first we calculate the complex permittivities of each layer:
  nz=n_elements(dz)
  eps=complexarr(nz)

  if n_elements(code) eq 1 then code=replicate(code, nz)

  ;ice:
  ind=where(type eq 3 or type eq 4, cnt)
  if cnt ne 0 then begin
    if n_elements(eps_func) ne 1 then begin
      case model of
        "vant": eps[ind]=eps_ice_vant(nu, t[ind], s[ind])
        "shokr": for i=0, cnt-1 do begin
	      	eps_ice_shokr, nu, t[ind[i]], s[ind[i]], code[ind[i]], eps1, no_rat=no_rat
		eps[ind[i]]=eps1
              endfor
      endcase
    endif else begin
      eps[ind]=call_function(eps_func, nu, t[ind], s[ind], code[ind])
    endelse
  endif
  ;snow:
  ind=where(type eq 1 or type eq 2, cnt)
  if cnt ne 0 then eps[ind]=complex(1+1.6*rho[ind]+1.86*rho[ind]^3, 0.0001)

  ;air and water:
  ind=where(type eq 5, cnt)
  if cnt ne 0 then for i=0, cnt-1 do $
		eps[ind[i]]=eps_oc2(nu, t[ind[i]], s[ind[i]])
  ;use a fixed value for ocean water:
  ;if cnt ne 0 then eps[ind]=complex(75., 45.)
  ind=where(type eq 0, cnt)
  if cnt ne 0 then eps[ind]=complex(1, 0)

  ;print, t
  ;print, s

  print, eps

  ;stop

  ppeps, nu, zangle, dz, t, eps, bound, tbv, tbh, theta=theta

end
  
;simply prepares the profile: (input taken from Rasmus's thermodynamic model)
pro plane_parallel, nu, zangle, dz1, t1, type1, rho1, s1, tbv1, tbh1, $
		theta=theta, allbound=allbound, eps=eps, tbsky=tbsky, $
		mix=model, ice_type=code, no_rat=no_rat, eps_calc=eps_func

;+*****************************************************************************
;		PLANE_PARALLEL
;******************************************************************************
;
; usage:	plane_parallel, nu, zangle, dz, t, type, rho, tbv, tbh $
;			[, mix=mix_model] [, /allbound] [, ice_type=ice_type]
;
; purpose:	Plane-parallel radiative-transfer model with discontinuous interfaces.
;		Given an ice vertical profile, calculates emitted brightness 
;		temperatures.  Main purpose is to prepare profile and calculate
;		complex permittivities.  Taylored to take output from Rasmus's
;		thermodynamic ice growth model.  Unlike above routine, top and
;		bottom layers (air and water) should not be included.
;
; input:	nu	Frequency in GHz
;
;		zangle	Zenith angle in degrees.
;
;		dz	Layer thickneses in m.  Should be from top to bottom.
;
;		t	Temperatures in K.
;
;		type	Ice or snow type.  Listed above.
;
;		rho	Ice or snow density.
;
; output:	tbv	Vertically polarized brightness temperature.  
;
; 		tbh	Horizontally polarized brightness temperature.  
;
; keywords:	allbound	Discontinuous boundary between all interfaces.
;			Default is strictly ice-air and water-ice.
;
;		mix	Mixture model.  Default is Vant model.
;			(Vant et al. (1978) J. Appl. Phys. 49 (3): 1264)
;			"Shokr" also supported (Shokr TGARS (1998) 36 (2): 463)
;
;		ice_type:	Ice classification.  Determines whether
;			to use spheres, oriented needles or random needles in
;			Shokr mixture model.  
;
; dependencies:	ice_model_util.pro, eps_calc.pro
;
; author: 	Peter Mills (pmills@iup.physik.uni-bremen.de)
;
; history:	Formally documented 2009-6-02
;
;-*****************************************************************************

  nz=n_elements(eps)
  if n_elements(tbsky) ne 1 then tbsky=5.
  if n_elements(model) ne 1 then model="vant"

  dz=[1, dz1, 1]
  t=[tbsky, t1, 271]
  type=[0, type1, 5]
  s=[0, s1, 35.]
  rho=[0.00124, rho1, 1.]
  nz=n_elements(dz)

  if keyword_set(allbound) then begin
    bound=lindgen(nz-1)
  endif else begin
    bound=where(type[0:nz-2] ne type[1:nz-1])
  endelse

  plane_parallel1, nu, zangle*!dtor, dz, t, type, rho, s, bound, tbv, tbh, $
		theta=theta, eps=eps, mix=model, ice_type=code, no_rat=no_rat, $
		eps_calc=eps_func

  tbv1=tbv[0]
  tbh1=tbh[0]

end
  

pro test_pp, plot=doplot
  nz=5
  dz=fltarr(nz)+0.2
  t=263+findgen(nz)*2
  t[0]=5.
  type=[0, 1, 3, 3, 5]
  rho=[0, 0.3, 0.9, 0.9, 1]
  s=[0, 3, 6, 6, 35.]
  bound=[0, 1, 3]
  ;bound=[0, 1, 2, 3]

  plane_parallel1, 1.4, !dtor*50., dz, t, type, rho, s, bound, tbv, tbh, $
		theta=theta
  theta=float(theta)

  if n_elements(doplot) ne 0 then begin
    ;grey=rgbtoindex(150, 150, 150)
    ps_start, doplot, /landscape
    grey=200

    z=-[0, total(dz, /cum)]
    thick=5
    plot, findgen(10), /nodata, xrange=[0, 1.5], yrange=[-1, 0], $
		xstyle=4, ystyle=4, thick=thick
    for i=1, nz-1 do begin
      oplot, [0.2, 1], [1, 1]*z[i], thick=thick
    endfor

    x=0.8-[0, total(dz*sin(theta), /cum)]

    oplot, x, z, thick=15, color=grey

    cs=2.5
    ct=5
    arrow, 0.3, z[2]+(z[3]-z[2])/3, 0.3, z[2], /data, thick=3
    arrow, 0.3, z[2]+2*(z[3]-z[2])/3, 0.3, z[3], /data, thick=3

    xyouts, 0.25, (z[2]+3*(z[3]-z[2])/5), textoidl("\Delta z_i"), $
		charsize=cs, charthick=ct
   
    xyouts, 1.05, z[2], textoidl("R_{i-1}"), charsize=cs, charthick=ct
    xyouts, 1.05, z[3], textoidl("R_i"), charsize=cs, charthick=ct
    xyouts, 0.8, (z[2]+z[3])/2, textoidl("\epsilon_i, T_i"), $
		charsize=cs, charthick=ct

    xyouts, 0.6, z[3]+dz[2]*0.76, textoidl("T_i\uparrow"), $
		charsize=cs, charthick=ct
    xyouts, 0.6, z[3]+dz[2]*0.1, textoidl("T_i\downarrow"), $
		charsize=cs, charthick=ct
    ps_end
  endif

  stop
end


