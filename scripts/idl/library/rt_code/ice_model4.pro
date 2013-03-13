@ice_scatter_util

;this version allows for variations in complex permittivity
;across the ice sheet

function error_error, x, params=params
  return, errorf(x)-params
end


;2009-02-25: instrument fov now represented by a Gaussian
;2009-08-08: corrected calculation of fov deviates
function ice_model, s, h, z, x0, y0, zangle, fov, err=err, $
		tsky=tsky, tice=tice, tw=tw, nsea=nsea, $
		niter=niter, plot=doplot, pol=pol, fail=fail, $
		n1=n1, n2=n2, freq=freq, ds=ds, offset=offset, $
		uniform_deviates=uniform_deviates, $
		xall=xall, yall=yall

common im_random_seed, seed

;error handler:
catch, err_code
if err_code ne 0 then begin
  message, !error_state.msg, /continue
  fail=err_code
  catch, /cancel
  ;stop
  return, -1
endif

if n_elements(tsky) eq 0 then tsky=5.		;sky temperature
if n_elements(tice) eq 0 then tice=269.		;ice temperature
if n_elements(tw) eq 0 then tw=271.		;water temperature
if n_elements(nsea) eq 0 then nsea=complex(9, 2.5)	;refractive index of sea water
;(since we never advance the ray further once it hits the water, 
;we can use a complex refractive index without worrying about its implications...)

if n_elements(offset) ne 1 then offset=0.
if n_elements(fov) ne 1 then fov=13.16*!dtor

if n_elements(ds) eq 0 then ds=0.1		;integration step

;special case: constant refractive index
if n_elements(n1) eq 1 then n1=replicate(n1, n_elements(s))
if n_elements(n2) eq 1 then n2=replicate(n2, n_elements(s))

zp=spl_init(s, z, double=doub)
hp=spl_init(s, h, double=doub)

dum1=where(finite(zp) ne 1, cnt1)
dum2=where(finite(zp) ne 1, cnt2)
;will introduce a selection bias, but right now we can't fix the
;the fitting problem...
if cnt1 ne 0 or cnt2 ne 0 then begin
  message, "Cubic spline fitting failed", /con
  fail=1
  return, -1
endif

fail=0

;we need to plot it to make sure things are working right:
if keyword_set(doplot) then begin
  s0=x0+tan(zangle)*y0
  xarr=s0+3-findgen(100)/10
  top=spl_interp(s, h, hp, xarr)
  bottom=spl_interp(s, z, zp, xarr)

  ;ps_start, "test_scatter.ps", /landscape

  plot, xarr, top, yrange=[-2, 1], thick=3
  oplot, xarr, bottom, thick=3
endif

if arg_present(xall) then begin
  ;for creating animations:
  xall=fltarr(maxn, niter)
  yall=fltarr(maxn, niter)
  xall[0, iter]=x0
  yall[0, iter]=y0
endif

;stop

if n_elements(niter) ne 1 then niter=1000
maxn=100

;approximate machine accuracy, 
;(!!!) not to be confused with complex permittivity (!!!)
eps=1e-7

bt=fltarr(niter)
nb=lonarr(niter)
sf=bytarr(niter)

maxh=max(h, min=minh)
hrange=[minh, maxh]
maxz=max(z, min=minz)
zrange=[minz, maxz]

maxht2=(maxh+0.1)

;extract real parts of refractive indices and calculate
;absorption coefficients:
n1r=float(n1)
n2r=float(n2)
c=0.3
abs1=4*!pi*imaginary(n1)*freq/c
abs2=4*!pi*imaginary(n2)*freq/c

t=fltarr(maxn)
r=fltarr(maxn)
stype=bytarr(maxn+1)

xarr=fltarr(maxn+1)
yarr=fltarr(maxn+1)

;minimum value for un-transformed deviate:
minrnum1=errorf(asin(offset)/sqrt(2)/fov)

for iter=0L, niter-1 do begin

  if keyword_set(uniform_deviates) then begin
    dangle=fov*(2*randomn(seed, 1)-1)
  endif else begin
    ;we find deviates by inverting the error function:
    rnum1=minrnum1+randomu(seed, 1)*(1-minrnum1)
    rnum=bisect("error_error", -10000., 10000., params=rnum1)*sqrt(2)*fov
    ;desc=0. > (rnum^2-offset^2)
    sgn=fix(randomu(seed, 1)*2)*2-1
    dangle=sgn*asin(sqrt((sin(rnum))^2-offset^2))
    ;print, dangle/!dtor
    ;stop
  end

  angle0=zangle+dangle
  ;print, angle0
  angle0=angle0[0]
  vx=sin(angle0)
  vy=-cos(angle0)

  ;to avoid unnecessary computation, we advance the ray
  ;until it is within striking distance of the ice surface:
  t0=(maxht2-y0)/vy
  y=maxht2
  x=x0+vx*t0

  xarr[0]=x
  yarr[0]=y

  absarr=fltarr(maxn)

  stype[0]=0		;0 for air, 1 for ice, 2 for sea...

  for i=0, maxn-1 do begin
    if stype[i] eq 0 then begin
      t[i]=spl_intercept(x, y, vx, vy, s, h, hp, $
		xint=xnew, yint=ynew, yrange=hrange)
      if t[i] lt 0 then break
      m=spl_slope(xnew, s, h, hp)
      nice=interpol(n1r, s, xnew)
      ;print, "n=", nice
      r[i]=m_refract(vx, vy, m, 1, nice, vxnew, vynew, p=pol)
      ;print, r[i]
      if randomu(seed, 1) gt r[i] then begin
        ;ray is refracted--now in ice...
        stype[i+1]=1
      endif else begin
        ;ray is reflected--still in air
        stype[i+1]=0
        m_reflect, vx, vy, m, vxnew, vynew
      endelse
      absarr[i]=interpol(abs1, s, xnew)
    endif else if stype[i] eq 1 then begin
      t1=spl_intercept(x, y, vx, vy, s, h, hp, $
		xint=xnew1, yint=ynew1, yrange=hrange)
      t2=spl_intercept(x, y, vx, vy, s, z, zp, $
		xint=xnew2, yint=ynew2, yrange=zrange)
      if t1 lt eps and t2 lt eps then begin
        t[i]=-1
        break
      endif
      if t1 lt eps then begin
        t[i]=t2
        interface=1
      endif else if t2 lt eps then begin
        t[i]=t1
        interface=0
      endif else if t1 lt t2 then begin
        t[i]=t1
        interface=0
      endif else begin
        t[i]=t2
        interface=1
      endelse

      if interface eq 0 then begin
        xnew=xnew1
        ynew=ynew1
        m=spl_slope(xnew, s, h, hp)
        nice=interpol(n1r, s, xnew)
        ;print, "n=", nice
        r[i]=m_refract(vx, vy, m, nice, 1, vxnew, vynew, p=pol)
        if randomu(seed, 1) gt r[i] then begin
          ;ray is refracted--now in air...
          stype[i+1]=0
        endif else begin
          ;ray is reflected--still in ice
          stype[i+1]=1
          m_reflect, vx, vy, m, vxnew, vynew
        endelse
        absarr[i]=interpol(abs1, s, xnew)
      endif else begin
        xnew=xnew2
        ynew=ynew2
        m=spl_slope(xnew, s, z, zp)
        nice=interpol(n2r, s, xnew)
        ;print, "n=", nice
        r[i]=m_refract(vx, vy, m, nice, nsea, vxnew, vynew, p=pol)
	;print, nsea
	;print, abs(r[i])
        if randomu(seed, 1) gt abs(r[i]) then begin
          ;ray is refracted--now in water...
          stype[i+1]=2
        endif else begin
          ;ray is reflected--still in ice
          stype[i+1]=1
          m_reflect, vx, vy, m, vxnew, vynew
        endelse
        absarr[i]=interpol(abs2, s, xnew)
      endelse
    endif else begin
      break
    endelse

    if keyword_set(doplot) then oplot, [x, xnew], [y, ynew]

    x=xnew
    y=ynew
    vx=vxnew
    vy=vynew

    if arg_present(xall) then begin
      xarr[i+1]=x
      yarr[i+1]=y
    endif
  endfor

  if keyword_set(doplot) then oplot, [x, x+vx], [y, y+vy]

  if arg_present(xall) then begin
    xall[i+1, iter]=x+vx
    yall[i+1, iter]=y+vy
  endif


  if i eq 0 then begin
    ;should fix the bug soon, but right now we want to get the thing running...
    message, "Premature ray trace exit", /con
    print, x, s[0], s[n_elements(s)-1]
    iter=iter-1
    ;stop
    continue
  endif

  nb[iter]=i

  if stype[nb[iter]] eq 1 then begin
    message, "Ray trace exit in ice", /con
    print, x, s[0], s[n_elements(s)-1]
  endif

  ;print, r
  ;print, t
  ;print, stype

  ;do the radiative transfer calculation:
  ;print, stype[0:nb[iter]]
  sf[iter]=stype[nb[iter]]
  if sf[iter] eq 0 then bt1=tsky else $
	if sf[iter] eq 1 then bt1=tice else $
	if sf[iter] eq 2 then bt1=tw
  for i=nb[iter]-1, 1, -1 do begin
    if stype[i] eq 1 then begin
      ;nstep=long(t[i]/ds+1)
      ;ds1=t[i]/nstep
      ;this is rather approximate:
      ;abs_coef=(absarr[i-1]-absarr[i])*(findgen(nstep)+0.5)/nstep+absarr[i]
      ;for j=0L, nstep-1 do begin
      ;  bt1=tice-(tice-bt1)*exp(-abs_coef[j]*ds1)
      ;endfor
      abs_coef=(absarr[i-1]+absarr[i])/2
      bt1=tice-(tice-bt1)*exp(-abs_coef*t[i])
      ;stop
    endif
  endfor
  bt[iter]=bt1

endfor

;print, bt

btave=total(bt)/niter
err=stddev(bt)

print, btave, "+/-", err

catch, /cancel

;stop

return, btave

;ps_end

end


