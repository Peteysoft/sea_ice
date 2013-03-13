@ice_scatter_util

function ice_model, s, h, z, x0, y0, zangle, fov, err=err, $
		tsky=tsky, tice=tice, tw=tw, rindex=nice, abs_coef=abs_coef, $
		niter=niter, plot=doplot, pol=pol, fail=fail, nsea=nsea, $
		animate=animate

common ice_model_ran_seed, seed

if n_elements(tsky) eq 0 then tsky=5.		;sky temperature
if n_elements(tice) eq 0 then tice=269.		;ice temperature
if n_elements(tw) eq 0 then tw=271.		;water temperature
if n_elements(nsea) eq 0 then nsea=complex(9, 2.5)	;refractive index of sea water
;(since we never advance the ray further once it hits the water, 
;we can use a complex refractive index without worrying about its implications...)

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
;  xarr=s
  top=spl_interp(s, h, hp, xarr)
  bottom=spl_interp(s, z, zp, xarr)

  ;ps_start, "test_scatter.ps", /landscape

  plot, xarr, top, yrange=[-2, 1], charsize=1.3, xtitle="x", ytitle="z"
  oplot, xarr, bottom
endif

;stop

if n_elements(niter) ne 1 then niter=1000
maxn=100
eps=1e-7

bt=fltarr(niter)
nb=lonarr(niter)
sf=bytarr(niter)

xall=fltarr(maxn, niter)
yall=fltarr(maxn, niter)

maxh=max(h, min=minh)
hrange=[minh, maxh]
maxz=max(z, min=minz)
zrange=[minz, maxz]

maxht2=(maxh+0.1)

frame=0L

if keyword_set(animate) then begin
  spawn, "rm -f mc_frames/mc_frame*"
  tvlct, [250, 100, 50], [0, 200, 250], [0, 250, 250], 1
  red=1
  blue_grey=2
  blue_green=3
  ray_color=red
  ice_air_col=blue_grey
  ice_water_col=blue_green
  surf_thick=5
  ray_thick=2
endif

for iter=0L, niter-1 do begin

  angle0=zangle+2*randomu(seed, 1)*fov-fov
  angle0=angle0[0]
  vx=sin(angle0)
  vy=-cos(angle0)

  ;to avoid unnecessary computation, we advance the ray
  ;until it is within striking distance of the ice surface:
  t0=(maxht2-y0)/vy
  y=maxht2
  x=x0+vx*t0
  ;y=y0
  ;x=x0

  t=fltarr(maxn)
  r=fltarr(maxn)
  stype=bytarr(maxn+1)

  stype[0]=0		;0 for air, 1 for ice, 2 for sea...

  xall[0, iter]=x0
  yall[0, iter]=y0

  for i=0, maxn-1 do begin
    if stype[i] eq 0 then begin
      t[i]=spl_intercept(x, y, vx, vy, s, h, hp, $
		xint=xnew, yint=ynew, yrange=hrange)
      if t[i] lt 0 then break
      m=spl_slope(xnew, s, h, hp)
      r[i]=m_refract(vx, vy, m, 1, nice, vxnew, vynew, p=pol)
      if randomu(seed, 1) gt r[i] then begin
        ;ray is refracted--now in ice...
        stype[i+1]=1
      endif else begin
        ;ray is reflected--still in air
        stype[i+1]=0
        m_reflect, vx, vy, m, vxnew, vynew
      endelse
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
        r[i]=m_refract(vx, vy, m, nice, 1, vxnew, vynew, p=pol)
        if randomu(seed, 1) gt r[i] then begin
          ;ray is refracted--now in air...
          stype[i+1]=0
        endif else begin
          ;ray is reflected--still in ice
          stype[i+1]=1
          m_reflect, vx, vy, m, vxnew, vynew
        endelse
      endif else begin
        xnew=xnew2
        ynew=ynew2
        m=spl_slope(xnew, s, z, zp)
        r[i]=m_refract(vx, vy, m, nice, nsea, vxnew, vynew, p=pol)
        if randomu(seed, 1) gt r[i] then begin
          ;ray is refracted--now in water...
          stype[i+1]=2
        endif else begin
          ;ray is reflected--still in ice
          stype[i+1]=1
          m_reflect, vx, vy, m, vxnew, vynew
        endelse
      endelse
    endif else begin
      break
    endelse

    if keyword_set(doplot) then oplot, [x, xnew], [y, ynew]

    x=xnew
    y=ynew
    vx=vxnew
    vy=vynew
    xall[i+1, iter]=x
    yall[i+1, iter]=y

    if keyword_set(animate) then begin
      s0=x0+tan(zangle)*y0
      xarr=s0+3-findgen(100)/10
      top=spl_interp(s, h, hp, xarr)
      bottom=spl_interp(s, z, zp, xarr)

      fname="mc_frames/mc_frame"+string(frame, format="(i4.4)")
      ps_start, fname+".ps", /landscape, /color

      plot, xarr, top, yrange=[-2, 1], charsize=1.3, xtitle="x", ytitle="z", $
	      		xstyle=5, ystyle=5, /nodata
      oplot, xarr, top, thick=surf_thick, color=ice_air_col
      oplot, xarr, bottom, thick=surf_thick, color=ice_water_col

      for j=0, iter-1 do oplot, xall[0:nb[j]+1, j], yall[0:nb[j]+1, j], $
	      	thick=ray_thick, color=ray_color
      oplot, xall[0:i+1, iter], yall[0:i+1, iter], thick=2, color=ray_color
      ps_end

      spawn, "convert "+fname+".ps -rotate -90 "+fname+".gif"
      frame=frame+1
    endif

  endfor
  nb[iter]=i
  xall[i+1, iter]=x+vx
  yall[i+1, iter]=y+vy

  if keyword_set(doplot) then oplot, [x, x+vx], [y, y+vy]
  if keyword_set(animate) then begin
    s0=x0+tan(zangle)*y0
    xarr=s0+3-findgen(100)/10
    top=spl_interp(s, h, hp, xarr)
    bottom=spl_interp(s, z, zp, xarr)

    fname="mc_frames/mc_frame"+string(frame, format="(i4.4)")
    ps_start, fname+".ps", /landscape, /color

    plot, xarr, top, yrange=[-2, 1], charsize=1.3, xtitle="x", ytitle="z", xstyle=5, ystyle=5, /nodata
    oplot, xarr, top, thick=surf_thick, color=ice_air_col
    oplot, xarr, bottom, thick=surf_thick, color=ice_water_col
    for j=0, iter do oplot, xall[0:nb[j]+1, j], yall[0:nb[j]+1, j], $
	    	color=ray_color, thick=ray_thick
    ps_end

    spawn, "convert "+fname+".ps -rotate -90 "+fname+".gif"
    frame=frame+1
  endif


  if i eq 0 then begin
    ;should fix the bug soon, but right now we want to get the thing running...
    message, "Premature ray trace exit", /con
    print, x, s[0], s[n_elements(s)-1]
    iter=iter-1
    ;stop
    continue
  endif

  r=r[0:nb[iter]-1]
  t=t[0:nb[iter]-1]
  stype=stype[0:nb[iter]]

  if stype[nb[iter]] eq 1 then begin
    message, "Ray trace exit in ice", /con
    print, x, s[0], s[n_elements(s)-1]
  endif

  ;print, r
  ;print, t
  ;print, stype

  ;do the radiative transfer calculation:
  sf[iter]=stype[nb[iter]]
  if sf[iter] eq 0 then bt[iter]=tsky else $
	if sf[iter] eq 1 then bt[iter]=tice else $
	if sf[iter] eq 2 then bt[iter]=tw
  for i=nb[iter]-1, 0, -1 do begin
    ;should be accounted for implicitly in the Monte-Carlo simulation?
    ;if stype[i+1] eq stype[i] then bt[iter]=bt[iter]*r[i] else bt[iter]=bt[iter]*(1-r[i])
    if stype[i] eq 1 then begin
      bt[iter]=tice-(tice-bt[iter])*exp(-abs_coef*t[i])
    endif
  endfor

endfor

if keyword_set(animate) then begin
  spawn, "convert -delay 2 mc_frames/mc_frame*.gif mc_animate.gif"
  stop
endif

btave=total(bt)/niter
err=stddev(bt)

print, btave, "+/-", err

;stop

return, btave

;ps_end

end


