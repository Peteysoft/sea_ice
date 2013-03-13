
;find the roots of a cubic:
function solve_cubic, a, b, c, d
  a2=a*a
  a3=a2*a
  q=(3*a*c-b^2)/9/a2
  r=(9*a*b*c-27*a2*d-2*b*b*b)/54/a3

  q3=q*q*q
  r2=r*r
  disc=q3+r2

  if disc ge 0 then begin
    i1=sqrt(q3+r2)
    s=(r+i1)^(1./3)
    t=(r-i1)^(1./3)
  endif else begin
    r1trd=r^(1./3)
    rtdiscover3=sqrt(-disc)/3
    s=complex(r1trd, rtdiscover3)
    t=complex(r1trd, -rtdiscover3)
  endelse

  i=complex(0, 1)
  x1=s+t-b/3/a
  x2=-(s+t)/2-b/3/a+sqrt(3)*(s-t)*i/2
  x3=-(s+t)/2-b/3/a-sqrt(3)*(s-t)*i/2

  return, [x1, x2, x3]

end

;find the slope of a cubic spline:
function spl_slope, x0, x, y, yp
  n=n_elements(x)
  ;xind=interpol(findgen(n), x, x0)
  lind=value_locate(x, x0)
  ;frac1=xind-lind
  h=x[lind+1]-x[lind]
  frac1=(x0-x[lind])/h
  frac2=1-frac1
  result=h*((frac1*frac1-1./3)*yp[lind+1]-(frac2*frac2-1./3)*yp[lind])/2 + $
		(y[lind+1]-y[lind])/h
  return, result
end

;find the nearest intercept of a line with a cubic spline:
function spl_intercept, x0_1, y0_1, vx, vy, x, y, yp, $
		xint=xint, yint=yint, yrange=yrange
  n=n_elements(x)

  if n_elements(yrange) ne 2 then begin
    miny=min(y)
    maxy=max(y)
  endif else begin
    miny=yrange[0]
    maxy=yrange[1]
  endelse

  ind=value_locate(x, x0_1)
  if ind eq n-1 then ind=n-2

  dir=vx/abs(vx)
  if dir lt 0 then iend=0 else iend=n-2

  eps=1e-6
  eps2=1e-10

  t=-1

  for i=ind, iend, dir do begin
    h=x[i+1]-x[i]
    if x0_1 lt x[i] or x0_1 gt x[i+1] then begin
      if dir lt 0 then begin
        t0=(x[i+1]-x0_1)/vx
        x0=h
        y0=y0_1+vy*t0
      endif else begin
        t0=(x[i]-x0_1)/vx
        x0=0
        y0=y0_1+vy*t0
      endelse
    endif else begin
      x0=x0_1-x[i]
      t0=0.
      y0=y0_1
    endelse
    ;check to see if intersections are even possible:
    ;t1=-x0/vx
    ;t2=(h-x0)/vx
    ;if t1 gt 0 then t=t1 else t=t2
    ;yt=vy*t+y0

    if (y0 gt maxy and vy gt 0) or (y0 lt miny and vy lt 0) then begin
      t=-1
      break
    endif

    ;generate polynomial coefficients:
    y1=y[i]
    y2=y[i+1]
    yp1=yp[i]
    yp2=yp[i+1]
    x0p2=x0*x0
    x0p3=x0p2*x0
    vx2=vx*vx
    vx3=vx2*vx
    h2=h*h

    ;this is infuriating...
    ;why don't I just solve the fucking thing directly??
    d=((x0p3-h2*x0)*yp2+(-x0p3+3*h*x0p2-2*h2*x0)*yp1+ $
		6*x0*y2+(6*h-6*x0)*y1-6*h*y0)/h/6

    c=((3*vx*x0p2-h2*vx)*yp2+(-3*vx*x0p2+6*h*vx*x0-2*h2*vx)*yp1+ $
		6*vx*y2-6*vx*y1-6*h*vy)/h/6

    b=vx2*(3*x0*yp2+(3*h-3*x0)*yp1)/h/6

    a=vx3*(yp2-yp1)/h/6

    ;print, "Coefficients:", a, b, c, d

    if abs(a) lt eps2 then begin
      if abs(b) eq 0 then begin
        t=-d/c
      endif else begin
        desc=c^2-4*b*d
        t=[-c+sqrt(complex(desc, 0)), -c-sqrt(complex(desc, 0))]/2/b
      endelse
    endif else begin
      t=fz_roots([d, c, b, a])
      ;t=solve_cubic(a, b, c, d)
    endelse
    ;print, "Roots:", t

    ;stop

    tf=float(t)
    if t0 eq 0. then begin
      gind=where(abs(imaginary(t)) lt eps and tf gt eps, cnt)
    endif else begin
      gind=where(abs(imaginary(t)) lt eps and tf ge 0., cnt)
    endelse

    if cnt eq 0 then begin
      t=-1		;continue may knock us out of the loop...
      continue
    endif

    t=min(tf[gind], k)
    xt=vx*t+x0
    ;print, h, xt
    if xt gt h or xt lt 0 then begin
      t=-1
      continue
    endif

    break
  endfor

  if t gt 0 then begin
    xint=x[i]+xt
    yint=spl_interp(x, y, yp, xint)
    t=t+t0
  endif

  ;print, ind
  ;stop

  return, t

end

;calculate the direction vector after refraction:
function m_refract, vx, vy, m, n1, n2, vxnew, vynew, pol=pol
  ;tangent vector:
  num=sqrt(1+m^2)
  tx=1./num
  ty=m/num
  
  ;project the slope onto the tangent vector:
  proj=vx*tx+vy*ty

  ;length of component along tangent:
  tanlen=proj*n1/n2

  ;if this is greater than one, we have total internal reflection:
  if abs(tanlen) gt 1 then begin
    m_reflect, vx, vy, m, vxnew, vynew
    return, 1
  endif

  ;length of component along normal:
  normlen1=-vx*ty+vy*tx
  normlen=-sqrt(1-tanlen^2)*normlen1/abs(normlen1)
  ;(this subroutine is highly redundant)

  vxnew=tanlen*tx+normlen*ty
  vynew=-normlen*tx+tanlen*ty

  if finite(vxnew) ne 1 then stop

  r=sqrt(vxnew*vxnew+vynew*vynew)
  ;print, r
  vxnew=vxnew/r
  vynew=vynew/r		;(did I mention that this function is quite redundant??)

  cost0=abs(normlen1)
  cost1=abs(normlen)

  if pol eq "h" then begin
    r=(n1*cost0-n2*cost1)/(n1*cost0+n2*cost1)
  endif else begin
    r=(n1*cost1-n2*cost0)/(n1*cost1+n2*cost0)
  endelse

  ;stop

  return, r^2

end

;calculate the direction vector after reflection:
pro m_reflect, vx, vy, m, vx2, vy2
  ;tangent vector:
  num=sqrt(1+m^2)
  tx=1./num
  ty=m/num

  proj1=vx*tx+vy*ty
  proj2=-vx*ty+vy*tx

  ;print, tx, ty
  ;print, proj1, proj2
  ;print

  vx2=proj1*tx+proj2*ty
  vy2=proj1*ty-proj2*tx

end

;calculate the direction vector after reflection:
pro m_reflect1, vx, vy, m, vx2, vy2
  ;normal vector:
  num=sqrt(1+m^2)
  nx=m/num
  ny=-1./num

  proj=vx*nx+vy*ny

  vx2=vx-2*proj*nx
  vy2=vy-2*proj*ny

end

pro test_reflect, vx, vy
  r=sqrt(vx*vx+vy*vy)
  vx=vx/r
  vy=vy/r
  for theta=0, 359 do begin
    m=tan(theta*!dtor)
    m_reflect, vx, vy, m, vx2, vy2
    plot, [0, -vx], [0, -vy], xrange=[-1, 1], yrange=[-1, 1]
    oplot, [0, cos(theta*!dtor)], [0, sin(theta*!dtor)]
    oplot, [0, vx2], [0, vy2]
    dum=get_kbrd(1)
  endfor

end

pro test_refract, vx, vy, n1, n2
  r=sqrt(vx*vx+vy*vy)
  vx=vx/r
  vy=vy/r
  for theta=0, 359 do begin
    m=tan(theta*!dtor)
    r=m_refract(vx, vy, m, n1, n2, vx2, vy2, pol="v")
    ;print, r
    plot, [0, -vx], [0, -vy], xrange=[-1, 1], yrange=[-1, 1]
    oplot, [0, cos(theta*!dtor)], [0, sin(theta*!dtor)]
    oplot, [0, vx2], [0, vy2]
    dum=get_kbrd(1)
  endfor

end

pro test_refract2
  nth=90
  theta=findgen(nth)*!dtor

  m=0.
  n=complex(2, 2)

  rv1=complexarr(nth)
  rh1=complexarr(nth)
  rv2=complexarr(nth)
  rh2=complexarr(nth)

  for i=0, nth-1 do begin
    thetai=theta[i]

    ;check against explicitly calculated from fresnel eqns:
    thetat=asin(sin(thetai)/n)
    print, thetat
    rv1[i]=(n*cos(thetai)-cos(thetat))/(n*cos(thetai)+cos(thetat))
    rh1[i]=(n*cos(thetat)-cos(thetai))/(n*cos(thetat)+cos(thetai))

    vx=-sin(thetai)
    vy=-cos(thetai)
    rv2[i]=m_refract(vx, vy, m, 1, n, vxnew, vynew, p="v")
    ;real only part of refraction vector vs. magnitude:
    rh2[i]=m_refract(vx, vy, m, 1, n, vxnew, vynew, p="h")
  endfor

  window, 0
  plot, theta, rv1^2
  oplot, theta, rv2
  oplot, theta, rh1^2
  oplot, theta, rh2

  window, 1
  plot, rv1^2, rv2, psym=2
  window, 2
  plot, rh1^2, rh2, psym=2

  stop

end

