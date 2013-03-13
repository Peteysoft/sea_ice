@time_class

;calculates relative (must be multiplied by the solar "constant")
;clear-sky short-wave flux based on latitude and date

function sw_flux, psi1, date1, avg=avg, lon=lon1

  ;Earth's tilt:
  gamma=23.44*!dtor
  ;if longitude is given, we assume UTC
  ;otherwise we assume local time
  if n_elements(lon1) eq 0 then lon=0. else lon=!dtor*lon1

  date=time_class_convert(date1)
  psi=!dtor*psi1

  ;angle representing time of year:
  ;date0={time_str, date.year-1, 12, 22, 0, 0, 0}	;my birthday
  ;alpha=2*!pi*day_value(date, date0)/365.4

  alpha=2*!pi*(date.dy-355.)/366.0

  cosalph=cos(alpha)
  singam=sin(gamma)
  cosgam=cos(gamma)

  denom=sqrt(cosalph*cosalph*singam*singam+cosgam*cosgam)

  if keyword_set(avg) then begin
    ;daily average:
    ;angle of sunrise and sunset:
    costheta0=tan(psi)*cosalph*tan(gamma)
    ind1=where(abs(costheta0) gt 1, n1, compl=ind2, ncomp=n2)
    q=costheta0
    if n1 gt 0 then q[ind1]=0.
    if n2 gt 0 then begin
      theta0=acos(costheta0[ind2])
      q[ind2]=(cos(psi[ind2])*cosgam[ind2]*sin(theta0)-sin(psi)*cosalph*singam*theta0)/denom[ind2]/!pi
    endif
  endif else begin
    ;angle representing time of day:
    ;theta=!pi*(date.hour+(date.minute+date.second/60.)/60.-12)/12
    ;theta=2*!pi*(date.dy-long(date.dy)-0.5)
    ;assume UTC (correct for local time):
    theta=2*!pi*(date.dy-long(date.dy)-0.5)+lon
    q=(cos(psi)*cos(theta)*cosgam-sin(psi)*cosalph*singam)/denom
  endelse

  ;stop

  ind=where(q gt 1, cnt)
  if cnt gt 0 then q[ind]=1.
  ind=where(q lt 0, cnt)
  if cnt gt 0 then q[ind]=0.

  return, q

end

