
;complex permittivity of brine:
function eps_br, nu, t, s
;+_____________________________________________________________________________
;
; usage:	eps=eps_br(nu, t, s)
;
; purpose:	Calculates complex relative permittivity of brine.
;		From Ulaby et al., Mircrowave Remote Sensing:
;		Active and Passive, Voume III (1986)
;
; input:	nu	Frequence in GHz
;
;		t	Temperature in degrees Celsius
;
;		s	Salinity in psu ~= ppt
;-_____________________________________________________________________________

  nu1=nu*1.e9
  eps0=8.85418782e-12

  ;normality of the solution:
  nb=0.9141*s*(1.707e-2 + s*(1.205e-5 + 4.058e-9*s))
  delta=25.-t

  ;real permittivity of pure water in the large frequency limit:
  eps_inf=4.9

  ;real permittivity of pure water in the stationary limit:
  epsw0=88.045 - t*(0.4147 + t*(6.295e-4 + t*1.075e-5))

  taub02pi=1.1109e-10+T*(-3.824e-12 + T*(6.938e-14 - T*5.096e-16))

  a1=1.0 + nb*(-0.255 + nb*(5.15e-2 - 6.89e-3*nb))
  b1=1.0 + nb*(0.146e-2*t - 4.89e-2 + nb*(-2.97e-2 + 5.64e-3*nb))
  c1=1. + delta*(-1.96e-2 + 8.08e-5*delta) - $
		nb*delta*(3.02e-5 + 3.92e-5*delta + $
		nb*(1.72e-5 - 6.58e-6*delta))


  ;real permittivity of brine in the stationary limit:
  epsb0=epsw0*a1

  tau2pi=taub02pi*b1

  ;conductivity:
  sigma25=nb*(10.39 + nb*(-2.378 + nb*(0.683 + nb*(-0.135 + 1.01e-2*nb))))
  sigma=sigma25*c1

  tmp1=(epsb0 - eps_inf)/(1+(tau2pi*nu1)^2)
  eps_re=eps_inf+tmp1
  eps_im=tau2pi*nu1*tmp1+sigma/2/!pi/nu1/eps0

  ;stop

  return, complex(eps_re, eps_im)

end

;copied more or less verbatim from Fortran SFT routines:
function eps_oc2, fnu, t, s
      TC=T-abs_zero()
      EP0=.8774D2-TC*(.40008-TC*(.9398D-3-.141D-5*TC))
      DELAM=FNU*(.111093-TC*(.38236D-2-TC*(.69375D-4-.50963D-6*TC)))
      IF S EQ 0.D0 then GOTO, lab1
      XN=S*(.1707D-1+S*(.12053D-4+.40575D-8*S))
      F=1.D0-XN*(.255059-XN*(.51514D-1-.68889D-2*XN))
      EP0=EP0*F
      F=1.D0-XN*(.48963D-1+XN*(.29667D-1-.56444D-2*XN))
      DELAM=DELAM*(F+.1463D-2*XN*TC)
lab1:     EP1=(EP0-.49D1)/(1.D0+DELAM*DELAM)
      B=DELAM*EP1
      IF S EQ 0.D then GOTO, lab2
      SIG=S*(.1825208-S*(.146192D-2-S*(.209324D-4-.128205D-6*S)))
      TAU=.25D2-TC
      A=.203318D-1+TAU*(.12664D-3+.24637D-5*TAU)-S*(.184911D-4-TAU*(0.25506D-6-.25506D-7*TAU))
      SIGMA=SIG*EXP(-TAU*A)
      B=B+SIGMA/(.55632485D-1*FNU)
lab2:     CKW=complex(.49D1+EP1,B)


  return, ckw

end


;salinity of brine as function of temperature:
function brine_salinity, t
  comp1=t gt -43.2
  comp2=t gt -36.8
  comp3=t gt -22.9
  comp4=t gt -8.2
  comp5=t gt -2.

  dim=size(t, /dim)
  if n_elements(dim) eq 1 and dim[0] eq 0 then begin
    sb=0.
  endif else begin
    sb=make_array(/float, size(t, /dim))-1
  endelse

  ind=where(comp5, cnt)
  if cnt ne 0 then begin
    s1=0.
    s2=37.6514
    t1=0.
    t2=-2.
    sb[ind]=s1+(t[ind]-t1)/(t2-t1)*(s2-s1)
  endif

  ind=where(comp1 and not comp2, cnt)
  if cnt ne 0 then begin
    t1=t[ind]
    sb[ind]=508.18+t1*(14.535 + t1*0.2018)
  endif

  ind=where(comp2 and not comp3, cnt)
  if cnt ne 0 then begin
    t1=t[ind]
    sb[ind]=242.94 + t1*(1.5299 + 0.0429*t1)
  endif

  ind=where(comp3 and not comp4, cnt)
  if cnt ne 0 then begin
    t1=t[ind]
    sb[ind]=57.041 - t1*(9.929 + t1*(0.16204 + 0.002396*t1))
  endif

  ind=where(comp4 and not comp5, cnt)
  if cnt ne 0 then begin
    t1=t[ind]
    sb[ind]=1.725 - t1*(18.756 + t1*0.3964)
  endif

  return, sb

end

function brine_volume_from_salinity, t, s, sb
  rhoi=0.917-1.403e-4*t
  sb=brine_salinity(t)
  rhob=1.+0.0008*sb
  vb=s*rhoi/(sb*rhob+s*(rhoi-rhob))
  ;stop
  return, vb
end

;brine volume as a function of temperature and salinity:
function brine_volume, t, s, basic=basic
  ;t=t-abs_zero()

  if keyword_set(basic) then return, s*(-49.185/t+0.532)/1000

  comp1=t lt -0.5
  comp2=t lt -2.06
  comp3=t lt -8.2
  comp4=t ge -22.9

  n=n_elements(t)
  vb=fltarr(n)-1

  ind=where(comp1 and comp2 ne 1, cnt)
  if cnt ne 0 then vb[ind]=s[ind]*(-52.56/t[ind]+2.28)
  ;if cnt ne 0 then vb[ind]=s[ind]*(-52.56/t[ind])
  ind=where(comp2 and comp3 ne 1, cnt)
  if cnt ne 0 then vb[ind]=s[ind]*(-45.917/t[ind]+0.930)
  ind=where(comp3 and comp4, cnt)
  if cnt ne 0 then vb[ind]=s[ind]*(-43.795/t[ind]+1.189)

  return, vb/1000.

end



