@ice_model_util

;mixture model for complex permittivity based on Shokr (1998)
;TGARS 36 (2): 463
pro eps_mix_shokr, eps_pi, eps_brn, v, type, eps_x, eps_z, $
		no_rationalize=no_rat


  ;if v lt 0.1 or keyword_set(no_rat) then begin
  if keyword_set(no_rat) then begin
    if type eq 1 or type eq 2 then begin
      eps_x=eps_pi*(1+3*v*(eps_brn-eps_pi)/ $
		(eps_brn+2*eps_pi))
      eps_z=eps_x
    endif else if type eq 3 then begin
      eps_x=eps_pi*(1+2*v*(eps_brn-eps_pi)/ $
		(eps_brn+eps_pi))
      eps_z=eps_pi+v*(eps_brn-eps_pi)
    endif else if type eq 4 or type eq 5 or type eq 6 then begin
      eps_x=eps_pi+v*(eps_brn-eps_pi)/ $
		(eps_brn+eps_pi)/3*(5*eps_pi+eps_brn)
      eps_z=eps_x
    endif else begin
      ;simplest possible mixture formula:
      eps_x=eps_pi+v*(eps_brn-eps_pi)
      eps_z=eps_x
    endelse

  endif else begin
    eps_re=float(eps_brn)
    eps_im=imaginary(eps_brn)

    if type eq 1 or type eq 2 then begin
      a=-2.
      b=2*eps_pi-eps_brn+3*v*(eps_brn-eps_pi)
      c=eps_pi*eps_brn
      d=sqrt(b^2+4*a*c)
      ;the real part can take on both positive and negative sign,
      ;the imaginary part cannot...
      d=complex(float(d), -imaginary(d))
      b=complex(-float(b), imaginary(b))
      ;seems an ugly hack...
      eps_x=(-b+d)/2/a
      ;stop
      eps_z=eps_x
    endif else if type eq 3 then begin
      a=-1.
      b=eps_brn-eps_pi+2*v*(eps_brn-eps_pi)
      c=eps_pi*eps_brn
      d=sqrt(b^2+4*a*c)
      d=complex(float(d), -imaginary(d))
      b=complex(-float(b), imaginary(b))
      eps_x=(b+d)/2/a
      ;stop
      eps_z=eps_pi+v*(eps_brn-eps_pi)
    endif else if type eq 4 or type eq 5 or type eq 6 then begin
      a=3.
      b=3*(eps_pi-eps_brn)+5*v*(eps_brn-eps_pi)
      c=3*eps_brn*eps_pi+v*(eps_brn-eps_pi)*eps_brn
      d=sqrt(b^2+4*a*c)
      d=complex(-float(d), -imaginary(d))
      b=complex(float(b), imaginary(b))
      eps_x=(b-d)/2/a
      ;stop
      eps_z=eps_x
    endif else begin
      ;simplest possible mixture formula:
      eps_x=eps_pi+vb*(eps_brn-eps_pi)
      eps_z=eps_x
    endelse
  endelse

end

;calculates the complex permittivity of sea ice using
;models by Shokr (1998):
;type codes correspond to those of Eicken and Lange (1989)
;JGR 94 (C6): 8193
;1=granular polygonal (spherical inclusions)
;2=granular orbicular (spherical inclusions)
;3=columnar (vertically oriented needles)
;4=intermediate (randomly oriented needles)
;5=mixed (randomly oriented needles)
;6=platelet (randomly oriented needles)

pro eps_ice_shokr, nu, t, s, type, eps_x, eps_z, $
	eps_brn=eps_brn, eps_pi=eps_pi, $
		no_rationalize=no_rat
  ;calculate salinity of brine:
  ;sb=brine_salinity(t)
  v=brine_volume_from_salinity(t, s, sb)

  ;calculate permittivities of brine and pure ice:
  if n_elements(eps_brn) eq 0 then eps_brn=eps_br(nu, t, sb)
  if n_elements(eps_pi) eq 0 then eps_pi=calc_eps_pi(nu, t)

  ;v=s/sb -- this is extremely approximate!!!!
  n=n_elements(v)
  ;print, v

  eps_mix_shokr, eps_pi, eps_brn, v, type, eps_x, eps_z, no_rat=no_rat

end

;empirical model for sea ice permittivity based on 
;Vant, Ramseier and Makios (1978) J. Appl. Phys. 49 (3): 1264
function eps_ice_vant, nu, t, s

  farr1=[0.1,    0.2,   0.4,  0.8,   1.,    2.,    4.]
  farr2=[0.1,    0.2,   0.4,  0.8,   1.,    2.,    4.,  7.5]
  b1arr=[3.22,   3.23,  3.26, 3.12,  3.12,  3.07,  3.05]
  a1arr=[20.6,  14.5,  12.3,  9.9,   9.0,   7.6,   7.2]
  b2arr=[0.161,  0.043, 0.043,0.048, 0.039, 0.034, 0.024, 0.032]
  a2arr=[13.24,  8.95,  7.15, 5.34,  5.04,  3.56,  3.29,  3.53]

  b1=interpol(b1arr, farr1, nu)
  a1=interpol(a1arr, farr1, nu)
  b2=interpol(b2arr, farr2, nu)
  a2=interpol(a2arr, farr2, nu)

  ;v=s*(-49.185/(t-273.2)+0.532)/1000.
  v=brine_volume_from_salinity(t, s)

  ;print, a2

  return, complex(a1*v+b1, a2*v+b2)

end

