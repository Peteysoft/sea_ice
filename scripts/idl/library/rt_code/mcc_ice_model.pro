@eps_parm

pro mcc_ice_model, s1, h1, z1, x0_1, y0_1, zangle1, fov1, $
		tbv, tbh, verr, herr, $
                tsky=tsky1, tice=tice1, tw=tw1, nsea=nsea, sw=sw1, $
                niter=niter1, plot=doplot, fail=fail, $
                eps=eps1, freq=freq1, offset=offset1, $
                xall=xall, yall=yall

  ;set defaults and convert to correct type:
  if n_elements(tsky1) ne 1 then tsky=5.D else tsky=double(tsky1)
  if n_elements(tice1) ne 1 then tice=272.D else tice=double(tice1)
  if n_elements(tw1) ne 1 then tw=272.D else tw=double(tw1)
  if n_elements(sw1) ne 1 then sw=5.D else sw=double(sw1)
  if n_elements(freq1) ne 1 then freq=1.4D else freq=double(freq1)
  if n_elements(fov1) ne 1 then fov=13.6D else fov=double(fov1)

  ;machine dependent parameters:
  if !version.arch eq "x86-64" then begin
    if n_elements(niter1) ne 1 then niter=long64(2000L) else niter=long64(niter1)
    nsurf=long64(n_elements(s1))
    ninit=long64(n_elements(x0_1))
  endif else begin
    if n_elements(niter1) ne 1 then niter=2000L else niter=long(niter1)
    nsurf=n_elements(s1)
    ninit=n_elements(x0_1)
  endelse

  if n_elements(offset1) eq 0 then offset1=0.D

  path="../../../bin"

  s=double(s1)
  h=double(h1)
  z=double(z1)

  if n_elements(eps1) eq 0 then eps1=complex(4., 0.1)
  if n_elements(eps1) eq 1 then eps1=replicate(eps1, nsurf)
  epsp=double(eps1)
  epspp=double(imaginary(eps1))

  if n_elements(h) ne nsurf or n_elements(z) ne nsurf or $
	  n_elements(eps1) ne nsurf then begin
    message, "Ice topography arrays must have the same number of elements"
    fail=1
    return
  endif

  x0=dblarr(ninit)
  x0[*]=double(x0_1)
  y0=dblarr(ninit)
  y0[*]=double(y0_1)
  zangle=dblarr(ninit)
  zangle[*]=double(zangle1)

  if n_elements(offset1) eq 1 then offset1=replicate(offset1, ninit)
  offset=double(offset1)

  if n_elements(y0) ne ninit or n_elements(zangle) ne ninit or $
	  n_elements(offset) ne ninit then begin
    message, "Initial conditions must have the same number of elements"
    fail=2
    return
  endif

  tbv=dblarr(ninit)
  tbh=dblarr(ninit)
  verr=dblarr(ninit)
  herr=dblarr(ninit)

  ;stop
  fail=call_external(path+"mcc_ice_idl_wrapper.so", "_Z15ice_idl_wrapperiPPv", $
		s, h, z, epsp, epspp, nsurf, x0, y0, zangle, offset, ninit, $
		niter, freq, fov, tice, tw, sw, tsky, tbv, tbh, verr, herr)

  ;stop

end

