;*** warning!!! there is a constant wrong here somewhere...
;complex permittivity of pure water:
function eps_pw, nu1, t1

  nu=nu1*1e9

  t=t1-abs_zero()

  tau=1.1109e-10 -t*(3.824e-12 + t*(6.938e-14 - 5.096e-16*t))

  eps0=88.045 - t*(0.4147 + t*(6.295e-4 + 1.075e-5*t))

  epsinf=4.9

  epsp=epsinf+(eps0-epsinf)/(1-(tau*nu)^2)

  epspp=(tau*nu*(eps0-epsinf))/(1+(tau*nu)^2)

  eps=epsinf+(eps0-epsinf)/complex(1, tau*nu)

  return, conj(eps)

end
