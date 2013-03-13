function calc_eps_pi, nu, t1
  t=t1-abs_zero()

  beta=(0.445+0.00211*t)*1e-4+0.585e-4/(1-t/29.1)^2

  theta=300./t1-1.
  alpha=(50.4+62.*theta)*1e-4*exp(-22.1*theta)

  eps_im=alpha/nu+beta*nu

  eps_re=3.1884+9.1e-4*t

  return, complex(eps_re, eps_im)

end

