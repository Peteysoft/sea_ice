
function eps_water_func, f, T,S ;empirical relationship for dielectric constant of water by Klein and Swift 
    lj=complex(0, 1)

    T=T-273.15
    omega=2*!pi*f
    eps_0=8.854e-12
    eps_inf=4.9

    eps_s_T=87.134-1.949e-1*T-1.276e-2*T^2+2.491e-4*T^3
    a_ST=1.+1.613e-5*S*T-3.656e-3*S+3.210e-5*S^2-4.232e-7*S^3
    eps_s=eps_s_T*a_ST

    tau_T0=1.768e-11-6.086e-13*T+1.104e-14*T^2-8.111e-17*T^3
    b_ST=1.+2.282e-5*S*T-7.638e-4*S-7.760e-6*S^2+1.105e-8*S^3
    tau=tau_T0*b_ST

    delta=25-T
    beta=2.0333e-2+1.266e-4*delta+2.464e-6*delta^2-S*(1.849e-5-2.551e-7*delta+2.551e-8*delta^2)
    sigma_25S=S*(0.182521-1.46192e-3*S+2.09324e-5*S^2-1.28205e-7*S^3)
    sigma=sigma_25S*exp(-delta*beta)

    eps_water=eps_inf+(eps_s-eps_inf)/(1+lj*omega*tau)-lj*sigma/(omega*eps_0)

    stop
    ;eps_water=77.9-81.5*1j
    return, eps_water

end

