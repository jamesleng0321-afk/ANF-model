function C = psBLIF_default_parameters()
    % pBLIF
    C.fs     = 1e6;
    C.tau    = 250e-6;
    C.varphi = 35e-6;

    C.mu     = 105e-6;
    C.sigma  = 4.6e-6;

    C.jit_a1 = 109e-6;
    C.jit_a2 = 3.24e-6;
    C.jit_a3 = 136e-6;

    C.lat_a1 = 106e-6;
    C.lat_a2 = 5.14e-6;
    C.lat_a3 = 368e-6;
    C.lat_a4 = 472e-6;

    % refractory
    C.use_refrac = true;
    C.arp = 0.3e-3;
    C.rrp = 1.5e-3;
    C.q = 0.76;
    C.r = 8.77e-3;

    % adaptation
    C.use_adap = true;
    C.t_a = 60e-3;
    C.m_a = 1.06;
    C.c_inc = 0.01;

    % facilitation
    C.use_facil = true;
    C.coeffs = [1.3e9, -2.42e6, 1.68e3, 0.51];

    % bookeeping
    C.max_path = 40;

    % else
    C.use_alifp_facil = 0;
    C.dont_carry_voltage = 0;

    C.precompute = false;
end
