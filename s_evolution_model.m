function s = s_evolution_model(s_prec, u, Tc, l_r, l_f)
    x0 = s_prec(1);
    y0 = s_prec(2);
    phi0 = s_prec(3);

    v = u(1);
    delta = u(2);
    beta = atan( l_r * tan(delta) / ( l_r+l_f ) ); 

    s(1) = x0 + Tc * v * cos( phi0 + beta);
    s(2) = y0 + Tc * v * sin( phi0 + beta);
    s(3) = phi0 + Tc * v * sin(beta)/l_r;

end