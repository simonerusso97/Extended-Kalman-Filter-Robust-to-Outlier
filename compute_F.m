function F = compute_F(Tc, v, phi, beta)
    F = [1 0 -Tc * v * sin(phi + beta);
         0 1  Tc * v * cos(phi + beta);
         0 0         1                 ];

end