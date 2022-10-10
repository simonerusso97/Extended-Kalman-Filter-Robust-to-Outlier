function z = predict_z(s, k, Tc, N_gps, N_compass)
    if mod( (k-1) , (Tc^-1)/N_gps ) == 0
        z(1) = s(1);
        z(2) = s(2);
    else
        z(1) = NaN;
        z(2) = NaN;
    end

    if mod( (k-1) , (Tc^-1)/N_compass ) == 0
        z(3) = s(3);
    else
        z(3) = NaN;
end