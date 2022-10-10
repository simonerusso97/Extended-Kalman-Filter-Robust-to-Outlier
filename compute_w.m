function w = compute_w(r, sigma)
    c = 2.5;
    for j = 1:length(r)
        if abs(r(j)/sigma) <= c
            w(j) = 1;
        else
            w(j) = 0;
        end
    end
end