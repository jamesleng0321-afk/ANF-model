function facil = aLIFP_facil(t, path, C)

    if isempty(path.facil)
        facil = ones(size(t));
        return
    end

    t_vec = t - path.facil;

    ff1 =  0.1e-3;
    ff2 = -1.4e-3;
    ff3 =  0.45;
    ff4 =  0.9e3;

    f1 = 1 - exp(- ff4 * (t_vec + ff1));
    f2 = 1 + ff3 * exp(-ff4 * (t_vec + ff2));
    facil = f1 .* f2;
end
