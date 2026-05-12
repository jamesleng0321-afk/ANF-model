function p_vec = psBLIF_per_pulse_prob(alifp_ret)
    p_vec = nan(size(alifp_ret));
    for k = 1:numel(alifp_ret)
        p_vec(k) = sum([alifp_ret{k}.path_prob]);
    end
end
