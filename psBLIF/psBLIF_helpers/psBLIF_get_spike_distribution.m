function [p_vec, t] = psBLIF_get_spike_distribution(alifp_ret, pulse_train, max_time, fs)
    t = 0:1/fs:max_time;
    p_vec = zeros(size(t));
    for k = 1:numel(alifp_ret)
        for m = 1:numel(alifp_ret{k})
            mu = pulse_train(k).pulse_onset+alifp_ret{k}(m).lat;
            sigma = alifp_ret{k}(m).jit;
            sigma = max(sigma, 1e-6); % enforce minimum
            p_vec = p_vec + normpdf(t, mu, sigma) ...
                             * alifp_ret{k}(m).path_prob;
        end
    end
    p_vec = p_vec/fs;
end
