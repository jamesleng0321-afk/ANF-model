function [p_vec, t] = psBLIF_get_spike_distribution_fast(alifp_ret, pulse_train, max_time, fs)
    t = 0:1/fs:max_time;
    max_sigma = 5;

    p_vec = zeros(size(t));
    for k = 1:numel(alifp_ret)
        for m = 1:numel(alifp_ret{k})
            mu = pulse_train(k).pulse_onset+alifp_ret{k}(m).lat;
            sigma = alifp_ret{k}(m).jit;
            sigma = max(sigma, 1e-6); % enforce minimum
            path_prob = alifp_ret{k}(m).path_prob;

            t_lower = mu-max_sigma*sigma;
            t_upper = mu+max_sigma*sigma;

            idx_lower = max(round(t_lower*fs)+1, 1);
            idx_upper = min(round(t_upper*fs)+1, length(t));

            p_temp = normpdf(t(idx_lower:idx_upper), mu, sigma);

            p_temp = p_temp*path_prob/sum(p_temp);

            p_vec(idx_lower:idx_upper) = p_vec(idx_lower:idx_upper) + ...
                p_temp;
        end
    end
end
