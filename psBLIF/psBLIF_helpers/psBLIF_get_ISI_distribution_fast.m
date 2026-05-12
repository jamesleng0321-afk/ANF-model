function [p_vec, t] = psBLIF_get_ISI_distribution_fast(isi_dist, fs)
    max_sigma = 5;
    t_max = max([isi_dist.mu] + max_sigma*[isi_dist.sigma]);
    t = 0:1/fs:t_max;

    p_vec = zeros(size(t));
    for k = 1:numel(isi_dist)
        mu = isi_dist(k).mu;
        sigma = isi_dist(k).sigma;
        sigma = max(sigma, 1e-6); % enforce minimum
        p = isi_dist(k).p;

        t_lower = mu-max_sigma*sigma;
        t_upper = mu+max_sigma*sigma;

        idx_lower = max(round(t_lower*fs)+1, 1);
        idx_upper = min(round(t_upper*fs)+1, length(t));

        p_temp = normpdf(t(idx_lower:idx_upper), mu, sigma);

        p_temp = p_temp*p/sum(p_temp);

        p_vec(idx_lower:idx_upper) = p_vec(idx_lower:idx_upper) + ...
            p_temp;
    end
end
