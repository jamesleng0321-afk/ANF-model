function [history, alifp_ret] = psBLIF(pulse_train, C)
% psBLIF  Path-based sequential BLIF model
%
% Inputs:
%   pulse_train : array of pulse structs
%   C           : parameter struct (extends pBLIF parameters)
%
% Outputs:
%   history     : cell array of path histories
%   alifp_ret   : spike paths per pulse

    % Global path index counter
    persistent IDX_CNT
    IDX_CNT = 0;

    history = {init_history()};
    alifp_ret = {};

    for k = 1:numel(pulse_train)
        pulse = pulse_train(k);

        this_pulse_history = [];
        this_pulse_alifp   = [];

        prev_hist = history{end};

        for i = 1:numel(prev_hist)
            path = prev_hist(i);

            % Cathodic phase times
            t_cath = get_positive_phase(pulse, C.fs);

            % Threshold modulation
            thr_mod = ones(size(t_cath));

            if C.use_facil
                thr_mod = thr_mod .* get_facil_at(t_cath, path, C);
            end
            if C.use_alifp_facil
                thr_mod = thr_mod .* aLIFP_facil(t_cath, path, C);
            end
            if C.use_refrac
                thr_mod = thr_mod .* get_refrac_at(t_cath, path, C);
            end
            if C.use_adap
                thr_mod = thr_mod .* get_adaptation_at(t_cath, path, C);
            end

            % Modified thresholds
            mu_new = C.mu * thr_mod;
            sigma_new = C.sigma * thr_mod;

            % Voltage at pulse start
            v0 = calc_voltage_at_start(pulse, path, C.tau);

            % Core pBLIF call
            res = pBLIF(pulse, mu_new, sigma_new, v0, C);

            % Append spike / no-spike paths
            [A, B, IDX_CNT] = append_to_hist(path, res, pulse, IDX_CNT, C);

            this_pulse_history = [this_pulse_history, A, B];
            this_pulse_alifp   = [this_pulse_alifp, A];
        end

        history{end+1} = cleanup_history(this_pulse_history, C.max_path);
        alifp_ret{end+1} = this_pulse_alifp;
    end

    history = history(2:end);
end

function ret = pBLIF(pulse_params, mu_new, sigma_new, v0, C)

    % Define numeric precision limit
    EPS = 1e-12;

    % Make waveform
    if C.precompute == true
        % Warning! this pulse does not have the correct amplitude set,
        % but it should not matter as we use
        % pulse_params.positive_amplitude
        pulse = C.precomputed_pulse;
    else
        pulse = psBLIF_synthesize_pulse(pulse_params, C.fs);
    end

    % Time vector
    dt = 1 / C.fs;
    t  = (0:numel(pulse)-1) * dt;

    % Cathodic phase
    cath_mask = pulse > 0;
    if sum(cath_mask) == 0 || pulse_params.positive_amplitude <= 0
        % If there is no pulse, there is nothing to do
        ret.prob = 0;
        ret.lat  = NaN;
        ret.jit  = NaN;
        return
    end

    % Membrane voltage
    v = at_end(v0, pulse_params.positive_amplitude, ...
               t(cath_mask) + dt, C.tau);

    % Spike probability
    prob = v2probvar(v, mu_new, sigma_new);
    prob = cummax(prob);
    v = prob2v(prob, C);

    % Eq. 13: find index of negative charge delivered
    if C.precompute == true
        cc_i = C.precomputed_canc_idx;
    else
        cc_i = psBLIF_find_canc_idx(pulse);
    end

    % Jitter for Eq. 16
    jit = vol2jit(v, C);
    jit_idx = min(cc_i, numel(v));
    jit_c = jit(jit_idx);

    % Eq. 16: always-cancelled spikes
    diff_t = (cc_i-1) * dt - t(cath_mask) - C.varphi;
    canceled_mask = diff_t < EPS;

    % Eq. 16: survival probability
    survival_prob = zeros(size(v));
    idx = ~canceled_mask;
    survival_prob(idx) = 1 - exp(-diff_t(idx) ./ jit_c(idx));

    % Eq. 14: combine probabilities
    diff_prob = diff([0, prob]);
    p_new = sum(diff_prob .* survival_prob);

    % Back to voltage
    v_new = prob2v(p_new, C);

    % Output
    ret.prob = p_new;
    ret.lat  = vol2lat(v_new, C);
    ret.jit  = vol2jit(v_new, C);
end


function p = v2probvar(u, mu, sigma)
%CDF with element-wise mu and sigma
    p = zeros(size(u));

    inf_idx = isinf(mu) | isinf(sigma);

    idx = ~inf_idx;
    p(idx) = 0.5 .* (1 + erf((u(idx) - mu(idx)) ./ ...
                     (sigma(idx) .* sqrt(2))));
end


function v = prob2v(p, C)
% PPF
    v = C.mu + C.sigma * sqrt(2) * erfinv(2*p - 1);
end


function jit = vol2jit(vol, C)
    jit = C.jit_a3 ./ (1 + exp((vol - C.jit_a1) ./ C.jit_a2));
end


function lat = vol2lat(vol, C)
    lat = C.lat_a3 ./ (1 + exp((vol - C.lat_a1) ./ C.lat_a2)) ...
        + C.lat_a4;
end


function positive_phase = get_positive_phase(pulse, fs)
    dt = 1/fs;
    pos_n = round(pulse.positive_duration * fs);
    positive_phase = (0:pos_n-1) *dt + pulse.pulse_onset;
end


%% ========================================================================
% Threshold functions
% ========================================================================

% Facilitation
% -------------------------------------------------------------------------
function facil = get_facil_at(t, path, C)

    if isempty(path.facil)
        facil = ones(size(t));
        return
    end

    dt = t - path.facil;
    facil = C.coeffs(1)*dt.^3 + C.coeffs(2)*dt.^2 + ...
            C.coeffs(3)*dt + C.coeffs(4);

    facil = min(facil, 1);
end

% Refractoriness
% -------------------------------------------------------------------------
function refrac = get_refrac_at(t, path, C)

    if isempty(path.last_spike)
        refrac = ones(size(t));
        return
    end

    t_since = t - path.last_spike;
    refrac  = zeros(size(t));

    inf_mask = t_since <= C.arp;
    refrac(inf_mask) = inf;

    idx = ~inf_mask;
    refrac(idx) = 1 ./ ...
        ((1 - exp((-t_since(idx) + C.arp) / (C.q * C.rrp))) .* ...
         (1 - C.r * exp((-t_since(idx) + C.arp) / C.rrp)));
end


% Long term adaptation
% -------------------------------------------------------------------------

function thr = get_adaptation_at(t, path, C)
    thr = 1 + (path.adap_c - 1) .* exp(-(t - path.adap_t) / C.t_a);
end

function new_c = get_adap_new_c(onset, prev, C)
    val = get_adaptation_at(onset, prev, C);
    new_c = min(val + C.c_inc, C.m_a);
end

%% ========================================================================
% Voltage helpers
% ========================================================================

% leaky integrator
% -------------------------------------------------------------------------
function [v_end, t_end] = calc_pulse_end_vol(pulse, path, tau)

    tpon  = pulse.pulse_onset;
    vpon  = at_end(path.pulse_end_vol, 0, tpon - path.pulse_end_t, tau);

    tcoff = tpon + pulse.positive_duration;
    vcoff = at_end(vpon, pulse.positive_amplitude, pulse.positive_duration, tau);

    tioff = tcoff + pulse.interphase_gap;
    vioff = at_end(vcoff, 0, pulse.interphase_gap, tau);

    teop  = tioff + pulse.negative_duration;
    veop  = at_end(vioff, pulse.negative_amplitude, pulse.negative_duration, tau);

    v_end = veop;
    t_end = teop;
end


function v0 = calc_voltage_at_start(pulse, path, tau)

    last_offset = path.pulse_end_t;
    onset = pulse.pulse_onset;
    v0 = at_end(path.pulse_end_vol, 0, onset-last_offset, tau);
end


function v = at_end(old, new, t, tau)
    if old == new
        v = old;
    else
        v = new + (old - new) * exp(-t / tau);
    end
end


%% ========================================================================
% Bookkeeping
% ========================================================================

function hist = init_history()
    hist.spike          = false;
    hist.prev_idx       = NaN;
    hist.prob           = 1;
    hist.lat            = NaN;
    hist.jit            = NaN;
    hist.path_idx       = 0;
    hist.path_prob      = 1;
    hist.pulse_end_vol  = 0;
    hist.pulse_end_t    = 0;
    hist.facil          = [];
    hist.last_spike     = [];
    hist.adap_c         = 1;
    hist.adap_t         = 0;
end


function [spiked, nospike, IDX_CNT] = ...
        append_to_hist(prev, res, pulse, IDX_CNT, C)

    % ---- Spike path ----
    IDX_CNT = IDX_CNT + 1;
    spiked.spike     = true;
    spiked.prev_idx  = prev.path_idx;
    spiked.prob      = res.prob;
    spiked.lat       = res.lat;
    spiked.jit       = res.jit;
    spiked.path_idx  = IDX_CNT;
    spiked.path_prob = res.prob * prev.path_prob;

    spiked.pulse_end_vol = 0;
    spiked.pulse_end_t   = 0;
    spiked.facil         = [];
    spiked.last_spike    = pulse.pulse_onset;
    spiked.adap_c        = get_adap_new_c(pulse.pulse_onset, prev, C);
    spiked.adap_t        = pulse.pulse_onset;

    % ---- No-spike path ----
    IDX_CNT = IDX_CNT + 1;
    nospike.spike     = false;
    nospike.prev_idx  = prev.path_idx;
    nospike.prob      = 1 - res.prob;
    nospike.lat       = NaN;
    nospike.jit       = NaN;
    nospike.path_idx  = IDX_CNT;
    nospike.path_prob = (1 - res.prob) * prev.path_prob;

    if C.dont_carry_voltage
        nospike.pulse_end_vol = 0;
        nospike.pulse_end_t = 0;
    else
        [nospike.pulse_end_vol, nospike.pulse_end_t] = ...
            calc_pulse_end_vol(pulse, prev, C.tau);
    end

    nospike.facil      = pulse.pulse_onset;
    nospike.last_spike = prev.last_spike;
    nospike.adap_c     = prev.adap_c;
    nospike.adap_t     = prev.adap_t;
end


function new_hist = cleanup_history(hist, max_path)

    if numel(hist) <= max_path
        new_hist = hist;
        return
    end

    probs = [hist.path_prob];
    [~, idx] = maxk(probs, max_path);

    prob_sum = sum(probs(idx));

    new_hist = hist(idx);
    for k = 1:numel(new_hist)
        new_hist(k).path_prob_raw = new_hist(k).path_prob;
        new_hist(k).path_prob     = new_hist(k).path_prob_raw / prob_sum;
    end
end
