%% Vector strength vs IPG (psBLIF version)

clear; clc;

%% ------------------------------------------------------------------------
% Parameters
% ------------------------------------------------------------------------

pulse_rate = 5000; % alternative: 2500

if pulse_rate == 5000
    ipg = [0, 10, 20, 30, 40, 50];     % µs
else
    ipg = [0, 10, 20, 30, 40, 50, 100];
end

amplitudes_dB = [1, 2];
duration_s = 0.3;

C = psBLIF_default_parameters();

phase_len_s = 40e-6;

pulse.pulse_onset        = 0;
pulse.positive_duration  = phase_len_s;
pulse.positive_amplitude = 1;

pulse.interphase_gap     = 0;

pulse.negative_duration  = phase_len_s;
pulse.negative_amplitude = -1;

%% ------------------------------------------------------------------------
% Preallocation
% ------------------------------------------------------------------------

vector_strength = nan(length(amplitudes_dB), length(ipg));

%% ========================================================================
% MAIN LOOP
% ========================================================================

for ipg_ind = 1:length(ipg)

    pulse.interphase_gap = ipg(ipg_ind) * 1e-6;

    % ---- Find I50 for this pulse ----
    eval_prob = @(amp) psBLIF_final_pulse_prob( ...
                        @(a) psBLIF_wrapper_scale_all(a, pulse, C), amp);

    I50 = psBLIF_find_threshold(eval_prob, 0.5);

    % ---- Build pulse train ----
    pulse_train = build_psBLIF_train(duration_s, pulse_rate, pulse);

    for a_ind = 1:length(amplitudes_dB)

        % ---- Scale amplitude (dB re I50) ----
        scale = 10^(amplitudes_dB(a_ind)/20);
        amp   = I50 * scale;

        % ---- Run model ----
        [~, out_ps] = psBLIF_wrapper_scale_all(amp, pulse_train, C);

        % ---- Get spike distribution ----
        t_max = pulse_train(end).pulse_onset + 2e-3;
        dist = psBLIF_get_spike_distribution_fast(out_ps, pulse_train, t_max, 1e6);

        % ---- Vector strength ----
        vector_strength(a_ind, ipg_ind) = ...
            aLIFP_calculate_vector_strength(dist, pulse_rate);

        fprintf("IPG: %d µs | Amp idx: %d\n", ipg(ipg_ind), a_ind);
    end
end

%% ========================================================================
% Plot
% ========================================================================

fig = figure("Name", "psBLIF_VectorStrength_IPG", ...
             "DefaultAxesFontSize",13);
fig.Position(3:4) = [650, 550];
hold on;

markers = ["o-", "--d", "-.x", ":s"];
labels  = strings(length(amplitudes_dB),1);

for a_ind = 1:length(amplitudes_dB)
    plot(ipg, vector_strength(a_ind,:), markers(a_ind), ...
        "LineWidth", 2);

    labels(a_ind) = sprintf("psBLIF %g dB re I_{50}", ...
                             amplitudes_dB(a_ind));
end

legend(labels, ...
    "Orientation","horizontal", ...
    "Location","northoutside");

ylabel("Vector Strength")
xlabel("IPG [µs]")
ylim([0 1])
grid on

%% ========================================================================
% Helper functions
% ========================================================================

function pulse = create_psBLIF_pulse(ipg_us)

    phase_len_s = 40e-6;

    pulse.pulse_onset        = 0;
    pulse.positive_duration  = phase_len_s;
    pulse.positive_amplitude = 1;

    pulse.interphase_gap     = ipg_us * 1e-6;

    pulse.negative_duration  = phase_len_s;
    pulse.negative_amplitude = -1;
end

function pulse_train = build_psBLIF_train(duration_s, rate, pulse)

    ipi = 1 / rate;
    n   = floor(duration_s / ipi);

    for k = 1:n
        pulse_train(k) = pulse;
        pulse_train(k).pulse_onset = (k-1) * ipi;
    end
end

function [history, out] = ...
        psBLIF_wrapper_scale_all(amplitude, pulse_train, C)

    scaled = pulse_train;

    for k = 1:numel(scaled)
        scaled(k).positive_amplitude = ...
            scaled(k).positive_amplitude * amplitude;

        scaled(k).negative_amplitude = ...
            scaled(k).negative_amplitude * amplitude;
    end

    [history, out] = psBLIF(scaled, C);
end

function p = psBLIF_final_pulse_prob(fun, amp)
    [~, out] = fun(amp);
    p = sum([out{end}.path_prob]);
end
