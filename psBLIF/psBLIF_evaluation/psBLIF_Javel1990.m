%% Javel1990_psBLIF_vs_aLIFP.m
% Reproduction of Javel 1990
% Mean firing probability vs stimulus level
% Comparison: psBLIF vs aLIFP

clear; clc;

%% ------------------------------------------------------------------------
% Load experimental data
% ------------------------------------------------------------------------

load('Javel1990.mat')

colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];

%% ------------------------------------------------------------------------
% Stimulation parameters
% ------------------------------------------------------------------------

pulse_rates = [100, 200, 400, 800];   % pps
duration_s  = 0.1;                    % 100 ms
amplitudes_a   = (1:0.1:3.2) * 1e-3;
amplitudes_ps  = (1:0.1:3.2) * 1e-3 *0.6;

%% ------------------------------------------------------------------------
% Pulse definition (100 µs biphasic)
% ------------------------------------------------------------------------

% ---- aLIFP ----
pulse_aLIFP = [-1*ones(50,1); ones(50,1)];

% ---- psBLIF ----
phase_len_s = 50e-6;

pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = phase_len_s;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0;
pulse_psBLIF.negative_duration  = phase_len_s;
pulse_psBLIF.negative_amplitude = -1;

C = psBLIF_default_parameters();


%% ------------------------------------------------------------------------
% Preallocate
% ------------------------------------------------------------------------

mean_prob_a  = zeros(length(pulse_rates), length(amplitudes_ps));
mean_prob_ps = zeros(length(pulse_rates), length(amplitudes_ps));

%% ========================================================================
% MAIN LOOP
% ========================================================================

for p_ind = 1:length(pulse_rates)

    rate = pulse_rates(p_ind);

    % ---- Build stimuli ----
    stimulus_a = aLIFP_get_pulse_train(duration_s*1e6, ...
                                       rate, pulse_aLIFP);

    pulse_train_ps = build_psBLIF_train(duration_s, ...
                                        rate, pulse_psBLIF);

    % ---- Amplitude loop ----
    for a_ind = 1:length(amplitudes_ps)

        % ---- aLIFP ----
        out_a = aLIFP(stimulus_a * amplitudes_a(a_ind));
        mean_prob_a(p_ind, a_ind) = ...
            mean([out_a.total_probability]);

        % ---- psBLIF ----
        [~, alifp_ret] = ...
            psBLIF_wrapper_scale_all(amplitudes_ps(a_ind), ...
                                     pulse_train_ps, C);

        per_pulse_prob = psBLIF_per_pulse_prob(alifp_ret);

        mean_prob_ps(p_ind, a_ind) = ...
            mean(per_pulse_prob);
    end
end

%% ------------------------------------------------------------------------
% Reference threshold (100 pps)
% ------------------------------------------------------------------------

[I50_a,  ~] = aLIFP_fit_gaussian(amplitudes_a, mean_prob_a(1,:));
[I50_ps, ~] = aLIFP_fit_gaussian(amplitudes_ps, mean_prob_ps(1,:));

%% ========================================================================
% Plot
% ========================================================================

fig = figure("Name","Javel1990_psBLIF_vs_aLIFP", ...
             "DefaultAxesFontSize",13);
fig.InnerPosition(3:4) = [700 550];
hold on

for p_ind = 1:length(pulse_rates)

    rate = pulse_rates(p_ind);

    [data_amp, data_prob] = ...
        extract_javel_data(data_Javel1990, rate);

    data_amp = 10.^(data_amp/20);

    if rate == 100
        [data_threshold,~] = ...
            aLIFP_fit_gaussian(data_amp, data_prob);
    end

    % ---- Plot experimental ----
    plot(20*log10(data_amp/data_threshold), ...
         data_prob, ...
         'x', ...
         'Color',colors(p_ind), ...
         'MarkerSize',10, ...
         'LineWidth',2);

    % ---- Plot aLIFP ----
    plot(20*log10(amplitudes_a/I50_a), ...
         mean_prob_a(p_ind,:), ...
         '--', ...
         'Color',colors(p_ind), ...
         'LineWidth',1.8);

    % ---- Plot psBLIF ----
    plot(20*log10(amplitudes_ps/I50_ps), ...
         mean_prob_ps(p_ind,:), ...
         '-', ...
         'Color',colors(p_ind), ...
         'LineWidth',1.8);
end

ylabel("Mean firing probability")
xlabel("Stimulus level [dB re 100 pps]")
xlim([-2 8])
grid on

legend( ...
    "data 100 pps", "aLIFP 100 pps", "psBLIF 100 pps", ...
    "data 200 pps", "aLIFP 200 pps", "psBLIF 200 pps", ...
    "data 400 pps", "aLIFP 400 pps", "psBLIF 400 pps", ...
    "data 800 pps", "aLIFP 800 pps", "psBLIF 800 pps", ...
    "Location","northoutside", ...
    "Orientation","horizontal", ...
    "NumColumns",3);

%% ========================================================================
% Helper Functions
% ========================================================================

function [history, alifp_ret] = ...
        psBLIF_wrapper_scale_all(amplitude, pulse_train, C)

    scaled = pulse_train;

    for k = 1:numel(scaled)
        scaled(k).positive_amplitude = ...
            scaled(k).positive_amplitude * amplitude;

        scaled(k).negative_amplitude = ...
            scaled(k).negative_amplitude * amplitude;
    end

    [history, alifp_ret] = psBLIF(scaled, C);
end

function pulse_train = build_psBLIF_train(duration_s, rate, pulse)

    ipi = 1/rate;
    n   = 1 + floor(duration_s/ipi);

    for k = 1:n
        pulse_train(k) = pulse;
        pulse_train(k).pulse_onset = (k-1)*ipi;
    end
end

function [amp, prob] = extract_javel_data(data_struct, rate)

    switch rate
        case 100
            amp  = data_struct.hundred(:,1);
            prob = data_struct.hundred(:,2);

        case 200
            amp  = data_struct.twoH(:,1);
            prob = data_struct.twoH(:,2);

        case 400
            amp  = data_struct.fourH(:,1);
            prob = data_struct.fourH(:,2);

        case 800
            amp  = data_struct.eightH(:,1);
            prob = data_struct.eightH(:,2);

        otherwise
            error("Unsupported pulse rate.")
    end
end
