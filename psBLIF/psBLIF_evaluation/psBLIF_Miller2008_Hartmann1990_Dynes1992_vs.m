%% VectorStrength_psBLIF_vs_aLIFP.m
% Reproduction of temporal coding (vector strength)
% Comparison: aLIFP vs psBLIF

clear; clc;

%% ------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------

load("vector_strength_data.mat")

ipgs = [0, 30];  % µs
pulse_rates = [50, 100, 200, 400, 800, 1000, 1250, 1600, 2500, 5000];

colors = ["#0072BD", "#D95319"];

C = psBLIF_default_parameters();

%% ------------------------------------------------------------------------
% Preallocate
% ------------------------------------------------------------------------

vs_a  = nan(length(ipgs), length(pulse_rates));
vs_ps = nan(length(ipgs), length(pulse_rates));

%% ========================================================================
% MAIN LOOP
% ========================================================================

aLIFP = false;

for i_ind = 1:length(ipgs)

    ipg = ipgs(i_ind);

    %% ---- Pulse definitions ----

    % aLIFP pulse (sample-based)
    pulse_a = [-1*ones(40,1); zeros(ipg,1); ones(40,1)];

    % psBLIF pulse (time-based)
    phase_len = 40e-6;
    ipg_s     = ipg * 1e-6;

    pulse_ps.pulse_onset        = 0;
    pulse_ps.positive_duration  = phase_len;
    pulse_ps.positive_amplitude = 1;
    pulse_ps.interphase_gap     = ipg_s;
    pulse_ps.negative_duration  = phase_len;
    pulse_ps.negative_amplitude = -1;

    %% ---- Threshold (90% firing) ----

    % aLIFP
    if aLIFP
        amp_a = aLIFP_get_threshold(pulse_a, (0.4:0.05:3)*1e-3, 0.9);
    end

    % psBLIF
    eval_prob_ps = @(amp) psBLIF_single_pulse_prob(pulse_ps, amp, C);
    amp_ps = psBLIF_find_threshold(eval_prob_ps, 0.9);

    %% ---- Pulse-rate loop ----
    for p_ind = 1:length(pulse_rates)

        rate = pulse_rates(p_ind);
        fprintf("IPG %d µs | Rate %d pps\n", ipg, rate);

        %% ---- aLIFP ----
        if aLIFP
            stim_a = aLIFP_get_pulse_train(0.3e6, rate, pulse_a);
            out_a  = aLIFP(stim_a * amp_a);

            dist_a = aLIFP_get_spike_distribution(out_a, 0.3);
            vs_a(i_ind, p_ind) = ...
                aLIFP_calculate_vector_strength(dist_a, rate);
        end

        %% ---- psBLIF ----
        stim_ps = build_psBLIF_train(0.3, rate, pulse_ps);

        [~, out_ps] = psBLIF_wrapper_scale_all(amp_ps, stim_ps, C);

        t_max = stim_ps(end).pulse_onset + 1e-2;
        dist_ps = psBLIF_get_spike_distribution_fast(out_ps, ...
            stim_ps, t_max, 1e6);
        vs_ps(i_ind, p_ind) = ...
            aLIFP_calculate_vector_strength(dist_ps, rate);
    end
end

%% ========================================================================
% Plot
% ========================================================================

fig = figure("Name","VectorStrength_psBLIF_vs_aLIFP");


legend_entries = strings(0);

for i_ind = 1:length(ipgs)

    if aLIFP
    % aLIFP
        semilogx(pulse_rates, vs_a(i_ind,:), ...
            '--o', ...
            "Color", colors(i_ind), ...
            "LineWidth", 2, ...
            "MarkerSize", 7);

        hold on;
        legend_entries(end+1) = sprintf("aLIFP %d µs", ipgs(i_ind));
    end

    % psBLIF
    semilogx(pulse_rates, vs_ps(i_ind,:), ...
        '-o', ...
        "Color", colors(i_ind), ...
        "LineWidth", 2, ...
        "MarkerSize", 7);
    hold on;
    legend_entries(end+1) = sprintf("psBLIF %d µs", ipgs(i_ind));
end

% ---- Experimental data ----

% Miller 2008
errorbar([250,1000,5000], ...
    temporal_coding_data.Miller.avg, ...
    temporal_coding_data.Miller.std, ...
    'xk','MarkerSize',8,'LineWidth',1.5);

% Hartmann & Klinke
p = temporal_coding_data.HartmannKlinke.cis(:,1) - ...
    temporal_coding_data.HartmannKlinke.averages;
n = temporal_coding_data.HartmannKlinke.averages - ...
    temporal_coding_data.HartmannKlinke.cis(:,2);

errorbar(temporal_coding_data.HartmannKlinke.frqs, ...
         temporal_coding_data.HartmannKlinke.averages, ...
         p, n, 'o', 'Color','#2d7e2e','LineWidth',1);

% Dynes & Delgutte
p = temporal_coding_data.DynesDelgutte.cis(:,1) - ...
    temporal_coding_data.DynesDelgutte.averages;
n = temporal_coding_data.DynesDelgutte.averages - ...
    temporal_coding_data.DynesDelgutte.cis(:,2);

errorbar(temporal_coding_data.DynesDelgutte.frq, ...
         temporal_coding_data.DynesDelgutte.averages, ...
         p, n, 'd', 'Color','#7e2d7d','LineWidth',1);

legend_entries(end+1:end+3) = ...
    ["Miller 2008","Hartmann & Klinke (1990)","Dynes & Delgutte (1992)"];

% ---- Formatting ----

xticks([50, 100, 500, 1000, 5000])
xlabel("Pulse rate [pps]")
ylabel("Vector Strength")
legend(legend_entries, "Location","southwest")
grid on
ylim([0, 1.1])

%% ========================================================================
% Helper Functions
% ========================================================================

function train = build_psBLIF_train(duration_s, rate, pulse)

    ipi = 1/rate;
    n   = floor(duration_s / ipi);

    for k = 1:n
        train(k) = pulse;
        train(k).pulse_onset = (k-1)*ipi;
    end
end

function p = psBLIF_single_pulse_prob(pulse, amp, C)

    scaled = pulse;
    scaled.positive_amplitude = scaled.positive_amplitude * amp;
    scaled.negative_amplitude = scaled.negative_amplitude * amp;

    [~, out] = psBLIF(scaled, C);
    p = sum(out{1}.path_prob);
end

function [history, out] = psBLIF_wrapper_scale_all(amp, train, C)

    scaled = train;
    for k = 1:numel(train)
        scaled(k).positive_amplitude = train(k).positive_amplitude * amp;
        scaled(k).negative_amplitude = train(k).negative_amplitude * amp;
    end

    [history, out] = psBLIF(scaled, C);
end
