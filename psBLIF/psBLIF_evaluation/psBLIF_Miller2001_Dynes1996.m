%% Miller2001_forward_masking_psBLIF_vs_aLIFP.m
% Forward masking (Miller 2001 Fig 7)
% Comparison: aLIFP vs psBLIF

clear; clc;

%% ------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------

measured_data = readtable("Miller_2001_fig7.csv");
load("Dynes1996_refraction.mat");

dynes_threshold_change = 10.^(data_Dynes1996.thrDynes / 20);
dynes_ipi = data_Dynes1996.ipiDynes;

%% ------------------------------------------------------------------------
% Model parameters
% ------------------------------------------------------------------------

phase_length_us = 100;
phase_length_s  = phase_length_us * 1e-6;

masker_probe_intervals = ...
    [300:100:2000, 2200:200:4000, 4500, 5000:1000:12000]';

C = psBLIF_default_parameters();

%% ------------------------------------------------------------------------
% Pulse definition
% ------------------------------------------------------------------------

% aLIFP pulse (monophasic)
pulse_aLIFP = [zeros(100,1); -1*ones(phase_length_us,1); zeros(100,1)];

% psBLIF pulse
pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = phase_length_s;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0;
pulse_psBLIF.negative_duration  = 0;
pulse_psBLIF.negative_amplitude = 0;

%% ------------------------------------------------------------------------
% Compute single pulse threshold (reference)
% ------------------------------------------------------------------------

% ---- psBLIF ----
eval_prob_psBLIF = @(amp) psBLIF_final_pulse_prob( ...
                        @(a) psBLIF_wrapper_scale_all(a, ...
                            pulse_psBLIF, C), amp);

psBLIF_I50   = psBLIF_find_threshold(eval_prob_psBLIF, 0.5);
psBLIF_lower = psBLIF_find_threshold(eval_prob_psBLIF, 0.16);
psBLIF_upper = psBLIF_find_threshold(eval_prob_psBLIF, 0.84);
psBLIF_rs    = (psBLIF_upper - psBLIF_lower) / (2 * psBLIF_I50);

% ---- aLIFP ----
eval_prob_aLIFP = @(amp) aLIFP(pulse_aLIFP * amp).total_probability;
aLIFP_I50   = psBLIF_find_threshold(eval_prob_aLIFP, 0.5);
aLIFP_lower = psBLIF_find_threshold(eval_prob_aLIFP, 0.16);
aLIFP_upper = psBLIF_find_threshold(eval_prob_aLIFP, 0.84);
aLIFP_rs    = (aLIFP_upper - aLIFP_lower) / (2 * aLIFP_I50);

%% ------------------------------------------------------------------------
% Preallocation
% ------------------------------------------------------------------------

models = ["aLIFP","psBLIF"];


psBLIF_I50_ratio = nan(length(masker_probe_intervals),1);
psBLIF_rs_ratio  = nan(length(masker_probe_intervals),1);
aLIFP_I50_ratio  = nan(length(masker_probe_intervals),1);
aLIFP_rs_ratio   = nan(length(masker_probe_intervals),1);


%% ========================================================================
% MAIN LOOP
% ========================================================================

for ind = 1:length(masker_probe_intervals)

    ipi = masker_probe_intervals(ind);

    % ---- psBLIF ----
    signal(1) = pulse_psBLIF;
    signal(1).positive_amplitude = ...
        signal(1).positive_amplitude*3*psBLIF_I50;

    signal(2) = pulse_psBLIF;
    signal(2).pulse_onset = ipi*1e-6;

    eval_prob_psBLIF = @(amp) psBLIF_final_pulse_prob( ...
                            @(a) psBLIF_wrapper_scale_probe(a, ...
                                signal, C), amp);

    lower = psBLIF_find_threshold(eval_prob_psBLIF, 0.16);
    I50   = psBLIF_find_threshold(eval_prob_psBLIF, 0.5);
    upper = psBLIF_find_threshold(eval_prob_psBLIF, 0.84);
    fprintf('%d %d\n', lower, upper)
    rs    = (upper - lower) / (2 * I50);

    psBLIF_I50_ratio(ind) = I50 / psBLIF_I50;
    psBLIF_rs_ratio(ind) = rs / psBLIF_rs;

    % ---- aLIFP ----

    eval_prob_aLIFP = @(amp) aLIFP_final_pulse_prob( ...
                             @(a) aLIFP(build_masked_stimulus_aLIFP( ...
                                a, aLIFP_I50, ipi, phase_length_us)), amp);

    lower = psBLIF_find_threshold(eval_prob_aLIFP, 0.16);
    I50   = psBLIF_find_threshold(eval_prob_aLIFP, 0.5);
    upper = psBLIF_find_threshold(eval_prob_aLIFP, 0.84);
    rs    = (upper - lower) / (2 * I50);

    aLIFP_I50_ratio(ind) = I50 / aLIFP_I50;
    aLIFP_rs_ratio(ind) = rs / aLIFP_rs;
end

%% ------------------------------------------------------------------------
% Convert units
% ------------------------------------------------------------------------

measured_data.ipi      = measured_data.ipi / 1e3;
masker_probe_intervals = masker_probe_intervals / 1e3;
dynes_ipi              = dynes_ipi / 1e3;

%% ========================================================================
% Plot
% ========================================================================

fig = figure("Name","Miller2001_psBLIF_vs_aLIFP", ...
             "DefaultAxesFontSize",13);
fig.InnerPosition(3:4) = [950 420];

tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

%% ---- Threshold shift ----
nexttile
hold on

plot(measured_data.ipi, measured_data.threshold_change, ...
     'x','Color','#000000','MarkerSize',10,'LineWidth',2);

plot(dynes_ipi, dynes_threshold_change, ...
     'o','Color','#7E2F8E','MarkerSize',9,'LineWidth',2);

plot(masker_probe_intervals, psBLIF_I50_ratio, ...
     '-','LineWidth',2,'Color','#0072BD');

plot(masker_probe_intervals, aLIFP_I50_ratio, ...
     '--','LineWidth',2,'Color','#D95319');

legend("Miller et al. (2001)", "Dynes (1996)", "psBLIF", "aLIFP");

ylabel("I_{50} / Unmasked I_{50}")
xlabel("Masker–probe Interval [ms]")
title("Threshold Shift")
grid on
ylim([0.8 3])
xlim([0 12])

%% ---- Relative spread ----
nexttile
hold on

plot(measured_data.ipi, measured_data.relative_spread_change, ...
     'x','Color','#000000','MarkerSize',10,'LineWidth',2);

plot(masker_probe_intervals, psBLIF_rs_ratio, ...
     '-','LineWidth',2,'Color','#0072BD');

plot(masker_probe_intervals, aLIFP_rs_ratio, ...
     '--','LineWidth',2,'Color','#D95319');

legend("Miller et al. (2001)", "psBLIF", "aLIFP");

ylabel("Relative Spread / Unmasked")
xlabel("Masker–probe Interval [ms]")
title("Relative Spread")
grid on
ylim([0.6 5])
xlim([0 12])

%% ========================================================================
% Helper Functions
% ========================================================================


function signal = build_masked_stimulus_aLIFP(amp, unmasked_I50, ...
                                        ipi, phase_us)

    signal = [  zeros(100,1); ...
                -1*ones(phase_us,1) * 3*unmasked_I50; ...
                   zeros(ipi - phase_us,1); ...
                -1*ones(phase_us,1) * amp; ...
                   zeros(100,1)];
end



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


function [history, alifp_ret] = ...
        psBLIF_wrapper_scale_probe(amplitude, pulse_train, C)

    scaled = pulse_train;

    % Only scale last pulse (probe)
    last = numel(scaled);

    scaled(last).positive_amplitude = ...
        scaled(last).positive_amplitude * amplitude;

    scaled(last).negative_amplitude = ...
        scaled(last).negative_amplitude * amplitude;

    [history, alifp_ret] = psBLIF(scaled, C);
end


function p = psBLIF_final_pulse_prob(fun, amp)
    [~, alifp_ret] = fun(amp);
    p = sum([alifp_ret{end}.path_prob]);
end

function p = aLIFP_final_pulse_prob(fun, amp)
    out = fun(amp);
    p = out(end).total_probability;
end
