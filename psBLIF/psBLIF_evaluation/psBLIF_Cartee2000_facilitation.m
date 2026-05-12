%% psBLIF_Cartee2000.m
% Reproduction of Cartee 2000 facilitation experiment using psBLIF

clear; clc;

%% ------------------------------------------------------------------------
% Load experimental data
% -------------------------------------------------------------------------

load("Cartee2000.mat")

data_ipis = [100, 200, 300]; % µs
ipis      = 100:1:300;       % µs sweep

%% ------------------------------------------------------------------------
% Model parameters
% -------------------------------------------------------------------------

C = psBLIF_default_parameters();

phase_dur  = 50e-6;               % 50 µs
phaseL     = 50;                  % samples for aLIFP
amplitudes = (0.01:0.01:2) * 1e-3;

I50_ratio_psBLIF = nan(size(ipis));
I50_ratio_aLIFP  = nan(size(ipis));

%% ------------------------------------------------------------------------
% Loop over IPIs
% -------------------------------------------------------------------------

for ind = 1:numel(ipis)

    ipi_us = ipis(ind);
    ipi    = ipi_us * 1e-6;

    %% ------------------- psBLIF -------------------

    single = build_pulse_train(1, ipi, phase_dur);
    dual   = build_pulse_train(2, ipi, phase_dur);

    fun_single = @(amp) psBLIF_wrapper(amp, single, C);
    fun_dual   = @(amp) psBLIF_wrapper(amp, dual, C);

    eval_single = @(amp) evaluate_probability_psBLIF(fun_single, amp);
    eval_dual   = @(amp) evaluate_probability_psBLIF(fun_dual, amp);

    I50_single = psBLIF_get_threshold(eval_single, amplitudes, 0.5);
    I50_dual   = psBLIF_get_threshold(eval_dual, amplitudes, 0.5);

    I50_ratio_psBLIF(ind) = I50_dual / I50_single;

    %% ------------------- aLIFP -------------------

    pulse = [-1 * ones(phaseL,1); ...
            ones(ipi_us - phaseL,1) ./ ((ipi_us - phaseL)/phaseL)];

    double_pulse = [zeros(100,1); pulse; pulse; zeros(100,1)];

    eval_single_aLIFP = @(amp) evaluate_probability_aLIFP(pulse, amp);
    eval_dual_aLIFP   = @(amp) evaluate_probability_aLIFP(double_pulse, amp);

    I50_single_aLIFP = psBLIF_get_threshold(eval_single_aLIFP, amplitudes, 0.5);
    I50_dual_aLIFP   = psBLIF_get_threshold(eval_dual_aLIFP, amplitudes, 0.5);

    I50_ratio_aLIFP(ind) = I50_dual_aLIFP / I50_single_aLIFP;

end

%% ------------------------------------------------------------------------
% Extract experimental statistics
% -------------------------------------------------------------------------

means(1) = mean(facilitation_Cartee2000.raw(1:15,2));
means(2) = mean(facilitation_Cartee2000.raw(16:30,2));
means(3) = mean(facilitation_Cartee2000.raw(31:43,2));

stds(1) = std(facilitation_Cartee2000.raw(1:15,2));
stds(2) = std(facilitation_Cartee2000.raw(16:30,2));
stds(3) = std(facilitation_Cartee2000.raw(31:43,2));

%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

fig = figure("Name","Cartee2000_comparison","DefaultAxesFontSize",13);
fig.OuterPosition(3:4) = [650 600];

errorbar(data_ipis, means, stds, 'xk', ...
    'MarkerSize',10,'LineWidth',2);
hold on;

plot(ipis, I50_ratio_psBLIF, ...
    'LineWidth',2, "Color","#0072BD");

plot(ipis, I50_ratio_aLIFP, ...
    'LineWidth',2, "Color","#D95319");

legend("data","psBLIF","aLIFP");
ylabel("I_{50} / Single pulse I_{50}");
xlabel("Masker pulse Interval [µs]");
xlim([100 300]);
grid on;

%% ========================================================================
% Local helper functions
% ========================================================================

function pulse_train = build_pulse_train(n_pulses, ipi, phase_dur)

    pulse_train = struct([]);

    for k = 1:n_pulses
        pulse_train(k).pulse_onset        = (k-1) * ipi;
        pulse_train(k).positive_duration  = phase_dur;
        pulse_train(k).positive_amplitude = 1;

        pulse_train(k).interphase_gap     = 0;
        pulse_train(k).negative_duration  = ipi - phase_dur;

        % charge balanced
        pulse_train(k).negative_amplitude = ...
            -phase_dur / (ipi - phase_dur);
    end
end


function [history, alifp_ret] = psBLIF_wrapper(amplitude, pulse_train, C)

    scaled = pulse_train;

    for k = 1:numel(scaled)
        scaled(k).positive_amplitude = ...
            scaled(k).positive_amplitude * amplitude;

        scaled(k).negative_amplitude = ...
            scaled(k).negative_amplitude * amplitude;
    end

    [history, alifp_ret] = psBLIF(scaled, C);

end

function p_total = evaluate_probability_psBLIF(the_function, amplitude)

    [~, alifp_ret] = the_function(amplitude);

    total_probability = nan(size(alifp_ret));

    for k = 1:numel(alifp_ret)
        total_probability(k) = sum([alifp_ret{k}.path_prob]);
    end

    p_total = 1-prod(1-total_probability);
end


function p_total = evaluate_probability_aLIFP(pulse, amplitude)
    [out_dist] = aLIFP(pulse * amplitude);
    p_total = 1-prod(1-[out_dist.total_probability]);
end
