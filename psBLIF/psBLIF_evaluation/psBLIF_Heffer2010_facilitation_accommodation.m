%% Heffer2010_psBLIF_vs_aLIFP.m
% Reproduction of Heffer 2010 rate-dependent spike probability changes
% Comparison between psBLIF and aLIFP
%
% INFO Heffer2010:
% Onset spike probability was calculated for each stimulation
% rate and current level for all recorded fibers. For a given current
% level, the probability of obtaining a spike within the onset
% period typically increased with stimulation rate. This increase
% in spike probability with stimulation rate may be due to the
% increased number of stimulus pulses within the 0 –2 ms onset
% period providing an increased opportunity for spiking. Alter-
% natively, it may also be caused by interactions between these
% pulses (i.e., facilitation). A simple approach was developed to
% examine the relative contribution of facilitation to the increased
% spike probability.
% A fiber’s response to a single, independent stimulus pulse
% was defined as the cumulative spike probability function,
% calculated over the onset period in 100 us intervals, in re-
% sponse to 200 pulse/s stimulation. The predicted ANF response
% to multiple stimulus pulses presented during the onset period,
% with no facilitation, was calculated as the **sum of multiple,
% independent response** functions. Importantly, as only a single
% spike was ever recorded during the onset period, the predicted
% spike probability was constrained to a maximum of 1. Facili-
% tation was then estimated by measuring the difference between
% the measured and predicted onset spike probability for each
% stimulation rate. For the example shown (Fig. 6), the single
% pulse response probability was ~0.04, the predicted onset
% probability at 5,000 pulse/s was ~0.28, and the measured onset
% probability at 5,000 pulse/s was 1. Thus for this representative
% fiber, facilitation was estimated to increase onset spike proba-
% bility at 5,000 pulse/s by ~0.72.
% .04 * 7 = .28
%
%

clear; clc;

%% ------------------------------------------------------------------------
% Load experimental data
% ------------------------------------------------------------------------

load("Heffer2010.mat")

stim_rates = [200, 1000, 2000, 5000];    % pulses per second
signal_length_s = 2e-3;                  % 2 ms

%% ------------------------------------------------------------------------
% Pulse definition
% ------------------------------------------------------------------------

phase_len_us = 25;
ipg_us       = 8;

pulse_aLIFP = [-1 * ones(phase_len_us,1); ...
                zeros(ipg_us,1); ...
                ones(phase_len_us,1)];

phase_len = phase_len_us*1e-6;
ipg = ipg_us*1e-6;

pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = phase_len;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = ipg;
pulse_psBLIF.negative_duration  = phase_len;
pulse_psBLIF.negative_amplitude = -1;

%% ------------------------------------------------------------------------
% Model parameters
% ------------------------------------------------------------------------

C = psBLIF_default_parameters();

% amplitudes = (0.8:0.1:4) * 1e-3;

%% ------------------------------------------------------------------------
% Target probability levels
% ------------------------------------------------------------------------

low_probs    = (0.02:0.02:0.18);
medium_probs = (0.3:0.05:0.7);
high_probs   = (0.75:0.025:0.95);

%% ------------------------------------------------------------------------
% Compute single pulse thresholds (psBLIF + aLIFP)
% ------------------------------------------------------------------------

% ---- psBLIF ----
single_ps(1) = pulse_psBLIF;

eval_single_ps = @(amp) psBLIF_final_pulse_prob( ...
                        @(a) psBLIF_wrapper_scale_probe(a, single_ps, C), amp);

% low_levels_ps    = arrayfun(@(p) psBLIF_get_threshold(eval_single_ps, amplitudes, p), low_probs);
% medium_levels_ps = arrayfun(@(p) psBLIF_get_threshold(eval_single_ps, amplitudes, p), medium_probs);
% high_levels_ps   = arrayfun(@(p) psBLIF_get_threshold(eval_single_ps, amplitudes, p), high_probs);

low_levels_ps    = arrayfun(@(p) find_threshold(eval_single_ps, p), low_probs);
medium_levels_ps = arrayfun(@(p) find_threshold(eval_single_ps, p), medium_probs);
high_levels_ps   = arrayfun(@(p) find_threshold(eval_single_ps, p), high_probs);


% ---- aLIFP ----
eval_single_a = @(amp) aLIFP_prob_last_pulse(pulse_aLIFP, amp);

% low_levels_a    = arrayfun(@(p) psBLIF_get_threshold(eval_single_a, amplitudes, p), low_probs);
% medium_levels_a = arrayfun(@(p) psBLIF_get_threshold(eval_single_a, amplitudes, p), medium_probs);
% high_levels_a   = arrayfun(@(p) psBLIF_get_threshold(eval_single_a, amplitudes, p), high_probs);

low_levels_a    = arrayfun(@(p) find_threshold(eval_single_a, p), low_probs);
medium_levels_a = arrayfun(@(p) find_threshold(eval_single_a, p), medium_probs);
high_levels_a   = arrayfun(@(p) find_threshold(eval_single_a, p), high_probs);


%% ------------------------------------------------------------------------
% Preallocate
% ------------------------------------------------------------------------

prob_change_ps.low    = nan(length(stim_rates), length(low_levels_ps));
prob_change_ps.medium = nan(length(stim_rates), length(medium_levels_ps));
prob_change_ps.high   = nan(length(stim_rates), length(high_levels_ps));

prob_change_a = prob_change_ps;

%% ------------------------------------------------------------------------
% Loop over stimulation rates
% ------------------------------------------------------------------------

for s_ind = 1:length(stim_rates)

    rate = stim_rates(s_ind);

    %% Build stimulus trains

    % ---- aLIFP stimulus ----
    [stimulus_a, num_pulses] = ...
        aLIFP_get_pulse_train(signal_length_s * 1e6, rate, pulse_aLIFP);

    % ---- psBLIF pulse train ----
    pulse_train_ps = build_rate_psBLIF_train(signal_length_s, rate, pulse_psBLIF);

    %% ---- LOW LEVELS ----
    prob_change_ps.low(s_ind,:) = ...
        compute_rate_change_psBLIF(pulse_train_ps, ...
                                   low_levels_ps, ...
                                   C);

    prob_change_a.low(s_ind,:) = ...
        compute_rate_change_aLIFP(stimulus_a, ...
                                  low_levels_a, ...
                                  num_pulses);

    %% ---- MEDIUM LEVELS ----
    prob_change_ps.medium(s_ind,:) = ...
        compute_rate_change_psBLIF(pulse_train_ps, ...
                                   medium_levels_ps, ...
                                   C);

    prob_change_a.medium(s_ind,:) = ...
        compute_rate_change_aLIFP(stimulus_a, ...
                                  medium_levels_a, ...
                                  num_pulses);

    %% ---- HIGH LEVELS ----
    prob_change_ps.high(s_ind,:) = ...
        compute_rate_change_psBLIF(pulse_train_ps, ...
                                   high_levels_ps, ...
                                   C);

    prob_change_a.high(s_ind,:) = ...
        compute_rate_change_aLIFP(stimulus_a, ...
                                  high_levels_a, ...
                                  num_pulses);
end

%% ------------------------------------------------------------------------
% Plot
% ------------------------------------------------------------------------

visual_off = 0.07;

fig = figure("Name","Heffer2010_comparison","DefaultAxesFontSize",13);
fig.OuterPosition(3:4) = [1500 600];
tiledlayout(1,3);

plot_panel(data_Heffer2010(3), prob_change_ps.high, ...
           prob_change_a.high, stim_rates, visual_off, "high");

plot_panel(data_Heffer2010(2), prob_change_ps.medium, ...
           prob_change_a.medium, stim_rates, visual_off, "medium");

plot_panel(data_Heffer2010(1), prob_change_ps.low, ...
           prob_change_a.low, stim_rates, visual_off, "low");

lg = legend("data", "psBLIF", "aLIFP", 'numcolumns', 3);
lg.Layout.Tile = "North";

%% ========================================================================
% Helper Functions
% ========================================================================

function pulse_train = build_rate_psBLIF_train(signal_len, rate, pulse)

    ipi = 1 / rate;
    n   = 1 + floor(signal_len / ipi);

    for k = 1:n
        pulse_train(k) = pulse;
        pulse_train(k).pulse_onset = (k-1) * ipi;
    end
end

function thr = find_threshold(eval_function, target_probability)
    thr = fzero(@(a) eval_function(a)-target_probability, [0,1]);
end

function p = psBLIF_final_pulse_prob(fun, amp)
    [~, alifp_ret] = fun(amp);
    p = sum([alifp_ret{end}.path_prob]);
end

function p = aLIFP_prob_last_pulse(pulse, amp)
    out = aLIFP(pulse * amp);
    p = sum([out(end).total_probability]);
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


function change = compute_rate_change_psBLIF(pulse_train, levels, C)

    change = zeros(1,length(levels));

    for i = 1:length(levels)

        level = levels(i);

        [~, ps_out] = psBLIF_wrapper_scale_all(level, pulse_train, C)

        per_pulse_prob = psBLIF_per_pulse_prob(ps_out);

        p_first = per_pulse_prob(1); % first pulse probability
        p_total = 1-prod(1-per_pulse_prob);

        change(i) = p_total - min(1, p_first * numel(pulse_train));
    end
end

function change = compute_rate_change_aLIFP(stimulus, levels, num_pulses)

    change = zeros(1,length(levels));

    for i = 1:length(levels)

        out = aLIFP(stimulus * levels(i));
        p_first = out(1).total_probability;
        p_onset = 1-prod(1-[out.total_probability]);          % correct imo
        % p_onset = min(1,sum([out.total_probability])); % wrong imo
        predicted = min(1, p_first(1) * num_pulses);

        change(i) = p_onset - predicted;
    end
end

function plot_panel(data_struct, ps_change, a_change, rates, off, label)

    nexttile

    lower_error = data_struct.change(1:4,2) - ...
                  data_struct.change(5:2:end,2);

    upper_error = data_struct.change(6:2:end,2) - ...
                  data_struct.change(1:4,2);

    errorbar((1:4)-off, data_struct.change(1:4,2), ...
             lower_error, upper_error, ...
             'square','linewidth',2,"color","black");
    hold on

    median_ps = median(ps_change,2);
    q25_ps = prctile(ps_change,25,2);
    q75_ps = prctile(ps_change,75,2);

    errorbar((1:4)+off, median_ps, ...
             median_ps-q25_ps, q75_ps-median_ps, ...
             'square','linewidth',2,"color","#0072BD");

    median_a = median(a_change,2);
    q25_a = prctile(a_change,25,2);
    q75_a = prctile(a_change,75,2);

    errorbar((1:4)+2.5*off, median_a, ...
            median_a-q25_a, q75_a-median_a, ...
            'square','linewidth',2,"color","#D95319");

    ylabel("Spike probability change")
    xticks(1:4)
    xticklabels(string(rates))
    xlabel("Pulse rate [pps]")
    title(label)
    ylim([-0.4 1])
    xlim([0.5 4.5])
    yline(0)
end
