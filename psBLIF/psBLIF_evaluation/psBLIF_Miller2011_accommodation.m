%% Miller2011_accommodation_psBLIF_vs_aLIFP.m
% Miller 2011 Fig 2 – Accommodation
% Probe recovery after high-rate masker
% Comparison: aLIFP vs psBLIF
%
% Unclear to me how to set the levels
% old Takanen script said 70% firing efficiency is targeted
% but is this for the first pulse or the whole train?
%
% Experimental maskers and probes were pulse
% trains, each typically having 300 ms train durations,
% unless noted otherwise. A post probe-train silent
% interval of 1,200 ms was used to facilitate neural
% recovery. Two masker pulse rates (250 and
% 5,000 pulse/s) were systematically varied, and masker
% current levels were varied so that probe responses
% were obtained for sub-threshold maskers (i.e., too low
% to elicit any spikes) and rates approaching maximal
% (saturation) rates. Probe pulses were presented at
% 100 pulse/s, a compromise between high rates for
% good temporal resolution and low rates that avoid
% cumulative stimulus-evoked effects
% Once a fiber was
% encountered, probe level was fixed to achieve a mid-
% range firing efficiency (FE), typically between 30%
% and 70%. Probe-alone FEs at or near 100% FE were
% avoided so that FE and probe-level changes would be
% correlated. Masker level was then varied to obtain a
% range of masking effects. Finally, when fiber contact
% time permitted, a series of post-stimulus-time histo-
% grams (PSTHs) were collected for multiple probe
% levels, so as to explore the effect of the magnitude of
% the probe response.
% ...
% Probe levels were selected so that
% responses were in the upper part of the modeled
% fiber’s dynamic range (i.e., FEs between 80% and
% 95%).
%
% PROBE
% level fixed (FE~90%)
% 250ms long (300?)
% 100pps
% Probe-train onset followed the masker offset by a
% delay equal to one masker interpulse interval
%
% MASK
% mask level varies (30-70%)
% dB rel Mask threshold (see below)
% 200ms long (300?)
% 250pps or 5000pps
%
% PULSE
% 40 μs/phase symmetric biphasic rectangular pulses NO IPG
%
% Group trends
% Group analyses of cat data (Fig. 2 A, B) show robust
% sub-threshold masking that is dependent on masker
% level and pulse rate. The abscissa (“Masker level”) is
% defined relative to the lowest masker level that evoked
% spikes for each ANF. The ordinate is the ratio of the
% number of spikes evoked by masked and unmasked
% probe trains
%
% Felsheim probe only levels WAY too low. 0.007! should be 0.7
% masker levels also too low!

clear; clc;

%% ------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------

load Miller_2011_fig2.mat

%% ------------------------------------------------------------------------
% Parameters
% ------------------------------------------------------------------------


target_probe_FE = 0.7;          % mid dynamic range
masker_rates    = [250 5000];   % pps
probe_rate      = 100;          % pps
threshold_mask_FE = 0.005;      % is this ok?


probe_duration  = 0.3;          % seconds
masker_duration = 0.3;          % seconds

levels_dB = -4:1:5;             % masker levels (relative to threshold)

% All electric stimuli were composed of 40 μs/phase
% symmetric biphasic rectangular pulses, with a leading
% cathodic phase
% no IPG but Takanen use 30us

phase_len = 40; % us
ipg       = 0;

pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = phase_len*1e-6;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = ipg*1e-6;
pulse_psBLIF.negative_duration  = phase_len*1e-6;
pulse_psBLIF.negative_amplitude = -1;


C = psBLIF_default_parameters();

%% ------------------------------------------------------------------------
% Find amplitudes for thresholds
% ------------------------------------------------------------------------

%% probe 100 PPS

probe_train = build_psBLIF_train(probe_duration, probe_rate, pulse_psBLIF);

eval_prob_psBLIF = @(amp) wrap_prob_psBLIF(@firing_efficiency_psBLIF, ...
    scale_train_psBLIF(probe_train, amp), C);

% probe_amp = psBLIF_find_threshold(eval_prob_psBLIF, target_probe_FE);
probe_amp = 1.2709*1e-3;

scaled_probe_train = scale_train_psBLIF(probe_train, probe_amp);

[~, probe_out] = psBLIF(scaled_probe_train, C);
probe_porbability = firing_efficiency_psBLIF(probe_out); % unmasked prob


%% masker 250 PPS

mask_250_train = build_psBLIF_train(masker_duration, masker_rates(1), pulse_psBLIF);

eval_prob_psBLIF = @(amp) wrap_prob_psBLIF(@firing_efficiency_psBLIF, ...
    scale_train_psBLIF(mask_250_train, amp), C);

% mask_250_amp = psBLIF_find_threshold(eval_prob_psBLIF, threshold_mask_FE);
mask_250_amp = 1.0601*1e-3;

scaled_mask_250_train = scale_train_psBLIF(mask_250_train, mask_250_amp);


%% masker 5000 PPS

mask_5000_train = build_psBLIF_train(masker_duration, masker_rates(2), pulse_psBLIF);

eval_prob_psBLIF = @(amp) wrap_prob_psBLIF(@firing_efficiency_psBLIF, ...
    scale_train_psBLIF(mask_5000_train, amp), C);

% mask_5000_amp = psBLIF_find_threshold(eval_prob_psBLIF, threshold_mask_FE);
mask_5000_amp = 1.1066*1e-3;

scaled_mask_5000_train = scale_train_psBLIF(mask_5000_train, mask_5000_amp);


%% ------------------------------------------------------------------------
% Preallocate
% ------------------------------------------------------------------------

n_levels = length(levels_dB);

recov250_ps  = zeros(n_levels,1);
recov5000_ps = zeros(n_levels,1);

%% ========================================================================
% MAIN LOOP
% ========================================================================

fprintf("Running accommodation simulation...\n")

for lvl_idx = 1:length(levels_dB)

    level = levels_dB(lvl_idx);
    fac = 10^(level/20);

    % 250 pps
    this_mask_250 = scale_train_psBLIF(scaled_mask_250_train, fac);
    full_train = concat_trains(this_mask_250, scaled_probe_train, 4e-3);
    [~, out250_ps]  = psBLIF(full_train, C);
    probe_response = out250_ps(length(this_mask_250):end);
    mask_prob_250 = firing_efficiency_psBLIF(probe_response);
    recov250_ps(lvl_idx) = mask_prob_250;

    % 5000 pps
    this_mask_5000 = scale_train_psBLIF(scaled_mask_5000_train, fac);
    full_train = concat_trains(this_mask_5000, scaled_probe_train, 0.2e-3);
    [~, out5000_ps]  = psBLIF(full_train, C);
    probe_response = out5000_ps(length(this_mask_5000):end);
    mask_prob_5000 = firing_efficiency_psBLIF(probe_response);
    recov5000_ps(lvl_idx) = mask_prob_5000;

end

%% Normalize

recov250_ps  = recov250_ps  / probe_porbability;
recov5000_ps = recov5000_ps / probe_porbability;

%% ========================================================================
% Plot
% ========================================================================

fig = figure("Name","Miller2011_psBLIF_vs_aLIFP", ...
             "DefaultAxesFontSize",13);
fig.Position(3:4) = [1000 420];

tiledlayout(1,2,"Padding","compact","TileSpacing","compact");

% ---- 250 pps ----
nexttile
hold on

scatter(accommodation_Miller2011.Masker250.Raw.MaskerLevel, ...
        accommodation_Miller2011.Masker250.Raw.Recovery, ...
        'o','MarkerEdgeColor','#5b802e','MarkerEdgeAlpha',0.4);

plot(accommodation_Miller2011.Masker250.Medians.MaskerLevel, ...
     accommodation_Miller2011.Masker250.Medians.Recovery, ...
     "kx-","LineWidth",2,"MarkerSize",8);

plot(levels_dB, recov250_ps, "o-","LineWidth",2);

xlabel("Masker level [dB re thr]")
ylabel("Probe recovery ratio")
title("250 pps")
grid on

% ---- 5000 pps ----
nexttile
hold on

scatter(accommodation_Miller2011.Masker5000.Raw.MaskerLevel, ...
        accommodation_Miller2011.Masker5000.Raw.Recovery, ...
        'o','MarkerEdgeColor','#5b802e','MarkerEdgeAlpha',0.4);

plot(accommodation_Miller2011.Masker5000.Medians.MaskerLevel, ...
     accommodation_Miller2011.Masker5000.Medians.Recovery, ...
     "kx-","LineWidth",2,"MarkerSize",8);

plot(levels_dB, recov5000_ps, "o-","LineWidth",2);

xlabel("Masker level [dB re thr]")
ylabel("Probe recovery ratio")
title("5000 pps")
grid on

lg = legend("raw data","median data","psBLIF", ...
            "Orientation","horizontal");
lg.Layout.Tile = "north";

%% ========================================================================
% Helper functions
% ========================================================================

function pulse_train = build_psBLIF_train(duration_s, rate, pulse)

    ipi = 1/rate;
    n   = 1 + floor(duration_s/ipi);

    for k = 1:n
        pulse_train(k) = pulse;
        pulse_train(k).pulse_onset = (k-1)*ipi;
    end
end

function scaled = scale_train_psBLIF(train, amp)

    scaled = train;
    for k = 1:length(train)
        scaled(k).positive_amplitude = train(k).positive_amplitude * amp;
        scaled(k).negative_amplitude = train(k).negative_amplitude * amp;
    end
end

function FE = firing_efficiency_psBLIF(out_struct)
    FE = mean(cellfun(@(x) sum([x.path_prob]), out_struct));
end

function p = wrap_prob_psBLIF(prob_calc_fun, stim, C)
    [~, alifp_ret] = psBLIF(stim, C);
    p = prob_calc_fun(alifp_ret);
end

function train_concat = concat_trains(first_train, second_train, gap)
    train_concat = first_train;
    last_pulse = train_concat(end).pulse_onset;
    for k = 1:length(second_train)
        pulse = second_train(k);
        pulse.pulse_onset = last_pulse + gap + pulse.pulse_onset;
        train_concat(end+1) = pulse;
    end
end
