%% Miller2011_accommodation_psBLIF_vs_aLIFP.m
% Miller 2011 Fig 1 – PSTH
% Probe recovery after high-rate masker
% PSTHs like in Takanen paper
% Mabye merge with the other Miller 2011 script.

clear; clc;

%% ------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------

% load Miller_2011_fig2.mat

%% ------------------------------------------------------------------------
% Parameters
% ------------------------------------------------------------------------

% In these cases, the maskers were 200-ms-long, 5,000-pulse/s trains,
% and the probes were 250-ms-long, 100-pulse/s pulses trains.
% probe onset occurring 0.2 ms after masker offset.
masker_duration = 0.2;          % seconds
masker_rate    = 5000;          % pps

probe_duration  = 0.25;          % seconds
probe_rate      = 100;          % pps

% target_probe_FE = 0.7;          % mid dynamic range
% threshold_mask_FE = 0.005;      % is this ok?
% --- lets use dB re I50 instead ---



% All electric stimuli were composed of 40 μs/phase
% symmetric biphasic rectangular pulses, with a leading
% cathodic phase
% no IPG but Takanen use 30us

pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = 40e-6;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0;
pulse_psBLIF.negative_duration  = 40e-6;
pulse_psBLIF.negative_amplitude = -1;


C = psBLIF_default_parameters();
C = psBLIF_precompute_canc_idx(pulse_psBLIF, C);

%% ------------------------------------------------------------------------
% Find amplitudes for thresholds
% ------------------------------------------------------------------------

%% Find I50
eval_prob = @(amp) psBLIF_single_pulse_prob(pulse_psBLIF, amp, C);
I50 = psBLIF_find_threshold(eval_prob, 0.5);

%% probe 100 PPS
probe_dB  = -.0;

probe_train = build_psBLIF_train(probe_duration-1e-9, ... % 25 pulses
    probe_rate, pulse_psBLIF);
amp_probe = I50*10^(probe_dB/20);
scaled_probe_train = scale_train_psBLIF(probe_train, amp_probe);

[~, out_probe]  = psBLIF(scaled_probe_train, C);

% firing efficiency
probe_only_fe = firing_efficiency_psBLIF(out_probe)
% Spike rate
probe_only_sr = probe_rate*probe_only_fe
% PSTH
t_max = probe_train(end).pulse_onset + ...
    max([out_probe{end}.lat] + 5*[out_probe{end}.jit]);

[p_probe_only, t_probe_only] = psBLIF_get_spike_distribution_fast(...
    out_probe, probe_train, t_max, 1e5);

figure;
bar(t_probe_only, p_probe_only);
ppp_probe_only = psBLIF_per_pulse_prob(out_probe);

%% masker 5000 PPS
mask_5000_train = build_psBLIF_train(masker_duration-1e-9, ...
    masker_rate, pulse_psBLIF);


%% ========================================================================
% MAIN LOOP
% ========================================================================

masker_dB = [.5, 0, -.4, -.5, -.7]; %

responses = [];
fprintf("Running accommodation simulation...\n")

for lvl_idx = 1:length(masker_dB)
    level = masker_dB(lvl_idx);
    amp_mask = I50*10^(level/20);

    this_scaled_mask = scale_train_psBLIF(mask_5000_train, amp_mask);
    full_train = concat_trains(this_scaled_mask, scaled_probe_train, 0.2e-3);
    [~, out5000_ps]  = psBLIF(full_train, C);

    responses{lvl_idx} = out5000_ps;
end

%% ========================================================================
% Statistics & Plot
% ========================================================================

fig = figure("Name","Miller2008_PSTH", ...
             "DefaultAxesFontSize",13);


tiledlayout(length(masker_dB)+1,1,"TileSpacing","compact","Padding","compact");

nexttile
plot(t_probe_only+full_train(length(mask_5000_train)+1).pulse_onset, ...
        p_probe_only, "LineWidth",1.5);
% pl.Color = [pl.Color, .6];
% bar([full_train.pulse_onset], ppp)
ylabel(sprintf('SR probe: %.1f \n@ %.1f dB re I_{50}', probe_only_sr, probe_dB));
xlim([0, 0.45]);
ylim([0, max(p_probe_only)]);
xticklabels({})
yticklabels({})

for lvl_idx = 1:length(masker_dB)
    resp = responses{lvl_idx};

    % firing efficiency

    mask_fe = firing_efficiency_psBLIF(resp(1:length(mask_5000_train)));
    probe_fe = firing_efficiency_psBLIF(resp(length(mask_5000_train)+1:end));

    % Spike rate
    mask_sr = masker_rate*mask_fe;
    probe_sr = probe_rate*probe_fe;
    % PSTH
    t_max = full_train(end).pulse_onset + ...
        max([resp{end}.lat] + 5*[resp{end}.jit]);

    [p, t] = psBLIF_get_spike_distribution_fast(...
    resp, full_train, t_max, 1e5);

    ppp = psBLIF_per_pulse_prob(resp);

    nexttile
    % figure;
    pl = plot(t_probe_only+full_train(length(mask_5000_train)+1).pulse_onset, ...
         p_probe_only, "LineWidth",1.5);

    % bar([probe_train.pulse_onset]+full_train(length(mask_5000_train)+1).pulse_onset...
        % , ppp_probe_only)
    hold on
    ll = plot(t, p, "LineWidth",1.5,"Color","#D95319");
    ll.Color = [ll.Color, .8];
    % bar([full_train.pulse_onset], ppp)
    ylabel(sprintf('mask %.1f dB \nSR mask: %.1f\nSR probe: %.1f', masker_dB(lvl_idx), mask_sr, probe_sr));
    xlim([0, 0.45]);
    ylim([0, max(p_probe_only)]);

    if lvl_idx ~= length(masker_dB)
        xticklabels({})
    end
    yticklabels({})
end



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

function p = psBLIF_single_pulse_prob(pulse, amp, C)

    scaled = pulse;
    scaled.positive_amplitude = scaled.positive_amplitude * amp;
    scaled.negative_amplitude = scaled.negative_amplitude * amp;

    [~, out] = psBLIF(scaled, C);
    p = sum(out{1}.path_prob);
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
