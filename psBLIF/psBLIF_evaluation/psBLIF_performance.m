%% psBLIF vs aLIFP performance


phase_len = 20; % us
ipg       = 10;

pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = phase_len*1e-6;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = ipg*1e-6;
pulse_psBLIF.negative_duration  = phase_len*1e-6;
pulse_psBLIF.negative_amplitude = -1;

pulse_aLIFP = [-ones(phase_len,1); zeros(ipg,1); ones(phase_len,1)];

C = psBLIF_default_parameters();


use_alifp = false;

psBLIF_procompute = 0;

if psBLIF_procompute == true
    C = psBLIF_precompute_canc_idx(pulse_psBLIF, C);
end


%% 0
% find threshold

target = 0.9;

% ---- psblif ----

eval_prob_psBLIF = @(amp) psBLIF_final_pulse_prob( ...
                        @(a) psBLIF_wrapper_scale_all(a, ...
                            pulse_psBLIF, C), amp);

psBLIF_thr = psBLIF_find_threshold(eval_prob_psBLIF, target);

pulse_psBLIF.positive_amplitude = pulse_psBLIF.positive_amplitude*psBLIF_thr;
pulse_psBLIF.negative_amplitude = pulse_psBLIF.negative_amplitude*psBLIF_thr;


% ---- aLIFP ----

if use_alifp
    eval_prob_aLIFP = @(amp) aLIFP(pulse_aLIFP * amp).total_probability;
    aLIFP_thr = psBLIF_find_threshold(eval_prob_aLIFP, target);
end

%% 1 - const duration, var rate

stim_dur_s = 500e-3;

rates = [100, 200, 400, 800, 2000];

timing_psblif_rate = zeros(size(rates));
timing_alifp_rate = zeros(size(rates));

for rate_idx = 1:length(rates)
    rate = rates(rate_idx);
    fprintf('rate %d ...', rate);

    % ---- psblif ----
    pulse_train = build_rate_psBLIF_train(stim_dur_s, rate, pulse_psBLIF);

    t1 = tic;
    [history, alifp_ret] = psBLIF(pulse_train, C);
    t2 = toc(t1);
    timing_psblif_rate(rate_idx) = t2;

    % ---- aLIFP ----
    if use_alifp
        [stimulus_a, num_pulses] = ...
            aLIFP_get_pulse_train(stim_dur_s * 1e6, rate, pulse_aLIFP);
    
        t1 = tic;
        alifp = aLIFP(stimulus_a*aLIFP_thr);
        t2 = toc(t1);
        timing_alifp_rate(rate_idx) = t2;
    end

    fprintf('done\n');
end


%% 2 - const n pulses, var duration

n_pulses = 50;

durations = [10, 50, 250, 500, 1000, 2000]*1e-3;

rates_comp = n_pulses./durations;

timing_psblif_dur = zeros(size(durations));
timing_alifp_dur = zeros(size(durations));

for dur_idx = 1:length(durations)
    rate = rates_comp(dur_idx);
    dur = durations(dur_idx);

    fprintf('dur %d ...', dur);

    % ---- psblif ----
    pulse_train = build_rate_psBLIF_train(dur, rate, pulse_psBLIF);

    t1 = tic;
    [history, alifp_ret] = psBLIF(pulse_train, C);
    t2 = toc(t1);
    timing_psblif_dur(dur_idx) = t2;

    % ---- aLIFP ----
    if use_alifp
        [stimulus_a, num_pulses] = ...
            aLIFP_get_pulse_train(dur * 1e6, rate, pulse_aLIFP);
    
        t1 = tic;
        alifp_ret = aLIFP(stimulus_a*aLIFP_thr);
        t2 = toc(t1);
        timing_alifp_dur(dur_idx) = t2;
    end

    fprintf('done\n');
end

%% plot
fig = figure();

tiledlayout(2,1, 'TileSpacing', 'compact')

nexttile
plot(rates, timing_psblif_rate, ...
    'x-', 'MarkerSize',10,'LineWidth',2,"Color","#D95319");
hold on
plot(rates, timing_alifp_rate, ...
    'x-', 'MarkerSize',10,'LineWidth',2,"Color","#0072BD");
yscale log
grid on
title('Runtime for constant stim duration (500ms)')
xlabel("pps")
ylabel("[s]");

nexttile
plot(durations, timing_psblif_dur, ...
    'x-', 'MarkerSize',10,'LineWidth',2,"Color","#D95319");
hold on
plot(durations, timing_alifp_dur, ...
    'x-', 'MarkerSize',10,'LineWidth',2,"Color","#0072BD");
yscale log
grid on
title('Runtime for constant number of pulses (n=50)')
xlabel("stim duration [s]")
ylabel("[s]");

leg = legend("psBLIF", "aLIFP", "orientation", "horizontal");
leg.Layout.Tile = 'north';

%% helpers

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

function p = psBLIF_final_pulse_prob(fun, amp)
    [~, alifp_ret] = fun(amp);
    p = sum([alifp_ret{end}.path_prob]);
end


function pulse_train = build_rate_psBLIF_train(signal_len, rate, pulse)

    ipi = 1 / rate;
    n   = 1 + floor(signal_len / ipi);

    for k = 1:n
        pulse_train(k) = pulse;
        pulse_train(k).pulse_onset = (k-1) * ipi;
    end
end
