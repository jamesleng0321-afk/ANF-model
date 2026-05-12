%% Miller1999_psBLIF_vs_aLIFP.m
% Reproduction of Miller 1999 Fig 2
% Spike latency and jitter vs amplitude
% Comparison: psBLIF vs aLIFP

clear; clc;

%% ------------------------------------------------------------------------
% Load experimental data
% ------------------------------------------------------------------------

data_latency = readtable("Miller_1999_fig2_latency_jitter.csv");
data_prob    = readtable("Miller_1999_fig2_firing_efficiency.csv");

data_amplitude = data_latency.x;
spiking_probability = data_prob.firingEffienciency / 100;

[data_I50,~] = ...
    aLIFP_fit_gaussian(data_amplitude, spiking_probability);

latency_data = data_latency.Latency;                  % ms
jitter_data  = (data_latency.Jitter - latency_data)*2; % ms

%% ------------------------------------------------------------------------
% Pulse definition
% ------------------------------------------------------------------------

% aLIFP pulse (monophasic)
pulse_aLIFP = [zeros(100,1); -1*ones(40,1); zeros(100,1)];

% psBLIF pulse
phase_len_s = 40e-6;

pulse_psBLIF.pulse_onset        = 100e-6;
pulse_psBLIF.positive_duration  = phase_len_s;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0;
pulse_psBLIF.negative_duration  = 0;
pulse_psBLIF.negative_amplitude = 0;

C = psBLIF_default_parameters();

%% ------------------------------------------------------------------------
% Compute single pulse threshold (reference)
% ------------------------------------------------------------------------

% ---- psBLIF ----

eval_prob_psBLIF = @(amp) psBLIF_final_pulse_prob( ...
                        @(a) psBLIF_wrapper_scale_all(a, ...
                            pulse_psBLIF, C), amp);

psBLIF_I50 = psBLIF_find_threshold(eval_prob_psBLIF, 0.5);

% ---- aLIFP ----

eval_prob_aLIFP = @(amp) aLIFP(pulse_aLIFP * amp).total_probability;
aLIFP_I50 = psBLIF_find_threshold(eval_prob_aLIFP, 0.5);

%% ------------------------------------------------------------------------
% Amplitude sweep
% ------------------------------------------------------------------------

amplitudes = (-3:0.1:3); % db

%% ------------------------------------------------------------------------
% Preallocate
% ------------------------------------------------------------------------

lat_a  = zeros(size(amplitudes));
jit_a  = zeros(size(amplitudes));
prob_a = zeros(size(amplitudes));

lat_ps  = zeros(size(amplitudes));
jit_ps  = zeros(size(amplitudes));
prob_ps = zeros(size(amplitudes));

pulse_start = 101 / 1e6;   % seconds

%% ========================================================================
% MAIN LOOP
% ========================================================================

for a = 1:length(amplitudes)

    scale = 10.^(amplitudes(a) / 20);

    % ---- aLIFP ----
    out_a = aLIFP(pulse_aLIFP * scale * aLIFP_I50);

    prob_a(a) = out_a.total_probability;
    lat_a(a)  = out_a(1).mu - pulse_start;
    jit_a(a)  = out_a(1).sigma;

    % ---- psBLIF ----
    scaled = pulse_psBLIF;
    scaled.positive_amplitude = scaled.positive_amplitude * ...
                                psBLIF_I50 * scale;

    [~, out_ps] = psBLIF(scaled, C);

    prob_ps(a) = sum(out_ps{1}.path_prob);
    lat_ps(a)  = out_ps{1}.lat;
    jit_ps(a)  = out_ps{1}.jit;
end

%% Convert units
lat_a  = lat_a  * 1e3;   % ms
lat_ps = lat_ps * 1e3;

jit_a  = jit_a  * 1e6;   % µs
jit_ps = jit_ps * 1e6;

%% ========================================================================
% Plot (single figure, two panels)
% ========================================================================

fig = figure("Name","Miller1999_psBLIF_vs_aLIFP", ...
             "DefaultAxesFontSize",13);
fig.InnerPosition(3:4) = [900 420];

tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

x_model = 10.^(amplitudes/20);          % Amplitude / I50 (model)
x_data  = data_amplitude/data_I50;      % Amplitude / I50 (data)

%% ------------------------------------------------------------------------
% Latency
% ------------------------------------------------------------------------

nexttile
plot(x_data, latency_data, ...
     'x','LineWidth',2,'MarkerSize',10);
hold on

plot(x_model, lat_a, ...
     '--','LineWidth',2,"Color","#0072BD");

plot(x_model, lat_ps, ...
     '-','LineWidth',2,"Color","#D95319");

ylabel("Spike latency (\mu_{st}-t_{start}) [ms]")
xlabel("Amplitude / I_{50}")
title("Latency")
grid on
xlim([0.89 1.21])
ylim([0.5 0.81])

%% ------------------------------------------------------------------------
% Jitter
% ------------------------------------------------------------------------

nexttile
plot(x_data, jitter_data*1e3, ...
     'x','LineWidth',2,'MarkerSize',10);
hold on

plot(x_model, jit_a, ...
     '--','LineWidth',2,"Color","#0072BD");

plot(x_model, jit_ps, ...
     '-','LineWidth',2,"Color","#D95319");

ylabel("Spike jitter (\sigma_{st}) [µs]")
xlabel("Amplitude / I_{50}")
title("Jitter")
grid on
xlim([0.89 1.21])
ylim([10 135])

%% ------------------------------------------------------------------------
% Shared legend
% ------------------------------------------------------------------------

lg = legend("data","aLIFP","psBLIF", ...
            "Orientation","horizontal");
lg.Layout.Tile = "north";
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

function p = psBLIF_final_pulse_prob(fun, amp)
    [~, alifp_ret] = fun(amp);
    p = sum([alifp_ret{end}.path_prob]);
end
