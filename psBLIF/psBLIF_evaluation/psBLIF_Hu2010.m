%% Hu2010_psBLIF_vs_aLIFP.m
% Reproduction of Hu et al. 2010
% Rate-dependent spike rate and vector strength
% Comparison: psBLIF vs aLIFP

% Felhheim: The stimulation amplitudes for
% the models were chosen such that they produce the mean
% spike values in the first time bin for the four groups

clear; clc;

%% ------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------

spike_rate_data      = readtable("Hu_etal_2010_fig5.csv");
vector_strength_data = readtable("Hu_etal_2010_fig6.csv");

%% ------------------------------------------------------------------------
% Stimulation parameters
% ------------------------------------------------------------------------

duration_s = 0.4;
rate       = 5000;      % pps
mod_freq   = 100;       % Hz
bin_width  = 0.05;      % 50 ms bins

time_edges = 0:bin_width:duration_s;
x_values   = time_edges(1:end-1) + bin_width/2;

%% ------------------------------------------------------------------------
% Pulse definitions
% ------------------------------------------------------------------------

% ---- aLIFP ----
pulse_aLIFP = [0; -1*ones(40,1); ones(40,1); 0];

% ---- psBLIF ----
phase_len_s = 40e-6;

pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = phase_len_s;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0;
pulse_psBLIF.negative_duration  = phase_len_s;
pulse_psBLIF.negative_amplitude = -1;

C = psBLIF_default_parameters();


%% ------------------------------------------------------------------------
% Modulation waveform
% ------------------------------------------------------------------------


mod_func = @(t) 1 + 0.1*sin(2*pi*mod_freq*t);

unmodulated_train = build_psBLIF_train(duration_s, rate, pulse_psBLIF);
modulated_train  = modulate_psBLIF_train(unmodulated_train, mod_func);


%% ------------------------------------------------------------------------
% Amplitudes (fit to match rate groups)
% ------------------------------------------------------------------------

amps.modulated   = [1.35, 1.6, 1.88, 2.24] * 1e-3;
amps.unmodulated = [1.43, 1.54, 1.72, 1.95] * 1e-3;

target_spikerates_mod = [];
target_spikerates_mod(1) = spike_rate_data.R1_modulated(1);
target_spikerates_mod(2) = spike_rate_data.R2_modulated(1);
target_spikerates_mod(3) = spike_rate_data.R3_modulated(1);
target_spikerates_mod(4) = spike_rate_data.R4_modulated(1);


target_spikerates_unmod = [];
target_spikerates_unmod(1) = spike_rate_data.R1_unmodulated(1);
target_spikerates_unmod(2) = spike_rate_data.R2_unmodulated(1);
target_spikerates_unmod(3) = spike_rate_data.R3_unmodulated(1);
target_spikerates_unmod(4) = spike_rate_data.R4_unmodulated(1);


%% fit amplitude modulated
amp_mod = [];
for a_ind = 1:length(target_spikerates_mod)
    target_rate = target_spikerates_mod(a_ind);
    n_elements = round(bin_width*rate);
    train = modulated_train(1:n_elements);
    amp = fzero(@(a) psBLIF_spike_rate(a, train, C, bin_width)-target_rate, [0+1e-9,1]);
    amp_mod(a_ind) = amp;
end

%%
amp_mod = [1.0776, 1.2166, 1.3265, 1.6752] * 1e-3;

%% fit amplitude unmodulated
amp_unmod = [];
for a_ind = 1:length(target_spikerates_unmod)
    target_rate = target_spikerates_unmod(a_ind);
    n_elements = round(bin_width*rate);
    train = unmodulated_train(1:n_elements);
    amp = fzero(@(a) psBLIF_spike_rate(a, train, C, bin_width)-target_rate, [0+1e-9,1]);
    amp_unmod(a_ind) = amp;
end


%%
amp_unmod = [1.1377, 1.17, 1.1917, 1.7666] * 1e-3;



%% ========================================================================
% MAIN LOOP (modulated vs unmodulated)
% ========================================================================

amplitudes = amp_unmod;
rate_unmod = {};
vector_strength_unmod = {};

for a_ind = 1:length(amplitudes)

    amp = amplitudes(a_ind);

    [history, alifp_ret] = ...
            psBLIF_wrapper_scale_all(amp, unmodulated_train, C);


    p_vec = psBLIF_per_pulse_prob(alifp_ret);

    %% split
    n_elements = round(bin_width*rate);
    last_element = round(duration_s*rate);
    split_p = reshape(p_vec(1:last_element),n_elements,[]);

    rate_unmod{a_ind} = sum(split_p)/bin_width;


    %% vector strength
    t_max = unmodulated_train(end).pulse_onset + 1e-3;
    [dist, t] = psBLIF_get_spike_distribution_fast(alifp_ret, unmodulated_train, t_max, 1e6);

    n_steps = round(time_edges(2)*1e6);
    split_p = reshape(dist(1:n_steps*8),n_steps,[]);

    vs = zeros(8, 1);
    for p_ind = 1:8
        vs(p_ind) = calculate_vector_strength(split_p(:,p_ind)', 5000);
    end
    vector_strength_unmod{a_ind} = vs;
end


amplitudes = amp_mod;
rate_mod = {};
vector_strength_mod = {};

for a_ind = 1:length(amplitudes)

    amp = amplitudes(a_ind);

    [history, alifp_ret] = ...
            psBLIF_wrapper_scale_all(amp, modulated_train, C);


    p_vec = psBLIF_per_pulse_prob(alifp_ret);

    %% split
    n_elements = round(bin_width*rate);
    last_element = round(duration_s*rate);
    split_p = reshape(p_vec(1:last_element),n_elements,[]);

    rate_mod{a_ind} = sum(split_p)/bin_width;


    %% vector strength
    t_max = modulated_train(end).pulse_onset + 1e-3;
    [dist, t] = psBLIF_get_spike_distribution_fast(alifp_ret, modulated_train, t_max, 1e6);

    n_steps = round(time_edges(2)*1e6);
    split_p = reshape(dist(1:n_steps*8),n_steps,[]);

    vs = zeros(8, 1);
    for p_ind = 1:8
        vs(p_ind) = calculate_vector_strength(split_p(:,p_ind)', 100);
    end
    vector_strength_mod{a_ind} = vs;
end



%% ------------------------------------------------------------------------
% Figure
% ------------------------------------------------------------------------

fig = figure("Name","Hu2010_psBLIF_vs_aLIFP", ...
             "DefaultAxesFontSize",13);
fig.Position(3:4) = [900 900];
tiledlayout(4,2,"TileSpacing","tight","Padding","tight");

for a_ind = 1:length(amplitudes)
    %% plot
    nexttile(a_ind*2 - 1);


    for modulated = [true, false]
        if modulated
            i_off = 0;
            color = "#0072BD";
            this_rate = rate_mod;
            this_vs = vector_strength_mod;
        else
            i_off = 1;
            color = "#D95319";
            this_rate = rate_unmod;
            this_vs = vector_strength_unmod;
        end

        plot(spike_rate_data.x, spike_rate_data{:, a_ind*2+i_off}, ...
            "x", "MarkerSize",8, "LineWidth",2, "Color",color);
        hold on

        plot(x_values*1e3, this_rate{a_ind}, ...
            "-o", "LineWidth",2, "Color",color);
    end

    ylabel("Spike rate [spike/s]")
    xlabel("Time [ms]")
    title(sprintf("Spike rate R%d", a_ind))
    grid on


    nexttile(a_ind * 2);

    for modulated = [true, false]
        if modulated
            i_off = 0;
            color = "#0072BD";
            this_rate = rate_mod;
            this_vs = vector_strength_mod;
        else
            i_off = 1;
            color = "#D95319";
            this_rate = rate_unmod;
            this_vs = vector_strength_unmod;
        end

        plot(vector_strength_data.x, vector_strength_data{:, a_ind*2+i_off}, ...
            "x", "MarkerSize",8, "LineWidth",2, "Color",color);
        hold on
        plot(x_values * 1e3, this_vs{a_ind}, ...
             "-o", "LineWidth", 2, "color", color);
    end


    ylabel("Vector Strength")
    xlabel("Time [ms]")
    title(sprintf("Vector Strength, R%d", a_ind))
    grid on;
end



lg = legend("Data modulated", "psBLIF modulated", ...
            "Data unmodulated", "psBLIF unmodulated", ...
            "Orientation","horizontal");
lg.Layout.Tile = "north";


%% ------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------


function pulse_train = build_psBLIF_train(duration_s, rate, pulse)

    ipi = 1/rate;
    n   = 1 + floor(duration_s/ipi);

    for k = 1:n
        pulse_train(k) = pulse;
        pulse_train(k).pulse_onset = (k-1)*ipi;
    end
end


function modulated = modulate_psBLIF_train(train, mod_func)

    modulated = train;

    for k = 1:numel(modulated)

        t = modulated(k).pulse_onset;
        amplitude = mod_func(t);

        modulated(k).positive_amplitude = ...
            modulated(k).positive_amplitude * amplitude;

        modulated(k).negative_amplitude = ...
            modulated(k).negative_amplitude * amplitude;
    end

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

function rate = psBLIF_spike_rate(amp, pulse_train, C, bin_width)
    [history, alifp_ret] = psBLIF_wrapper_scale_all(amp, pulse_train, C);
    p_vec = psBLIF_per_pulse_prob(alifp_ret);
    rate = sum(p_vec)/bin_width;
end

function vs = calculate_vector_strength(distribution, pulse_rate)

    T = 1/pulse_rate * 1e6;

    t = 1:length(distribution);
    theta = 2*pi*mod(t,T)/T;

    x = cos(theta).*distribution;
    y = sin(theta).*distribution;

    denom = sum(distribution);

    if denom == 0
        vs = 0;
    else
        vs = sqrt(sum(x).^2 + sum(y).^2)/denom;
    end
end

calculate_vector_strength_fun = @calculate_vector_strength;
