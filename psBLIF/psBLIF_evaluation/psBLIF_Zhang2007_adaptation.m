%% Zhang2007_psBLIF.m
% Spike rate over time bins
clear; clc;

%% ------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------

data = readtable("Zhang_2007_fig2.csv");


%% ------------------------------------------------------------------------
% Pulse definition (biphasic)
% ------------------------------------------------------------------------

pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = 40e-6;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0;
pulse_psBLIF.negative_duration  = 40e-6;
pulse_psBLIF.negative_amplitude = -1;

C = psBLIF_default_parameters();


%% ------------------------------------------------------------------------
% Amplitudes
% ------------------------------------------------------------------------

%% fit I50

eval_prob = @(amp) psBLIF_final_pulse_prob( ...
                    @(a) psBLIF_wrapper_scale_all(a, pulse_psBLIF, C), amp);

psBLIF_I50 = psBLIF_find_threshold(eval_prob, 0.5);

%% felsheim
% felsheim amplitude calculations are strange to me.

% base_dB = 0;
% lower_offset_dB = -1.7;
% upper_offset_dB = 0.9;
% amplitudes_dB = [base_dB + lower_offset_dB, base_dB + lower_offset_dB, base_dB + lower_offset_dB; ...
                %  base_dB,                   base_dB,                   base_dB + upper_offset_dB];



%% fitted amplitudes
% I fitted the the amplitudes to the spikerates in the first bin,
% not sure if this is the right way to go.
amplitudes_dB = [58.4555,     0.496576,      3.07247,
                 58.4555,     6.76132,       5.96307,];

amplitudes_dB(1:2,1) = [.5, 1.3]; % more reasonable

amplitudes = psBLIF_I50 * 10.^(amplitudes_dB / 20);


%% ------------------------------------------------------------------------
% Time bins
% ------------------------------------------------------------------------

x_bins = [0;4;12;24;36;48;100;200;300] * 1e-3; % seconds
x_values = x_bins(2:end) - diff(x_bins)/2;

pulse_rates = [250, 1000, 5000];

mean_spiking_probability = nan(length(pulse_rates), 2, length(x_bins)-1);

%% ========================================================================
% MAIN LOOP
% ========================================================================

data_inds = [4, 7, 9; 2, 5, 8]; % Data indexing (unchanged!)

for p_ind = 1:length(pulse_rates)

    rate = pulse_rates(p_ind);

    % ---- Build pulse train ----
    pulse_train = build_psBLIF_train(0.3, rate, pulse_psBLIF);

    for amp_ind = 1:2

        amp = amplitudes(amp_ind, p_ind);

        % ---- fit amp to spike rate ----
        % target_rate = data{1, data_inds(amp_ind, p_ind)};
        % bin_width = x_bins(2);
        % train = build_psBLIF_train(bin_width, rate, pulse_psBLIF);
        % amp = fzero(@(a) psBLIF_spike_rate(a, train, C, bin_width)-target_rate, [0+1e-9,1]);

        % amplitudes(amp_ind, p_ind) = amp;
        % amplitudes_dB(amp_ind, p_ind) = 20*log10(amp/psBLIF_I50);

        % ---- Run model ----
        [~, out_ps] = psBLIF_wrapper_scale_all(amp, pulse_train, C);

        % ---- Spike distribution ----
        t_max = pulse_train(end).pulse_onset + 1e-2;
        dist = psBLIF_get_spike_distribution_fast(out_ps, pulse_train, ...
                                                    t_max, 1e6);

        % ---- Bin into spike rate ----
        for ind = 1:(length(x_bins)-1)

            i1 = round(x_bins(ind)     * 1e6) + 1;
            i2 = round(x_bins(ind + 1) * 1e6);

            curr_dist = dist(i1:i2);

            mean_spiking_probability(p_ind, amp_ind, ind) = ...
                sum(curr_dist) / (x_bins(ind+1) - x_bins(ind));
        end
    end
end

%% ========================================================================
% Plot
% ========================================================================





fig = figure("Name","Zhang2007_psBLIF", ...
             "DefaultAxesFontSize",13);
fig.OuterPosition(3:4) = [1200, 700];

t = tiledlayout(2,3,"Padding","none");
t.TileIndexing = 'columnmajor';

for p_ind = 1:length(pulse_rates)
    for amp_ind = 1:2

        ax = nexttile;
        hold on

        % ---- psBLIF ----
        plot(x_values, ...
            squeeze(mean_spiking_probability(p_ind, amp_ind, :)), ...
            'o-', ...
            'LineWidth',2, ...
            "Color","#D95319");

        % ---- data ----
        plot(x_values, ...
            data{:, data_inds(amp_ind, p_ind)}, ...
            'xk', ...
            'LineWidth',2, ...
            'MarkerSize',8);

        xlabel("Time [s]")
        ylabel("Spike rate [spike/s]")

        title(sprintf("%d pps, %.1f dB re I_{50}", ...
            pulse_rates(p_ind), amplitudes_dB(amp_ind, p_ind)))

        grid on

        % ---- dynamic ylim ----
        m = max([ ...
            squeeze(mean_spiking_probability(p_ind, amp_ind, :)); ...
            data{:, data_inds(amp_ind, p_ind)}], [], 'all');

        ylim([0, ceil(m/100)*100]);
    end
end

leg = legend(ax, "psBLIF", "data", ...
    "Orientation","horizontal");
leg.Layout.Tile = 'north';

%% ========================================================================
% Helper functions
% ========================================================================

function pulse_train = build_psBLIF_train(duration_s, rate, pulse)

    ipi = 1/rate;
    n   = floor(duration_s / ipi);

    pulse_train = repmat(pulse, 1, n);

    for k = 1:n
        pulse_train(k).pulse_onset = (k-1)*ipi;
    end
end

function [history, out] = ...
        psBLIF_wrapper_scale_all(amplitude, pulse_train, C)

    scaled = pulse_train;

    for k = 1:numel(scaled)
        scaled(k).positive_amplitude = ...
            scaled(k).positive_amplitude * amplitude;

        scaled(k).negative_amplitude = ...
            scaled(k).negative_amplitude * amplitude;
    end

    [history, out] = psBLIF(scaled, C);
end

function p = psBLIF_final_pulse_prob(fun, amp)
    [~, out] = fun(amp);
    p = sum([out{end}.path_prob]);
end

function rate = psBLIF_spike_rate(amp, pulse_train, C, bin_width)
    [history, alifp_ret] = psBLIF_wrapper_scale_all(amp, pulse_train, C);
    p_vec = psBLIF_per_pulse_prob(alifp_ret);
    rate = sum(p_vec)/bin_width;
end
