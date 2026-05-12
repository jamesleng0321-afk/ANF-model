%% Miller2008_spike_rate_adaptation_psBLIF.m
% ISI distribution vs pulse rate
clear; clc;

%% ------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------

miller_fig1 = readtable("Miller_2008_fig1.csv");

%% ------------------------------------------------------------------------
% Parameters
% ------------------------------------------------------------------------

pulse_rates = [250, 1000, 5000];
duration_s  = 500e-3;   % as in Miller 2008

C = psBLIF_default_parameters();

% stimulus levels (taken from original script)
lvl = [1100 1200 1100]*1e-6; % shitty Takanen values.
onset_bins = [0 12;
              4 50;
              200 300] * 1e-3; % s


%% ------------------------------------------------------------------------
% Pulse definition
% ------------------------------------------------------------------------

% stimuli consisted of single, 40-μs/phase biphasic pulses
pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = 40e-6;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0; % Takanen use 10us for whatever reason
pulse_psBLIF.negative_duration  = 40e-6;
pulse_psBLIF.negative_amplitude = -1;

% find threshold I50

eval_prob_ps = @(amp) psBLIF_single_pulse_prob(pulse_psBLIF, amp, C);
I50 = psBLIF_find_threshold(eval_prob_ps, .5);

db_offset = [.8, .8, .8];

%% ------------------------------------------------------------------------
% Histogram bins
% ------------------------------------------------------------------------

xBins = (0:0.05:30)*1e-3; % seconds

ydata = zeros(length(xBins), length(pulse_rates));

%% ========================================================================
% MAIN LOOP
% ========================================================================

rate_isi_structs = {};

for rate_ind = 1:length(pulse_rates)

    rate = pulse_rates(rate_ind);

    % ---- Build pulse train ----
    pulse_train = build_psBLIF_train(duration_s, rate, pulse_psBLIF);

    % ---- Scale stimulus ----
    % amp = lvl(rate_ind);
    amp = I50 * 10^(db_offset(rate_ind)/20);

    [history, out_ps] = psBLIF_wrapper_scale_all(amp, pulse_train, C);

    % ---- Extract ISI distribution ----
    isi_struct = psBLIF_get_ISI_struct(history, out_ps, pulse_train);

    rate_isi_structs{rate_ind} = isi_struct;
end

% % normalize like original
% ydata(1:length(counts), rate_ind) = counts / max(counts);

%% split into bins

isi_results_p = {};
isi_results_t = {};

for rate_ind = 1:length(pulse_rates)

    rate = pulse_rates(rate_ind);

    for bin_idx = 1:length(onset_bins)
        early = onset_bins(bin_idx,1);
        late  = onset_bins(bin_idx,2);

        % filter ISIs based on spike timing
        this_isi_struct = rate_isi_structs{rate_ind};
        mask = (early <= [this_isi_struct.onset]) & ...
               ([this_isi_struct.onset] <= late);
        binned_struct = rate_isi_structs{rate_ind}(mask);

        [p_vec, t] = psBLIF_get_ISI_distribution_fast(binned_struct, 1e6);

        isi_results_p{rate_ind, bin_idx} = p_vec;
        isi_results_t{rate_ind, bin_idx} = t;
    end
end


%% ========================================================================
% Plot (3x3 grid)
% ========================================================================

fig = figure("Name","Miller2008_psBLIF", ...
             "DefaultAxesFontSize",13);
fig.Position(3:4) = [1000 900];

tiledlayout(3,3,"TileSpacing","compact","Padding","compact");

rate_labels = ["250 pps","1000 pps","5000 pps"];
bin_labels  = ["Early","Mid","Late"];

for rate_ind = 1:3
    for bin_idx = 1:3

        nexttile
        hold on

        t = isi_results_t{rate_ind, bin_idx};
        p = isi_results_p{rate_ind, bin_idx};

        l1 = plot(t*1e3, p / max(p), 'LineWidth',2, 'Color',"#D95319");

        % ---- Add Miller data ONLY for middle column ----
        if bin_idx == 2
            switch pulse_rates(rate_ind)
                case 250
                    data_vec = miller_fig1.x250_c2_r3;
                case 1000
                    data_vec = miller_fig1.x1000_c2_r3;
                case 5000
                    data_vec = miller_fig1.x5000_c2_r4;
            end

            data_vec(isnan(data_vec)) = 0;
            data_vec = data_vec / max(data_vec);

            lnh = plot(xBins*1e3, data_vec, 'Color',"black", 'LineWidth',1);
            lnh.Color = [lnh.Color, 0.5];
            uistack(lnh,'top')
            lgh = [l1, lnh];
        end
        xlim([0 22])

        if rate_ind == 1
            title(sprintf('%d-%d ms', onset_bins(bin_idx,:)*1e3))
        end

        if rate_ind == 3
            xlabel("ISI [ms]")
        end

        if bin_idx == 1
            ylabel(sprintf("%s\n%.1f dB re I50\n", ...
                            rate_labels(rate_ind), db_offset(rate_ind)));
        end

        grid on
    end
end

lg = legend(lgh, ["psBLIF", "Miller 2008"], ...
            "Orientation","horizontal");
lg.Layout.Tile = "north";

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

function p = psBLIF_single_pulse_prob(pulse, amp, C)

    scaled = pulse;
    scaled.positive_amplitude = scaled.positive_amplitude * amp;
    scaled.negative_amplitude = scaled.negative_amplitude * amp;

    [~, out] = psBLIF(scaled, C);
    p = sum(out{1}.path_prob);
end
