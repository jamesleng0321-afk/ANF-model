%% ========================================================================
% Miller2008_psBLIF_vs_ANF.m
% [EN] ISI distribution vs pulse rate (psBLIF and ANF models comparison)
% [ZH] ISI分布与脉冲频率的关系 (psBLIF 与 ANF 模型对比)
% ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
% Load data / 加载数据
% ------------------------------------------------------------------------
miller_fig1 = readtable("Miller_2008_fig1.csv");

%% ------------------------------------------------------------------------
% Parameters / 参数设置
% ------------------------------------------------------------------------
pulse_rates = [250, 1000, 5000];
duration_s  = 500e-3;   % 500 ms, as in Miller 2008
C = psBLIF_default_parameters();

onset_bins = [0   12;
              4   50;
              200 300] * 1e-3; % [Early; Mid; Late] bins in seconds

db_offset = [0.8, 0.8, 0.8]; % dB relative to threshold
nTrials_ANF = 50;            

%% ------------------------------------------------------------------------
% Pulse definition & Thresholding / 脉冲定义与阈值计算
% ------------------------------------------------------------------------
% ---- psBLIF Pulse & Threshold ----
pulse_psBLIF.pulse_onset        = 0;
pulse_psBLIF.positive_duration  = 40e-6;
pulse_psBLIF.positive_amplitude = 1;
pulse_psBLIF.interphase_gap     = 0; 
pulse_psBLIF.negative_duration  = 40e-6;
pulse_psBLIF.negative_amplitude = -1;

disp('Calculating psBLIF threshold...');
eval_prob_ps = @(amp) psBLIF_single_pulse_prob(pulse_psBLIF, amp, C);
psBLIF_I50 = psBLIF_find_threshold(eval_prob_ps, .5);

% ---- ANF Pulse & Threshold ----
Fs = 1e6;
NoiseAlpha = 0.8;
% 40us positive, 0us IPG, 40us negative
SinglePulse_ANF = [0, -1*ones(1,40), 1*ones(1,40), 0]; 

disp('Calculating ANF dynamic range & threshold...');
[Level, Prob] = Library.FindThreshold([SinglePulse_ANF, zeros(1,2000)], Fs, NoiseAlpha, 0.0001e-6, @Model_SinglePulse, 1000);
[muSingle, sigmaSingle] = Library.FitNeuronDynamicRange(Level', Prob);
ANF_I50 = muSingle;

%% ------------------------------------------------------------------------
% Histogram bins / 直方图分箱
% ------------------------------------------------------------------------
xBins = (0:0.05:30)*1e-3; % seconds for plotting
xBins_ANF_edges = xBins;  % edges for ANF histogram
xBins_ANF_centers = xBins(1:end-1) + diff(xBins)/2;

%% ========================================================================
% MAIN LOOP / 主循环
% ========================================================================
isi_results_p_ps = cell(length(pulse_rates), length(onset_bins));
isi_results_t_ps = cell(length(pulse_rates), length(onset_bins));
isi_results_p_ANF = cell(length(pulse_rates), length(onset_bins));

for rate_ind = 1:length(pulse_rates)
    rate = pulse_rates(rate_ind);
    fprintf('\nProcessing Rate: %d pps...\n', rate);
    
    % ---- 1. psBLIF Model Processing ----
    pulse_train = build_psBLIF_train(duration_s, rate, pulse_psBLIF);
    amp_psBLIF = psBLIF_I50 * 10^(db_offset(rate_ind)/20);
    
    [history, out_ps] = psBLIF_wrapper_scale_all(amp_psBLIF, pulse_train, C);
    isi_struct = psBLIF_get_ISI_struct(history, out_ps, pulse_train);
    
    for bin_idx = 1:size(onset_bins, 1)
        early = onset_bins(bin_idx,1);
        late  = onset_bins(bin_idx,2);
        % Filter psBLIF ISIs
        mask = (early <= [isi_struct.onset]) & ([isi_struct.onset] <= late);
        binned_struct = isi_struct(mask);
        [p_vec, t] = psBLIF_get_ISI_distribution_fast(binned_struct, 1e6);
        isi_results_p_ps{rate_ind, bin_idx} = p_vec;
        isi_results_t_ps{rate_ind, bin_idx} = t;
    end
    
% ---- 2. ANF Model Processing (Monte Carlo) ----
    amp_ANF = ANF_I50 * 10^(db_offset(rate_ind)/20);
    Istim = Experiment.stim_PulseTrain(SinglePulse_ANF, rate, 100, 0, duration_s, Fs);
    input = Istim * amp_ANF;
    
    % Prepare containers for raw ISIs in each time bin
    raw_isis_bin = cell(1, size(onset_bins, 1));
    
    for tr = 1:nTrials_ANF
        p_noise = Library.oneonfnoise(length(input), NoiseAlpha);
        c_noise = Library.oneonfnoise(length(input), NoiseAlpha);
        [~, SpTimes, ~, ~] = Model_PulseTrain(input, p_noise, c_noise, Fs);
        
        % Calculate ISIs if at least 2 spikes occurred
        if length(SpTimes) > 1
            % ==========================================
            % 【修改这里】强制将 SpTimes 转换为列向量
            % Ensure SpTimes is a column vector to prevent vertcat errors
            % ==========================================
            SpTimes = SpTimes(:); 
            
            ISIs = diff(SpTimes);
            onsets = SpTimes(1:end-1); % The time the ISI started
            
            % Assign ISIs to appropriate bins
            for bin_idx = 1:size(onset_bins, 1)
                early = onset_bins(bin_idx,1);
                late  = onset_bins(bin_idx,2);
                mask = (onsets >= early) & (onsets <= late);
                current_ISIs = ISIs(mask);
                raw_isis_bin{bin_idx} = [raw_isis_bin{bin_idx}; current_ISIs(:)];
            end
        end
    end

    % Compute histograms for ANF ISIs and normalize
    for bin_idx = 1:size(onset_bins, 1)
        counts = histcounts(raw_isis_bin{bin_idx}, xBins_ANF_edges);
        if max(counts) > 0
            isi_results_p_ANF{rate_ind, bin_idx} = counts / max(counts); % Normalize to 1
        else
            isi_results_p_ANF{rate_ind, bin_idx} = zeros(1, length(counts));
        end
    end
end

%% ========================================================================
% Plot (3x3 grid) / 绘图
% ========================================================================
fig = figure("Name","Miller2008_psBLIF_vs_ANF", "DefaultAxesFontSize",13);
fig.Position(3:4) = [1200 900];
tiledlayout(3, 3, "TileSpacing","compact", "Padding","compact");

rate_labels = ["250 pps", "1000 pps", "5000 pps"];
bin_labels  = ["Early", "Mid", "Late"];

for rate_ind = 1:3
    for bin_idx = 1:3
        nexttile; hold on;
        
        % ---- Plot psBLIF (Blue) ----
        t_ps = isi_results_t_ps{rate_ind, bin_idx};
        p_ps = isi_results_p_ps{rate_ind, bin_idx};
        if ~isempty(p_ps) && max(p_ps) > 0
            l_ps = plot(t_ps*1e3, p_ps / max(p_ps), 'LineWidth',2, 'Color',"#0072BD");
        else
            l_ps = plot(NaN, NaN, 'LineWidth',2, 'Color',"#0072BD"); % Empty placeholder
        end
        
        % ---- Plot ANF (Orange) ----
        p_ANF = isi_results_p_ANF{rate_ind, bin_idx};
        l_anf = plot(xBins_ANF_centers*1e3, p_ANF, '-', 'LineWidth',2, 'Color',"#D95319");
        
        % ---- Add Miller Data (Black) ONLY for middle column ----
        lgh = [l_ps, l_anf];
        if bin_idx == 2
            switch pulse_rates(rate_ind)
                case 250,  data_vec = miller_fig1.x250_c2_r3;
                case 1000, data_vec = miller_fig1.x1000_c2_r3;
                case 5000, data_vec = miller_fig1.x5000_c2_r4;
            end
            data_vec(isnan(data_vec)) = 0;
            if max(data_vec) > 0
                data_vec = data_vec / max(data_vec);
            end
            l_data = plot(xBins*1e3, data_vec, 'Color',"black", 'LineWidth',1.5);
            l_data.Color = [l_data.Color, 0.6]; % Transparency
            uistack(l_data, 'top');
            lgh = [l_ps, l_anf, l_data];
        end
        
        xlim([0 22]);
        ylim([0 1.05]);
        
        % Annotations
        if rate_ind == 1
            title(sprintf('%d - %d ms (Onset)', onset_bins(bin_idx,:)*1e3));
        end
        if rate_ind == 3
            xlabel("ISI [ms]");
        end
        if bin_idx == 1
            ylabel(sprintf("%s\n%.1f dB re I50\nNorm. Prob.", ...
                            rate_labels(rate_ind), db_offset(rate_ind)));
        end
        grid on; box off;
    end
end

% Construct Legend
if length(lgh) == 3
    lg_labels = ["psBLIF", "ANF", "Miller 2008"];
else
    lg_labels = ["psBLIF", "ANF"];
end
lg = legend(lgh, lg_labels, "Orientation", "horizontal");
lg.Layout.Tile = "north";

%% ========================================================================
% Helper functions / 辅助函数
% ========================================================================
function pulse_train = build_psBLIF_train(duration_s, rate, pulse)
    ipi = 1/rate;
    n   = floor(duration_s / ipi);
    pulse_train = repmat(pulse, 1, n);
    for k = 1:n
        pulse_train(k).pulse_onset = (k-1)*ipi;
    end
end

function [history, out] = psBLIF_wrapper_scale_all(amplitude, pulse_train, C)
    scaled = pulse_train;
    for k = 1:numel(scaled)
        scaled(k).positive_amplitude = scaled(k).positive_amplitude * amplitude;
        scaled(k).negative_amplitude = scaled(k).negative_amplitude * amplitude;
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