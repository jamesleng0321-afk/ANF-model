clear; clc; close
% this is an example code to find out a threshold for a pulse train

%%%%%%%%%%%% Model params %%%%%%%%%%%%
Fs = 1e6;
NoiseAlpha = 0.8;

%%%%%%%%%%% Stimulus params %%%%%%%%%%%%
stim_rate = [100];
stim_duration = 0.3;
stim_PPD = 40e-6;
stim_IPG = 8e-6;
stim_leadingPol =-1;
stim_laggingPol =1;

% Threshold for a single pulse
SinglePulse = [0, stim_leadingPol*ones(1,stim_PPD*Fs),...
    zeros(1,stim_IPG*Fs),...
    stim_laggingPol*ones(1,stim_PPD*Fs), 0]; % make a single sample pulse

% Find out threshold for a single pulse
[Level,Probability]=Library.FindThreshold([SinglePulse, zeros(1,2000)],Fs,NoiseAlpha,0.0001e-6,@Model_SinglePulse,1000);
[muSingle,sigmaSingle]=Library.FitNeuronDynamicRange(Level',Probability);

% Make a pulse train and adjust level
Istim=Experiment.stim_PulseTrain(SinglePulse,stim_rate,100,0,stim_duration,Fs);
input = Istim * muSingle;

% make memebrane noise waveforms
p_noise = Library.oneonfnoise(length(input),NoiseAlpha);
c_noise = Library.oneonfnoise(length(input),NoiseAlpha);

% Run the model
[nspikes,SpTimes,pSpikes,cSpikes] = Model_PulseTrain(input,p_noise,c_noise,Fs);


% Increase the stimulus intensity
% Linear_Level = Threshold * 10^(dB/20)
dB_above_threshold = 0; 
amplitude = muSingle * 10^(dB_above_threshold / 20);

% Generate the unit unit sequence and multiply it by the new amplitude
Istim = Experiment.stim_PulseTrain(SinglePulse,stim_rate,100,0,stim_duration,Fs);
input = Istim * amplitude; 

% Increase the number of repeated trials
nTrials = 10; 
allSpikeTimes = [];

fprintf(' %d  trials...\n', nTrials);
tic();
for i = 1:nTrials

    p_noise = Library.oneonfnoise(length(input),NoiseAlpha);
    c_noise = Library.oneonfnoise(length(input),NoiseAlpha);

    [nspikes, SpTimes, pSpikes, cSpikes] = Model_PulseTrain(input, p_noise, c_noise, Fs);

    allSpikeTimes = [allSpikeTimes; SpTimes(:)]; 
end
fprintf('Plotting PSTH...\n');
timing = toc();

% Plotting PSTH
bin_width = 0.001; 
edges = 0 : bin_width : stim_duration; 

figure;
histogram(allSpikeTimes, edges, 'FaceColor', [0.5 0.5 0.5]); 
xlabel('Time after pulse train onset (s)');
ylabel('# of Spikes');
title(sprintf('PSTH (Rate: %d pps, Level: +%d dB, Trials: %d)', stim_rate, dB_above_threshold, nTrials));
xlim([0 stim_duration]);


max_time = 0.3; 
valid_spikes = allSpikeTimes(allSpikeTimes <= max_time); 

% Plotting gray normal PSTH (bin-width = 1 ms)
bin_width = 0.001; 
edges = 0 : bin_width : max_time;
counts = histcounts(valid_spikes, edges);

% Count to Spikes/second 
rate_PSTH = counts / (nTrials * bin_width); 
t_PSTH = edges(1:end-1) + bin_width/2; 

figure;
hold on;

bar(t_PSTH * 1000, rate_PSTH, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'BarWidth', 1);

% Calculate and plot aPSTH 
aPSTH_edges_ms = [0, 4, 12, 24, 48, 100, 200, 300]; 
aPSTH_edges_s = aPSTH_edges_ms / 1000; 

aPSTH_rates = zeros(1, length(aPSTH_edges_s)-1);
aPSTH_t_centers = zeros(1, length(aPSTH_edges_s)-1);

for i = 1:(length(aPSTH_edges_s)-1)
    t_start = aPSTH_edges_s(i);
    t_end = aPSTH_edges_s(i+1);

    spike_count_in_win = sum(valid_spikes >= t_start & valid_spikes < t_end);

    win_width = t_end - t_start;
    aPSTH_rates(i) = spike_count_in_win / (nTrials * win_width);

    aPSTH_t_centers(i) = aPSTH_edges_ms(i) + 1; 
end

plot(aPSTH_t_centers, aPSTH_rates, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');

xlabel('Time after pulse train onset (ms)');
ylabel('Spikes / second');
title(sprintf('Model PSTH & aPSTH (Rate: %d pps)', stim_rate));
xlim([0 300]);
legend('PSTH', 'aPSTH', 'Location', 'northeast');
box off;
set(gca, 'TickDir', 'out');
hold off;