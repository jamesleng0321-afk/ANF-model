% Demo script for the path-based sequential BLIF (psBLIF) model
clear; clc;

%% ------------------------------------------------------------------------
% Model parameters
% -------------------------------------------------------------------------

C = psBLIF_default_parameters();

%% ------------------------------------------------------------------------
% Build pulse train
% -------------------------------------------------------------------------

n_pulses = 12;
ipi      = 1e-3;   % interpulse interval (1 ms)

pulse_train = struct([]);

for k = 1:n_pulses
    pulse_train(k).pulse_onset        = (k-1) * ipi;
    pulse_train(k).positive_duration  = 100e-6;
    pulse_train(k).positive_amplitude = 385e-6;
    pulse_train(k).interphase_gap     = 5e-6;
    pulse_train(k).negative_duration  = 100e-6;
    pulse_train(k).negative_amplitude = -385e-6;
end

%% ------------------------------------------------------------------------
% Run psBLIF
% -------------------------------------------------------------------------

[history, alifp_ret] = psBLIF(pulse_train, C);

%% ------------------------------------------------------------------------
% Extract spike probabilities per pulse
% -------------------------------------------------------------------------

p_spike = zeros(1, n_pulses);

for k = 1:n_pulses
    paths = alifp_ret{k};
    if isempty(paths)
        p_spike(k) = 0;
    else
        p_spike(k) = sum([paths.path_prob]);
    end
end

%% ------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------

figure;
bar(1:n_pulses, p_spike);
xlabel('Pulse index');
ylabel('Spike probability');
title('psBLIF spike probability per pulse');
grid on;

%% ------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------

figure;
t_max = pulse_train(end).pulse_onset + 1e-3;
[p, t] = psBLIF_get_spike_distribution(alifp_ret, pulse_train, t_max, 1e6);
plot(t, p);
grid on;

hold on;
for k = 1:numel(alifp_ret)
    for m = 1:numel(alifp_ret{k})
        mu = pulse_train(k).pulse_onset+alifp_ret{k}(m).lat;
        sigma = alifp_ret{k}(m).jit;
        p_vec = normpdf(t, mu, sigma) ...
                         * alifp_ret{k}(m).path_prob;
        p_vec = p_vec/1e6;
        plot(t, p_vec);
    end
end

%% ------------------------------------------------------------------------
% Inter-spike interval ISI
% -------------------------------------------------------------------------

isi_dist = psBLIF_get_ISI_struct(history, alifp_ret, pulse_train);

[p_vec, t] = psBLIF_get_ISI_distribution_fast(isi_dist, 1e6);

figure;
plot(t*1e3, p_vec);
grid on;
