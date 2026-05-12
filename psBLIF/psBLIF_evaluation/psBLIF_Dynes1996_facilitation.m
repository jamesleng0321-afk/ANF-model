%% Dynes1996_psBLIF_vs_aLIFP.m
% Reproduction of Dynes 1996 Fig 4-1
% Comparison between psBLIF and aLIFP
%
% INFO:
% Dynes1996 use monophasic. Felsheim use biphasic.
% In the model, monophasic makes significant difference to biphasic
% as the voltage carried over is negative for biphasic pulses therfore requiring
% much larger amps for the subsequent pulses.

clear; clc;

%% ------------------------------------------------------------------------
% Load experimental data
% ------------------------------------------------------------------------

measured_data = readtable("Dynes_1996_fig4-1.csv");

masker_counts = [1 4 40];        % conditions in paper
phase_len     = 100e-6;          % 100 us
% In the paper, pulses are 100us cathodic monophasic pulses
% as of aLIPF we use 100us biphasic pulses

amplitudes = [0.01:0.01:2] * 1e-3;

%% ------------------------------------------------------------------------
% Model parameters
% ------------------------------------------------------------------------

C = psBLIF_default_parameters();
C.use_facil = 1;
C.use_alifp_facil = 0;
% C.tau = C.tau/4
C.dont_carry_voltage = 0;

%% ------------------------------------------------------------------------
% Compute single pulse threshold (reference)
% ------------------------------------------------------------------------

% ---- psBLIF ----
single_pulse = build_psBLIF_pulse_train(1, phase_len);

eval_prob_single = @(amp) psBLIF_final_pulse_prob( ...
                        @(a) psBLIF_wrapper_scale_probe(a, single_pulse, C), amp);

I50_single_pulse = psBLIF_get_threshold(eval_prob_single, amplitudes, 0.5);

% ---- aLIFP ----
% pulse_aLIFP = -1 * ones(phaseL_us,1);

% eval_single_a = @(amp) aLIFP_prob_last_pulse(pulse_aLIFP, amp);

% I50_single_a = psBLIF_get_threshold(eval_single_a, amplitudes, 0.5);

%% ------------------------------------------------------------------------
% Masker amplitude (89% of single pulse threshold)
% ------------------------------------------------------------------------
% as of paper, maskers are 1dB below resting threshold

masker_amp = 0.89 * I50_single_pulse;
% masker_amp = 0.01 * I50_single_pulse;

% masker_amp_a  = 0.89 * I50_single_a;

%% ------------------------------------------------------------------------
% Loop over masker conditions
% ------------------------------------------------------------------------

for m_ind = 1:numel(masker_counts)

    num_masker = masker_counts(m_ind);

    masker_probe_interval = measured_data.(sprintf('x%d', num_masker)) * 1e-3; % data is in ms -> s
    masker_probe_ration_db = measured_data.(sprintf('y%d', num_masker));

    % remove nan
    not_nan = ~isnan(masker_probe_interval);
    masker_probe_interval = masker_probe_interval(not_nan);
    masker_probe_ration_db = masker_probe_ration_db(not_nan);

    masker_probe_ratio = 10.^(masker_probe_ration_db / 20);  % convert dB → linear ratio

    %% preallocate for output
    I50_ratio = nan(size(masker_probe_interval));
    % I50_ratio_a  = nan(size(masker_probe_interval));

    %% Build fixed masker blocks

    % ---- psBLIF masker train ----
    % masker_ps = build_psBLIF_pulse_train(num_masker, phase_len);

    % ---- aLIFP masker signal ----
    % masker_block_a = [];
    % for k = 1:num_masker
    %     masker_block_a = [masker_block_a; ...
    %                       zeros(1000 - phaseL_us,1); ...
    %                       pulse_aLIFP * masker_amp_a];
    % end

    %% Loop over IPIs

    for ind = 1:numel(masker_probe_interval)

        mpi  = masker_probe_interval(ind);

        % ipi_us = round(masker_probe_interval(ind));

        %% ---------------- psBLIF ----------------

        % Build full pulse train: maskers + probe
        probe_train = build_masker_probe_psBLIF( ...
                        num_masker, ...
                        masker_amp, ...
                        mpi, ...
                        phase_len);

        eval_probe_ps = @(amp) psBLIF_final_pulse_prob( ...
                            @(a) psBLIF_wrapper_scale_probe(a, probe_train, C), amp);

        I50_probe_ps = psBLIF_get_threshold(eval_probe_ps, amplitudes, 0.5);

        I50_ratio(ind) = I50_probe_ps / I50_single_pulse;

        %% ---------------- aLIFP ----------------

        % eval_probe_a = @(amp) aLIFP_prob_last_pulse( ...
        %                       build_probe_aLIFP(masker_block_a, ...
        %                                         pulse_aLIFP, ...
        %                                         ipi_us, ...
        %                                         amp), 1);

        % I50_probe_a = psBLIF_get_threshold(eval_probe_a, amplitudes, 0.5);

        % I50_ratio_a(ind) = I50_probe_a / I50_single_a;

    end

    %% --------------------------------------------------------------------
    % Plot
    %% --------------------------------------------------------------------

    fig = figure("Name", sprintf("Dynes1996_%dMaskers", num_masker), ...
                 "DefaultAxesFontSize",13);
    fig.OuterPosition(3:4) = [650 600];

    hold on;

    plot(masker_probe_interval*1e3, masker_probe_ratio, 'x', ...
        'MarkerSize',10,'LineWidth',2);

    plot(masker_probe_interval*1e3, I50_ratio, ...
        'LineWidth',2, "Color","#0072BD");

    % plot(masker_probe_interval*1e3, I50_ratio_a, ...
    %     '--','LineWidth',2, "Color","#D95319");

    legend("data","psBLIF","aLIFP","Location","southeast");
    ylabel("I_{50} / Single pulse I_{50}");
    xlabel("Masker probe Interval [ms]");
    ylim([0.2 1.3]);
    grid on;

end

%% ========================================================================
% Helper functions
% ========================================================================

function pulse_train = build_psBLIF_pulse_train(n, phase_dur)
    % 1 ms ipi
    % amp must be scaled later

    pulse_train = struct([]);

    for k = 1:n
        pulse_train(k).pulse_onset        = (k-1)*1e-3; % dummy spacing
        pulse_train(k).positive_duration  = phase_dur;
        pulse_train(k).positive_amplitude = 1;
        pulse_train(k).interphase_gap     = 0;
        pulse_train(k).negative_duration  = phase_dur;
        pulse_train(k).negative_amplitude = -1;
    end
end

function p = psBLIF_final_pulse_prob(fun, amp)
    [~, alifp_ret] = fun(amp);
    p = sum([alifp_ret{end}.path_prob]);
    % fprintf("end %d\n", p);
    % fprintf("fst %d %d\n", sum([alifp_ret{1}.path_prob]), amp);
end

% function p = aLIFP_prob_last_pulse(signal, amp)
%     out = aLIFP(signal * amp);
%     p = sum([out(end).total_probability]);
% end

% function sig = build_probe_aLIFP(masker_block, pulse, ipi_us, amp)
%     sig = [masker_block; ...
%            zeros(ipi_us - numel(pulse),1); ...
%            pulse * amp; ...
%            zeros(100,1)];
% end


function pulse_train = build_masker_probe_psBLIF( ...
                        num_masker, masker_amp, mpi, phase_dur)
    % amp of masker is fixed,
    % amp of probe must be scaled later

    pulse_train = struct([]);

    % ---- Maskers ----
    for k = 1:num_masker
        pulse_train(k).pulse_onset        = (k-1) * 1e-3; % 1 ms spacing
        pulse_train(k).positive_duration  = phase_dur;
        pulse_train(k).positive_amplitude = masker_amp;
        pulse_train(k).interphase_gap     = 0;
        pulse_train(k).negative_duration  = phase_dur;
        pulse_train(k).negative_amplitude = -masker_amp;
    end

    % ---- Probe ----
    pulse_train(num_masker+1).pulse_onset        = ...
        pulse_train(num_masker).pulse_onset + mpi;

    pulse_train(num_masker+1).positive_duration  = phase_dur;
    pulse_train(num_masker+1).positive_amplitude = 1;   % scaled later
    pulse_train(num_masker+1).interphase_gap     = 0;
    pulse_train(num_masker+1).negative_duration  = phase_dur;
    pulse_train(num_masker+1).negative_amplitude = -1;
end

function [history, alifp_ret] = ...
        psBLIF_wrapper_scale_probe(amplitude, pulse_train, C)

    scaled = pulse_train;

    % Only scale last pulse (probe)
    last = numel(scaled);

    scaled(last).positive_amplitude = ...
        scaled(last).positive_amplitude * amplitude;

    scaled(last).negative_amplitude = ...
        scaled(last).negative_amplitude * amplitude;

    [history, alifp_ret] = psBLIF(scaled, C);
end
