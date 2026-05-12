function [aLIFP_train, t] = psBLIF_to_aLIFP_train(psBLIF_train, fs)

    last_pulse = psBLIF_train(end);

    t_max = last_pulse.pulse_onset + ...
                last_pulse.positive_duration + ...
                last_pulse.interphase_gap + ...
                last_pulse.negative_duration;

    t = 0:1/fs:t_max;

    aLIFP_train = zeros(size(t));

    for k = 1:numel(psBLIF_train)
        curr_pulse = psBLIF_train(k);
        pulse_shape = psBLIF_synthesize_pulse(curr_pulse, fs)';

        onset_idx = round(curr_pulse.pulse_onset * fs) + 1;

        aLIFP_train(onset_idx:onset_idx+length(pulse_shape)-1) ...
            = pulse_shape;
    end

end
