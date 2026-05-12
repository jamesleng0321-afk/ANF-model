function C = psBLIF_precompute_canc_idx(pulse_params, C)
    % If the pulse shape is always the same, we can precompute
    pulse = psBLIF_synthesize_pulse(pulse_params, C.fs);
    canc_idx = psBLIF_find_canc_idx(pulse);

    % store in parameters
    C.precompute = true;
    C.precomputed_canc_idx = canc_idx;
    C.precomputed_pulse = pulse;
end
