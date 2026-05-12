function pulse = psBLIF_synthesize_pulse(params, fs)

    pos_n = round(params.positive_duration * fs);
    ipg_n = round(params.interphase_gap    * fs);
    neg_n = round(params.negative_duration * fs);

    pulse = zeros(pos_n + ipg_n + neg_n, 1);
    pulse(1:pos_n) = params.positive_amplitude;

    if neg_n > 0
        pulse(end-neg_n+1:end) = params.negative_amplitude;
    end
end
