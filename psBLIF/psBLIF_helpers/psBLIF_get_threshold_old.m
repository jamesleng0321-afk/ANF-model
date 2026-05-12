% -------------------------------------------------------------------------------------------------
% GPL3
% -------------------------------------------------------------------------------------------------
function [threshold] = psBLIF_get_threshold_old(the_function, amplitudes, target_probability)

    last_probability = nan;
    probabilities = nan(size(amplitudes));
    threshold = nan;

    for a_ind = 1:length(amplitudes)

        [~, alifp_ret] = the_function(amplitudes(a_ind));

        % sum is fine for single thresholds, but Cartee00 treats two low IPI pulses
        % as single pulse with single response, so sum is fine there as well.
        % -> spiking in response to first or second pulse is fine. For longer IPI
        % sum is not correct anymore.
        curr_probability = 0;
        for k = 1:numel(alifp_ret)
            curr_probability = curr_probability + sum([alifp_ret{k}.path_prob]);
        end

        if curr_probability >= target_probability
            if a_ind == 1
                error("Too low amplitudes were given. Please start with a lower values.")
            end
            m = (amplitudes(a_ind) - amplitudes(a_ind - 1)) / (curr_probability - last_probability);
            threshold = m * (target_probability - last_probability) + amplitudes(a_ind - 1);
            break;
        end

        last_probability = curr_probability;
        probabilities(a_ind) = curr_probability;
        if curr_probability >= 1
            break;
        end
    end
end
