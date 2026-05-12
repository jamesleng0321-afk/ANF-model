function isi_dist = psBLIF_get_ISI_struct(history, alifp_ret, pulse_train)
    % Should just store "last spike ancestor" during simulation instead of
    % doing this path searching.

    isi_dist = struct('mu', {}, 'sigma', {}, 'p', {}, ...
                      'onset', {});
                      % save onset time to split into time bins later

    % loop backwards through pulses
    for late_pulse_idx = length(pulse_train):-1:2

        late_pulse_ret = alifp_ret{late_pulse_idx};

        this_onset = pulse_train(late_pulse_idx).pulse_onset;

        % check every spike distribution for current pulse
        for inner_idx = 1:length(late_pulse_ret)
            late_node = late_pulse_ret(inner_idx);

            prev = late_node.prev_idx;
            found_prev_spike = false;
            first_spike = false;

            % find previous spike
            % walk pulses backward
            for early_pulse_idx = late_pulse_idx-1:-1:1
                if found_prev_spike
                    break % exit early_pulse_idx loop
                end
                current_hist = history{early_pulse_idx};

                % check every path node for ancestor
                for inner_prev_idx = 1:length(current_hist)
                    early_node = current_hist(inner_prev_idx);
                    if early_node.path_idx == prev
                        if early_node.spike == true
                            % fprintf('found prev spike\n');
                            found_prev_spike = true;
                            break
                            % continue next inner_idx
                        else
                            % fprintf('found prev node\n');
                            prev = early_node.prev_idx;
                            if prev == 0
                                first_spike = true;
                            end
                            break
                            % continue next early_pulse_idx
                        end
                    end
                end


            end

            % ---- compute ISI ----
            if found_prev_spike
                dist_late_lat = early_node.lat;
                dist_late_jit = early_node.jit;
                dist_late_onset = this_onset;

                dist_early_lat = late_node.lat;
                dist_early_jit = late_node.jit;
                dist_early_onset = pulse_train(early_pulse_idx).pulse_onset;

                this_isi.mu = dist_late_onset - dist_early_onset + ...
                              dist_late_lat - dist_early_lat;
                this_isi.sigma = sqrt(dist_late_jit^2 + dist_early_jit^2);

                this_isi.p = late_node.path_prob;
                this_isi.onset = dist_early_onset;

                isi_dist(end+1) = this_isi;
                if this_isi.mu<0
                    fprintf('mu<0!\n');
                end
            end

            if ~found_prev_spike & ~first_spike
                fprintf('WTFF?\n');
            end
        end
    end
end
