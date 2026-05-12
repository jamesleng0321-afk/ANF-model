function ret = psBLIF_find_canc_idx(i_sig)
% Find cancellation index

    idx_a = find(i_sig > 0, 1, 'last');

    idx_b = idx_a + 1;
    ret = inf(1, idx_a);

    cumm  = 0;
    state = 0; % 0 = add_left, 1 = add_right

    while true
        if state == 0  % add_left
            if idx_a < 1
                break
            end
            cumm = cumm + i_sig(idx_a);
            if cumm >= 0
                state = 1; % add_right
            else
                ret(idx_a) = idx_b;
                idx_a = idx_a - 1;
            end
        end

        if state == 1  % add_right
            if idx_b > numel(i_sig)
                break
            end
            cumm = cumm + i_sig(idx_b);
            if cumm < 0
                ret(idx_a) = idx_b;
                idx_b = idx_b + 1;
                idx_a = idx_a - 1;
                state = 0; % add_left
            else
                idx_b = idx_b + 1;
            end
        end
    end
end
