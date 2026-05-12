function thr = psBLIF_find_threshold(eval_function, target_probability)
    thr = fzero(@(a) eval_function(a)-target_probability, [0+1e-9,1]);
end
