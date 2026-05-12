function [threshold, info] = psBLIF_get_threshold( ...
    eval_function, amplitudes, target_probability, varargin)
% psBLIF_get_threshold
%
% Robust threshold finder using bracketing + bisection.
%
% Inputs:
%   eval_function        function handle @(amp) -> [~, alifp_ret]
%   amplitudes          vector of test amplitudes (must be sorted)
%   target_probability  desired spike probability (0 < p < 1)
%
% Optional name-value:
%   'Tolerance'         default: 1e-4
%   'MaxIter'           default: 50
%
% Outputs:
%   threshold           interpolated threshold
%   info                struct with diagnostic info

    %% ------------------------------------------------------------
    % Parse optional arguments
    %% ------------------------------------------------------------

    p = inputParser;
    addParameter(p, 'Tolerance', 1e-4, @(x) x > 0);
    addParameter(p, 'MaxIter', 50, @(x) x > 0);
    parse(p, varargin{:});

    tol     = p.Results.Tolerance;
    maxIter = p.Results.MaxIter;

    %% ------------------------------------------------------------
    % Validate inputs
    %% ------------------------------------------------------------

    if ~isa(eval_function, 'function_handle')
        error('eval_function must be a function handle.');
    end

    if any(diff(amplitudes) <= 0)
        error('amplitudes must be strictly increasing.');
    end

    if target_probability <= 0 || target_probability >= 1
        error('target_probability must be between 0 and 1.');
    end

    %% ------------------------------------------------------------
    % Evaluate coarse sweep and bracket solution
    %% ------------------------------------------------------------

    probs = nan(size(amplitudes));

    for i = 1:numel(amplitudes)
        probs(i) = eval_function(amplitudes(i));

        if i > 1 && probs(i) < probs(i-1)
            warning('Probability is not monotonic increasing.');
        end

        if probs(i) >= target_probability
            break
        end
    end

    if probs(1) >= target_probability
        error('Target probability already exceeded at lowest amplitude.');
    end

    if probs(i) < target_probability
        error('Target probability not reached. Increase amplitude range.');
    end

    % Bracketing interval
    a_low  = amplitudes(i-1);
    a_high = amplitudes(i);
    p_low  = probs(i-1);
    p_high = probs(i);

    %% ------------------------------------------------------------
    % Bisection refinement
    %% ------------------------------------------------------------

    for iter = 1:maxIter

        a_mid = 0.5 * (a_low + a_high);
        p_mid = eval_function(a_mid);

        if abs(p_mid - target_probability) < tol
            threshold = a_mid;
            break
        end

        if p_mid < target_probability
            a_low = a_mid;
            p_low = p_mid;
        else
            a_high = a_mid;
            p_high = p_mid;
        end

        if abs(a_high - a_low) < tol
            threshold = 0.5 * (a_low + a_high);
            break
        end

        if iter == maxIter
            warning('Maximum iterations reached in bisection.');
            threshold = 0.5 * (a_low + a_high);
        end
    end

    %% ------------------------------------------------------------
    % Diagnostics
    %% ------------------------------------------------------------

    info.probabilities = probs;
    info.bracket = [a_low a_high];
    info.iterations = iter;

end
