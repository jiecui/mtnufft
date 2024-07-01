function t = get_time_points(this, num_points, method, options)
    % MTNUSP.BRONEZGPSS.GET_TIME_POINTS generates time points for the simulation
    %
    % Syntax:
    %
    % Input(s):
    %
    % Output(s):
    %
    % Example:
    %
    % Note:
    %   DRS         - Deterministic Regular Sampling
    %   DAS         - Deterministic Arithmetic Sampling
    %   RJS         - Random Jittered Sampling
    %   RBUS        - Random Bernoulli Unit Sampling
    %   RBFS        - Random Bernoulli Fractional Sampling
    %   RPS         - Random Poisson Sampling
    %   RUS         - Random Uniform Sampling
    %   TB1         - Regular sampling (DRS) at the Nyquist rate (sampling
    %                 method 1 in Bronez, 1988)
    %   TB2         - Regular over-sampling with some missing samples (RBFS)
    %                 (Sampling method 3 in Bronez, 1988)
    %   TB3         - Arithmetic sampling (DAS) (Sampling method 2 in Bronez,
    %                 1988)
    %   MD_2Hz
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Wed 12/13/2023 10:57:11.432 PM
    % $Revision: 0.3 $  $Date: Tue 06/18/2024 07:36:57.286 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) mtnusp.BronezGPSS
        num_points (1, 1) double {mustBeInteger} % number of time points to generate
        method (1, 1) string {mustBeMember(method, ...
                                  ["DRS", "DAS", "RJS", "RBUS", "RBFS", "RPS", "RUS", "TB1", "TB2", "TB3", "MD_2Hz"])} ...
            = "TB1"
    end % positional

    arguments
        options.ControlProbability (1, 1) double {mustBeNonnegative, mustBeLessThanOrEqual(options.ControlProbability, 1)} = 0.5
        options.ExponentialMu (1, 1) double {mustBePositive} = 1 % mu for mean of exponential probability distribution
        options.JitterMu (1, 1) double {mustBeNonnegative} = 0 % mu for mean of jittered sampling (RJS only)
        options.JitterSigma (1, 1) double {mustBePositive} = .01 % sigma for standard deviation of jittered sampling (RJS only)
        options.MaxNumPoints (1, 1) double {mustBeInteger, mustBePositive} = 1000 % maximum number of time points to generate
        options.Seed (1, 1) double {mustBeInteger, mustBePositive} = 42 % random seed for reproducibility
        options.TimeStart (1, 1) double = 0 % start time (seconds)
        options.TimeEnd (1, 1) double = 1 % end time (seconds)
    end % optional

    jit_mu = options.JitterMu;
    jit_sigma = options.JitterSigma;
    max_num_points = options.MaxNumPoints;
    t_start = options.TimeStart;
    t_end = options.TimeEnd;

    if isempty(t_end)
        T = this.T;
        t_end = t_start + T;
    end % if

    ctrl_p = options.ControlProbability;
    rng_seed = options.Seed;
    mu = options.ExponentialMu;

    % ======================================================================
    % main
    % ======================================================================
    % generate random stream
    % ----------------------
    s = RandStream('mt19937ar', 'Seed', rng_seed);

    % generate time points
    % --------------------
    switch method
        case "DRS"
            t = linspace(t_start, t_end, num_points);
        case "DAS"
            assert(num_points > 2, 'DAS requires number of points > 2')
            ctrl_r = 1 / ctrl_p;
            a = 2 / (1 + ctrl_r);
            b = a * (ctrl_r - 1) / (num_points - 2);
            n = 1:num_points;
            t = t_start + a * (n - 1) + b * (n - 1) .* (n - 2) / 2;
        case "RJS"
            t = this.get_time_points(num_points, "DRS", TimeStart = t_start, ...
                TimeEnd = t_end);
            t = t + normrnd(jit_mu, jit_sigma, size(t));
            t = sort(t);
        case "RBUS"
            t = linspace(t_start, t_end, num_points);

            for k = 2:num_points

                if rand(s, 1) < ctrl_p
                    t(k) = nan;
                end % if

            end % for

            t(isnan(t)) = []; % remove nan values
        case "RBFS"
            n1 = round(.5 + num_points / (1 - ctrl_p));
            tau = num_points / n1;
            t = tau * linspace(t_start, t_end, num_points);

            for k = 2:num_points

                if rand(s, 1) < ctrl_p
                    t(k) = nan;
                end % if

            end % for

            t(isnan(t)) = []; % remove nan values

        case "RPS"
            pd = makedist('Exponential', 'mu', mu);
            t_intv = pd.random(1, max_num_points);
            t = cumsum([t_start, t_intv]);
            t = t(t <= t_end);

        case "RUS"
            t = t_start + (t_end - t_start) .* rand(s, 1, num_points);
        case "TB1"
            t = this.get_time_points(num_points, "DRS", TimeStart = 1, ...
                TimeEnd = num_points);
        case "TB2"
            n = 1:60;
            t = n * 5/6;
            t([1, 5, 17, 18, 19, 23, 27, 32, 53, 56]) = [];
        case "TB3"
            t = this.get_time_points(num_points, "DAS", ...
                TimeStart = t_start, ...
                TimeEnd = t_end, ...
                ControlProbability = 1/4);
        case "MD_2Hz"
            n = 1:60;
            t = n * 5/12;
            t([1, 5, 17, 18, 19, 23, 27, 32, 53, 56]) = [];
    end % switch-case

end % function get_time_points

% [EOF]
