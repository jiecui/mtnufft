function s_table = mtnu_speed_analysis(rng_seed, options)
    % ANALYSIS.MTNU_SPEED_ANALYSIS compare the speed of multitaper nonuniform methods
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
    %
    % References:
    %
    % See also .

    % Copyright 2024 Richard J. Cui. Created: Mon 07/08/2024 16:57:45.065 PM
    % $Revision: 0.1 $  $Date: Mon 07/08/2024 16:57:45.065 PM $
    %
    % Mayo Clinic Foundation
    % Rochester, MN 55902, USA
    %
    % Email: cui.jie@maya.edu

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        rng_seed (1, 1) double {mustBeInteger, mustBePositive} = 42 % random number generator seed
    end % positional

    arguments
        options.BandHalfWidth (1, 1) double {mustBePositive} = .05 % half width of the band (Hz)
        options.Duration (1, 1) double {mustBePositive} = 50 % duration of the signal (s)
        options.FcMax (1, 1) double = nan % maximum center frequency of analysis bands (Hz)
        options.FcMin (1, 1) double = nan % minimum center frequency of analysis bands (Hz)
        options.GWNPower (1, 1) double = 0 % power of Gaussian white noise (dB)
        options.MaxFrequency (1, 1) double = 1/2 % maximum frequency of the signal (Hz)
        options.MinFrequency (1, 1) double = 0 % minimum frequency of the signal (Hz
        options.NumTimePoints (1, 1) double {mustBeInteger, mustBePositive} = 50 % number of time points
        options.NumTrials (1, 1) double {mustBeInteger, mustBePositive} = 1000 % number of trials
        options.NumFreqPoints (1, 1) double {mustBeInteger, mustBePositive} = 200 % number of frequency points
        options.TimePointsMethod (1, :) string ...
            {mustBeMember(options.TimePointsMethod, ...
             ["Uniform", "MissingPoints", "ArithmeticSampling", "Jittering"])} ...
            = ["Uniform", "MissingPoints", "ArithmeticSampling", "Jittering"] % method to generate time points
    end % optional

    f_w = options.BandHalfWidth;
    T = options.Duration;
    fc_max = options.FcMax;
    fc_min = options.FcMin;
    fmax = options.MaxFrequency;
    fmin = options.MinFrequency;
    num_points = options.NumTimePoints;
    num_trials = options.NumTrials;
    num_fc = options.NumFreqPoints;
    tp_method = options.TimePointsMethod;
    gw_std = 10 ^ (options.GWNPower / 20);

    % ======================================================================
    % main
    % ======================================================================
    % get time points and process samples
    % -----------------------------------
    s_old = rng(rng_seed); % save random number generator
    s_table = table();

    for tp_method_k = tp_method
        fprintf("Estimating speed for %s time points ...", tp_method_k)
        s_table_k = speed_ananlysis_tp(tp_method_k, rng_seed, num_points, ...
            num_trials, num_fc, fmax, T, f_w, gw_std);
        s_table = cat(1, s_table, s_table_k);
        fprintf(" done.\n")
    end % for

    % restore random number generator
    % -------------------------------
    rng(s_old)
end % function mtnu_speed_analysis

% ==========================================================================
% subroutines
% ==========================================================================
function s_table = speed_ananlysis_tp(tp_method, rng_seed, num_points, num_trials, num_fc, fc_min, fc_max, fmin, fmax, T, f_w, gw_std)

    arguments
        tp_method (1, :) string
        rng_seed (1, 1) double
        num_points (1, 1) double {mustBeInteger, mustBePositive}
        num_trials (1, 1) double {mustBeInteger, mustBePositive}
        num_fc (1, 1) double {mustBeInteger, mustBePositive}
        fc_min (1, 1) double
        fc_max (1, 1) double
        fmin (1, 1) double
        fmax (1, 1) double {mustBePositive}
        T (1, 1) double {mustBePositive}
        f_w (1, 1) double {mustBePositive}
        gw_std (1, 1) double {mustBePositive}
    end % positional

    % pramaeters
    % ----------
    if isnan(fc_min)
        fc_min = fmin + f_w;
    end % if

    if isnan(fc_max)
        fc_max = fmax - f_w;
    end % if

    % * sampling points of time instants
    bg = mtnusp.BronezGPSS();

    switch tp_method
        case "Uniform"
            t = bg.get_time_points(num_points, "DRS", ...
                TimeStart = 1, ...
                TimeEnd = num_points);
        case "MissingPoints"
            t = bg.get_time_points(0, "TB2");
            num_points = length(t);
        case "ArithmeticSampling"
            t = bg.get_time_points(num_points, "TB3", ...
                TimeStart = 1, ...
                TimeEnd = num_points);
        case "Jittering"
            t = bg.get_time_points(num_points, "RJS", ...
                JitterMu = 0, ...
                JitterSigma = 0.1, ...
                TimeStart = 1, ...
                TimeEnd = num_points);
    end % switch

    % convert to column vectors
    t = t(:);

    % * frequency points of interest
    fc = linspace(fc_min, fc_max, num_fc)';
    TW = T * f_w;
    K = 2 * TW - 1;

    % error analysis by Monte Carlo simulation
    % ----------------------------------------
    rng(rng_seed)

    % * process sample points
    RB = bg.gpssmat(0, fmax, t);
    L = chol(RB, 'lower');

    % * parameters for BronezGPSS
    A = [fc, ones(num_fc, 1) * .05]; % analysis bands
    B = [0, fmax]; % signal bands

    tc_mtnu = zeros(num_trials, 1); % time cost of MTNUFFT
    tc_mtls = zeros(num_trials, 1); % time cost of MTLS
    tc_bgfx = zeros(num_trials, 1); % time cost of BGFixed
    tc_bgad = zeros(num_trials, 1); % time cost of BGAdaptive

    parfor k = 1:num_trials
        % generate samples
        x_k = L * randn(num_points, 1) * gw_std; % zero mean with std gw_std
        nus_k = mtnusp.NUContinuous(x_k, t - t(1));

        % * MTNUFFT
        t0 = tic;
        nus_k.mtnuspectrum('QuerryFrequencies', fc, ...
            'Timehalfbandwidth', TW, ...
            'NumberTapers', K);
        t_int = toc(t0);
        tc_mtnu(k) = t_int; % time cost of MTNUFFT

        % * MTLS
        t0 = tic;
        nus_k.mtlsspectrum('Timehalfbandwidth', TW, 'InputFrequencies', fc);
        t_int = toc(t0);
        tc_mtls(k) = t_int; % time cost of MTLS

        % * BGFixed
        bgfx_k = mtnusp.BronezGPSS(t, x_k, ...
            T = T, ...
            SignalBand = B, ...
            AnalysisBand = A, ...
            SelectionMethod = 'auto', ...
            NumTapers = [K, 2 * K], ...
            lambdaFactor = K);
        t0 = tic;
        bgfx_k.spectrumgpss();
        t_int = toc(t0);
        tc_bgfx(k) = t_int; % time cost of BGFixed

        % * BGAdaptive
        bgad_k = mtnusp.BronezGPSS(t, x_k, ...
            T = T, ...
            SignalBand = B, ...
            AnalysisBand = A, ...
            SelectionMethod = 'adaptive', ...
            NumTapers = [K, 2 * K], ...
            lambdaFactor = -30);
        t0 = tic;
        bgad_k.spectrumgpss();
        t_int = toc(t0);
        tc_bgad(k) = t_int; % time cost of BGAdaptive
    end % for

    % build speed table
    % -----------------
    avg_mtnu = mean(1 ./ tc_mtnu); % mean number of trials per second
    std_mtnu = std(1 ./ tc_mtnu); % standard deviation of number of trials per second
    sem_mtnu = std_mtnu / sqrt(num_trials); % standard error of number of trials per second

    avg_mtls = mean(1 ./ tc_mtls); % mean number of trials per second
    std_mtls = std(1 ./ tc_mtls); % standard deviation of number of trials per second
    sem_mtls = std_mtls / sqrt(num_trials); % standard error of number of trials per second

    avg_bgfx = mean(1 ./ tc_bgfx); % mean number of trials per second
    std_bgfx = std(1 ./ tc_bgfx); % standard deviation of number of trials per second
    sem_bgfx = std_bgfx / sqrt(num_trials); % standard error of number of trials per second

    avg_bgad = mean(1 ./ tc_bgad); % mean number of trials per second
    std_bgad = std(1 ./ tc_bgad); % standard deviation of number of trials per second
    sem_bgad = std_bgad / sqrt(num_trials); % standard error of number of trials per second

    s_table = table(tp_method, avg_mtnu, std_mtnu, sem_mtnu, ...
        avg_mtls, std_mtls, sem_mtls, ...
        avg_bgfx, std_bgfx, sem_bgfx, ...
        avg_bgad, std_bgad, sem_bgad, ...
        'VariableNames', ["TimePointMethod", ...
           "MTNUFFT_mean", "MTNUFFT_std", "MTNUFFT_sem", ...
           "MTLS_mean", "MTLS_std", "MTLS_sem", ...
           "BGFixed_mean", "BGFixed_std", "BGFixed_sem", ...
           "BGAdaptive_mean", "BGAdaptive_std", "BGAdaptive_sem"]);

end % function

% [EOF]
