function e_table = mtnu_error_analysis(rng_seed, options)
    % MTNU_ERROR_ANALYSIS error analysis of MTNUFFT method
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

    % Copyright 2024 Richard J. Cui. Created: Tue 07/02/2024 10:37:03.476 AM
    % $Revision: 0.2 $  $Date: Sun 07/07/2024 12:15:57.004 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

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

    fw = options.BandHalfWidth;
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

    for tp_method_k = tp_method
        fprintf("Processing %s time points ...", tp_method_k)
        e_table.(tp_method_k) = err_ananlysis_tp(tp_method_k, rng_seed, num_points, ...
            num_trials, num_fc, fc_max, fc_min, fmin, fmax, T, fw, gw_std);
        fprintf(" done.\n")
    end % for

    % restore random number generator
    % -------------------------------
    rng(s_old)
end % function mtnu_error_analysis

% ==========================================================================
% subroutines
% ==========================================================================
function e_table = err_ananlysis_tp(tp_method, rng_seed, num_points, num_trials, num_fc, fc_min, fc_max, fmin, fmax, T, fw, gw_std)

    arguments
        tp_method (1, :) string
        rng_seed (1, 1) double
        num_points (1, 1) double {mustBeInteger, mustBePositive}
        num_trials (1, 1) double {mustBeInteger, mustBePositive}
        num_fc (1, 1) double {mustBeInteger, mustBePositive}
        fc_min (1, 1) double
        fc_max (1, 1) double
        fmin (1, 1) double
        fmax (1, 1) double
        T (1, 1) double {mustBePositive}
        fw (1, 1) double {mustBePositive}
        gw_std (1, 1) double {mustBePositive}
    end % positional

    % pramaeters
    % ----------
    if isnan(fc_min)
        fc_min = fmin + fw;
    end % if

    if isnan(fc_max)
        fc_max = fmax - fw;
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
    TW = T * fw;
    K = 2 * TW - 1;

    % * grand truth of spectrum
    S_fc = gw_std ^ 2 * ones(num_fc, 1); % GWN spectrum

    % error analysis by Monte Carlo simulation
    % ----------------------------------------
    rng(rng_seed)

    % * process sample points
    RB = bg.gpssmat(fmin, fmax, t);
    L = chol(RB, 'lower');

    % * parameters for BronezGPSS
    A = [fc, ones(num_fc, 1) * .05]; % analysis bands
    B = [fmin, fmax]; % signal bands

    err_mtnu = zeros(num_fc, num_trials); % error of MTNUFFT
    err_mtls = zeros(num_fc, num_trials); % error of MTLS
    err_bgfx = zeros(num_fc, num_trials); % error of BGFixed
    err_bgad = zeros(num_fc, num_trials); % error of BGAdaptive

    parfor k = 1:num_trials
        % generate samples
        x_k = L * randn(num_points, 1) * gw_std; % zero mean with std gw_std
        nus_k = mtnusp.NUContinuous(x_k, t - t(1));

        % * MTNUFFT
        pxx_mtnu_k = nus_k.mtnuspectrum('QuerryFrequencies', fc, ...
            'Timehalfbandwidth', TW, ...
            'NumberTapers', K);
        err_mtnu(:, k) = abs(pow2db(pxx_mtnu_k ./ S_fc)) .^ 2; % squared error

        % * MTLS
        pxx_mtls_k = nus_k.pmtlomb('Timehalfbandwidth', TW, 'InputFrequencies', fc);
        pxx_mtls_k = pxx_mtls_k * num_points; % scale to FFT definition of periodogram
        pxx_mtls_k = mean(pxx_mtls_k(:, :, 1), 2); % average over the number of tapers
        err_mtls(:, k) = abs(pow2db(pxx_mtls_k ./ S_fc)) .^ 2; % squared error

        % * BGFixed
        bgfx_k = mtnusp.BronezGPSS(t, x_k, ...
            T = T, ...
            SignalBand = B, ...
            AnalysisBand = A, ...
            SelectionMethod = 'auto', ...
            NumTapers = [K, 2 * K], ...
            lambdaFactor = K);
        pxx_bgfx_k = bgfx_k.spectrumgpss();
        err_bgfx(:, k) = abs(pow2db(pxx_bgfx_k(:) ./ S_fc)) .^ 2; % squared error

        % * BGAdaptive
        bgad_k = mtnusp.BronezGPSS(t, x_k, ...
            T = T, ...
            SignalBand = B, ...
            AnalysisBand = A, ...
            SelectionMethod = 'adaptive', ...
            NumTapers = [K, 2 * K], ...
            lambdaFactor = -30);
        pxx_bgad_k = bgad_k.spectrumgpss();
        err_bgad(:, k) = abs(pow2db(pxx_bgad_k(:) ./ S_fc)) .^ 2; % squared error
    end % for

    % build error table
    % -----------------
    avg_mtnu = mean(err_mtnu, 2); % mean squared error
    std_mtnu = std(err_mtnu, 0, 2); % standard deviation
    sem_mtnu = std_mtnu / sqrt(num_trials); % standard error of the mean

    avg_mtls = mean(err_mtls, 2); % mean squared error
    std_mtls = std(err_mtls, 0, 2); % standard deviation
    sem_mtls = std_mtls / sqrt(num_trials); % standard error of the mean

    avg_bgfx = mean(err_bgfx, 2); % mean squared error
    std_bgfx = std(err_bgfx, 0, 2); % standard deviation
    sem_bgfx = std_bgfx / sqrt(num_trials); % standard error of the mean

    avg_bgad = mean(err_bgad, 2); % mean squared error
    std_bgad = std(err_bgad, 0, 2); % standard deviation
    sem_bgad = std_bgad / sqrt(num_trials); % standard error of the mean

    e_table = table(fc, avg_mtnu, std_mtnu, sem_mtnu, ...
        avg_mtls, std_mtls, sem_mtls, ...
        avg_bgfx, std_bgfx, sem_bgfx, ...
        avg_bgad, std_bgad, sem_bgad, ...
        'VariableNames', ["Frequency", ...
           "MTNUFFT_mean", "MTNUFFT_std", "MTNUFFT_sem", ...
           "MTLS_mean", "MTLS_std", "MTLS_sem", ...
           "BGFixed_mean", "BGFixed_std", "BGFixed_sem", ...
           "BGAdaptive_mean", "BGAdaptive_std", "BGAdaptive_sem"]);

end % function

% [EOF]
