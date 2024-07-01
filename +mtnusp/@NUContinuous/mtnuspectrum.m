function [S, f] = mtnuspectrum(this, x, t, options)
    % NUCONTINUOUS.MTNUSPECTRUM spectral estimation using MTNUFFT method
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

    % Copyright 2023 Richard J. Cui. Created: Wed 10/11/2023 10:27:08.208 PM
    % $Revision: 0.4 $  $Date: Wed 06/26/2024 12:21:02.181 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) mtnusp.NUContinuous
        x (:, :) double = double.empty(1, 0) % input signal either vector or matrix n x channels
        t (1, :) double = double.empty(1, 0) % time vector of sampling instants
    end % positional

    arguments
        options.Halfbandwidth (1, 1) double = .05 % half bandwidth of analysis bands (Hz)
        options.MaxFrequency (1, 1) double = 1/2 % maximum frequency of signal x (Hz)
        options.NormMethod (1, 1) string ...
            {mustBeMember(options.NormMethod, {'L2Norm', 'BGNorm'})} = 'L2Norm' % normalization method
        options.NumberTapers (1, 1) double = nan % number of tapers
        options.TimeHalfbandwidth (1, 1) double = 0 % TW = time x halfbandwidth
        options.Tapers (:, :) double = double.empty(1, 0) % tapers to use for multitaper method
        options.QuerryFrequencies (1, :) double = double.empty(1, 0) % frequencies of evaluation
    end % optional

    f_w = options.Halfbandwidth;
    fmax = options.MaxFrequency;
    assert(fmax > 0, 'NUContinuous.mtnufft: fmax must be positive')
    norm_method = options.NormMethod;

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    if isempty(x)
        x = this.Signal;
        assert(~isempty(x), 'NUContinuous.mtnufft: x is empty')
    end % if

    if isvector(x)
        x = x(:);
    end % if

    if isempty(t)
        t = this.SamplingTimes;
        assert(~isempty(t), 'NUContinuous.mtnufft: t is empty')
    end % if

    % normalize t
    t = t - t(1);

    assert(size(x, 1) == numel(t), 'NUContinuous.mtnufft: x and t must have the same length')

    % * options
    TW = options.TimeHalfbandwidth;
    fvec = options.QuerryFrequencies;
    tapers = options.Tapers;
    num_taper = options.NumberTapers;

    % * average sampling interval
    N = length(t);
    D = t(end) - t(1);
    avg_d_t = D / (N - 1); % consistent with uniformly sampled signal
    T = avg_d_t * N;

    % * frequency points of interest
    if isempty(fvec)
        fvec = (0:N - 1) / T;
    end % if

    f = fvec;

    % * calculate uniform Slepian (DPSS)
    if TW == 0
        TW = T;
    end % if

    if isnan(num_taper)
        K = floor(2 * TW - 1);
    else
        K = num_taper;
    end % if

    % estimate the power spectrum
    % ---------------------------
    J = this.mtnufft('QuerryFrequencies', f, ...
        'MaxFrequency', fmax, ...
        'NormMethod', norm_method, ...
        'Halfbandwidth', f_w, ...
        'TimeHalfbandwidth', TW, ...
        'Tapers', tapers, ...
        'NumberTapers', K); % J = f x tapers x channels
    S = mean(abs(J) .^ 2, 2); % average over tapers and channels
    S = squeeze(S); % remove singleton dimension

end % function mtnuspectrum

% [EOF]
