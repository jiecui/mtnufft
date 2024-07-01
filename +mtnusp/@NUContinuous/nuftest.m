function [Fval, A, f, sig, sd] = nuftest(this, x, t, options)
    % MTNUSP.NUCONTINUOUS.NUFTEST compute F-test of nonuniformly sampled signals using Thomson's multitaper method
    %
    % Syntax:
    %
    % Input(s):
    %
    % Output(s):
    %   Fval        - (F-statistic in frequency x channels/trials format)
    %   A           - (Line amplitude for X in frequency x channels/trials form)
    %   f           - (frequencies of evaluation)
    %   sig         - (F distribution (1-p)% confidence level)
    %   sd          - (standard deviation of the amplitude A)
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    %   Percival, D. B., & Walden, A. T. (1993). Spectral analysis for
    %       physical applications : multitaper and conventional univeriate
    %       techniques. Cambridge; New York, NY, USA: Cambridge University
    %       Press.
    %
    % See also ftestc (chronux).

    % 2022 Richard J. Cui. Created: Sun 06/19/2022  5:00:37.191 PM
    % $Revision: 0.6 $  $Date: Thu 10/12/2023 10:01:53.716 AM $
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
        options.TimeHalfbandwidth (1, 1) double = 0 % TW = time x halfbandwidth
        options.QuerryFrequencies (1, :) double = double.empty(1, 0) % frequencies of evaluation
        options.Tapers (:, :) double = double.empty(1, 0) % tapers to use for multitaper method
        options.NumberTapers (1, 1) double = nan % number of tapers (cannot be set with Tapers together)
        options.PValue (1, 1) double = 0.05
    end % optional

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    if isempty(x)
        x = this.Signal;
        assert(~isempty(x), 'NUContinuous.nuftest: x is empty')
    end % if

    if isvector(x)
        x = x(:);
    end % if

    if isempty(t)
        t = this.SamplingTimes;
        assert(~isempty(t), 'NUContinuous.nuftest: t is empty')
    end % if

    % normalize t
    t = t - t(1);

    assert(size(x, 1) == numel(t), 'NUContinuous.nuftest: x and t must have the same length')

    % * options
    TW = options.TimeHalfbandwidth;
    fvec = options.QuerryFrequencies;
    tapers = options.Tapers;
    K = options.NumberTapers;
    p = options.PValue;

    % compute taper sequence (default: DPSS)
    % --------------------------------------
    % * average sampling interval
    [N, C] = size(x);
    D = t(end) - t(1);
    avg_d_t = D / (N - 1); % consistent with uniformly sampled signal
    T = avg_d_t * N;
    assert(length(t) == N, ...
    'NUContinuous.nuftest: Taper must have the same length as t')

    % * calculate uniform Slepian (DPSS)
    if TW == 0
        TW = T;
    end % if

    if isnan(K)
        K = floor(2 * TW - 1);
    end % if

    if isempty(tapers)
        tapers = dpss(N, TW, K); % DPSS - t x K
    elseif isrow(tapers) == true
        tapers = tapers.';
        K = 1;
    else
        K = size(tapers, 2);
    end % if

    if isempty(fvec) == true
        fvec = (0:N - 1) / T;
    end % if

    % calculate F-test
    % ----------------
    Kodd = 1:2:K; % even indexed taper (start from 0, see Percival's book)
    Keven = 2:2:K;
    [J, f] = this.mtnufft(x, t, QuerryFrequencies = fvec, ...
        TimeHalfbandwidth = TW, ...
        NumberTapers = K, ...
        Tapers = tapers);
    Jp = J(:, Kodd, :); % drop the even nuffts
    tapers = tapers(:, :, ones(1, C)); % add channel indices to the tapers - t x K x C
    H0 = squeeze(sum(tapers(:, Kodd, :), 1)); % calculate sum of tapers for even prolates - K x C

    if C == 1
        H0 = H0';
    end

    Nf = length(fvec); % number of frequencies (no selection yet)
    H0 = H0(:, :, ones(1, Nf)); % add frequency indices to H0 - K x C x f
    H0 = permute(H0, [3 1 2]); % permute H0 to get dimensions to match those of Jp - f x K x C
    H0sq = sum(H0 .* H0, 2); % sum of squares of H0^2 across taper indices - f x C
    JpH0 = sum(Jp .* squeeze(H0), 2); % sum of the product of Jp and H0 across taper indices - f x C
    A = squeeze(JpH0 ./ H0sq); % amplitudes for all frequencies and channels
    Kp = size(Jp, 2); % number of even prolates
    Ap = A(:, :, ones(1, Kp)); % add the taper index to C
    Ap = permute(Ap, [1 3 2]); % permute indices to match those of H0
    Jhat = Ap .* H0; % fitted value for the fft

    % F-statistic
    num = (K - 1) .* (abs(A) .^ 2) .* squeeze(H0sq); %numerator for F-statistic
    den = squeeze(sum(abs(Jp - Jhat) .^ 2, 2) + sum(abs(J(:, Keven, :)) .^ 2, 2)); % denominator for F-statistic
    Fval = num ./ den; % F-statisitic
    sig = finv(1 - p, 2, 2 * K - 2); % F-distribution (1-p) percential
    var = den ./ (K * squeeze(H0sq)); % variance of amplitude
    sd = sqrt(var); % standard deviation of amplitude

    % A = A*Fs;
end % function nuftest

% [EOF]
