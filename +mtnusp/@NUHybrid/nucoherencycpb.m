function [C, phi, S12, S1, S2, f, zerosp, confC, phistd, Cerr] = nucoherencycpb(this, y, x, t, options)
    % NUSPECTRUM.NUCOHERENCYCPB multitaper coherency between nonuniformly sampled continuous and binned point process
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
    %   Need Chronux software toolbox.
    %
    % References:
    %
    % See also .

    % Richard J. Cui. Created: Wed 07/06/2022 10:31:30.356 PM
    % $Revision: 0.2 $  $Date: Thu 07/14/2022 11:17:36.027 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) NUSpectrum
        y (:, :) double % binned point process - samples x channels/trials
        x (:, :) double = [] % input nonuniform signal either vector or matrix - samples x channels/trials
        t (1, :) double = [] % time vector of nonuniformly sampling instants
    end % positional

    arguments
        options.TimeHalfbandwidth (1, 1) double = 0 % TW = time x halfbandwidth
        options.QuerryFrequencies (1, :) double = double.empty(1, 0) % frequencies of evaluation
        options.Tapers (:, :) double = double.empty(1, 0) % tapers to use for multitaper method
        options.NumberTapers (1, 1) double = nan % number of tapers
        options.TrialAverage (1, 1) logical = false % average across trials/channels
        options.FiniteSizeCorrection (1, 1) logical = false % wheter use finite size correction
        options.ErrorType (1, :) double = 0 % error calculation  - optional (Default 0).
        % [1 p] - Theoretical error bars;
        % [2 p] - Jackknife error bars;
        % [0 p] or 0 - no error bras
        % p     - p-value for estimation (e.g., 0.05)
    end % optional

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    % * required
    if isvector(y)
        y = y(:);
    end % if

    if isempty(x)
        x = this.Signal;
        assert(~isempty(x), 'NUSpectrum.nucoherencycpb: x is empty')
    end % if

    if isvector(x)
        x = x(:);
    end % if

    if isempty(t)
        t = this.SamplingTimes;
        assert(~isempty(t), 'NUSpectrum.nucoherencycpb: t is empty')
    end % if

    % normalize t
    t = t - t(1);

    assert(size(x, 1) == numel(t), 'NUSpectrum.nucoherencycpb: x and t must have the same length')
    assert(size(x, 1) == size(y, 1), 'NUSpectrum.nucoherencycpb: x and y must have the same number of samples')
    assert(size(x, 2) == size(y, 2), 'NUSpectrum.nucoherencycpb: x and y must have the same number of trials/channels')

    % * options
    TW = options.TimeHalfbandwidth;
    fvec = options.QuerryFrequencies;
    tapers = options.Tapers;
    num_taper = options.NumberTapers;
    trialave = options.TrialAverage;
    fscorr = options.FiniteSizeCorrection;
    err = options.ErrorType;

    % compute taper sequence (default: DPSS)
    % --------------------------------------
    % * average sampling interval
    [N, Ch] = size(y);
    D = t(end) - t(1);
    avg_d_t = D / (N - 1); % consistent with uniformly sampled signal
    T = avg_d_t * N;

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

    if isempty(tapers)
        tapers = dpss(N, TW, K);
    elseif isrow(tapers) == true
        tapers = tapers.';
    end % if

    assert(size(tapers, 1) == N, ...
    'NUSpectrum.nucoherencycpb: Taper must have the same length as t')

    % compute multitaper coherency of nonuniformly sampled and binned point process data
    % ----------------------------------------------------------------------------------
    % * Fourier transform of nonuniformly sampled continous signal
    J1 = this.mtnufft(x, t, QuerryFrequencies = fvec, TimeHalfbandwidth = TW, ...
    Tapers = tapers);

    % * FFT of nonuniformly binned point process
    [J2, ~, ~, Npt] = this.mtnufftpb(y, QuerryFrequencies = fvec, TimeHalfbandwidth = TW, ...
    Tapers = tapers);

    zerosp = zeros(1, Ch);
    zerosp(Npt == 0) = true; % set the channles/trials where no points were found to be TRUE

    % * power spectra
    S1 = squeeze(mean(conj(J1) .* J1, 2));
    S2 = squeeze(mean(conj(J2) .* J2, 2));
    S12 = squeeze(mean(conj(J1) .* J2, 2));

    if trialave == true
        S1 = squeeze(mean(S1, 2));
        S2 = squeeze(mean(S2, 2));
        S12 = squeeze(mean(S12, 2));
    end

    % * coherency
    C12 = S12 ./ sqrt(S1 .* S2);
    C = abs(C12);
    phi = angle(C12);

    % outpu
    % -----
    if nargout == 10

        if fscorr == true
            [confC, phistd, Cerr] = coherr(C, J1, J2, err, trialave, [], Npt);
        else
            [confC, phistd, Cerr] = coherr(C, J1, J2, err, trialave);
        end

    elseif nargout == 9

        if fscorr == true
            [confC, phistd] = coherr(C, J1, J2, err, trialave, [], Npt);
        else
            [confC, phistd] = coherr(C, J1, J2, err, trialave);
        end

    end

end % function nucoherencycpb

% [EOF]
