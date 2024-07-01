function [J, f, Mpt, Npt] = mtnufftpb(this, y, t, options)
    % NUPOINTBINNED.MTNUFFTPB calculate NUFFT of binned point process
    %
    % Syntax:
    %
    % Input(s):
    %
    % Output(s):
    % 	J 	        - [array] FFT of nonuniformly sampled data in form f x K x C,
    %		          where 'f' is the frequency index, 'K' is the number of
    % 		          tapers and 'C' is the number of channels/trials.
    %   f 	        - [array] frequencies of interest
    %   Mpt         - [array] mean point rate per sample (i.e. number of points
    %                 per sample) in each channel - 1 x C
    %   Npt 	    - [array] number of points in each channel - 1 x C
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Richard J. Cui. Created: Mon 07/11/2022 11:22:09.279 PM
    % $Revision: 0.2 $  $Date: Sun 07/17/2022  2:38:22.499 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) mtnusp.NUPointbinned
        y (:, :) double % binned point process - samples x channels/trials
        t (1, :) double = double.empty(1, 0) % time vector of sampling instants
    end % positional

    arguments
        options.TimeHalfbandwidth (1, 1) double = 0 % TW = time x halfbandwidth
        options.QuerryFrequencies (1, :) double = double.empty(1, 0) % frequencies of evaluation
        options.Tapers (:, :) double = double.empty(1, 0) % tapers to use for multitaper method
        options.NumberTapers (1, 1) double = nan % number of tapers
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

    if isempty(t)
        t = this.SamplingTimes;
        assert(~isempty(t), 'NUPointbinned.mtnufftpb: t is empty')
    end % if

    % normalize t
    t = t - t(1);

    assert(size(y, 1) == numel(t), 'NUSpectrum.mtnufftpb: y and t must have the same length')

    % * options
    TW = options.TimeHalfbandwidth;
    fvec = options.QuerryFrequencies;
    tapers = options.Tapers;
    num_taper = options.NumberTapers;

    % point process
    % -------------
    [N, C] = size(y); % size of binned data
    Npt = sum(y, 1); % size of data - 1 x C
    Mpt = Npt / N; % mean rate for each channel - 1 x C

    % compute taper sequence (default: DPSS)
    % --------------------------------------
    % * average sampling interval
    D = t(end) - t(1);
    avg_d_t = D / (N - 1); % consistent with uniformly sampled signal
    T = avg_d_t * N;

    % * frequency points of interest
    if isempty(fvec)
        fvec = (0:N - 1) / T;
    end % if

    f = fvec;
    num_fvec = length(fvec);

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
    'NUSpectrum.mtnufftpb: Taper must have the same length as samples')

    % FFT of tapers
    tapers = tapers(:, :, ones(1, C)); % add channel/trial indices
    H = nufft(tapers, (0:N - 1) * avg_d_t, fvec); % - f x K x C

    % * interpolate Slepian (DPSS) to time instants of sampling
    v = (0:N - 1) * avg_d_t; % time vector of sampling instants of DPSS

    dps_x = zeros(N, K);

    for k = 1:K
        dps_k = interp1(v, tapers(:, k), t, 'spline');
        dps_x(:, k) = dps_k / norm(dps_k, 2); % normalize DPSS
    end % for

    dps_x = dps_x(:, :, ones(1, C)); % add channel/trial indices

    % nonuniform Fourier transform of binned point process
    % ----------------------------------------------------
    y = y(:, :, ones(1, K)); % add taper indices to the data
    y = permute(y, [1 3 2]); % permute the data to be of the same dimensions as dps_x - sample x K x C
    y_proj = y .* dps_x; % multiply data by the tapers
    % * nonuniform FFT of the projected data
    J = nufft(y_proj, t, fvec);
    % * crucial to subtract the mean rate
    M = Mpt(:); % column vector
    meanpt = M(:, ones(1, K), ones(1, num_fvec)); % add taper and frequency indices
    meanpt = permute(meanpt, [3 2 1]); % permute to get the same dimensions as H - f x K x C
    J = J - H .* meanpt; % subtract the dc

end % function mtnufftpb

% [EOF]
