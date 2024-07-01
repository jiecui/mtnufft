function [pxx, f, x, t] = pmtlomb(this, x, t, options)
    % MTNUSP.NUCONTINUOUS.PMTLOMB calculates the Lomb-Scargle periodogram using multitaper method
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
    %   If both InputFrequencies and OversamplingFactor are specified, the
    %   InputFrequencies will be used.
    %
    % References:
    %
    %   Babu, P., & Stoica, P. (2010). Spectral analysis of nonuniformly
    %        sampled data – a review. Digital Signal Processing, 20(2),
    %        359-378. doi:10.1016/j.dsp.2009.06.019
    %
    %   Springford, A., Eadie, G. M., & Thomson, D. J. (2020). Improving the
    %        Lomb–Scargle Periodogram with the Thomson Multitaper. The
    %        Astronomical Journal, 159(5), 205. doi:10.3847/1538-3881/ab7fa1
    %
    % See also .

    % 2022 Richard J. Cui. Created: Sun 04/03/2022  2:00:53.607 AM
    % $Revision: 0.7 $  $Date: Sun 12/24/2023 10:35:42.120 AM $
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
        options.InputFrequencies (1, :) double = nan % input frequencies to be considered
        options.MaxFrequency (1, 1) double = 0 % maximum frequency to be considered
        options.OversamplingFactor (1, 1) double = 0 % oversampling factor
        options.TimeHalfbandwidth (1, 1) double = 0 % TW = time x halfbandwidth
        options.Tapers (:, :) double = double.empty(1, 0) % tapers to use for multitaper method
    end % optional

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    if isempty(x)
        x = this.Signal;
        assert(~isempty(x), 'NUContinuous.pmtlomb: x is empty')
    end % if

    if isempty(t)
        t = this.SamplingTimes;
        assert(~isempty(t), 'NUContinuous.pmtlomb: t is empty')
    end % if

    assert(numel(x) == numel(t), 'NUContinuous.pmtlomb: x and t must have the same length')

    % * options
    fvec = options.InputFrequencies;
    TW = options.TimeHalfbandwidth;
    taper_seq = options.Tapers;
    fmax = options.MaxFrequency;
    ofact = options.OversamplingFactor;

    % compute Slepian (DPSS)
    % ----------------------
    % * average sampling interval
    N = length(t);
    T = t(end) - t(1);
    avg_d_t = T / N;

    % * calculate uniform Slepian (DPSS)
    if TW == 0
        TW = T;
    end % if

    num_seq = floor(2 * TW - 1);

    if isempty(taper_seq)
        taper_seq = dpss(N, TW, num_seq);
    elseif isrow(taper_seq) == true
        taper_seq = taper_seq.';
    end % if

    assert(size(taper_seq, 1) == N, ...
    'NUContinuous.pmtlomb: Taper must have the same length as t')

    % * interpolate Slepian (DPSS) to time instants of sampling
    v = (0:N - 1) * avg_d_t; % time vector of sampling instants of DPSS

    dps_x = zeros(size(taper_seq));

    for k = 1:num_seq
        dps_k = interp1(v, taper_seq(:, k), t, 'spline');
        dps_x(:, k) = dps_k / norm(dps_k, 2); % normalize DPSS
    end % for

    % compute multitaper lomb-scargle periodogram
    % -------------------------------------------
    if isvector(x) == true
        x = x(:);
    end % if

    num_chans = size(x, 2);

    % get frequency vector to be returned
    if isnan(fvec) == true
        [~, f] = plomb(x(:, 1), t, fmax, ofact);
    else
        f = fvec;
    end % if

    num_f = length(f);
    pxx = zeros(num_f, num_seq, num_chans);
    % channel by channel
    for k = 1:num_chans
        x_k = x(:, k);
        X_k = x_k * ones(1, num_seq) .* dps_x; % multiply DPSS with input signal

        % compute LS periodogram of channel k
        if isnan(fvec) == true
            pxx_k = plomb(X_k, t, fmax, ofact); % compute periodogram of each tapered sequence
        else
            pxx_k = plomb(X_k, t, fvec); % compute periodogram of each tapered sequence
        end % if

        pxx(:, :, k) = pxx_k; % store average spectrum
    end % for

    if isreal(x) == true
        pxx = pxx / 2; % for real signal only half of the power is in the positive frequencies
    end % if

end % function pmtlomb

% [EOF]
