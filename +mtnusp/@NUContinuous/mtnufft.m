function [J, f, D] = mtnufft(this, x, t, options)
    % MTNUSP.NUCONTINUOUS.MTNUFFT Compute nonuniform FFT using multitaper method
    %
    % Syntax:
    %   [J, f, D] = mtnufft(this)
    %   [J, f, D] = mtnufft(this, x)
    %   [J, f, D] = mtnufft(this, x, t)
    %   [J, f, D] = mtnufft(this, x, t, Name, Value)
    %
    % Input(s):
    %   this        - mtnusp.NUContinuous object
    %   x           - input signal either vector or matrix n x channels
    %   t           - time vector of sampling instants
    %
    % Name-Value Pair(s):
    %   'Halfbandwidth'     - half bandwidth of analysis bands (Hz)
    %                         Default: 0.05
    %   'MaxFrequency'      - maximum frequency of signal x (Hz)
    %                         Default: 1/2
    %   'NumberTapers'      - number of tapers
    %                         Default: nan (use 2*TW-1)
    %   'NormMethod'        - normalization method
    %                         'L2Norm'  - normalize taper to unit L2 norm (default)
    %                         'BGNorm'  - normalize taper to be consistent with Bronez GPSS
    %   'TimeHalfbandwidth' - TW = time x halfbandwidth
    %                         Default: 0 (use T)
    %   'Tapers'            - tapers to use for multitaper method
    %                         Default: [] (use DPSS)
    %   'QuerryFrequencies' - frequencies of evaluation
    %                         Default: [] (use (0:N-1)/T)
    %
    % Output(s):
    %   J           - nonuniform FFT of signal x at frequencies f
    %   f           - frequencies of evaluation
    %   D           - taper sequence at time instants t
    %
    % Example:
    %   % generate nonuniformly sampled signal
    %   fs = 1; % sampling frequency (Hz)
    %   T = 100; % duration (seconds)
    %   N = fs * T; % number of samples
    %   t = sort(T * rand(N, 1)); % random sampling instants
    %   x = sin(2 * pi * .1 * t) + 0.5 * randn(N, 1); % signal samples
    %   nus = mtnusp.NUContinuous(x, t);
    %
    % Note:
    %
    % References:
    %
    %   Dutt, A., & Rokhlin, V. (1993). Fast Fourier Transforms for
    %        Nonequispaced Data. SIAM Journal on Scientific Computing,14(6),
    %        1368-1393. doi:10.1137/0914081
    %
    %   Springford, A., Eadie, G. M., & Thomson, D. J. (2020). Improving the
    %        Lomb–Scargle Periodogram with the Thomson Multitaper. The
    %        Astronomical Journal, 159(5), 205. doi:10.3847/1538-3881/ab7fa1
    %
    % See also .

    % 2022-2025 Richard J. Cui. Created: Sat 06/18/2022 10:27:31.613 PM
    % $Revision: 1.4 $  $Date: Tue 08/26/2025 11:59:41.015 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.jie.cui@gmail.com

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
        options.NumberTapers (1, 1) double = nan % number of tapers
        options.NormMethod (1, 1) string ...
            {mustBeMember(options.NormMethod, {'L2Norm', 'BGNorm'})} = 'L2Norm' % normalization method
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

    % compute taper sequence (default: DPSS)
    % --------------------------------------
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

    if isempty(tapers)
        tapers = dpss(N, TW, K);
        % ***********************************************************
        % * interpolate Slepian (DPSS) to time instants of sampling *
        % ***********************************************************
        % This part is the core idea. Need error analysis.
        v = (0:N - 1) * avg_d_t; % time vector of sampling instants of DPSS

        dps_x = zeros(size(tapers));

        for k = 1:K
            dps_k = interp1(v, tapers(:, k), t, 'spline');
            dps_x(:, k) = dps_k;
        end % for

        % normalization
        switch norm_method
            case 'L2Norm'

                for k = 1:K
                    dps_k = dps_x(:, k);
                    dps_x(:, k) = dps_k / norm(dps_k, 2); % normalize DPSS
                end % for

                dps_x = dps_x * sqrt(1 / avg_d_t); % make it consistent with Percival book (Chronux)

            case 'BGNorm'
                RB = mtnusp.BronezGPSS.gpssmat(0, fmax, t);
                dps_x = normalize_eigenvectors(dps_x, RB, f_w);
        end % switch

    else

        if isrow(tapers) == true
            dps_x = tapers.';
        else
            dps_x = tapers;
        end % if

    end % if

    assert(size(dps_x, 1) == N, ...
    'NUContinuous.mtnufft: Taper must have the same length as t')
    assert(size(dps_x, 2) == K, ...
    'NUContinuous.mtnufft: Taper must have the same number of tapers as K')

    % compute multitaper spectrum of nonuniformly sampled data
    % --------------------------------------------------------
    num_chans = size(x, 2);
    D = dps_x(:, :, ones(1, num_chans));
    Y = x(:, :, ones(1, K));
    Y = permute(Y, [1, 3, 2]);
    X = Y .* D;
    J = nufft(X, t, fvec);

    switch norm_method
        case 'L2Norm'
            J = J * avg_d_t; % make consistent with Percival book (Chronux)
        case 'BGNorm'
            J = J / sqrt(2 * f_w);
    end % switch

end % function mtnufft

% ==========================================================================
% subroutines
% ==========================================================================
function V_nm = normalize_eigenvectors(V, R, fw)

    arguments
        V (:, :) double
        R (:, :) double
        fw (1, 1) double % analysis hald-band withd (hz)
    end % positional

    % calculate normalization matrix
    % ------------------------------
    V_nm = zeros(size(V));
    num_vec = size(V, 2);

    for k = 1:num_vec
        w_k = V(:, k);
        r_k = w_k' * R * w_k;
        s_k = sqrt(2 * fw / r_k);
        V_nm(:, k) = s_k * w_k;
    end % for k

end % function

% [EOF]
