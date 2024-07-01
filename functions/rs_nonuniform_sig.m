function [y, ty, b] = rs_nonuniform_sig(x, tx, fs, p, q, b, options)
    % RS_NONUNIFORM_SIG resampling nonuniform signals
    %
    % Syntax:
    %   [y, ty, b] = rs_nonuniform_sig(x, tx, fs, p, q, b, options) 
    %
    % Input(s):
    %   x               - [array] input signal: samples x channels
    %   tx              - [array] timestamps of x
    %   fs              - [double] (optional) desired sampling frequency 
    %                     (samples/unit time, time unit consistent with tx)
    %   p             
    %
    % Output(s):
    %   y               - [double] resampled signal
    %   ty              - [double] time stamps of resampled signal
    %   b               - [double] FIR filter coefficients
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % 2021 Richard J. Cui. Created: Thu 12/23/2021 11:20:21.260 PM
    % $Revision: 0.3 $  $Date: Wed 10/19/2022 12:58:59.465 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        x (:, :) double % input signal
        tx (1, :) double % time instants
        fs (1, 1) double {mustBeNonnegative(fs)} = 4 % desired frequency (unit time)
        p (1, 1) double {mustBeNonnegative(p)} = 0 % numerator of resampling factor
        q (1, 1) double {mustBeNonnegative(q)} = 0 % donominator of resampling factor
        b (1, :) double = [] % FIR filter coefficients
    end % positional

    arguments
        options.AssumedNyquist (1, 1) {mustBePositive(options.AssumedNyquist)} = 1 % cycles per unit time
        options.NorminalSampleRate (1, 1) {mustBeNonnegative(options.NorminalSampleRate)} = 0
        options.Alpha (1, 1) double {mustBeNonnegative} = 4 % fMN >= alpha * fs
        options.UseLPFilter (1, 1) logical = false % whether to use low pass filter
        options.CutoffRatio (1, 1) double {mustBeInRange(options.CutoffRatio, 0, 1)} = .25
        options.Beta (1, 1) double {mustBePositive(options.Beta)} = 5 % LFP parameter
        options.Method (1, 1) string {mustBeMember(options.Method, ["linear", "pchip", "spline"])} = "linear"
    end % optional

    % ======================================================================
    % main
    % ======================================================================
    fN = options.AssumedNyquist;
    fNM = options.NorminalSampleRate;
    alpha = options.Alpha;
    use_lpf = options.UseLPFilter;
    cutoff_ratio = options.CutoffRatio;
    beta = options.Beta;
    method = options.Method;

    tx = tx(:);

    % set parameters
    % --------------
    % estimate sample average sample rate
    Tx = mean(diff(tx));
    fx = 1 / Tx;

    if fNM == 0
        fNM = fx;
    end % if

    % check fN
    if fNM < 2 * fN
        cprintf([1, .5, 0], 'Norminal frequency %.4f may be too low for Nyquist frequency %.4f\n', fNM, fN)
    end % if

    % check alpha
    if fNM < alpha * fs
        cprintf([1, .5, 0], 'Desired frequency %.4f may be too high\n', fs)
    end % if

    % set resampling factor
    % ---------------------
    if p * q == 0
        [p, q] = rat(fs / fNM);
    end % if

    % set LPF
    % -------
    if use_lpf == true

        if isempty(b) == true
            % ensure an odd length filter
            n = round(beta * max(p, q));

            if rem(n, 2) == 0
                n = n + 1;
            end % if

            % use .25 of Nyquist range of desired sample rate
            % construct LPF
            b = p * fir1(n, cutoff_ratio / q);
        end % if

    else
        b = [];
    end % if

    % resample signals one by one
    % ---------------------------
    if isvector(x)
        x = x(:);
    end % if

    num_signals = size(x, 2);
    y = [];

    for k = 1:num_signals
        [y_k, ty] = resample_one_signal(x(:, k), tx, fs, p, q, b, method);
        y = cat(2, y, y_k);
    end % for

end % function rs_nonuniform_sig

% ==========================================================================
% subroutines
% ==========================================================================
function [y, ty] = resample_one_signal(x, tx, fs, p, q, b, method)
    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        x (:, 1) double % input signal
        tx (:, 1) double % time instants
        fs (1, 1) double {mustBeNonnegative(fs)} = 0 % desired frequency (unit time)
        p (1, 1) double {mustBeNonnegative(p)} = 0 % numerator of resampling factor
        q (1, 1) double {mustBeNonnegative(q)} = 0 % donominator of resampling factor
        b (1, :) double = [] % FIR filter coefficients
        method (1, 1) string = "linear";
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    % detrend the signal
    % ------------------
    % compute slope and offset (y = a1 x + a2)
    a(1) = (x(end) - x(1)) / (tx(end) - tx(1));
    a(2) = x(1);

    xdetrend = x - polyval(a, tx);

    % resample the signal
    % -------------------
    if isempty(b)
        [ydetrend, ty] = resample(xdetrend, tx, fs, p, q, method);
    else
        [ydetrend, ty] = resample(xdetrend, tx, fs, p, q, b, method);
    end % if

    % add back the trend
    % ------------------
    y = ydetrend + polyval(a, ty);

end % function

% [EOF]
