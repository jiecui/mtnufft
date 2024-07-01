function [S, f, phi_square] = jitter_spectrum_model(len, lambda, sigma, Fs, options)
    % MTNUSP.NUPOINT.JITTER_SPECTRUM_MODEL power spectrum model for jittering
    %
    % Syntax:
    %   [S, f, phi_square] = jitter_spectrum_model(len, lambda, sigma, Fs, options)
    %
    % Input(s):
    %
    % Output(s):
    %
    % Example:
    %
    % Note:
    %   Currently, we assume that the jittering is a Gaussian process.
    %
    % References:
    %   [1] Brémaud, P. (2014). Fourier analysis of stochastic processes. In
    %       Fourier analysis and stochastic processes (pp. 119-179). Cham:
    %       Springer International Publishing (Jittering pp. 311-313).
    %
    % See also .

    % Copyright 2024 Richard J. Cui. Created: Tue 02/13/2024  1:24:54.548 PM
    % $Revision: 0.1 $  $Date: Tue 02/13/2024  1:24:54.573 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        len (1, 1) double {mustBeNumeric, mustBeNonempty} % signal length
        lambda (1, 1) double {mustBeNumeric, mustBeNonempty} = 1 % point rate (number / second)
        sigma (1, 1) double = .1 % standard deviation of jittering
        Fs (1, 1) double {mustBeNumeric, mustBeNonempty} = 1 % sampling frequency (Hz)
    end % positional

    arguments
        options.FPass (1, 2) double = [0, Fs / 2] % frequency range
        options.Frequency (1, :) double = double.empty(1, 0) % frequency of interest
        options.NumberTapers (1, 1) double = 6
        options.TimeHalfbandwidth (1, 1) double = 3.5
    end % optional

    fpass = options.FPass;
    fmax = fpass(2);
    f = options.Frequency;

    if isempty(f)
        f = 0:lambda / Fs:fmax;
    end % if

    TW = options.TimeHalfbandwidth;
    K = options.NumberTapers;

    % ======================================================================
    % main
    % ======================================================================
    % estimate the scaling factor for chronux spectrum function
    % ---------------------------------------------------------
    ts = 0:1 / lambda:(len - 1) / lambda; % Time grid for regular grid points
    taper = [TW, K];

    [mintime, maxtime] = minmaxsptimes(ts);
    dt = 1 / Fs; % sampling time
    t = mintime - dt:dt:maxtime + dt; % time grid for prolates
    tapers = dpsschk(taper, length(t), Fs); % check tapers
    data_proj = interp1(t', tapers, ts);
    z = sum(data_proj) .^ 2; %
    mz = mean(z); % TODO: need to make it clear

    % Define the power of the characteristic function of a Gaussian process
    % ---------------------------------------------------------------------
    phi2 = @(x) exp(- sigma ^ 2 * (x / mz * lambda * len) .^ 2);

    % spectrum of regular grid
    % ------------------------
    S_N = zeros(1, length(f));
    S_N(ismember(f, lambda:lambda:fmax)) = mz; % TODO: need to make it clear

    % characteristic function power of jittering function
    phi_square = phi2(f);
    S = phi_square .* (S_N - lambda) + lambda; % spectrum of jitterred signal

end % function jitter_spectrum_model

% [EOF]
