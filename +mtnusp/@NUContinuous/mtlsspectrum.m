function [S, f, S_eig, x, t] = mtlsspectrum(this, x, t, options)
    % NUCONTINUOUS.MTLSSPECTRUM spectral esitmation using PMTLOMB method
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

    % Copyright 2023 Richard J. Cui. Created: Mon 12/25/2023 12:49:32.806 PM
    % $Revision: 0.1 $  $Date: Mon 12/25/2023 12:49:32.826 PM $
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

    fc = options.InputFrequencies;
    fmax = options.MaxFrequency;
    olfc = options.OversamplingFactor;
    fw = options.TimeHalfbandwidth;
    tapers = options.Tapers;

    % ======================================================================
    % main
    % ======================================================================
    % calculate the eigenspectrum
    % ---------------------------
    [S_eig, f, x, t] = this.pmtlomb(x, t, ...
        InputFrequencies = fc, ...
        MaxFrequency = fmax, ...
        OversamplingFactor = olfc, ...
        TimeHalfbandwidth = fw, ...
        Tapers = tapers);

    % calculate the multitaper spectrum
    % ---------------------------------
    num_points = length(x);
    S_scaled = S_eig * num_points; % scale to FFT definition of periodogram
    S = mean(S_scaled(:, :, 1), 2); % average over tapers

end % function mtlsspectrum

% [EOF]
