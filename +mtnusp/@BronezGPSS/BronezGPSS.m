classdef BronezGPSS < handle
    % Class MTNUSP.BRONEZGPSS spectrum estimation using generalized prolate spheroidal sequences
    %
    % Syntax:
    %
    % Note:
    %   signal, x, can be real/complex
    %
    %   Method of selecting the dominant eigenvectors:
    %   * adaptive
    %   * auto
    %   * leakage
    %
    % Reference:
    %
    %   [1] Bronez 1985
    %   [2] Bronez 1988
    %
    % See also: .

    % Copyright 2023-2024 Richard J. Cui. Created: Tue 12/12/2023 12:10:09.917 AM
    % $Revision: 0.9 $  $Date: Sun 06/23/2024 10:34:24.556 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    properties
        AnalysisBand (:, 2) double % analysis band = [center, half-width] (Hz)
        BandIncrement (1, 1) double % band increment (Hz)
        DoParallel (1, 1) logical % do parallel computation
        LambdaFactor (1, 1) double % factor for choosing eigenvectors
        MaxFrequency (1, 1) double % maximum frequency of the signal (Hz)
        MaxHalfWidth (:, 1) double % maximum half width at each frequency band (Hz)
        NumTapers (1, 2) int8 % [min, max] number of tapers allowed
        SelectionMethod (1, 1) string % method for selecting dominant eigenvectors
        SelectedDomainEigenvectors (1, :) struct % selected domain eigenvectors
        SignalSamples (:, 1) double % set of signal samples
        SignalBand (1, 2) double % signal band = [center, half-width] (Hz)
        TimePoints (:, 1) double % set of time points (seconds)
        T (1, 1) double % duration
        Verbose (1, 1) logical % verbose mode
    end % properties

    % ======================================================================
    % constructor
    % ======================================================================
    methods

        function this = BronezGPSS(t, x, options)

            % parse input parameters
            % ----------------------
            arguments
                t (:, 1) double {mustBeNumeric} = [] % time sampling points (seconds)
                x (:, 1) double {mustBeNumeric} = [] % signal process samples
            end % positional

            arguments
                options.AnalysisBand (:, 2) double = ...
                    [linspace(0, .5 -1/40, 20).', ones(20, 1) * 1/40] % [center, half-width] (Hz)
                options.BandIncrement (1, 1) double = 1/100 % band increment (Hz)
                options.DoParallel (1, 1) logical = false % do parallel computation
                options.LambdaFactor (1, 1) double = 1e-2 % see Note above
                options.MaxFrequency (1, 1) double = 1/2 % maximum frequency (Hz)
                options.MaxHalfWidth (:, 1) double = 1/2 % maximum half-width (Hz)
                options.NumTapers (1, 2) int8 = [1, 100] % [min, max] number of tapers allowed
                options.SignalBand (1, 2) double = [0, .5] % [center, half-width] (Hz)
                options.SelectionMethod (1, 1) string ...
                    {mustBeMember(options.SelectionMethod, ["auto", "adaptive", "leakage"])} = "auto" % method for selecting dominant eigenvectors
                options.T (1, 1) double = 1 % duration (seconds)
                options.Verbose (1, 1) logical = false % verbose mode
            end % optional

            A = options.AnalysisBand;
            B = options.SignalBand;
            do_par = options.DoParallel;
            l_factor = options.LambdaFactor;
            fmax = options.MaxFrequency;
            max_fw = options.MaxHalfWidth;
            method = options.SelectionMethod;
            num_tapers = options.NumTapers;
            dt_fw = options.BandIncrement;
            T = options.T;
            verbose = options.Verbose;

            % initialize supper class
            % -----------------------

            % initialize properties
            % ---------------------
            % * time points, sampling points and duration
            if isempty(t)
                this.TimePoints = [];
                this.SignalSamples = [];
                this.T = T;
            else
                assert(length(t) == length(x), 'MTNUSP:BronezGPSS:constructor', ...
                't and x must have the same length.')
                [this.TimePoints, idx] = sort(t);
                this.SignalSamples = x(idx);

                if isempty(T)
                    this.T = max(this.TimePoints) - min(this.TimePoints);
                else
                    this.T = T;
                end % if

            end % if

            % * other parameters
            this.DoParallel = do_par;
            this.MaxFrequency = fmax;
            this.MaxHalfWidth = max_fw;
            this.NumTapers = num_tapers;
            this.SignalBand = B;
            this.AnalysisBand = A;
            this.SelectionMethod = method;
            this.SelectedDomainEigenvectors = struct([]);
            this.BandIncrement = dt_fw;
            this.LambdaFactor = l_factor;
            this.Verbose = verbose;

        end % constructor

    end % constructor

    % ======================================================================
    % other methods
    % ======================================================================
    methods (Static)
        varargout = gpssmat(varargin) % construct GPSS matrix
        varargout = leakage(varargin) % estimate leakage of weight sequence
    end % static methods

    methods
        varargout = est_int_spectrum(this, varargin) % estimate integrated spectrum
        varargout = get_time_points(this, varargin) % get time points
        varargout = spectrumgpss(this, varargin) % estimate spectrum using gpss
        varargout = sel_dominant_eig(this, varargin) % dominant eigenvalues and eigenvectors
    end % methods

end % classdef

% [EOF]
