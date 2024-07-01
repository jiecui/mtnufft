classdef NUContinuous < matlab.mixin.SetGet & dynamicprops
    % Class NUCONTINUOUS spectral analysis of non-uniformly sampled continuous signals

    % 2022 Richard J. Cui. Created: Sun 04/03/2022  1:40:12.575 AM
    % $Revision: 0.8 $  $Date: Mon 12/25/2023 12:46:02.730 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % properties of the class
    % ======================================================================
    properties
        Signal
        SamplingTimes
    end % properties

    % ======================================================================
    % methods
    % ======================================================================
    % constructor
    % -----------
    methods

        function this = NUContinuous(x, t)
            % parse inputs
            % ------------
            arguments
                x (:, :) double = double.empty(1, 0) % signal either vector or matrix n x channels
                t (1, :) double = double.empty(1, 0) % time instances of sampling signal
            end % arguments

            % superclass constructor (if any)
            % -------------------------------

            % initialize class properties
            % ---------------------------
            if isvector(x)
                x = x(:);
            end % if

            set(this, 'Signal', x);
            set(this, 'SamplingTimes', t);

        end % constructor

    end % methods

    % other methods
    % -------------
    methods
        varargout = mtnufft(this, varargin) % NUFFT using Thomson's multitaper method
        varargout = mtlsspectrum(this, varargin) % multitaper spectrum using Lomb-Scargle method
        varargout = mtnuspectrum(this, varargin) % nonuniform sampling spectrum using Thomson's multitaper method
        varargout = nuftest(this, varargin) % F-test nonuniformly sampled signals using Thomson's multitaper method
        varargout = pmtlomb(this, varargin) % Lomb-Scargle periodogram using Thomson's multitaper method
    end % methods

end % classdef

% [EOF]
