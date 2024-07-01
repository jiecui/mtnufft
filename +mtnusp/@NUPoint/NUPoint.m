classdef NUPoint < matlab.mixin.SetGet & dynamicprops
    % Class MTNUSP.NUPOINT.NUPOINT multi-taper NUFFT of point sources
    %
    % Syntax:
    %
    % Input(s):
    %
    % Output(s):
    %
    % See also .
    %
    % Note:
    %   Require Chronux toolbox.

    % Copyright 2024 Richard J. Cui. Created: Tue 02/13/2024 12:42:54.695 PM
    % $Revision: 0.1 $  $Date: Tue 02/13/2024 12:42:54.726 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % properties of the class
    % ======================================================================
    properties

    end % properties

    % ======================================================================
    % methods
    % ======================================================================
    % constructor
    % -----------
    methods

        function this = NUPoint()

        end

    end % methods

    % other methods
    % -------------
    methods

    end % methods

    % static methods
    methods (Static)
        varargout = jitter_spectrum_model(varargin) % power spectrum of jittering model
    end % static method

end % classdef

% [EOF]
