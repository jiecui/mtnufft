classdef NUPointbinned < matlab.mixin.SetGet & dynamicprops
    % Class MTNUSP.NUPOINTBINNED.NUPOINTBINNED multi-taper NUFFT of binned point process

    % Copyright 2022 Richard J. Cui. Created: Sun 07/17/2022  2:28:34.804 PM
    % $Revision: 0.2 $  $Date: Tue 02/13/2024 12:46:22.715 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    properties

    end % properties

    methods

        function this = NUPointbinned()

        end

    end % methods

    % other methods
    % -------------
    methods
        varargout = mtnufftpb(this, varargin) % NUFFT of binned point process using Thomson's multitaper method
    end % methods

end % classdef

% [EOF]
