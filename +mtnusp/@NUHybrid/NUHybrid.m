classdef NUHybrid < matlab.mixin.SetGet & dynamicprops
    % Class NUHYBRID (summary)

    % Copyright 2022 Richard J. Cui. Created: Sun 07/17/2022  2:31:32.066 PM
    % $Revision: 0.1 $  $Date: Sun 07/17/2022  2:31:32.108 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    properties

    end % properties

    methods

        function this = NUHybrid()

        end

    end % methods

    % other methods
    % -------------
    methods
        varargout = nucoherencycpb(this, varargin) % coherency between nonuniformly sampled continous signal and binned point process
    end % methods
end % classdef

% [EOF]
