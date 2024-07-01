function gamma = leakage(lambda)
    % MTNUSP.BRONEZGPSS.LEAKAGE estimate the leakage, or the sidelobe energy, of weight sequences
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

    % Richard J. Cui. Created: Tue 12/19/2023  9:07:47.643 AM
    % $Revision: 0.1 $  $Date: Tue 12/19/2023  9:07:47.806 AM $
    %
    % Mayo Clinic Foundation
    % Rochester, MN 55901, USA
    %
    % Email: Cui.Jie@mayo.edu

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        lambda (1, :) double {mustBeReal, mustBeFinite, mustBeNonnegative} % eigenvalues
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    gamma = 10 * log10(1 - lambda); % sidelobe energy (dB)

end % function leakage

% [EOF]
