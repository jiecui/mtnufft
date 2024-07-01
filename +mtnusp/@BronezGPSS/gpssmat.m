function R = gpssmat(fc, fw, t)
    % MTNUSP.BRONEZGPSS.GPSSMAT construct GPSS matrix
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

    % Richard J. Cui. Created: Tue 12/12/2023  5:46:36.143 PM
    % $Revision: 0.2 $  $Date: Fri 12/22/2023 10:09:51.109 AM $
    %
    % Mayo Clinic Foundation
    % Rochester, MN 55901, USA
    %
    % Email: Cui.Jie@mayo.edu

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        fc (1, 1) double % center frequency of frequency band of interest (Hz)
        fw (1, 1) double % half-band width (Hz)
        t (:, 1) double = [] % time points
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    % get the time points
    if isempty(t)
        t = this.TimePoints;
    end % if

    t = sort(t);

    % calculate the GPSS Matrix
    % -------------------------
    N = length(t);
    R = zeros(N, N);

    for i = 1:N

        for j = i:N

            if i == j
                R(i, i) = 2 * fw;
            else
                tau_ij = t(i) - t(j);
                R(i, j) = sin(2 * pi * fw * tau_ij) / (pi * tau_ij) ...
                    * exp(2 * pi * 1j * fc * tau_ij);
                R(j, i) = conj(R(i, j));
            end % for j

        end % for i

    end % function gpssmat

    % R = real(R);

    % [EOF]
