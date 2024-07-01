function Px = est_int_spectrum(this, domin_ev, x)
    % MTNUSP.BRONEZGPSS.EST_INT_SPECTRUM estimate the integrated spectrum
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
    % See also BronezGPSS.cal_generalized_eig.

    % Richard J. Cui. Created: Mon 12/18/2023  4:16:23.634 PM
    % $Revision: 0.6 $  $Date: Sun 06/23/2024 10:39:27.108 PM $
    %
    % Mayo Clinic Foundation
    % Rochester, MN 55901, USA
    %
    % Email: Cui.Jie@mayo.edu

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) mtnusp.BronezGPSS
        domin_ev (1, :) struct = struct([]) % the domain eigenvalues and eigenvectors (see cal_generalized_eig)
        x (:, 1) double = [] % signal samples
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    verbose = this.Verbose;

    % get domain eigenvectors
    % -----------------------
    if isempty(domin_ev)
        domin_ev = this.SelectedDomainEigenvectors;
    end % if

    % get signal
    % ----------
    if isempty(x)
        x = this.SignalSamples;
    end % if

    % get analysis bands
    % ------------------
    % estimate the integrated spectrum
    % --------------------------------
    num_bands = length(domin_ev);
    Px = table('Size', [num_bands, 4], 'VariableTypes', {'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'BandCenter_Hz', 'Lower_Hz', 'Upper_Hz', 'IntPower'});

    if verbose
        pb = CmdLineProgressBar('Estimating the integrated spectrum...');
    end % if

    for k = 1:num_bands

        if verbose
            pb.print(k, num_bands);
        end % if

        if isempty(domin_ev(k).Eigenvectors)
            continue;
        end % if

        % get the analysis band
        fc_k = domin_ev(k).CenterFrequency;
        fw_k = domin_ev(k).BandHalfWidth;
        lb_k = fc_k - fw_k;
        ub_k = fc_k + fw_k;
        Px.BandCenter_Hz(k) = fc_k;
        Px.Lower_Hz(k) = lb_k;
        Px.Upper_Hz(k) = ub_k;

        % get the domain eigenvectors
        V_k = domin_ev(k).Eigenvectors;
        K_k = size(V_k, 2); % number of eigenvectors
        px_k = zeros(1, K_k);

        for m = 1:K_k
            % get the eigenvector
            w_m = V_k(:, m);

            % get the projection
            px_k(m) = conj(w_m' * x) * (w_m' * x);
        end % for

        % average
        Px.IntPower(k) = mean(px_k);

    end % for

    clear pb

end % function est_int_spectrum

% [EOF]
