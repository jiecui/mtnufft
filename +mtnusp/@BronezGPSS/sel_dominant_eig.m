function domin_ev = sel_dominant_eig(this, options)
    % MTNUSP.BRONEZGPSS.SEL_DOMINANT_EIG select dominant eigenvalues and eigenvectors
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

    % Richard J. Cui. Created: Mon 12/18/2023 12:39:55.293 PM
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
    end % positional

    arguments
        options.DoParallel (1, 1) ...
            {mustBeA(options.DoParallel, {'logical', 'double'})} = nan % use parallel computing
        options.Verbose (1, 1) ...
            {mustBeA(options.Verbose, {'logical', 'double'})} = nan % show progress
    end % optional

    do_par = options.DoParallel;

    if isnan(do_par)
        do_par = this.DoParallel;
    end % if

    verbose = options.Verbose;

    if isnan(verbose)
        verbose = this.Verbose;
    end % if

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    t = this.TimePoints;
    B = this.SignalBand; % [center, half-width] (Hz)
    A = this.AnalysisBand; % [center, half-width] (Hz)
    method = this.SelectionMethod; % method to select eigenvectors
    max_fw = this.MaxHalfWidth; % maximum analysis half-band width (Hz)
    num_tapers = this.NumTapers; % [min, max] number of tapers allowed in a band
    dt_fw = this.BandIncrement; % time step for analysis half-band width (Hz)
    l_factor = this.LambdaFactor; % see Note in constructor

    c = parcluster('local');

    if do_par
        num_workers = c.NumWorkers;
    else
        num_workers = 0;
    end % if

    % calculate dominant eigenvalues and eigenvectors for each analysis band
    % ----------------------------------------------------------------------
    num_bands = size(A, 1);
    domin_ev = struct('Eigenvalues', [], 'Eigenvectors', []);
    RB = this.gpssmat(B(1), B(2), t); % gpss matrix for signal band
    fc = A(:, 1); % center frequency of analysis bands
    fw = A(:, 2); % half-band width of analysis bands

    parfor (k = 1:num_bands, num_workers)

        if verbose
            fprintf('Calculating dominant eigs: %d out of %d... ', k, num_bands);
        end % if

        fc_k = fc(k); % center frequency of k-th analysis band
        fw_k = fw(k); % half-band width of k-th analysis band
        lambda_k = []; % eigenvalues of k-th analysis band
        V_k = []; % eigenvectors of k-th analysis band

        % select K dominant eigs according to method
        switch method
            case "auto"
                [lambda_k, V_k, fw_k] = auto_domi_eigs(this, t, fc_k, fw_k, RB, l_factor);
            case "adaptive"
                [lambda_k, V_k, fw_k] = adaptive_domi_eigs(this, t, fc_k, fw_k, ...
                    max_fw, dt_fw, RB, l_factor, num_tapers);
            case "leakage"
                [lambda_k, V_k, fw_k] = leakage_domi_eigs(this, t, fc_k, fw_k, RB, l_factor);
        end % switch

        % save eigenvectors and eigenvalues
        domin_ev(k).CenterFrequency = fc_k;
        domin_ev(k).BandHalfWidth = fw_k;
        domin_ev(k).Eigenvalues = lambda_k;
        domin_ev(k).Eigenvectors = V_k;

        if verbose
            fprintf('done!\n');
        end % if

    end % for k

    this.SelectedDomainEigenvectors = domin_ev;

end % function cal_generalized_eig

% ==========================================================================
% subroutines
% ==========================================================================
function [d_lambda, d_V, d_fw] = adaptive_domi_eigs(this, t, fc, fw0, max_fw, dt_fw, RB, ul_leakage, num_tapers)
    % find the dominant eigenvalues and eigenvectors adaptively
    % first search the maximum number of eigs allowed in a band. If not found,
    % increase the analysis half-band width and search again until the upper
    % limit of leakage is reached.

    arguments
        this (1, 1) mtnusp.BronezGPSS
        t (:, 1) double % time points (s)
        fc (1, 1) double % center frequency (Hz)
        fw0 (1, 1) double % initial analysis half-band width (Hz)
        max_fw (1, 1) double % maximum analysis half-band width (Hz)
        dt_fw (1, 1) double % time step for analysis half-band width (Hz)
        RB (:, :) double % gpss matrix for signal band
        ul_leakage (1, 1) double % upper limit of leakage
        num_tapers (1, 2) double % number of tapers allowed in a band
    end % positional

    % search the band width
    for fw_k = fw0:dt_fw:max_fw
        % gpss matrix of analysis band
        RA_k = this.gpssmat(fc, fw_k, t);

        % calculate generalized eigenvalues and eigenvectors
        [V_k, D_k] = eig(RA_k, RB);

        % sort eigenvalues in descending order
        [lambda_k, idx_k] = sort(diag(D_k), 'descend');
        V_k = V_k(:, idx_k);

        % adaptive selection
        % ------------------
        % search the maximum number of eigs allowed in a band
        for K_j = num_tapers(2):-1:num_tapers(1) % search from maximum to minimum
            lambda_kj = lambda_k(K_j);
            gamma_kj = this.leakage(abs(lambda_kj)); % deal with numeric error

            if gamma_kj <= ul_leakage
                found_flag = true;
                break
            else
                found_flag = false;
            end % if

        end % for

        if found_flag
            break
        end % if

    end % while

    % get the dominant eigenvalues and eigenvectors
    if found_flag
        d_fw = fw_k;
        d_lambda = lambda_k(1:K_j);
        V = V_k(:, 1:K_j);
        % normalize eigenvectors
        d_V = normalize_eigenvectors(V, RB, d_fw);

    else
        d_fw = [];
        d_lambda = [];
        d_V = [];
        fprintf("No dominant eig found for band centered at %.2f!\n", fc);
    end % if

end % function

function [d_lambda, d_V, fw] = leakage_domi_eigs(this, t, fc, fw, RB, l_factor)

    arguments
        this (1, 1) mtnusp.BronezGPSS
        t (:, 1) double % time points (s)
        fc (1, 1) double % center frequency (Hz)
        fw (1, 1) double % analysis half-band width (Hz)
        RB (:, :) double % gpss matrix for signal band
        l_factor (1, 1) double % see lambda factor for selection (Note in constructor)
    end % positional

    % gpss matrix of analysis band
    RA = this.gpssmat(fc, fw, t);

    % calculate generalized eigenvalues and eigenvectors
    [V, D] = eig(RA, RB);

    % sort eigenvalues in descending order
    [lambda, idx] = sort(diag(D), 'descend');
    V = V(:, idx);

    % find dominant eigenvalues and eigenvectors
    [d_lambda, V] = find_dominant_eigenvectors(this, lambda, V, "leakage", l_factor);

    % normalize eigenvectors
    d_V = normalize_eigenvectors(V, RB, fw);

end % function

function [d_lambda, d_V, fw] = auto_domi_eigs(this, t, fc, fw, RB, l_factor)

    arguments
        this (1, 1) mtnusp.BronezGPSS
        t (:, 1) double % time points (s)
        fc (1, 1) double % center frequency (Hz)
        fw (1, 1) double % analysis half-band width (Hz)
        RB (:, :) double % gpss matrix for signal band
        l_factor (1, 1) double % see lambda factor for selection (Note in constructor)
    end % positional

    % gpss matrix of analysis band
    RA = this.gpssmat(fc, fw, t);

    % calculate generalized eigenvalues and eigenvectors
    [V, D] = eig(RA, RB);

    % sort eigenvalues in descending order
    [lambda, idx] = sort(diag(D), 'descend');
    V = V(:, idx);

    % find dominant eigenvalues and eigenvectors
    [d_lambda, V] = find_dominant_eigenvectors(this, lambda, V, "auto", l_factor);

    % normalize eigenvectors
    d_V = normalize_eigenvectors(V, RB, fw);

end % function

function [d_lambda, d_V, K] = find_dominant_eigenvectors(this, lambda, V, method, l_factor)

    arguments
        this (1, 1) mtnusp.BronezGPSS
        lambda (:, 1) double % eigenvalues
        V (:, :) double % eigenvectors
        method (1, 1) string % method to select eigenvectors
        l_factor (1, 1) double % see Note in constructor
    end % positional

    switch method
        case 'auto'

            if l_factor < 1
                idx = lambda >= lambda(1) * (1 - l_factor);
            else
                idx = 1:round(l_factor);
            end % if

        case 'leakage'
            gamma = this.leakage(abs(lambda)); % deal with numeric error
            idx = gamma <= l_factor;

    end % switch

    d_lambda = lambda(idx);
    d_V = V(:, idx);
    K = length(d_lambda);

end % function

function V_nm = normalize_eigenvectors(V, R, fw)

    arguments
        V (:, :) double
        R (:, :) double
        fw (1, 1) double % analysis hald-band withd (hz)
    end % positional

    % calculate normalization matrix
    % ------------------------------
    V_nm = zeros(size(V));
    num_vec = size(V, 2);

    for k = 1:num_vec
        w_k = V(:, k);
        r_k = w_k' * R * w_k;
        s_k = sqrt(2 * fw / r_k);
        V_nm(:, k) = s_k * w_k;
    end % for k

end % function

% [EOF]
