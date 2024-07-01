function [pxx, px] = spectrumgpss(this)
    % MTNUSP.BRONEZGPSS.GPSSSPECTRUM power spectrum estimation using GPSS tapering
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

    % Richard J. Cui. Created: Fri 12/22/2023  5:46:06.170 PM
    % $Revision: 0.4 $  $Date: Sun 06/23/2024 10:39:27.108 PM $
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

    % ======================================================================
    % main
    % ======================================================================
    % select dominant eigs
    % --------------------
    d_ev = this.SelectedDomainEigenvectors;

    if isempty(d_ev)
        this.sel_dominant_eig();
    end % if

    % estmate the integrated power
    % ----------------------------
    px = this.est_int_spectrum();

    % convert to power density
    % ------------------------
    bw = px.Upper_Hz - px.Lower_Hz;
    pxx = px.IntPower ./ bw; % accumulated power

end % function gpssspectrum

% [EOF]
