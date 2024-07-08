function plot_error_analysis(e_table, options)
    % PLOT_ERROR_ANALYSIS plot the error analysis of the NUSP object
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

    % Copyright 2024 Richard J. Cui. Created: Sun 07/07/2024 01:03:08.462 AM
    % $Revision: 0.1 $  $Date: Sun 07/07/2024 01:03:08.463 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        e_table (1, 1) struct
    end % positional

    arguments
        options.TimePointsMethod (1, :) string ...
            {mustBeMember(options.TimePointsMethod, ...
             ["Uniform", "MissingPoints", "ArithmeticSampling", "Jittering"])} ...
            = ["Uniform", "MissingPoints", "ArithmeticSampling", "Jittering"] % method to generate time points
        options.XLim (1, 2) double {mustBeReal, mustBeFinite} = [0, 0.5] % x axis limit (Hz)
        options.YLim (1, 2) double {mustBeReal, mustBeFinite} = [2, 10] % y axis limit (dB)
    end % optional

    tp_method = options.TimePointsMethod;
    x_lim = options.XLim;
    y_lim = options.YLim;

    % ======================================================================
    % main
    % ======================================================================
    % plot error analysis
    % -------------------
    for tp_k = tp_method
        et_k = e_table.(tp_k);

        ah_k = show_error_analysis(et_k);
        xlim(ah_k, x_lim)
        ylim(ah_k, y_lim)

        xlabel(ah_k, 'Normalized frequency (Hz)')
        ylabel(ah_k, 'Squared error (dB)')
        title(ah_k, sprintf('Time point method: %s', tp_k))
    end % for

end % function plot_error_analysis

% ==========================================================================
% subroutines
% ==========================================================================
function [fc_o, avg_o, err_o] = remove_inf_nan(fc_i, avg_i, err_i)

    arguments
        fc_i (:, 1) double
        avg_i (:, 1) double
        err_i (:, 1) double
    end % positional

    % remove inf
    % -----------
    inf_idx = isinf(avg_i);
    fc_o = fc_i(~inf_idx);
    avg_o = avg_i(~inf_idx);
    err_o = err_i(~inf_idx);

    % remove nan
    % -----------
    nan_idx = isnan(avg_o);
    fc_o = fc_o(~nan_idx);
    avg_o = avg_o(~nan_idx);
    err_o = err_o(~nan_idx);

end % function

function ah = show_error_analysis(e_table)

    arguments
        e_table (:, :) table
    end % positional

    co = colororder;
    fc = e_table.Frequency;

    figure
    hold on
    % mtls
    [mtls_fc, mtls_avg, mtls_sem] = remove_inf_nan(fc, ...
        e_table.MTLS_mean, e_table.MTLS_sem);
    h_mtls = plot_mean_error(mtls_fc, mtls_avg, mtls_sem, ...
        CurrentAxes = gca, ...
        Color = co(2, :));
    % bgfixed
    [bgfixed_fc, bgfixed_avg, bgfixed_sem] = remove_inf_nan(fc, ...
        e_table.BGFixed_mean, e_table.BGFixed_sem);
    h_bgfx = plot_mean_error(bgfixed_fc, bgfixed_avg, bgfixed_sem, ...
        CurrentAxes = gca, ...
        Color = co(3, :));
    % bgadaptive
    [bgadaptive_fc, bgadaptive_avg, bgadaptive_sem] = remove_inf_nan(fc, ...
        e_table.BGAdaptive_mean, e_table.BGAdaptive_sem);
    h_bgad = plot_mean_error(bgadaptive_fc, bgadaptive_avg, bgadaptive_sem, ...
        CurrentAxes = gca, ...
        Color = co(4, :));
    % mtnufft
    [mtnufft_fc, mtnufft_avg, mtnufft_sem] = remove_inf_nan(fc, ...
        e_table.MTNUFFT_mean, e_table.MTNUFFT_sem);
    h_mtnf = plot_mean_error(mtnufft_fc, mtnufft_avg, mtnufft_sem, ...
        CurrentAxes = gca, ...
        Color = co(1, :));

    % add legend
    legend([h_mtls.XHandle, h_bgfx.XHandle, h_bgad.XHandle, h_mtnf.XHandle], ...
        "MTLS", "BGFixed", "BGAdaptive", "MTNUFFT", ...
        Location = "best")

    ah = gca;

end % function

% [EOF]
