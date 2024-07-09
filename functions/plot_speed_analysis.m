function plot_speed_analysis(s_table, options)
    % PLOT.PLOT_SPEED_ANALYSIS plot the speed analysis of the NUSP object
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

    % Richard J. Cui. Created: Mon 07/08/2024 11:43:03.811 PM
    % $Revision: 0.1 $  $Date: Mon 07/08/2024 11:43:03.811 PM $
    %
    % Mayo Clinic Foundation
    % Rochester, MN 55901, USA
    %
    % Email: Cui.Jie@mayo.edu

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        s_table (:, :) table
    end % positional

    arguments
        options.TimePointsMethod (1, :) string ...
            {mustBeMember(options.TimePointsMethod, ...
             ["Uniform", "Jittering", "MissingPoints", "ArithmeticSampling"])} ...
            = ["Uniform", "Jittering", "MissingPoints", "ArithmeticSampling"] % method to generate time points
        options.YLim (1, 2) double {mustBeReal, mustBeFinite} = [1, 2000] % y axis limit (dB)
    end % optional

    tp_method = options.TimePointsMethod;
    y_lim = options.YLim;

    % ======================================================================
    % main
    % ======================================================================
    % plot speed analysis
    % -------------------
    d_str = ["MTNUFFT_mean", "MTLS_mean", "BGFixed_mean", "BGAdaptive_mean"];
    e_str = ["MTNUFFT_std", "MTLS_std", "BGFixed_std", "BGAdaptive_std"];

    % * build array
    data = table2array(s_table(ismember(s_table.TimePointMethod, tp_method), d_str));
    err = table2array(s_table(ismember(s_table.TimePointMethod, tp_method), e_str));

    % * create grouped bar plot
    figure
    bh = bar(gca, data);
    hold on

    % * add error bars
    num_bars = size(data, 1);
    num_groups = size(data, 2);
    group_width = min(.8, num_bars / (num_bars + 1.5));

    for k = 1:num_bars
        % Calculate center of each bar
        x_k = (1:num_groups) - group_width / 2 + (2 * k - 1) * group_width / (2 * num_bars);
        errorbar(x_k, data(:, k), err(:, k), 'k', 'linestyle', 'none');
    end

    hold off

    % * add other info
    set(gca, 'XTickLabel', ["Uniform", "Jittering", "Missing", "Arithematic"])
    set(gca, 'YScale', 'log')
    ylim(y_lim)
    legend(bh([2:4, 1]), ["MTLS", "BGFiexed", "BGAdaptive", "MTNUFFT"])
    xlabel('Sampling scheme')
    ylabel('Number of spectra (/s)')
    title('Speed of spectral estimation')

end % function plot_speed_analysis

% [EOF]
