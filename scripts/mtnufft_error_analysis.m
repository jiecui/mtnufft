% MTNUFFT_ERROR_ANALYSIS error analysis of MTNUFFT method

% do Monte Carlo simulation
% -------------------------
e_table = mtnu_error_analysis(NumTrials = 1000);

% plot results
% ------------
plot_error_analysis(e_table)

% Copyright 2024 Richard J. Cui. Created: Mon 07/08/2024 16:36:27.301 PM
% $Revision: 0.2 $  $Date: Wed 07/10/2024 09:59:21.776 AM $
%
% Mayo Clinic Foundation
% Rochester, MN 55901, USA
%
% Email: Cui.Jie@mayo.edu
%
% [EOF]
