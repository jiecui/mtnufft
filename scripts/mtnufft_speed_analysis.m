% MTNUFFT_SPEED_ANALYSIS speed analysis of MTNUFFT method

% do Monte Carlo simulation
% -------------------------
s_table = mtnu_speed_analysis(NumTrials = 1000);

% plot results
% ------------
plot_speed_analysis(s_table);

% Copyright 2024 Richard J. Cui. Created: Wed 07/10/2024 09:59:21.776 AM
% $Revision: 0.1 $  $Date: Wed 07/10/2024 09:59:21.776 AM $
%
% Mayo Clinic Foundation
% Rochester, MN 55901, USA
%
% Email: Cui.Jie@mayo.edu
%
% [EOF]
