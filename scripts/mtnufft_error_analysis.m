% MTNUFFT_ERROR_ANALYSIS error analysis of MTNUFFT method

% do Monte Carlo simulation
% -------------------------
e_table = mtnu_error_analysis(NumTrials = 1000);

% plot results
% ------------
show_error_analysis(e_table);

% Copyright 2024 Richard J. Cui. Created: Mon 07/08/2024 16:36:27.301 PM
% $Revision: 0.1 $  $Date: Mon 07/08/2024 16:36:27.305 PM $
%
% Mayo Clinic Foundation
% Rochester, MN 55901, USA
%
% Email: Cui.Jie@mayo.edu
%
% [EOF]
