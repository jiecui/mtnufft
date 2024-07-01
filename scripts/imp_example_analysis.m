% IMP_EXAMPLE_ANALYSIS analysis of example impedance signals

%% hyper-parameters
% -----------------
co = colororder;
ver_bose = true; % show more detials
do_par = true; % do parallel computation

% * interval of extracting the signal
num_bands = 240; % number of analysis bands
fmax = 46; % frequency limit of the signal
intv = [150, 160]; % [start, end] in days
fs = 24; % desired sampling frequency (cycles/day)
fc_min = 0;
fc_max = fs / 2; % maximum frequency of interest

% * number of Raleigh frequency
Nr = fs * diff(intv);

% * multitaper parameters
TW = 3.5; % time-bandwidth product
K = 6; % number of tapers
params.tapers = [TW, K];
params.Fs = fs;
params.fpass = [fc_min, fc_max];
params.pad = 2;

%% read in sample data
% -------------------
z = readtable('impedance_example.csv');
t = z.Time; % time from implantaion (days)
imp = z.Impedance; % (Ohm)

% * extract the signal
intv1 = [intv(1) - 1, intv(2) + 1];
tx_idx = t >= intv1(1) & t <= intv1(2);
tx = t(tx_idx);
x = imp(tx_idx);

idx = tx >= intv(1) & tx <= intv(2);
txx = tx(idx);
xx = x(idx);

% * resampling to the desired sampling rate
[y, ty] = rs_nonuniform_sig(x, tx, fs);

% * remove the boundary effect
idx = ty >= intv(1) & ty <= intv(2);
ty = ty(idx);
y = y(idx);

% * plot
figure
hr = plot(ty, y, LineWidth = 1);
hold on
hi = scatter(txx, xx, 10, 'o', 'filled');
grid on
grid minor
set(gca, 'XTick', intv(1):1:intv(2))
legend([hr, hi], 'Resampled signal', 'Original signal')
xlabel('Time from implantation (days)')
ylabel('Impedance (\Omega)')
title('Resampled signal of M1 EL_00', 'Interpreter', 'none')

%% estimate spectrum using uniformly re-sampled signal
% ----------------------------------------------------
yc = detrend(y, 'constant');
% * power spectrum estimation with CHRONUX multitaper method
[pxx_mtc, f] = mtspectrumc(yc, params);
figure
plot(f, pow2db(pxx_mtc), 'LineWidth', 1)
xlim([fc_min fc_max])
ylim([-20, 40])
grid on
grid minor
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('MTSPECTRUMC of resampled signal')

% * power spectrum estimation with MATLAB periodogram
pxx_p = periodogram(yc, [], f, fs);
figure
plot(f, pow2db(pxx_p), 'LineWidth', 1)
xlim([fc_min fc_max])
ylim([-20, 40])
grid on
grid minor
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('Periodogram of resampled signal')

% * power spectrum estimation with MATLAB pwelch method
pxx_wl = pwelch(yc, 4 * fs, 3 * fs, f, fs);
figure
plot(f, pow2db(pxx_wl), 'LineWidth', 1)
xlim([fc_min fc_max])
ylim([-20, 40])
grid on
grid minor
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('PWELCH of resampled signal')

% * power spectrum estimation with MATLAB multitaper method
pxx_mtm = pmtm(yc, {TW, K}, f, fs, "Tapers", "slepian", "unit");
figure
plot(f, pow2db(pxx_mtm), 'LineWidth', 1)
xlim([fc_min fc_max])
ylim([-20, 40])
grid on
grid minor
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('PMTM of resampled signal')

%% Spectrum of sampling times
% ---------------------------
fs_hat = length(txx) / (txx(end) - txx(1)); % estimated sampling frequency
par_p = params;
par_p.Fs = 40 * fs_hat;
par_p.fpass = [0, par_p.Fs / 2];
par_p.err = [1, 0.05]; % 95 % confidence interval
[pxx_p, f_p, R_p, err_p] = mtspectrumpt(txx - txx(1), par_p);

figure
plot(f_p, pow2db(err_p(1, :)), 'g', 'LineWidth', .5)
hold on
plot(f_p, pow2db(err_p(2, :)), 'g', 'LineWidth', .5)
plot(f_p, pow2db(pxx_p), 'LineWidth', 1)
yline(pow2db(R_p), 'r', 'LineWidth', .5)
xlim tight
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('Spectrum of sampling times')

%% try to figure out a jittering model fit to the spectrum
% --------------------------------------------------------
len0 = length(txx);
Fs0 = 1350 * 2;
TW0 = TW;
K0 = K;
lambda0 = 4 * hours(days(1));
sigma0 = 20 / seconds(days(1));

[S0, f0] = mtnusp.NUPoint.jitter_spectrum_model(len0, lambda0, sigma0, Fs0, ...
    TimeHalfbandwidth = TW0, ...
    NumberTapers = K0);

% plot
figure
plot(f_p, pow2db(pxx_p))
hold on
plot(f0, pow2db(S0), 'LineWidth', 1)
yline(pow2db(lambda0), 'LineWidth', 1, 'Color', co(5, :))
ylim([0 45])
xlim([0 Fs0 / 2])
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('Spectrum model fit to sampling times')

%% F-test of uniformly sampled signal
% -----------------------------------
f_u = f;
Fval_u = ftestc(yc, params);

p1 = 0.01;
p2 = 0.001;
rot = 1 / Nr; % rule of thumb (rot)
sig1_u = finv(1 - p1, 2, 2 * K - 2); % F-inverse cumulative distribution
sig2_u = finv(1 - p2, 2, 2 * K - 2);
sig_rot_u = finv(1 - rot, 2, 2 * K - 2);

%% estimate spectrum using multitaper Lomb-Scargle (MTLS) periodogram
% -------------------------------------------------------------------
% * MATLAB LS periodogram
s = detrend(xx, 0);
ts = txx - txx(1);
pxx_ls = plomb(s, ts, f);
% plot
figure
plot(f, pow2db(pxx_ls), 'LineWidth', 1)
xlim([fc_min fc_max])
ylim([-20, 40])
grid on
grid minor
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('LS periodogram of nonuniform signal')

% * MTLS periodogram
N = length(ts);
ofact = 10;
nus = mtnusp.NUContinuous(s, ts);
[pxx, f] = nus.pmtlomb('TimeHalfbandwidth', TW, 'MaxFrequency', fs / 2, 'OverSamplingfactor', 10);
pxx = pxx * N; % scale to FFT definition of periodogram
pxx = mean(pxx(:, :, 1), 2); % average over the number of tapers

% plot
figure
plot(f, pow2db(pxx), 'LineWidth', 1)
xlim([fc_min fc_max])
ylim([-20, 40])
grid on
grid minor
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('MTLS power spectrum of nonuniform signal')

%% estimate spectrum using multitaper nonuniform Fourier transform (MTNUFFT)
% --------------------------------------------------------------------------
fc = linspace(fc_min, fc_max, num_bands)'; % frequency grid
s = detrend(xx, 'Constant');
ts = txx - txx(1);
T = ceil(ts(end)); % signal length (days)
f_w = TW / T; % half bandwidth of analysis bands (1/day)
nus = mtnusp.NUContinuous(s, ts);
pxx_mtnu = nus.mtnuspectrum('QuerryFrequencies', fc, ...
    'MaxFrequency', fmax, ...
    'NormMethod', 'BGNorm', ...
    'Halfbandwidth', f_w, ...
    'TimeHalfbandwidth', TW, ...
    'NumberTapers', K); % average over tapers; one channel only
pow_mtnu = pow2db(pxx_mtnu); % power spectrum (dB)

% plot
figure
plot(fc, pow_mtnu, 'LineWidth', 1)
xlim([fc_min fc_max])
ylim([-20, 40])
grid on
grid minor
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('MTNUFFT power spectrum of nonuniform signal')

%% nonuniform F-test
% ------------------
N = length(s);
Fval = nus.nuftest(QuerryFrequencies = fc, ...
    TimeHalfbandwidth = TW, ...
    NumberTapers = K);

p1 = 0.01;
p2 = 0.001;
rot = 1 / Nr; % rule of thumb (rot)
sig1 = finv(1 - p1, 2, 2 * K - 2); % F-inverse cumulative distribution
sig2 = finv(1 - p2, 2, 2 * K - 2);
sig_rot = finv(1 - rot, 2, 2 * K - 2);

%% estimate spectrum using Bronez-fixed multitaper method
% -------------------------------------------------------
num_bands = length(fc); % number of analysis bands
fw = f_w * ones(num_bands, 1); % half bandwidth
A = [fc, fw]; % analysis bands
B = [0, fmax]; % signal band

bgfx = mtnusp.BronezGPSS(ts, s, ...
    T = T, ...
    SignalBand = B, ...
    AnalysisBand = A, ...
    MaxFrequency = fmax, ...
    SelectionMethod = 'auto', ...
    NumTapers = [K, 2 * K], ...
    LambdaFactor = K, ...
    DoParallel = do_par, ...
    Verbose = ver_bose);
pxx_bgfx = bgfx.spectrumgpss();
pow_bgfx = pow2db(pxx_bgfx); %

%% MTNUFF0 - use precomputed eignevectors at fc = 0
% -------------------------------------------------
d_ev = bgfx.SelectedDomainEigenvectors;
pxx_mtnu0 = nus.mtnuspectrum('QuerryFrequencies', fc, ...
    'Timehalfbandwidth', TW, ...
    'NormMethod', 'BGNorm', ...
    'Halfbandwidth', f_w, ...
    'Tapers', d_ev(fc == 0).Eigenvectors, ... % factor 10 to make it consistent with BronezGPSS
    'NumberTapers', K);
pow_mtnu0 = pow2db(pxx_mtnu0); % power spectrum (dB)

%% suboptimality of the MTNUFFT method
% ------------------------------------
lambda_0 = d_ev(fc == 0).Eigenvalues;
Delta_B = zeros(1, num_bands);

for i = 1:num_bands
    lambda_i = d_ev(i).Eigenvalues;
    db_i = abs(d_ev(i).Eigenvalues - lambda_0);
    Delta_B(i) = sum(db_i) / K;
end % for

% plot
figure
plot(fc, Delta_B, 'LineWidth', 1)
xlim([fc_min, fc_max])
xlabel('Center Frequency (cycles/day)')
ylabel('Suboptimality')
title('Approximate suboptimality of MTNUFFT method')

%% plot spectrum of Bronez-fixed multitaper and MTNUFFT0
% ------------------------------------------------------
figure
h1 = plot(fc, pow_bgfx, 'LineWidth', 1);
hold on
h2 = plot(fc, pow_mtnu0, 'LineWidth', 1);
xlim([fc_min, fc_max])
ylim([-10, 30])
grid on
grid minor
legend([h1, h2], {'BGFixed', 'MTNUFFT0'})
xlabel('Frequency (cycles/day)')
ylabel('Power density (dB \times day/cycle)')
title('Bronez-fixed multitaper power spectrum of nonuniform signal')

%% plot spectrum and F-test of uniformly sampled signal
% -----------------------------------------------------
figure
yyaxis left
plot(f_u, Fval_u)
xlim([fc_min fc_max])
hold on
plot([0 fs / 2], [sig1_u sig1_u], '--', 'LineWidth', .5)
plot([0 fs / 2], [sig2_u sig2_u], '--', 'LineWidth', .5)
plot([0 fs / 2], [sig_rot_u sig_rot_u], 'LineWidth', 1)
hold off
xlabel('Frequency (cycles/day)')
ylabel('F-statistic')

yyaxis right
h1 = plot(fc, pow_bgfx, 'g-', 'LineWidth', 1);
hold on
h2 = plot(f_u, pow2db(pxx_mtc), 'LineStyle', '-', 'LineWidth', 1);
xlim([fc_min fc_max])
ylim([-20, 30])
grid on
grid minor
legend([h1, h2], {'BGFixed', 'MTSPECTRUMC'})
ylabel('Power density (dB \times day/cycle)')
title('F-test and MTSPECTRUMC of resampled signal')

%% plot spectrum and F-test of nonuniformly sampled signal
% --------------------------------------------------------
figure
yyaxis left
plot(fc, Fval)
xlim([fc_min fc_max])
hold on
plot([0 fs / 2], [sig1 sig1], '--', 'LineWidth', .5)
plot([0 fs / 2], [sig2 sig2], '--', 'LineWidth', .5)
plot([0 fs / 2], [sig_rot sig_rot], 'LineWidth', 1)
hold off
xlabel('Frequency (cycles/day)')
ylabel('F-statistic')

yyaxis right
h1 = plot(fc, pow_bgfx, 'g-', 'LineWidth', 1);
hold on
h2 = plot(fc, pow2db(pxx_mtnu), 'LineStyle', '-', 'LineWidth', 1);
xlim([fc_min fc_max])
ylim([-10, 30])
grid on
grid minor
legend([h1, h2], {'BGFixed', 'MTNUFFT'})
title('F-test and MTNUFFT of nonuniformly sampled signal')

%% 2024 Richard J. Cui
% Created: Mon 07/01/2024 12:28:10.331 PM
% Revision: 0.1  Date: Mon 07/01/2024 12:28:10.331 PM
%
% Rocky Creek Dr. NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% [EOF]
