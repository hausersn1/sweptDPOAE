% Analyze swept tone DPOAE data using least-squares fit of chirp model

% Abdala et al., 2015: Optimizing swept-tone protocols for recording
% distortion-product otoacoustic emissions in adults and newborns

windowdur = 0.25;
offsetwin = 0.0; % not finding additional delay
npoints = 512;

%% Set variables from the stim
phi1_inst = 2 * pi * stim.phi1_inst;
phi2_inst = 2 * pi * stim.phi2_inst;
phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
t = stim.t;
rdp = 2 / stim.ratio - 1; % f_dp = f2 * rdp

if stim.speed < 0 % downsweep
    f_start = stim.fmax;
    f_end = stim.fmin;
else
    f_start = stim.fmin;
    f_end = stim.fmax;
end


%% Artifact Rejection
trials = stim.resp;
energy = squeeze(sum(trials.^2, 2));
good = energy < median(energy) + 2*mad(energy);

count = 0;
trials_clean = zeros(sum(good), size(trials, 2));
for y = 1:size(trials, 1) % just get trials w/o artifact ("good" trials)
    if good(y) == 1
        count = count +1;
        trials_clean(count, :) = trials(y,:);
    end
end
DPOAE = mean(trials_clean, 1);

count_2x = floor(count/2)*2; % need even number of trials
noise = zeros(count_2x, size(trials, 2));
count = 0;
for x = 1:2:count_2x
    count = count + 1;
    noise(count,:) = (trials_clean(x,:) - trials_clean(x+1,:)) / 2;
end
NOISE = mean(noise,1);


%% Set up for analysis
% set freq we're testing and the timepoints when they happen.


if stim.speed < 20
    freq_f2 = 2 .^ linspace(log2(f_start), log2(f_end), npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = log2(freq_f2/f_start)/stim.speed + stim.buffdur;
else
    freq_f2 = linspace(f_start, f_end, npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = (freq_f2-f_start)/stim.speed + stim.buffdur;
end



% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
coeffs_n = zeros(npoints, 2);
tau_dp = zeros(npoints, 1); % delay if offset > 0

%% Least Squares fit of Chirp model
for k = 1:npoints
    fprintf(1, 'Running window %d / %d\n', k, npoints);
    
    win = find( (t > (t_freq(k) - windowdur/2)) & ...
        (t < (t_freq(k) + windowdur/2)));
    taper = hanning(numel(win))';
    
    model_dp = [cos(phi_dp_inst(win)) .* taper;
        -sin(phi_dp_inst(win)) .* taper];
    
    model_f1 = [cos(phi1_inst(win)) .* taper;
        -sin(phi1_inst(win)) .* taper];
    
    model_f2 = [cos(phi2_inst(win)) .* taper;
        -sin(phi2_inst(win)) .* taper];
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 6);
    coeff_n = zeros(maxoffset, 6);
    resid = zeros(maxoffset, 3);
    
    for offset = 0:maxoffset
        resp = DPOAE(win+offset) .* taper;
        resp_n = NOISE(win+offset) .* taper;
        
        % for model_dp
        coeff(offset + 1, 1:2) = model_dp' \ resp';
        coeff_n(offset + 1, 1:2) = model_dp' \ resp_n';
        resid(offset + 1, 1) = sum( (resp  - coeff(offset + 1, 1:2) * model_dp).^2);
        
    end
    resp = DPOAE(win) .* taper;
    resp_n = NOISE(win) .* taper;
    
    % for model_f1
    coeff(1, 3:4) = model_f1' \ resp';
    coeff_n(1, 3:4) = model_f1' \ resp_n';
    resid(1, 2) = sum( (resp  - coeff(1, 3:4) * model_f1).^2);
    
    % for model_f2
    coeff(1, 5:6) = model_f2' \ resp';
    coeff_n(1, 5:6) = model_f2' \ resp_n';
    resid(1, 3) = sum( (resp  - coeff(1, 5:6) * model_f2).^2);
    
    [~, ind] = min(resid(:,1));
    coeffs(k, 1:2) = coeff(ind, 1:2);
    coeffs_n(k, 1:2) = coeff_n(ind, 1:2);
    coeffs(k, 3:6) = coeff(1,3:6);
    
    tau_dp(k) = (ind(1) - 1) * 1/stim.Fs; % delay in sec
end

%% Amplitude and Delay calculations
a_dp = coeffs(:, 1);
b_dp = coeffs(:, 2);
a_f1 = coeffs(:, 3);
b_f1 = coeffs(:, 4);
a_f2 = coeffs(:, 5);
b_f2 = coeffs(:, 6);
a_n = coeffs_n(:, 1);
b_n = coeffs_n(:, 2);

% for f1
mag_f1 = abs(complex(a_f1, b_f1)) .* stim.VoltageToPascal .*stim.PascalToLinearSPL; %dBSPL
theta_f1 = unwrap(angle(complex(a_f1,b_f1)))/(2*pi); % cycles;
tau_pg_f1 = -(diff(theta_f1)./diff(freq_f1))/1000; % ms

% for f2
mag_f2 = abs(complex(a_f2, b_f2)) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_f2 = unwrap(angle(complex(a_f2,b_f2)))/(2*pi);
tau_pg_f2 = -(diff(theta_f2)./diff(freq_f2))/1000;

% for dp
phi_dp = tau_dp.*freq_dp'; % cycles (from delay/offset)
phasor_dp = exp(-1j * phi_dp * 2 * pi);
mag_dp = abs(complex(a_dp, b_dp).*phasor_dp) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_dp = unwrap(angle(complex(a_dp,b_dp).*phasor_dp))/(2*pi); % cycles from angle
theta_dp = theta_dp-max(theta_dp); % starting at zero
tau_pg_dp = -(diff(theta_dp)./diff(freq_dp')).*1000; % ms
f_x_dp = (freq_dp(2:end) + freq_dp(1:end-1))/2;

% for nf
mag_nf = abs(complex(a_n, b_n).*phasor_dp) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;

%% Separating D and R components by IFFT

complex_dp = complex(a_dp, b_dp).*phasor_dp * stim.VoltageToPascal * ...
    stim.PascalToLinearSPL;

fs = stim.Fs; % Exact value doesn't matter, but simplest to match orig

% Reconstruct time-domain response to 100 ms. Then the first 50 ms can be
% used, with the recognition that part of the right half is negative time.
% 50 ms should still be plenty for all OAE components to have come back
% and the impulse response to have decayed to noise floor.
timedur = 100e-3;
% Check if 1/(frequency resolution) is long enough
if 1/mean(abs(diff(freq_dp))) < timedur/2
    warning('Too few DPOAE points, aliasing  will likely  occur');
end
N = roundodd(timedur * fs);
f = (0:(N-1))*fs/N; % FFT bin frequencies

% Window the frequency domain response to avoid sharp edges while filling
% in zeros for empty bins
rampsamps = 8;
ramp = hanning(2*rampsamps);
complex_dp_ramp  = complex_dp;
complex_dp_ramp(1:rampsamps) = complex_dp_ramp(1:rampsamps)...
    .*ramp(1:rampsamps);
complex_dp_ramp((end-rampsamps+1):end) = ...
    complex_dp_ramp((end-rampsamps+1):end).* ramp((end-rampsamps+1):end);

% Fill in FFT bins from dp data
FFT_dp =  interp1(freq_dp, complex_dp_ramp, f);
FFT_dp(isnan(FFT_dp)) = 0;

% First bin is f=0, then you have even number of bins
% Need to take first half and conjugate mirror to fill other.
Nhalf = (N-1)/2;
nonzeroHalf =  FFT_dp(2:(2+Nhalf-1));
FFT_dp((Nhalf+2):end) = conj(nonzeroHalf(end:-1:1));

impulse_dp = ifftshift(ifft(FFT_dp));

% Make time vector: t=0 will be at the center owing to ifftshift
t = (0:(N-1))/fs - timedur/2; % in milleseconds

% Start a bit negative because D component has close to 0 delay and go to
% half of the reconstructed duration (50 ms) as planned
t_min = -4e-3;
t_max = timedur*1e3/2;
inds_valid = t > t_min & t < t_max;
impulse_dp = impulse_dp(inds_valid);
t_dp = t(inds_valid);

% Plot envelope of impulse response to see if there are two peaks, with the
% notch between peaks being somewhere in the 1-5 ms range.
figure(40);
impulse_dp_env = abs(hilbert(impulse_dp));
plot(t_dp*1e3,  impulse_dp_env, 'linew', 2);
xlim([t_min*1e3, 20]);
xlabel('Time (ms)', 'FontSize', 16);
ylabel('Impulse Response Envelope', 'FontSize', 16);
set(gca, 'FontSize', 16);
title('IFFT method',  'FontSize', 16);

% Location of one of the notches..
D_only_dur = 3.5e-3;

% Do windowing of signals
% Start with  hard windows (box) and then smooth edges by 0.5 ms
smoothing_kernel = blackman(ceil(0.5e-3*fs));
smoothing_kernel  = smoothing_kernel / sum(smoothing_kernel);

% Error using conv2
% First and second arguments must be single or double.
% original code (which works elsewhere?)
% win_D_only = conv(t_dp < D_only_dur, smoothing_kernel, 'same');
% win_R_only = conv(t_dp > D_only_dur, smoothing_kernel, 'same');
% trying the following to solve problem:
t_dp_D_only = t_dp.*(t_dp < D_only_dur);
t_dp_R_only = t_dp.*(t_dp > D_only_dur);
win_D_only = conv(t_dp_D_only, smoothing_kernel, 'same');
win_R_only = conv(t_dp_R_only, smoothing_kernel, 'same');

impulse_dp_D_only = impulse_dp .* win_D_only;
impulse_dp_R_only = impulse_dp .* win_R_only;

complex_dp_D_only_allbins = fft(impulse_dp_D_only);
complex_dp_R_only_allbins = fft(impulse_dp_R_only);
f_complex_dp_D_only_allbins = (0:(numel(impulse_dp_D_only)-1)) ...
    * fs/numel(impulse_dp_D_only);
phasor_tmin_correction = ...
    exp(1j*2*pi*f_complex_dp_D_only_allbins*abs(t_min));
complex_dp_D_only_allbins_corrected = complex_dp_D_only_allbins ...
    .* phasor_tmin_correction;
complex_dp_R_only_allbins_corrected = complex_dp_R_only_allbins ...
    .* phasor_tmin_correction;
complex_dp_D_IFFT = interp1(f_complex_dp_D_only_allbins,...
    complex_dp_D_only_allbins_corrected, freq_dp);
complex_dp_R_IFFT = interp1(f_complex_dp_D_only_allbins,...
    complex_dp_R_only_allbins_corrected, freq_dp);

%% Plot ifft Method
figure(41);
hold on;
plot(freq_dp, db(abs(complex_dp)), 'linew', 2);
hold on;
semilogx(freq_dp, db(abs(complex_dp_D_IFFT)), 'linew', 2);
hold on;
semilogx(freq_dp, db(abs(complex_dp_R_IFFT)), 'linew', 2);
xlabel('DPOAE Frequency (Hz)', 'FontSize', 16);
ylabel('DPOAE level (dB SPL)', 'FontSize', 16);
set(gca, 'FontSize', 16);
legend('Total', 'D only', 'R only');
ylim([-60, 80]);
title('Separating by IFFT',  'FontSize', 16);


%% Plot figures
figure(10);
hold on;
semilogx(freq_f2, db(mag_dp), 'linew', 2);
hold on;
semilogx(freq_f1, db(mag_f1), 'linew', 2);
hold on;
semilogx(freq_f2, db(mag_f2), 'linew', 2);
hold on;
semilogx(freq_f2, db(mag_nf),'--', 'linew', 2);
xlim([stim.fmin * rdp, stim.fmax * rdp]);
xlabel('Frequency (Hz)', 'FontSize', 16);
ylabel('DPOAE level (dB SPL)', 'FontSize', 16);
set(gca, 'FontSize', 16);
legend('DP', 'F1', 'F2', 'NF');
title('Total DPOAE',  'FontSize', 16);
xticks([500, 1000, 2000, 4000, 8000, 16000])

figure(30);
semilogx(freq_dp, theta_dp, 'o')
set(gca, 'FontSize', 16);
title('Theta, angle(a,b)', 'FontSize', 16)
xlabel('Probe Frequency (Hz)', 'FontSize', 16);
ylabel('Phase (cycles)', 'FontSize', 16);
xticks([500, 1000, 2000, 4000, 8000, 16000])

figure(20);
semilogx(f_x_dp, tau_pg_dp, 's');
title('Group Delay (-d\theta/df)', 'FontSize', 16);
xlabel('2F_1-F_2 Frequency (Hz)', 'FontSize', 16);
ylabel('Group Delay (ms)', 'FontSize', 16);
set(gca, 'FontSize', 16);
ylim([-10, 30])
xlim([stim.fmin*rdp, stim.fmax*rdp])
xticks([500, 1000, 2000, 4000, 8000, 16000])

figure(21);
loglog(f_x_dp, tau_pg_dp/1000.*f_x_dp', 's');
title('N - for calculating Qerb', 'FontSize', 16);
xlabel('2F_1-F_2 Frequency (Hz)', 'FontSize', 16);
ylabel('N_{DPOAE}', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlim([stim.fmin*rdp, stim.fmax*rdp])
xticks([500, 1000, 2000, 4000, 8000, 16000])

