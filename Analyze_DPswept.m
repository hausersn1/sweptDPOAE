% Analyze swept tone DPOAE data using least-squares fit of chirp model

% Abdala et al., 2014: Distortion-product otoacoustic emission 
% reflection-component delays and cochlear tuning: Estimates from 
% across the human lifespan

windowdur = 0.031; 
offsetwin = 0.02; 
npoints = floor(max(stim.t)/(windowdur/20)); 


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

trials = stim.resp; 

%% Artifact Rejection 
energy = squeeze(sum(trials.^2, 2)); % same cut off for both trial types
good = energy < median(energy) + 2*mad(energy);

count = 0; 
trials_clean = zeros(sum(good), size(trials, 2)); 
for y = 1:size(trials, 1) 
    if good(y) == 1 
        count = count +1; 
        trials_clean(count, :) = trials(y,:); 
    end
end
DPOAE = mean(trials_clean, 1); 

count_2x = floor(count/2)*2; 
noise = zeros(count_2x, size(trials, 2)); 
count = 0; 
for x = 1:2:count_2x
    count = count + 1; 
    noise(count,:) = (trials_clean(x,:) - trials_clean(x+1,:)) / 2; 
end
NOISE = mean(noise,1); 

%% Set up for analysis 
% set freq we're testing and the timepoints when they happen.    
freq_f2 = linspace(f_start, f_end, npoints); 
freq_f1 = freq_f2 ./ stim.ratio; 
freq_dp = 2.*freq_f1 - freq_f2; 
t_freq = (freq_f2-f_start)/stim.speed + stim.buffdur; 

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2); 
coeffs_n = zeros(npoints, 2); 
tau_dp = zeros(npoints, 1); 
tau_f1 = zeros(npoints, 1); 
tau_f2 = zeros(npoints, 1); 

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
    
        % for model_f1
        coeff(offset + 1, 3:4) = model_f1' \ resp';
        coeff_n(offset + 1, 3:4) = model_f1' \ resp_n';
        resid(offset + 1, 2) = sum( (resp  - coeff(offset + 1, 3:4) * model_f1).^2);   
    
        % for model_f2
        coeff(offset + 1, 5:6) = model_f2' \ resp';
        coeff_n(offset + 1, 5:6) = model_f2' \ resp_n';
        resid(offset + 1, 3) = sum( (resp  - coeff(offset + 1, 5:6) * model_f2).^2);   
    
    end
    
    ind = zeros(3,1); 
    for i = 1:3
        [~, ind(i)] = min(resid(:,i));
        coeffs(k, (i*2-1):(i*2)) = coeff(ind(i), (i*2-1):(i*2));
        coeffs_n(k,(i*2-1):(i*2)) = coeff_n(ind(i),(i*2-1):(i*2)); 
    end
    tau_dp(k) = (ind(1) - 1) * 1/stim.Fs; % delay in sec
    tau_f1(k) = (ind(2) - 1) * 1/stim.Fs; % delay in sec
    tau_f2(k) = (ind(3) - 1) * 1/stim.Fs; % delay in sec
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

f_x_dp = (freq_dp(2:end) + freq_dp(1:end-1))/2; 
f_x_f1 = (freq_f1(2:end) + freq_f1(1:end-1))/2; 
f_x_f2 = (freq_f2(2:end) + freq_f2(1:end-1))/2; 

% for f1
phi_f1 = tau_f1.*freq_f1'; % cycles
phasor_f1 = exp(-1j * phi_f1 * 2 * pi); 
mag_f1 = abs(complex(a_f1, b_f1).*phasor_f1) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_f1 = unwrap(angle(complex(a_f1,b_f1).*phasor_f1))/(2*pi); % cycles; 
tau_pg_f1 = -diff(theta_f1)./diff(freq_f1./1000); 

% for f2
phi_f2 = tau_f2.*freq_f2'; % cycles
phasor_f2 = exp(-1j * phi_f2 * 2 * pi); 
mag_f2 = abs(complex(a_f2, b_f2).*phasor_f2) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_f2 = unwrap(angle(complex(a_f2,b_f2).*phasor_f2))/(2*pi); % cycles; 
tau_pg_f2 = -diff(theta_f2)./diff(freq_f2./1000); 

phi_correction = 2.*phi_f1 - phi_f2; 
% for dp
phi_dp = tau_dp.*freq_dp'; %- phi_correction; % cycles
phasor_dp = exp(-1j * phi_dp * 2 * pi); 
mag_dp = abs(complex(a_dp, b_dp).*phasor_dp) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_dp = unwrap(angle(complex(a_dp,b_dp).*phasor_dp))/(2*pi); % cycles; 
tau_pg_dp = -diff(theta_dp)./diff(freq_dp./1000); %millisec

% for nf
mag_nf = abs(complex(a_n, b_n).*phasor_dp) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;

%% Separating D and R components


% maxfreq = max(freq_dp); 
% minfreq = min(freq_dp); 
% testfreq = linspace(floor(minfreq), floor(maxfreq), floor(((maxfreq)-minfreq)/25)); 
% % 
% for m = 1:size(testfreq, 2)-1
%     
%     f_win = find( testfreq(m) < freq_dp & ...
%         freq_dp < testfreq(m)+300);
%     f_taper = hanning(f_win)'; 
% 
%     dp = mag_dp(f_win); 
%     
%     N = 2.*nextpow2(2*size(dp, 1)); 
%     timewave = ifft(mag_dp, N); 
%     t = 1000.*linspace(0,N/stim.Fs, N); 
%     figure(10); plot(t, abs(timewave))
% 
% end

% % %
% % X = complex(a_f2, b_f2).*phasor_f2; 
% % N=2.^16; 
% % timewave = ifft(X, N); 
% % t = 1000.*linspace(0,N/stim.Fs, N); 
% % plot(abs(timewave))
% % 
% % window = find(t < 5, 1, 'last'); 
% % winend = find(t<10,1, 'last'); 
% % dist = abs(fft(timewave(1:window))); 
% % figure(11); plot(dist)
% % refl = abs(fft(timewave(window+1:winend))); 
% % figure(11); hold on; plot(refl)


%% Plot figures
figure(1); 
hold on; 
semilogx(freq_dp, db(mag_dp), 'linew', 2); 
hold on; 
semilogx(freq_f1, db(mag_f1), 'linew', 2); 
hold on; 
semilogx(freq_f2, db(mag_f2), 'linew', 2); 
hold on; 
semilogx(freq_dp, db(mag_nf),'--', 'linew', 2); 
xlim([stim.fmin * rdp, stim.fmax * rdp]);
xlabel('DPOAE Frequency (Hz)', 'FontSize', 16);
ylabel('DPOAE level (dB SPL)', 'FontSize', 16);
set(gca, 'FontSize', 16);
legend('dp', 'f1', 'f2', 'nf'); 

figure(2); 
hold on; 
semilogx(freq_dp, theta_dp-max(theta_dp), 'o')
hold on; 
semilogx(freq_dp, theta_f1-max(theta_f1),'o');
hold on; 
semilogx(freq_dp, theta_f2-max(theta_f2), 'o');
title('Theta, angle(a,b)')
xlabel('Probe Frequency (Hz)', 'FontSize', 16); 
ylabel('Phase (cycles)', 'FontSize', 16); 

% 
% figure(3); 
% semilogx(f_x, tau_pg_dp, 's', testfreq, (tau_dp.*1000), 'or');
% title('Group Delay, pg from slope, tau from offset');
% xlabel('Frequency (Hz)', 'FontSize', 16); 
% ylabel('Group Delay (ms)', 'FontSize', 16); 
% ylim([-10, 30])
% xlim([stim.fmin, stim.fmax])
% legend('tau_p_g (-d\theta/df)', 'tau_n (offset)')
