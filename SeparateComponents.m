%% Run analysis to separate D and R components of DPOAE

% requires load of data file with struct called stim. 

%% Run Total DPOAE
% for the total DPOAE // short window
windowdur = 0.04; % seconds
offsetwin = 0.010;
[freq_f2, coeffs_tot, coeffs_n_tot, phi_dp_tot] = Analyze_DPcomponents(stim, windowdur, offsetwin); 

freq_f1 = freq_f2 ./ stim.ratio;
freq_dp = 2.*freq_f1 - freq_f2;

a_dp_tot = coeffs_tot(:, 1);
b_dp_tot = coeffs_tot(:, 2);
a_f1_tot = coeffs_tot(:, 3);
b_f1_tot = coeffs_tot(:, 4);
a_f2_tot = coeffs_tot(:, 5);
b_f2_tot = coeffs_tot(:, 6);
a_n_tot = coeffs_n_tot(:, 1);
b_n_tot = coeffs_n_tot(:, 2);

% Magnitudes of total
% for stimuli
mag_f1 = abs(complex(a_f1_tot, b_f1_tot)) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
mag_f2 = abs(complex(a_f2_tot, b_f2_tot)) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;

% for dp
phasor_dp_tot = exp(-1j * phi_dp_tot * 2 * pi);
mag_dp_tot = abs(complex(a_dp_tot, b_dp_tot).*phasor_dp_tot) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_dp_tot = unwrap(angle(complex(a_dp_tot,b_dp_tot).*phasor_dp_tot))/(2*pi); % cycles;

% for nf
mag_nf = abs(complex(a_n_tot, b_n_tot).*phasor_dp_tot) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;

%% Distortion component
windowdur = 0.2; 
offsetwin = 0.0;

% For the distortion component // long window
[~, coeffs_D, ~, ~] = Analyze_DPcomponents(stim, windowdur, offsetwin); 

a_dp_D = coeffs_D(:, 1);
b_dp_D = coeffs_D(:, 2);

% Magnitudes of D only
% for dp
mag_dp_D = abs(complex(a_dp_D, b_dp_D)) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_dp_D = unwrap(angle(complex(a_dp_D,b_dp_D)))/(2*pi); % cycles;

%% Calculate Reflection by subtracting the Distortion
a_dp_R = a_dp_tot - a_dp_D; 
b_dp_R = b_dp_tot - b_dp_D; 

% Magnitudes of R only
% for dp
mag_dp_R = abs(complex(a_dp_R, b_dp_R).*phasor_dp_tot) .* stim.VoltageToPascal .*stim.PascalToLinearSPL;
theta_dp_R = unwrap(angle(complex(a_dp_R,b_dp_R).*phasor_dp_tot))/(2*pi); % cycles;
delay = -diff(theta_dp_R)./diff(freq_dp); 

% %% Saving to Struct
% 
% res.f2 = freq_f2; 
% res.f1 = freq_f1; 
% res.dp = freq_dp; 
% res.tot.dp.mag = mag_dp_tot; 
% res.mag.dp.R = mag_dp_R; 
% res.mag.dp.D

%% Plot figures
figure(1);
hold on;
semilogx(freq_f2, db(mag_dp_tot), 'linew', 2);
hold on; 
semilogx(freq_f2, db(mag_dp_D), 'linew', 2); 
hold on;
semilogx(freq_f1, db(mag_f1), 'linew', 2);
hold on;
semilogx(freq_f2, db(mag_f2), 'linew', 2);
hold on;
semilogx(freq_f2, db(mag_nf),'--', 'linew', 2);
xlim([stim.fmin, stim.fmax]);
xlabel('DPOAE Frequency (Hz)', 'FontSize', 16);
ylabel('DPOAE level (dB SPL)', 'FontSize', 16);
set(gca, 'FontSize', 16);
legend('Total DP', 'D only', 'f1', 'f2', 'nf');

figure(2); 
hold on; 
semilogx(freq_dp, db(mag_dp_R), 'linew', 2); 
hold on; 
semilogx(freq_dp, db(mag_nf)+3, 'linew', 2); 
title('Reflection Component');
legend('R component', 'NF')

figure(3);
hold on; 
semilogx(freq_dp, theta_dp_tot-max(theta_dp_tot),'o');
hold on;
semilogx(freq_dp, theta_dp_D-max(theta_dp_D), 's')
hold on;
semilogx(freq_dp, theta_dp_R-max(theta_dp_R), '+')
hold on; 
title('Theta, angle(a,b)')
xlabel('Probe Frequency (Hz)', 'FontSize', 16);
ylabel('Phase (cycles)', 'FontSize', 16);
legend('Total', 'Distortion', 'Reflection'); 

figure(4); 
loglog(freq_dp(1:end-1), delay.*freq_dp(1:end-1)', 'o'); 

