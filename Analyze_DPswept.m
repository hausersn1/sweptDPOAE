% DPOAE swept Analysis
% Author: Samantha Hauser
% Created: May 2023
% Last Updated: Sept 16, 2023
% Purpose:
% Helpful info: Ensure a datafile is loaded

% For quick visualization after data collection.

%%%%%%%%% Set these parameters %%%%%%%%%%%%%%%%%%

windowdur = 0.25; %0.25;
offsetwin = 0.0; % not finding additional delay
npoints = 512;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to be compatible with old data
if exist('data', 'var')
    stim = data.stim;
    subj = data.info.subj.ID;
    ear = data.info.subj.ear;
    DPOAEtrials = data.resp.AllBuffs;
    trials = data.resp.trialsCollected;
else
    subj = stim.subj;
    ear = stim.ear;
    DPOAEtrials = stim.resp;
    trials = size(stim.resp,1);
end

%% Set variables from the stim
phi1_inst = 2 * pi * stim.phi1_inst;
phi2_inst = 2 * pi * stim.phi2_inst;
phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
rdp = 2 / stim.ratio - 1;    % f_dp = f2 * rdp

t = stim.t;

if stim.speed < 0 % downsweep
    f_start = stim.fmax;
    f_end = stim.fmin;
else
    f_start = stim.fmin;
    f_end = stim.fmax;
end

% set freq we're testing and the timepoints when they happen.
if strcmp(stim.scale, 'log')        % in octave scaling
    freq_f2 = 2 .^ linspace(log2(f_start), log2(f_end), npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = log2(freq_f2/f_start)/stim.speed + stim.buffdur;
else                            % otherwise linear scaling
    freq_f2 = linspace(f_start, f_end, npoints);
    freq_f1 = freq_f2 ./ stim.ratio;
    freq_dp = 2.*freq_f1 - freq_f2;
    t_freq = (freq_f2-f_start)/stim.speed + stim.buffdur;
end

nfreqs = stim.nearfreqs;

%% Artifact Rejection
% Set empty matricies for next steps
coeffs = zeros(npoints, 2);
a_temp = zeros(trials, npoints);
b_temp = zeros(trials, npoints);

% Least Squares fit of DP after artifact rejection
for x = 1:trials
    DPOAE = DPOAEtrials(x, :);
    fprintf(1, 'Checking trial %d / %d for artifact\n', x, trials);
    
    for k = 1:npoints
        win = find( (t > (t_freq(k) - windowdur/2)) & ...
            (t < (t_freq(k) + windowdur/2)));
        taper = hanning(numel(win))';
        
        model_dp = [cos(phi_dp_inst(win)) .* taper;
            -sin(phi_dp_inst(win)) .* taper];
        
        resp = DPOAE(win) .* taper;
        
        coeffs(k, 1:2) = model_dp' \ resp';
    end
    
    a_temp(x,:) = coeffs(:, 1);
    b_temp(x,:) = coeffs(:, 2);
end

oae = abs(complex(a_temp, b_temp));
median_oae = median(oae);
std_oae = std(oae);

resp_AR = DPOAEtrials;
for j = 1:trials
    for k = 1:npoints
        if oae(j,k) > median_oae(1,k) + 3*std_oae(1,k)
            win = find( (t > (t_freq(k) - windowdur.*.1)) & ...
                (t < (t_freq(k) + windowdur.*.1)));
            resp_AR(j,win) = NaN;
        end
    end
end

DPOAE = mean(resp_AR, 1, 'omitNaN');

%% LSF analysis

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2);
tau_dp = zeros(npoints, 1); % delay if offset > 0
coeffs_noise = zeros(npoints,8);

% if duration changes with frequency
%durs = -.5*(2.^(0.003*t_freq)-1)/ (0.003*log(2)) + 0.5;

% Least Squares fit of Chirp model (stimuli, DP, noise)
for k = 1:npoints
    
    fprintf(1, 'Running window %d / %d\n', k, npoints);
    
    % if using durs: windowdur = durs(k);
    win = find( (t > (t_freq(k) - windowdur/2)) & ...
        (t < (t_freq(k) + windowdur/2)));
    taper = hanning(numel(win))';
    
    % set the response
    resp = DPOAE(win) .* taper;
    
    % DP Coeffs with variable delay calculation
    model_dp = [cos(phi_dp_inst(win)) .* taper;
        -sin(phi_dp_inst(win)) .* taper];
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 6);
    coeff_n = zeros(maxoffset, 6);
    resid = zeros(maxoffset, 3);
    for offset = 0:maxoffset
        resp = DPOAE(win+offset) .* taper;
        coeff(offset + 1, 1:2) = model_dp' \ resp';
        resid(offset + 1, 1) = sum( (resp  - coeff(offset + 1, 1:2) * model_dp).^2);
    end
    [~, ind] = min(resid(:,1));
    coeffs(k, 1:2) = coeff(ind, 1:2);
    % Calculate delay
    tau_dp(k) = (ind(1) - 1) * 1/stim.Fs; % delay in sec

    % F1 Coeffs
    model_f1 = [cos(phi1_inst(win)) .* taper;
        -sin(phi1_inst(win)) .* taper];
    coeffs(k, 3:4) = model_f1' \ resp';
    
    % F2 Coeffs
    model_f2 = [cos(phi2_inst(win)) .* taper;
        -sin(phi2_inst(win)) .* taper];
    coeffs(k, 5:6) = model_f2' \ resp';
    
    % Noise Coeffs
    model_noise = ...
        [cos(nfreqs(1)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(1)*phi_dp_inst(win)) .* taper;
        cos(nfreqs(2)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(2)*phi_dp_inst(win)) .* taper;
        cos(nfreqs(3)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(3)*phi_dp_inst(win)) .* taper;
        cos(nfreqs(4)*phi_dp_inst(win)) .* taper;
        -sin(nfreqs(4)*phi_dp_inst(win)) .* taper];
    coeffs_noise(k,:) = model_noise' \ resp';
    
end

%% Amplitude and Delay calculations
a_dp = coeffs(:, 1);
b_dp = coeffs(:, 2);
a_f1 = coeffs(:, 3);
b_f1 = coeffs(:, 4);
a_f2 = coeffs(:, 5);
b_f2 = coeffs(:, 6);

% complex DPOAE
oae_complex = complex(a_dp, b_dp);

% complex average noise
noise = zeros(npoints,4);
for i = 1:2:8
    noise(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
end
noise_complex = mean(noise,2);

% delay
phi_dp = tau_dp.*freq_dp'; % cycles (from delay/offset)
phasor_dp = exp(-1j * phi_dp * 2 * pi);

VtoSPL = stim.VoltageToPascal .* stim.PascalToLinearSPL;

%% Plot Results Figure
figure;
plot(freq_f2/1000, db(abs(oae_complex).*VtoSPL), 'linew', 2, 'Color', [0 0.4470 0.7410]);
hold on;
plot(freq_f2/1000, db(abs(noise_complex).*VtoSPL), '--', 'linew', 2, 'Color', [0.6350 0.0780 0.1840]);
plot(freq_f2/1000, db(abs(complex(a_f2,b_f2)).*VtoSPL), 'linew', 2, 'Color', [0.4940 0.1840 0.5560]);
plot(freq_f1/1000, db(abs(complex(a_f1, b_f1)).*VtoSPL), 'linew', 2, 'Color', [0.9290 0.6940 0.1250]);
title('DPOAE', 'FontSize', 14)
set(gca, 'XScale', 'log', 'FontSize', 14)
xlim([.5, 16])
ylim([-50, 90])
xticks([.5, 1, 2, 4, 8, 16])
ylabel('Amplitude (dB SPL)', 'FontWeight', 'bold')
xlabel('F2 Frequency (kHz)', 'FontWeight', 'bold')
legend('OAE', 'NF', 'F2', 'F1')
drawnow;

%% TODO: Add creation of summary points and conversion to EPL