function [freq_f2, coeffs, coeffs_n, phi_dp] = Analyze_DPcomponents(stim, windowdur, offsetwin, npoints)

% Analyze swept tone DPOAE data using least-squares fit of chirp model

% Abdala et al., 2015: Optimizing swept-tone protocols for recording
% distortion-product otoacoustic emissions in adults and newborns

npoints = 512;

%% Set variables from the stim
phi1_inst = 2 * pi * stim.phi1_inst;
phi2_inst = 2 * pi * stim.phi2_inst;
phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
t = stim.t;

if stim.speed < 0 % downsweep
    f_start = stim.fmax;
    f_end = stim.fmin;
else
    f_start = stim.fmin;
    f_end = stim.fmax;
end


trials = stim.resp;
DPOAE = stim.DPOAE; 
NOISE = stim.NOISE; 

% %% Artifact Rejection
% energy = squeeze(sum(trials.^2, 2)); % same cut off for both trial types
% good = energy < median(energy) + 2*mad(energy);
% 
% count = 0;
% trials_clean = zeros(sum(good), size(trials, 2));
% for y = 1:size(trials, 1)
%     if good(y) == 1
%         count = count +1;
%         trials_clean(count, :) = trials(y,:);
%     end
% end
% DPOAE = mean(trials_clean, 1);
% 
% count_2x = floor(count/2)*2;
% noise = zeros(count_2x, size(trials, 2));
% count = 0;
% for x = 1:2:count_2x
%     count = count + 1;
%     noise(count,:) = (trials_clean(x,:) - trials_clean(x+1,:)) / 2;
% end
% NOISE = mean(noise,1);

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
    phi_dp = tau_dp.*freq_dp'; % cycles
end
