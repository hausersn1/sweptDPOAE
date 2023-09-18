%% Swept DPOAEs for Humans (ER-10X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Samantha Hauser, AuD
% Modified from: Hari Bharadwaj, PhD (SNAP Lab)
% Created: November 2021
% Last revision: 6-Sep-2023 (added artifact checking)
%
% References:
%   SNR based endpoint: Abdala, C., Luo, P. & Shera, C.A. Characterizing
%       the Relationship Between Reflection and Distortion Otoacoustic 
%       Emissions in Normal-Hearing Adults. JARO 23, 647-664 (2022). 
%       https://doi.org/10.1007/s10162-022-00857-z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up data storage and subject info

% Measure-General info
info.measure = 'DPOAEswept';
info.version = 'Auto_v02';

% Visit info
if exist('C:\Experiments\Sam\current_visit.mat','file')
    load('C:\Experiments\Sam\current_visit.mat', 'visit')
    ask = questdlg(sprintf('Is this subject %s?', visit.subj.ID), ...
        'Check Subject', 'Yes', 'No', 'No');
else
    ask = 'No';
end

if strcmp(ask, 'No')
    cd ..
    startVisit
    cd(info.measure)
end

subj = visit.subj;
info.room = visit.room;
info.univ = visit.univ;
info.researcher = visit.researcher;

% Get ear info
subj.ear = questdlg('Which ear?', 'Ear', 'L', 'R', 'R');

% Get date/time
datetag = datestr(clock);
info.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '-';

% Make directory to save results
paraDir = 'C:\Experiments\Sam\DPOAEswept\Results\';
respDir = strcat(paraDir,filesep,visit.subj.ID,filesep);
addpath(genpath(paraDir));
if(~exist(respDir,'dir'))
    mkdir(respDir);
end

fname = strcat(respDir, info.measure, '_', ...
    subj.ID, '_', subj.ear, '_', datetag, '.mat');

%% Run Test
tic;

try
    
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    initializeER10X;
    
    % Initializing TDT and specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    % Get stimulus structure
    stim = Make_DPswept;
    
    % Set stims in buffdata:
    buffdata = [stim.y1; stim.y2];
    
    % Check for clipping and load to buffer
    if(any(abs(buffdata(:)) > 1))
        error('What did you do!? Sound is clipping!! Cannot Continue!!\n');
    end
    
    % Set attenuation
    drop_f1 = stim.drop_f1;
    drop_f2 = stim.drop_f2;
    delayComp = 1; % Always
    
    % Live analysis parameters
    windowdur = stim.windowdur;
    SNRcriterion = stim.SNRcriterion;
    maxTrials = stim.maxTrials;
    minTrials = stim.minTrials;
    
    phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi; %dp
    t = stim.t;
    npoints = stim.npoints; 
    nearfreqs = stim.nearfreqs;
    VtoSPL = stim.VoltageToPascal .* stim.PascalToLinearSPL;
    
    edges = 2 .^ linspace(log2(stim.fmin), log2(stim.fmax), 21);
    bandEdges = edges(2:2:end-1);
    centerFreqs = edges(3:2:end-2);
    
    if stim.speed < 0
        f_start = stim.fmax;
        f_end = stim.fmin;
    else
        f_start = stim.fmin;
        f_end = stim.fmax;
    end
    
    testfreq = 2 .^ linspace(log2(f_start), log2(f_end), npoints); 
        
    if strcmp(stim.scale, 'log')
        t_freq = log2(testfreq/f_start)/stim.speed + stim.buffdur;
    else
        t_freq = (testfreq-f_start)/stim.speed + stim.buffdur;
    end
    
    k = 0;
    doneWithTrials = 0;
    figure;
    
    % set matrix for storing response
    resp = zeros(maxTrials, size(buffdata,2));
    
    while doneWithTrials == 0
        k = k + 1;
        k_kept = k - stim.ThrowAway;
        
        %Start playing from the buffer:
        vins = playCapture2(buffdata, card, 1, 0,...
            drop_f1, drop_f2, delayComp);
        
        % save the response
        if k > stim.ThrowAway
            resp(k - stim.ThrowAway,  :) = vins;
        end
        
        WaitSecs(0.15);
        
        if k > stim.ThrowAway
            % Set empty matricies for next steps
            coeffs_resp = zeros(npoints, 2);
            coeffs_noise = zeros(npoints, 8);
            
            for p = 1:npoints
                win = find( (t > (t_freq(p) - windowdur/2)) & ...
                    (t < (t_freq(p) + windowdur/2)));
                taper = hanning(numel(win))';
                
                resp_trial = vins(win) .* taper; % just current trial
                
                model_dp = [cos(phi_dp_inst(win)) .* taper;
                    -sin(phi_dp_inst(win)) .* taper];
                model_noise = ...
                    [cos(nearfreqs(1)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(1)*phi_dp_inst(win)) .* taper;
                    cos(nearfreqs(2)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(2)*phi_dp_inst(win)) .* taper;
                    cos(nearfreqs(3)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(3)*phi_dp_inst(win)) .* taper;
                    cos(nearfreqs(4)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(4)*phi_dp_inst(win)) .* taper];
                
                coeffs_resp(p,:) = model_dp' \ resp_trial';
                coeffs_noise(p,:) = model_noise' \ resp_trial';
            end
            
            % calculate amplitudes
            oae_trials(k_kept,:) = abs(complex(coeffs_resp(:, 1),  coeffs_resp(:, 2)));
            median_oae = median(oae_trials,1);
            dpoae_full = db(median_oae.*VtoSPL); 
            
            noise_trial = zeros(npoints,4);
            for i = 1:2:8
                noise_trial(:,ceil(i/2)) = complex(coeffs_noise(:,i), coeffs_noise(:,i+1));
            end
            noise_trials(k_kept,:) = abs(mean(noise_trial, 2));
            median_noise = median(noise_trials,1); 
            dpnf_full = db(median_noise.*VtoSPL); 
        
            dpoae = zeros(length(centerFreqs),1);
            dpnf = zeros(length(centerFreqs),1);
            dpoae_w = zeros(length(centerFreqs),1);
            dpnf_w = zeros(length(centerFreqs),1);
            
            % weighted average around 9 center frequencies
            for z = 1:length(centerFreqs)
                band = find( testfreq >= bandEdges(z) & testfreq < bandEdges(z+1));
          
                % TO DO: NF from which SNR was calculated included median of 7 points
                % nearest the target frequency.
                SNR = dpoae_full(band) - dpnf_full(band);
                weight = (10.^(SNR./10)).^2;
                
                dpoae(z, 1) = mean(dpoae_full(band));
                dpnf(z,1) = mean(dpnf_full(band));
                
                dpoae_w(z,1) = sum(weight.*dpoae_full(band))/sum(weight);
                dpnf_w(z,1) = sum(weight.*dpnf_full(band))/sum(weight);
                
            end
            
            % median SNR
            SNR_temp = dpoae_w - dpnf_w;
            
            noisy_trials = 0;
            % artifact check
            if k_kept >= stim.minTrials
                std_oae = std(oae_trials,1);
                for r = 1:k_kept
                    for q = 1:npoints
                        if oae_trials(r,q) > median_oae(1,q) + 3*std_oae(1,q)
                            noisy_trials = noisy_trials+1;
                            break;
                        end
                    end
                end
            end
            
            % if SNR is good enough and we've hit the minimum number of
            % trials, then stop.
            if SNR_temp(1:8) >= SNRcriterion
                if k_kept >= minTrials + noisy_trials
                    doneWithTrials = 1;
                end
            elseif k == maxTrials
                doneWithTrials = 1;
            end
            
            pass = (SNR_temp>=SNRcriterion);
            oae_pass = dpoae_w;
            oae_fail = dpoae_w;
            oae_pass(~pass) = NaN;
            oae_fail(pass) = NaN;
            
            % Plot amplitudes from live analysis
            hold off;
            plot(centerFreqs./1000,oae_pass, 'o', 'linew', 2, 'color', [0 0.4470 0.7410]);
            hold on;
            plot(centerFreqs/1000,oae_fail, 'o', 'linew', 2, 'color', 'k'),
            plot(centerFreqs./1000,dpnf_w, 'x', 'linew', 2, 'color', [0.6350 0.0780 0.1840]);
            hold off;
            legend('DPOAE', '', 'NOISE', 'location', 'northeast');
            title('DPOAE')
            xlabel('Frequency (Hz)')
            ylabel('Median Amplitude (dB SPL)')
            set(gca, 'XScale', 'log', 'FontSize', 14)
            xlim([0.5, 16]);
            xticks([.5, 1, 2, 4, 8, 16]);
            ylim([-45, 45]);
            yticks((-45:15:45))
            grid on;
            drawnow;
            
            fprintf(1, 'Trials run: %d / Noisy Trials: %d \n', k_kept, noisy_trials);
        end
    end
    
    %% Full data struct
    data.info = info;
    data.stim = stim;
    data.info.subj = subj;
    data.resp.trialsCollected = k_kept;
    data.resp.AllBuffs = resp(1:k_kept,:);
    data.resp.testDur_s = toc;
    data.resp.noisyTrials = noisy_trials; 
    
    save(fname,'data');
    
    %% Close TDT, ER-10X connections etc. and cleanup
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    
    %% Quick Analysis
    analyze = questdlg('Do you want to see the analysis?', 'Analysis?', 'Yes', 'No', 'No'); 
    if strcmp(analyze, 'Yes')
        Analyze_DPswept; 
    end
    
catch me
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    rethrow(me);
end
