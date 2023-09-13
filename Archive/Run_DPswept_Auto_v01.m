%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Swept DPOAEs for Humans (ER-10X)
%
%
% Samantha Hauser
% Created: November 2021
% Last revision: 5-Sep-2023
%
% References:
%   SNR based endpoint:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up data storage and subject info

% Measure-General info
data.info.name = 'DPOAEswept';
data.info.version = 'Auto_v01';

% Visit info
if exist('C:\Experiments\Sam\current_visit.mat','file')
    load('C:\Experiments\Sam\current_visit.mat', 'visit')
    ask = questdlg(sprintf('Is this subject %s?', visit.subj), 'Check Subject', 'Yes', 'No', 'No');
else
    ask = 'No';
end

if strcmp(ask, 'No')
    cd ..
    startVisit
    cd(data.info.name)
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

fname = strcat(respDir, info.name, '_', ...
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
    
    t = stim.t;
    testfreq = stim.testfreq;
    
    if stim.speed < 0
        f_start = stim.fmax;
        f_end = stim.fmin;
    else
        f_start = stim.fmin;
        f_end = stim.fmax;
    end
    
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
        
        %Start playing from the buffer:
        vins = playCapture2(buffdata, card, 1, 0,...
            drop_f1, drop_f2, delayComp);
        
        % save the response
        if k > stim.ThrowAway
            resp(k - stim.ThrowAway,  :) = vins;
        end
        
        WaitSecs(0.25);
        fprintf(1, 'Trials run: %d \n', k);
        
        % test OAE
        if k > stim.ThrowAway
            OAEtrials = resp(1:k-stim.ThrowAway, :);
            OAE = median(OAEtrials,1);
            coeffs_temp = zeros(length(testfreq), 2);
            coeffs_noise = zeros(length(testfreq), 8);
            for m = 1:length(testfreq)
                win = find( (t > (t_freq(m)-windowdur/2)) & ...
                    (t < (t_freq(m)+windowdur/2)));
                taper = hanning(numel(win))';
                
                oae_win = OAE(win) .* taper;
                
                phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi; %dp
                
                model_dp = [cos(phi_dp_inst(win)) .* taper;
                    -sin(phi_dp_inst(win)) .* taper];
                
                nearfreqs = stim.nearfreqs;
                
                model_noise = ...
                    [cos(nearfreqs(1)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(1)*phi_dp_inst(win)) .* taper;
                    cos(nearfreqs(2)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(2)*phi_dp_inst(win)) .* taper;
                    cos(nearfreqs(3)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(3)*phi_dp_inst(win)) .* taper;
                    cos(nearfreqs(4)*phi_dp_inst(win)) .* taper;
                    -sin(nearfreqs(4)*phi_dp_inst(win)) .* taper];
                
                coeffs_temp(m,:) = model_dp' \ oae_win';
                coeffs_noise(m,:) = model_noise' \ oae_win';
            end
            
            % for noise
            noise = zeros(length(testfreq),4);
            for i = 1:2:8
                noise(:,ceil(i/2)) = abs(complex(coeffs_noise(:,i), coeffs_noise(:,i+1)));
            end
            noise = mean(noise, 2);
            
            oae = abs(complex(coeffs_temp(:,1), coeffs_temp(:,2)));
            
            SNR_temp = db(oae) - db(noise);
            
            VtoSPL = stim.VoltageToPascal .* stim.PascalToLinearSPL;
            
            % Plot amplitudes from live analysis
            hold off;
            plot(testfreq./1000,db(oae.*VtoSPL), 'o', 'linew', 2);
            hold on;
            plot(testfreq./1000,db(noise.*VtoSPL), 'x', 'linew', 2);
            hold off; 
            legend('DPOAE', 'NOISE', 'location', 'southwest');
            title('DPOAE')
            xlabel('Frequency (Hz)')
            ylabel('Median Amplitude (dB SPL)')
            set(gca, 'XScale', 'log', 'FontSize', 14)
            xlim([0.5, 16]);
            xticks([.5, 1, 2, 4, 8, 16]);
            ylim([-30, 30]);
            yticks((-30:6:30))
            grid on;
            
            % Plot SNRs
            xstart = 0.6;
            xend = .9;
            ystart = 0.6;
            yend = .9;
            
            axes('Position',[xstart ystart xend-xstart yend-ystart], 'XScale', 'log')
            box on
            hold on
            bar(testfreq./1e3, SNR_temp,'Color', 'k');
            plot([0, 16], [stim.SNRcriterion, stim.SNRcriterion], '--', 'Color', [0.5, 0.5, 0.5]); 
            xlim([0.5,16]);
            ylim([-6,12]);
            yticks((-6:3:12))
            xlabel('Frequency (Hz)');
            ylabel('SNR')
            hold off
            set(gcf,'Position',[1557 538 560 420])
            drawnow;
            
            % if SNR is good enough and we've hit the minimum number of
            % trials, then stop.
            if SNR_temp(1:8) > SNRcriterion
                if k - stim.ThrowAway >= minTrials
                    doneWithTrials = 1;
                end
            elseif k == maxTrials
                doneWithTrials = 1;
            end
        end
    end

    %% Full data struct
    data.info = info; 
    data.stim = stim;
    data.info.subj = subj;
    data.resp.trialsCollected = k-stim.ThrowAway;
    data.resp.AllBuffs = resp(1:k-stim.ThrowAway,:);
    data.info.testdur_s = toc; 
    
    save(fname,'data');
    
    %% Close TDT, ER-10X connections etc. and cleanup
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    
catch me
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    rethrow(me);
end
