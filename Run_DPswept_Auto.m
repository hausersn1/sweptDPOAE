% Swept DPOAEs for Humans (ER-10X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other m-files required: Make_DPswept_log.m, etc.
% Other files required: none
% MAT-files required: none
%
% References: 
%   SNR based endpoint: 
% 
% Samantha Hauser
% November 2021; Last revision: 6-June-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    initializeER10X;
    
    % Initializing TDT
    % Specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    % Get stimulus structure
    stim = Make_DPswept_log;
    
    % Get subject and ear info
    subj = input('Please subject ID:', 's');
    earflag = 1;
    while earflag == 1
        ear = input('Please enter which year (L or R):', 's');
        switch ear
            case {'L', 'R', 'l', 'r', 'Left', 'Right', 'left', 'right', 'LEFT', 'RIGHT'}
                earname = strcat(ear, 'Ear');
                earflag = 0;
            otherwise
                fprintf(2, 'Unrecognized ear type! Try again!');
        end
    end
    
    stim.subj = subj;
    stim.ear = ear;
    
    % Make directory to save results
    paraDir = 'C:\Experiments\Sam\DPOAEswept\Results\';
    
    addpath(genpath(paraDir));
    if(~exist(strcat(paraDir,'\',subj),'dir'))
        mkdir(strcat(paraDir,'\',subj));
    end
    respDir = strcat(paraDir,'\',subj,'\');
    
    % Set stims in buffdata:
    buffdata = [stim.y1; stim.y2];
    % Check for clipping and load to buffer
    if(any(abs(buffdata(:)) > 1))
        error('What did you do!? Sound is clipping!! Cannot Continue!!\n');
    end
    
    button = input('Do you want the subject to press a button to proceed? (Y or N):', 's');
    switch button
        case {'Y', 'y', 'yes', 'Yes', 'YES'}
            getResponse(RZ);
            fprintf(1, '\nSubject pressed a button...\nStarting Stimulation...\n');
        otherwise
            fprintf(1, '\nStarting Stimulation...\n');
    end
    
    %% Set attenuation and play
    drop_f1 = stim.drop_f1;
    drop_f2 = stim.drop_f2;
    delayComp = 1; % Always
    
    windowdur = 0.5;
    SNRcriterion = stim.SNRcriterion;
    maxTrials = stim.maxTrials;
    minTrials = stim.minTrials;
    doneWithTrials = 0;
    figure;
    
    %% Add useful info to structure
    mic_sens = 50e-3; % mV/Pa
    mic_gain = db2mag(gain + 6); % +6 for balanced cable
    P_ref = 20e-6;
    DR_onesided = 1;
    stim.VoltageToPascal = 1 / (DR_onesided * mic_gain * mic_sens);
    stim.PascalToLinearSPL = 1 /  P_ref;
    
    resp = zeros(maxTrials, size(buffdata,2));
    
    %% Loop for presenting stimuli
    % variable for live analysis
    k = 0;
    t = stim.t;
    testfreq = [.75, 1, 1.5, 2, 3, 4, 6, 8, 12].* 1000;
    
    if stim.speed < 0
        f1 = stim.fmax;
        f2 = stim.fmin;
    else
        f1 = stim.fmin;
        f2 = stim.fmax;
    end
    
    if stim.speed < 20
        t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
    else
        t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
    end
    
    
    while doneWithTrials == 0
        k = k + 1;
        
        %Start playing from the buffer:
        vins = playCapture2(buffdata, card, 1, 0,...
            drop_f1, drop_f2, delayComp);
        
        if k > stim.ThrowAway
            resp(k - stim.ThrowAway,  :) = vins;
        end
        
        WaitSecs(0.25);
        
        fprintf(1, 'Done with # %d trials \n', k);
        
        % test OAE
        
        OAEtrials = resp(1:k-stim.ThrowAway, :);
        OAE = median(OAEtrials,1);
        coeffs_temp = zeros(length(testfreq), 2);
        coeffs_noise = zeros(length(testfreq), 8);
        for m = 1:length(testfreq)
            win = find( (t > (t_freq(m)-windowdur/2)) & ...
                (t < (t_freq(m)+windowdur/2)));
            taper = hanning(numel(win))';
            
            oae_win = OAE(win) .* taper;
            
            phiProbe_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
            
            model_dp = [cos(phiProbe_inst(win)) .* taper;
                -sin(phiProbe_inst(win)) .* taper];
            
        if stim.speed < 0 % fixed
            nearfreqs = [1.10, 1.12, 1.14, 1.16];
        else
            nearfreqs = [.90, .88, .86, .84];
        end
        
        model_noise = ...
            [cos(nearfreqs(1)*phiProbe_inst(win)) .* taper;
            -sin(nearfreqs(1)*phiProbe_inst(win)) .* taper;
            cos(nearfreqs(2)*phiProbe_inst(win)) .* taper;
            -sin(nearfreqs(2)*phiProbe_inst(win)) .* taper;
            cos(nearfreqs(3)*phiProbe_inst(win)) .* taper;
            -sin(nearfreqs(3)*phiProbe_inst(win)) .* taper;
            cos(nearfreqs(4)*phiProbe_inst(win)) .* taper;
            -sin(nearfreqs(4)*phiProbe_inst(win)) .* taper];
            
            coeffs_temp(m,:) = model_dp' \ oae_win';
            coeffs_noise(m,:) = model_noise' \ oae_win';
        end
        
        % for noise
        noise2 = zeros(length(testfreq),4);
        for i = 1:2:8
            noise2(:,ceil(i/2)) = abs(complex(coeffs_noise(:,i), coeffs_noise(:,i+1)));
        end
        noise = mean(noise2, 2);
        
        oae = abs(complex(coeffs_temp(:,1), coeffs_temp(:,2)));
        
        SNR_temp = db(oae) - db(noise);
        
        mult = stim.VoltageToPascal .* stim.PascalToLinearSPL;
        hold off;
        plot(testfreq./1000,db(oae.*mult), 'o', 'linew', 2);
        hold on;
        plot(testfreq./1000,db(noise.*mult), 'x', 'linew', 2);
        legend('DPOAE', 'NOISE');
        xlabel('Frequency (Hz)')
        ylabel('Median Amplitude dB')
        set(gca, 'XScale', 'log', 'FontSize', 14)
%         xticks([.5, 1, 2, 4, 8, 16])
        xlim([0.5, 16])
        ylim([-30, 30]); 
        yticks((-30:6:30))
        grid on; 
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
    
    stim.resp = resp(1:k-stim.ThrowAway,:);
    
    
    %% Save results
    datetag = datestr(clock);
    stim.date = datetag;
    datetag(strfind(datetag,' ')) = '_';
    datetag(strfind(datetag,':')) = '_';
    fname = strcat(respDir, 'DPOAEswept_', stim.subj, '_', stim.ear,...
        '_', datetag, '.mat');
    save(fname,'stim');
    
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
