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
    paraDir = 'C:\Experiments\Sam\sweptDPOAE-main\Results\';
    
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

    resp = zeros(stim.Averages, size(buffdata,2)); 
    
    for k = 1:(stim.ThrowAways + stim.Averages)

        %Start playing from the buffer:
        vins = playCapture2(buffdata, card, 1, 0,...
            drop_f1, drop_f2, delayComp);
        
        if k > stim.ThrowAways
            resp(k - stim.ThrowAways,  :) = vins;
        end
        
        WaitSecs(0.25); 

        fprintf(1, 'Done with # %d / %d trials \n', k, (stim.Averages + stim.ThrowAways));
    end
    
    stim.resp = resp; 

    %% Add useful info to structure
    mic_sens = 50e-3; % mV/Pa
    mic_gain = db2mag(gain + 6); % +6 for balanced cable
    P_ref = 20e-6;
    DR_onesided = 1;
    stim.VoltageToPascal = 1 / (DR_onesided * mic_gain * mic_sens);
    stim.PascalToLinearSPL = 1 /  P_ref;
    
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