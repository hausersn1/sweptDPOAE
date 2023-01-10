%% Run swept DPOAE using NEL (but outside it)

% Create subject
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

% Make directory to save results
paraDir = 'C:\Users\Heinz Lab - NEL2\Desktop\sweptDPOAE\Results';
addpath(genpath(paraDir));
if(~exist(strcat(paraDir,'\',subj),'dir'))
    mkdir(strcat(paraDir,'\',subj));
end
respDir = strcat(paraDir,'\',subj,'\');

% Tell user to make sure ER-10B+ gain is set correctly
uiwait(warndlg('Set ER-10B+ GAIN to 40 dB','SET ER-10B+ GAIN WARNING','modal'));

%% Start w/ Delay if needed
% Get stimulus structure
stim = Make_DPswept_Chin;
stim.subj = subj;
stim.ear = ear;

% Set attenuation and stims:
buffdata = [stim.y1; stim.y2];
drop_f1 = stim.drop_f1;
drop_f2 = stim.drop_f2;

% Make arrays to store measured mic outputs
resp = zeros(stim.Averages, size(buffdata,2));

% Ask if we want a delay (for running yourself)
button = input('Do you want a 10 second delay? (Y or N):', 's');
switch button
    case {'Y', 'y', 'yes', 'Yes', 'YES'}
        DELAY_sec=10;
        fprintf(1, '\n%.f seconds until START...\n',DELAY_sec);
        pause(DELAY_sec)
        fprintf(1, '\nWe waited %.f seconds ...\nStarting Stimulation...\n',DELAY_sec);
    otherwise
        fprintf(1, '\nStarting Stimulation...\n');
end

% Initializing TDT
fig_num=99;
GB_ch=1;
FS_tag = 3;
Fs = 48828.125;
[f1RZ,RZ,~]=load_play_circuit_Nel2(FS_tag,fig_num,GB_ch);

%% Loop for presenting stimuli
for k = 1:(stim.ThrowAway + stim.Averages)
    
    % Load the 2ch variable data into the RZ6:
    invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', buffdata(1, :));
    invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', buffdata(2, :));
    % Set the delay of the sound
    invoke(RZ, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
    playrecTrigger = 1;
    % Set attenuations
    rc = PAset([0, 0, drop_f1, drop_f2]);
    % Set total length of sample
    RZ6ADdelay = 97; % Samples
    resplength = size(buffdata,2) + RZ6ADdelay; % How many samples to read from OAE buffer
    invoke(RZ, 'SetTagVal', 'nsamps', resplength);
    
    %Start playing from the buffer:
    invoke(RZ, 'SoftTrg', playrecTrigger);
    currindex = invoke(RZ, 'GetTagVal', 'indexin');
    
    while(currindex < resplength)
        currindex=invoke(RZ, 'GetTagVal', 'indexin');
    end
    
    vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
        'F32','F64',1);
    
    % Save data
    if k > stim.ThrowAway
        resp(k - stim.ThrowAway,  :) = vin((RZ6ADdelay + 1):end);
    end
    
    % Get ready for next trial
    invoke(RZ, 'SoftTrg', 8); % Stop and clear "OAE" buffer
    %Reset the play index to zero:
    invoke(RZ, 'SoftTrg', 5); %Reset Trigger
    
    pause(0.05);
    
    fprintf(1, 'Done with trial %d / %d\n', k,...
        (stim.ThrowAway + stim.Averages));
    
end
stim.resp = resp;

%% Add useful info to structure
mic_sens = 0.05; % mV / Pa
mic_gain = db2mag(40);

P_ref = 20e-6; % * sqrt(2);
DR_onesided = 1;

stim.VoltageToPascal = 1 / (DR_onesided * mic_gain * mic_sens);
stim.PascalToLinearSPL = 1 /  P_ref;
%% Save Measurements
datetag = datestr(clock);
click.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '_';
fname = strcat(respDir,'DPOAEswept_',subj,'_', earname,'_',datetag, '.mat');
save(fname,'stim');
fprintf(1, 'Saved!');
%% Close TDT, ER-10X connections etc. and cleanup
close_play_circuit(f1RZ, RZ);
Analyze_DPswept;
