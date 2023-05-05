%% Run swept DPOAE using NEL (but outside it)
mainDir = pwd;

% Create subject
subj = input('Please subject ID:', 's');

date = datetime('now');
date.Format = 'yyyy_MM_dd';
datetag = string(date);
stim.date = datetag;

% Get CalibFile from NEL DATA storage
ExpDataDir = 'C:\NEL\ExpData\';
cd(ExpDataDir);
findDir = dir(sprintf('*%s*%s*', datetag, subj));
calibflag = 1;
noCalibFlag = 0;

while calibflag == 1
    calibNum = input('What NEL calib #? ', 's');
    if isempty(findDir)
        fprintf(2,'You need to run a NEL calibration first!\n');
        return
    elseif length(findDir)~=1
        fprintf(2,'Multiple Directories. I am confused. \n')
        return
    else
        datadir = [ExpDataDir findDir.name];
        cd(datadir);
        calib_file = dir(sprintf('coef_00%s_calib.mat', calibNum));
        if isempty(calib_file)
            if calibNum == '999'
                noCalibFlag = 1;
                stim.b = 0;
                calibflag = 0;
            else
                fprintf(2, 'No calib of that number. Try again! \n')
            end
        else
            load(calib_file.name, 'b');
            stim.b = b;
            calibflag = 0;
        end
    end
end

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
cd(mainDir)

% Make directory to save results
paraDir = 'C:\Users\Heinz Lab - NEL2\Desktop\OAEs\sweptDPOAE\Results';
addpath(genpath(paraDir));
if(~exist(strcat(paraDir,'\',subj),'dir'))
    mkdir(strcat(paraDir,'\',subj));
end
respDir = strcat(paraDir,'\',subj,'\');

% Tell user to make sure ER-10B+ gain is set correctly
uiwait(warndlg('Set ER-10B+ GAIN to 40 dB','SET ER-10B+ GAIN WARNING','modal'));

%% Start w/ Delay if needed
% Get stimulus structure
stim = Make_DPswept_NEL;
stim.subj = subj;
stim.ear = ear;

windowdur = 0.5;
SNRcriterion = stim.SNRcriterion;
maxTrials = stim.maxTrials;
minTrials = stim.minTrials;

doneWithTrials = 0;
figure;

% Set attenuation and stims:
buffdata = [stim.y1; stim.y2];
drop_f1 = stim.drop_f1;
drop_f2 = stim.drop_f2;

% filter data
if noCalibFlag ~= 1
    buffdata = filter(b, 1, buffdata, [], 1);
end

% Make arrays to store measured mic outputs
resp = zeros(maxTrials, size(buffdata,2));

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
        
        phi_dp_inst = (2.*stim.phi1_inst - stim.phi2_inst) * 2 * pi;
        phiProbe_inst = phi_dp_inst;
        model_dp = [cos(phiProbe_inst(win)) .* taper;
            -sin(phiProbe_inst(win)) .* taper];
        
        model_noise = ...
            [cos(0.9*phiProbe_inst(win)) .* taper;
            -sin(0.9*phiProbe_inst(win)) .* taper;
            cos(0.88*phiProbe_inst(win)) .* taper;
            -sin(0.88*phiProbe_inst(win)) .* taper;
            cos(0.86*phiProbe_inst(win)) .* taper;
            -sin(0.86*phiProbe_inst(win)) .* taper;
            cos(0.84*phiProbe_inst(win)) .* taper;
            -sin(0.84*phiProbe_inst(win)) .* taper];
        
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
    
    hold off;
    plot(testfreq./1000,db(oae.*10000), 'o', 'linew', 2);
    hold on;
    plot(testfreq./1000,db(noise.*10000), 'x', 'linew', 2);
    legend('DPOAE', 'NOISE');
    xlabel('Frequency (Hz)')
    ylabel('Median Amplitude dB')
    set(gca, 'XScale', 'log', 'FontSize', 14)
    xticks([.5, 1, 2, 4, 8, 16])
    xlim([0.5, 16])
    drawnow;
    
    if SNR_temp(1:8) > SNRcriterion
        if k - stim.ThrowAway >= minTrials
            doneWithTrials = 1;
        end
    elseif k == maxTrials
        doneWithTrials = 1;
    end
    
    fprintf(1, 'Done with trial %d\n', k);
    
end
stim.resp = resp(1:k-stim.ThrowAway,:);

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
fprintf(1, 'Saved!\n');
%% Close TDT, ER-10X connections etc. and cleanup
close_play_circuit(f1RZ, RZ);

