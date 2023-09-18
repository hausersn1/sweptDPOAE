function stim = Make_DPswept()

% Stimulus Parameters
stim.fmin = 500; %Hz of F2
stim.fmax = 16000; % Hz
stim.speed = 1; % oct/sec
stim.drop_f1 = 29; %19;  % levels %53
stim.drop_f2 = stim.drop_f1 + 10 ; %f2 is 10 dB softer (more attn)
stim.ratio = 1.22;
stim.scale = 'log';
stim.buffdur = 0.25;% buffer duration
stim.Fs = 48828.125;

% Trial Parameters
stim.SNRcriterion = 6;
stim.maxTrials = 50;
stim.minTrials = 12;
stim.ThrowAway = 0;

% Add useful info to structure
gain = 30; 
stim.mic_sens = 50e-3; % mV/Pa
stim.mic_gain = db2mag(gain + 6); % +6 for balanced cable
stim.P_ref = 20e-6;
stim.DR_onesided = 1;
stim.VoltageToPascal = 1 / (stim.DR_onesided * stim.mic_gain * stim.mic_sens);
stim.PascalToLinearSPL = 1 /  stim.P_ref;

%% For live analysis
stim.windowdur = 0.25;
%stim.testfreq = [.75, 1, 1.5, 2, 3, 4, 6, 8, 12].* 1000;
stim.npoints = 128; 

if stim.speed < 0 % fixed
    stim.nearfreqs = [1.10, 1.12, 1.14, 1.16];
else
    stim.nearfreqs = [.90, .88, .86, .84];
end

%% Create stimulus
buffdur = stim.buffdur;
Fs = stim.Fs;

if stim.speed < 0 % downsweep
    f_start = stim.fmax;
    f_end = stim.fmin;
else % upsweep
    f_start = stim.fmin;
    f_end = stim.fmax;
end

if strcmp(stim.scale, 'log')
    dur = log2(stim.fmax/stim.fmin) / abs(stim.speed) + (2*buffdur);
else
    dur = abs(f_start - f_end) / abs(stim.speed) + (2*buffdur);
end

t = 0: (1/Fs): (dur - 1/Fs);
stim.t = t;

buffinst1 = find(t < buffdur, 1, 'last');
buffinst2 = find(t > (dur - buffdur), 1, 'first');
buffdur_exact = t(buffinst1);

if strcmp(stim.scale, 'log')
    % in Phase
    start_2 = f_start*t(1:buffinst1);
    phi2_inst = f_start*(2.^( (t-buffdur_exact) * stim.speed) - 1) / (stim.speed * log(2)) + start_2(end); % Cycles
    end_2 = f_end*t(1:(length(t)-buffinst2+1)) + phi2_inst(buffinst2);
    phi2_inst(1:buffinst1) = start_2;
    phi2_inst(buffinst2:end) = end_2;
else % linear sweep
    % in Phase
    start_2 = f_start*t(1:buffinst1);
    buffdur_exact = t(buffinst1);
    phi2_inst = f_start*(t-buffdur_exact) + stim.speed*((t-buffdur_exact).^2)/2 + start_2(end); % Cycles
    end_2 = f_end*t(1:(length(t)-buffinst2+1)) + phi2_inst(buffinst2);
    phi2_inst(1:buffinst1) = start_2;
    phi2_inst(buffinst2:end) = end_2;
end

phi1_inst = phi2_inst / stim.ratio;
stim.y1 = scaleSound(rampsound(cos(2 * pi * phi1_inst), stim.Fs, 0.005));
stim.y2 = scaleSound(rampsound(cos(2 * pi * phi2_inst), stim.Fs, 0.005));
stim.phi1_inst = phi1_inst;
stim.phi2_inst = phi2_inst;
