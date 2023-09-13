function stim = Make_DPswept_log()

stim.fmin = 500;
stim.fmax = 16000;
stim.speed = 1; % oct/sec
stim.drop_f1 = 40; %19;  % levels %53
stim.drop_f2 = stim.drop_f1 + 10 ; %f2 is 10 dB softer (more attn)
stim.ratio = 1.22;

stim.buffdur = 0.25;% buffer duration
stim.Fs = 48828.125;

stim.ThrowAways = 0; 
stim.Averages = 48; 

stim.SNRcriterion = 6; 
stim.maxTrials = 50; 
stim.minTrials = 12; 
stim.ThrowAway = 1;

% Create stimulus
buffdur = stim.buffdur;
Fs = stim.Fs; 

if stim.speed < 0 % downsweep
    f1 = stim.fmax; 
    f2 = stim.fmin; 
else % upsweep
    f1 = stim.fmin; 
    f2 = stim.fmax; 
end 

dur = log2(stim.fmax/stim.fmin) / abs(stim.speed) + (2*buffdur);
t = 0: (1/Fs): (dur - 1/Fs);
stim.t = t;

buffinst1 = find(t < buffdur, 1, 'last');
buffinst2 = find(t > (dur - buffdur), 1, 'first');

% in Phase
start_2 = f1*t(1:buffinst1);
buffdur_exact = t(buffinst1);
phi2_inst = f1*(2.^( (t-buffdur_exact) * stim.speed) - 1) / (stim.speed * log(2)) + start_2(end); % Cycles 
end_2 = f2*t(1:(length(t)-buffinst2+1)) + phi2_inst(buffinst2); 
phi2_inst(1:buffinst1) = start_2;
phi2_inst(buffinst2:end) = end_2;

phi1_inst = phi2_inst / stim.ratio;
stim.y1 = scaleSound(rampsound(cos(2 * pi * phi1_inst), stim.Fs, 0.005));
stim.y2 = scaleSound(rampsound(cos(2 * pi * phi2_inst), stim.Fs, 0.005));
stim.phi1_inst = phi1_inst;
stim.phi2_inst = phi2_inst;
