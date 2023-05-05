function stim = Make_DPswept_log()
% rawstim (structure) should contain fields fmin, fmax, speed, Fs, ratio,
% VtoPforH

stim.fmin = 500;
stim.fmax = 16000;
stim.speed = -1; %change to linear, Hz/sec
stim.drop_f1 = 53;  % levels
stim.drop_f2 = stim.drop_f1 + 10 ; %f2 is 10 dB softer (more attn)
stim.buffdur = 0.25;% buffer duration
stim.ThrowAway = 1;
%stim.Averages = 75; % num of trials
stim.Averages = 20;
stim.SNRcriterion = 6; 
stim.minTrials = 12; 
stim.maxTrials = 50; 
stim.ratio = 1.22;
stim.Fs = 48828.125;

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

% in Frequency
% buffdur_exact = t(buffinst1);
% f2_inst = f1 + stim.speed*(t-buffdur_exact);  
% f2_inst(1:buffinst1) = f1;
% f2_inst(buffinst2:end) = f2;

% in Phase
start_2 = f1*t(1:buffinst1);
buffdur_exact = t(buffinst1);
phi2_inst = f1*(2.^( (t-buffdur_exact) * stim.speed) - 1) / (stim.speed * log(2)) + start_2(end); % Cycles 
end_2 = f2*t(1:(length(t)-buffinst2+1)) + phi2_inst(buffinst2); 
phi2_inst(1:buffinst1) = start_2;
phi2_inst(buffinst2:end) = end_2;

% f2_inst = fmax * 2.^( (t - buffdur) * stim.speed);
% f2_inst(1:buffinst1) = fmax;
% f2_inst(buffinst2:end) = fmin;
% phi2_inst = fmax * (2.^( (t - buffdur) * stim.speed) - 1) / (stim.speed * log(2));
% phi2_inst(1:buffinst1) = fmax .* (t(1:buffinst1) - buffdur);
% phi2_inst(buffinst2:end) = fmin .* (t(buffinst2:end) - t(buffinst2)) + phi2_inst(buffinst2);


%f1_inst = f2_inst / stim.ratio;
phi1_inst = phi2_inst / stim.ratio;
stim.y1 = scaleSound(rampsound(cos(2 * pi * phi1_inst), stim.Fs, 0.005));
stim.y2 = scaleSound(rampsound(cos(2 * pi * phi2_inst), stim.Fs, 0.005));
%stim.f1_inst = f1_inst;
%stim.f2_inst = f2_inst;
stim.phi1_inst = phi1_inst;
stim.phi2_inst = phi2_inst;
