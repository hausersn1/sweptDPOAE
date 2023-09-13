function stim = Make_DPswept()
% rawstim (structure) should contain fields fmin, fmax, speed, Fs, ratio,
% VtoPforH

stim.fmin = 500;
stim.fmax = 16000;
stim.speed = 1; % octaves/sec
stim.drop_f1 = 19;  % levels
stim.drop_f2 = stim.drop_f1 + 10 ; %f2 is 10 dB softer (more attn)
stim.buffdur = 0.25;% buffer duration
stim.ThrowAways = 1;
stim.Averages = 75; % num of trials
stim.ratio = 1.22;
stim.Fs = 48828.125;

% Create stimulus
buffdur = stim.buffdur;
Fs = stim.Fs; 

if stim.speed < 0 % downsweep
    f_start = stim.fmax; 
    f_end = stim.fmin; 
else % upsweep
    f_start = stim.fmin; 
    f_end = stim.fmax; 
end 

dur = abs(f_start - f_end) / abs(stim.speed) + (2*buffdur);
t = 0: (1/Fs): (dur - 1/Fs);
stim.t = t;

buffinst1 = find(t < buffdur, 1, 'last');
buffinst2 = find(t > (dur - buffdur), 1, 'first');

% in Frequency
buffdur_exact = t(buffinst1);
f2_inst = f_start + stim.speed*(t-buffdur_exact);  
f2_inst(1:buffinst1) = f_start;
f2_inst(buffinst2:end) = f_end;

% in Phase
start_2 = f_start*t(1:buffinst1);
buffdur_exact = t(buffinst1);
phi2_inst = f_start*(t-buffdur_exact) + stim.speed*((t-buffdur_exact).^2)/2 + start_2(end); % Cycles 
end_2 = f_end*t(1:(length(t)-buffinst2+1)) + phi2_inst(buffinst2); 
phi2_inst(1:buffinst1) = start_2;
phi2_inst(buffinst2:end) = end_2;

% f2_inst = fmax * 2.^( (t - buffdur) * stim.speed);
% f2_inst(1:buffinst1) = fmax;
% f2_inst(buffinst2:end) = fmin;
% phi2_inst = fmax * (2.^( (t - buffdur) * stim.speed) - 1) / (stim.speed * log(2));
% phi2_inst(1:buffinst1) = fmax .* (t(1:buffinst1) - buffdur);
% phi2_inst(buffinst2:end) = fmin .* (t(buffinst2:end) - t(buffinst2)) + phi2_inst(buffinst2);


f1_inst = f2_inst / stim.ratio;
phi1_inst = phi2_inst / stim.ratio;
stim.y1 = scaleSound(rampsound(cos(2 * pi * phi1_inst), stim.Fs, 0.005));
stim.y2 = scaleSound(rampsound(cos(2 * pi * phi2_inst), stim.Fs, 0.005));
stim.f1_inst = f1_inst;
stim.f2_inst = f2_inst;
stim.phi1_inst = phi1_inst;
stim.phi2_inst = phi2_inst;
