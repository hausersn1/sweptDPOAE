function y = rampsound(x,fs,risetime)
% Function to ramp a sound file using a dpss ramp
% USAGE:
% y = rampsound(x,fs,risetime)
%
% risetime in seconds, fs in Hz
% Hari Bharadwaj

Nramp = ceil(fs*risetime*2)+1;
w = dpss(Nramp,1,1);

w = w - w(1);
w = w/max(w);
sz = size(x);
half = ceil(Nramp/2);
wbig = [w(1:half); ones(numel(x)- 2*half,1); w((end-half+1):end)];

if(sz(1)== numel(x))
    y = x.*wbig;
else
    y = x.*wbig';
end

