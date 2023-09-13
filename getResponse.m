function resp = getResponse(RP,dur)
% Queries the TDT button response variables for 'dur' seconds and returns
% sequence
% USAGE:
%  resp = getResponse(RP,dur);
%  OR
% resp = getResponse(RP);
%
%  If dur is given (in seconds), it will check for responses for 'dur'
%  seconds...if no argument is given, it will wait for  nresp presses
%
% nresp is hardcoded in the first line of the code

nresp = 1;
if(~exist('dur','var'))
    dur = [];
    waitflag = 1;
else
    waitflag = 0;
end
rate = 0.05;


if(isempty(dur))
    dur = 0.5;
end

resp = [];

carry = 0;
while(numel(resp) < nresp)
    nquery = floor(dur/rate);
    resptemp  = zeros(1,nquery);
    for j = 1:nquery
        resptemp(j)= RP.GetTagVal('respMax');
        RP.SetTagVal('respRst',1);
        pause(rate);
    end
    
    resp = [resp, 1+log2(resptemp(diff([carry resptemp])>0))]; %#ok<AGROW>
    carry = resptemp(end);
    RP.SetTagVal('respRst',1);
    if(~waitflag)
        break;
    end
end