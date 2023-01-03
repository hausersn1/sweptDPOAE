function rc = PAset(attn)
%

% AF 8/27/01

global PA

if ((length(attn) > 1) && (length(PA) ~= length(attn)))
   nelerror('PAset: incompatable number of PA devices and attenuations');
   rc = -1;
   return;
end
rc = 1;
for i = 1:length(PA)
   attnval = attn(min(i,end));
   if (attnval ~= PA(i).attn)
      PA(i).attn = attnval;
      lrc = invoke(PA(i).activeX,'SetAtten',attnval);
      if (lrc==0)
         nelerror(['PAset: Can''t set attenuation on PA #' int2str(PA(i).serial_no)]);
         rc = 0;
      end
   end
end
