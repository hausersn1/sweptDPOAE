
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f1,RP,FS]=load_play_circuit_Nel2(FS_tag,fig_num,GB_ch)
% Loads the TDT circuit and makes actx links necessary
% The TDT matlab syntax has now changed to look more
% like OOPS but this old style is still supported.
%
%------------
% Hari Bharadwaj, September 6, 2010
%------------
warning('off'); 

CIR_PATH='BasicPlay_OAE_Nel2.rcx'; %The *.rco circuit used to play the files

%Generate the actx control window in a specified figure:
%-------------------------------------------------------
f1=figure(fig_num);
set(f1,'Position',[5 5 30 30],'Visible','off'); %places and hides the ActX Window
RP=actxcontrol('RPco.x',[5 5 30 30],f1); %loads the actx control for using rco circuits
invoke(RP,'ConnectRX8','GB',GB_ch); 
% The rco circuit can be run at the following set of discrete sampling
% frequencies (in Hz): 0=6k, 1=12k, 2=25k, 3=50k, 4=100k, 5=200k.
% Use the tags listed above to specify below:
%--------------------------------------------------------------------------
invoke(RP,'LoadCOFsf',CIR_PATH,FS_tag); %loads the circuit using the specified sampling freq.
FS_vec=[6 12 24.4140625 48.828125 100 200]*1e3; %NOT EXACT!!!!!!!!!!!!!!!!!!!
FS=FS_vec(FS_tag+1);

invoke(RP,'Run'); %start running the circuit

Status = double(invoke(RP,'GetStatus'));
if bitget(double(Status),1)==0
    error('Error connecting to RP2');
elseif bitget(double(Status),2)==0
    error('Error loading circuit');
elseif bitget(double(Status),3)==0
    error('Error running circuit');
end

% End of load_play_circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%