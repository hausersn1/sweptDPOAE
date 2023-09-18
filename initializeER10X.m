% Initialize ER-10X if used for any measurement
p10x = genpath('C:\Experiments\ER-10X\MATLABAPI\Matlab\');
addpath(p10x);
loaded = ER10XConnectLoadLib('C:\Experiments\ER-10X\MATLABAPI\');
[err, ER10XHandle] = er10x_open();
fprintf(1, 'Result of ER10_X_OPEN: %s\n', err);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Could not open ER10X!');
end
err = er10x_connect(ER10XHandle);
fprintf(1, 'Result of ER10_X_CONNECT: %s\n', err);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Could not connect to ER10X!');
end

probeIndex = 1; % 0=A/1=B (changed to B 9/15/23)
gain = 30; % dB
temperatureF = 0; % Turn off heater
er10x_set_gain(ER10XHandle,probeIndex,gain);
er10x_set_output_limiter(ER10XHandle, 1); % Enable output limiter

%er10x_output_limiter_confirm_disable(ER10XHandle); % Confirm output limiter off
er10x_set_mic_response(ER10XHandle, probeIndex, 1); % Flat frequency response

er10x_set_heater_temperature(ER10XHandle,probeIndex,temperatureF);