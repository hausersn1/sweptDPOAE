temperatureF = 86; % 30C
probeIndex = 1; % 0=A/1=B (changed to B 9/15/23)
er10x_set_heater_temperature(ER10XHandle,probeIndex,temperatureF);
err = er10x_disconnect(ER10XHandle);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Could not close ER10X!');
end
[err, ER10XHandle] = er10x_close(ER10XHandle);
if strcmp(err, 'ER10X_ERR_OK')
    fprintf('Continuing...\n');
else
    error('Something wrong! Could not close ER10X!');
end
ER10XCloseAll();
ER10XConnectUnloadLib();
rmpath(p10x);