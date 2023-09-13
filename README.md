# sweptDPOAE

`Make_DPswept.m` generates a sweep of a probe tone and a suppressor, as well as sets some other stimulus parameters. This saves a `stim` struct so which is read by `Run_DPswept`

`Run_DPswept.m` generates the stimuli and plays them with the ER10X. The `stim` struct saves all the raw the data and the stimulus parameters. 

`Analyze_DPswept` will analyze a loaded stimuli structure that contains the data and plot the amplitude and phase. It will also plot the total DPOAE as well as the separated distortion and reflection components. 

`rampsound.m` and `scalesound.m` are needed to make the stimuli. The ER10x files are needed for running this in the SNAP/ARDC lab. 

`simpleAnalysis.m` just does an fft of the response. 
