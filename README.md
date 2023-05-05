# sweptDPOAE

`Make_DPswept_log.m` generates a sweep of a two tones, as well as sets some other stimulus parameters. This saves a `stim` struct so which is read by `Run_DPswept_Auto.m`. Versions with _Auto will stop collecting data once a set SNR criterion is reached. 

`Run_DPswept.m` generates the stimuli and plays them with the ER10X. The `stim` struct saves all the raw the data and the stimulus parameters. 

`Analyze_DPswept` will analyze a loaded stimuli structure that contains the data and plot the amplitude and phase. It will also plot the total DPOAE as well as the separated distortion and reflection components. 

`rampsound.m` and `scalesound.m` are needed to make the stimuli. The ER10x files are needed for running this in the SNAP/ARDC lab. 

## Running for Chins using Heinz Lab NEL set up
Equivalent scripts exist within this repository for running on either the ER10X or the NEL set up. Just run `Run_DPswept_NEL_Auto.m` instead. 
