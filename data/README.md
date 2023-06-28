# RSG2prior data

## Bayesian Computation through Cortical Latent Dynamics, Sohn et al., Neuron (2019)

This folder has neural data, spike times, for each session.
```
|-- RSG2prior/data
|   |-- behavior/
|       |-- H_RSGprior_DMFC.mat
|       |-- G_RSGprior_DMFC.mat
|   |-- neurons/
|       |-- raw_spike_data
|       |-- processed_data
|       |-- PSTH
|       |-- raster
```
### DETAILS 
```
|-- RSG2prior/data/behavior
```
 Each mat file contains behavioral data from all recording sessions, separately for animal H(Haydn) and G(Gershwin). Important variables are as below:
	- T: sample time interval (t_s) [# trials x 1]
	- t: production time interval (t_p) [# trials x 1]
	- idShortTrial: index for each trials' prior condition (1 if 'short', 0 for 'long') [# trials x 1]
	- idHandEye: index for each trials' response modality (1 if 'eye', 0 for 'hand') [# trials x 1]
	- theta: target location (0 for right, 1 for left) [# trials x 1]
	- idOut: index to indicate each trial is outlier or not (e.g., aborting trials). 1 if a trial is outlier, 0 otherwise [# trials x 1] 
	- sessId: index for session [# trials x 1]

```
|-- RSG2prior/data/neurons/raw_spike_data
```
 All data in 2016 (16XXXX.mat) is from animal H(Haydn), those in 2017 (17XXXX) is from G(Gershwin). G's data collected in 2017/8 included noSet trials.

 - sp [# spikes x 1]: spike times in clock of the blackrock recording system. Given its sample frequency of 30000Hz, this can be converted into spike times in [sec] by dividing it with 30000 (to be compatible with tInfo's time.
 - idClust [# spikes x 1]: for each spike, this index indicates which neuron generated the spike. Its element is 4 digit and the first digit indicates different v probes from which that neurons is recorded (so 1st digits can be only 1XXX, 2XXX, 3XXX).
 - idSU [# neurons x 2]: index indicating whether each neuron is single unit (1 in 2nd column) or not (0).
 - tInfo [# events x 12 (or 14 in G for no-set experiment)]: information of all events and behavior in all trials
	% columns: 
	1) trialId("movingBall_trials"), 
	2) stageCode 
		% 1. Fix On
		% 2. Fixation complete
		% 3. target On
		% 4. Ready
		% 5. Set
		% 6. Production 
		% 7. target acquired 
		% 8. reward delay
		% 9. reward
		% 10. trial end
		% 0. ITI
		% 8. incorrect
		% 0. Bad
	3) time in blackrock clock [sec], 
	4) idShortTrial, 
	5) idHandEye, 
	6) theta, 
	7) T, 
	8) t, 
	9) fixTimeDur, 
	10) targetTimeDur, 
	11) iti, 
	12) reward, 
	13) idNoSet (1 if noSet trials; only for G), 
	14) idSuccessNoSet (1 if not aborting noSet trials)
```
|-- RSG2prior/data/neurons/processed_data
```
- PSTH_periSet_G_full_SUMU.mat
- PSTH_periSet_G_full_SUMU.mat
```
|-- RSG2prior/data/neurons/PSTH
```
 This folder has Peri-Stimulus Time Histogram (PSTH) of a representative cell (YYMMDD_ID). Due to github storage limit, all cells' plot cannot be uploaded. Inside the plot, each panel shows firing rates around task events, sorted by behaviorally relevant trial types.
```
|-- RSG2prior/data/neurons/raster
```
 This folder has raster plot of a representative cell (YYMMDD_ID). Due to github storage limit, all cells' plot cannot be uploaded. Inside the plot, each panel shows spikes around task events. Note that raster is not sorted by behaviorally relevant trial types.
