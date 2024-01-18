# RSG2prior data/code

### ["Bayesian Computation through Cortical Latent Dynamics", Sohn H*, Narain D*, Meirhaeghe N*, Jazayeri M , Neuron (2019)](https://www.sciencedirect.com/science/article/pii/S0896627319305628)

### Folder structure is as follows:
```
|-- figures
|   |-- plot_scripts/
|   |-- Movies/
|   |-- PDF Figures/
|   |-- PNG Figures/
|-- data
|   |-- behavior/
|   |-- neurons/
|       |-- raw_spike_data
|       |-- PSTH
|       |-- raster
|       |-- processed_data
|-- task
|-- analysis
|   |-- utils/
```
### SET UP GUIDELINE 

1. Change base_directory into your 'RSG2prior' directory in initRSG2prior.m
    (e.g., '/Users/hansem/Documents/RSG2prior')
2. add all subfolders of 'RSG2prior' directory to MATLAB path, particularly /analysis/utils
3. open and run MATLAB live scripts in figure/plot_scripts 
    - in many codes, processed data is used for figures
    - note that the codes are not fully cleaned up
4. Start from raw_spike_data to reproduce results!


### DETAILS


#### RSG2prior data
```
|-- data
|   |-- behavior/
|       |-- H_RSGprior_DMFC.mat
|       |-- G_RSGprior_DMFC.mat
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
|       |-- H_20190315_RSGadapt.mat: animal H's generalization test data (Figure S1F)
|       |-- G_20190322_RSGadapt.mat: animal G's generalization test data (Figure S1F)
```
```
|-- data
|   |-- neurons/
|       |-- raw_spike_data
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
Many .mat files here are used to run plot_scripts. Some partial information about what these files are:
    - 'trajKS' refers to neural trajectory after spike sorting by kilosort. It's mostly in DataHigh format.
    - 'PSTH' is a trial-averaged, smoothed neural firing rate
    - relevant task epoch are prior support ('prior'), measurement ('ts'), production ('tp'), or period from ready flash to go ('periSet')
    - in many cases, data are processed in a condition- and animal-specific manner (e.g., ER: eye right, HL: hand left, _H_: animal H, _G_: animal G)
    - 'avgAttB4PCA' means that data are averaged with attrition before PCA
```
|-- RSG2prior/data/neurons/PSTH
```
 This folder has Peri-Stimulus Time Histogram (PSTH) of a representative cell (YYMMDD_ID). Due to github storage limit, all cells' plot cannot be uploaded. Inside the plot, each panel shows firing rates around task events, sorted by behaviorally relevant trial types.
```
|-- RSG2prior/data/neurons/raster
```
 This folder has raster plot of a representative cell (YYMMDD_ID). Due to github storage limit, all cells' plot cannot be uploaded. Inside the plot, each panel shows spikes around task events. Note that raster is not sorted by behaviorally relevant trial types.
```
|-- task
```
RSG_twoPrior_handEye.xml is an old MWorks code; caution: it is not cleaned up and largely not human-readable without MWorks editor, which is not supported anymore (as of 2023 Jun 17).
```
|-- analysis
```
: analysis code used to process data and plot figures; caution: it is not cleaned up.
|--
#### GPFA (Gaussian-Process Factor Analysis; Yu et al., J Neurophysiol, 2009)

GPFA was run for representative sessions (animal H's 161218; animal G's 170511) using [DataHigh](https://github.com/BenjoCowley/DataHigh) toolbox. Used data is in RSG2prior/data/neurons/processed_data/[PSTH_161218_ReadySetGo_single_trial.mat](https://github.com/hansem/RSG2prior/blob/main/data/neurons/processed_data/PSTH_161218_ReadySetGo_single_trial.mat) and in the format compatible with DataHigh toolbox. Briefly, the MATLAB file has a structural array, 'D', where each element is from single trials and its data field has a 2D matrix, 53 neurons x variable-size time bins (1ms). Each bin has 1 if there is a spike, 0 otherwise. EpochStarts field indicates when Ready and Set stimulus is given and condition field is each trials' condition identifier, indicating prior (short/long), ts (from 480 to 1200), modality (eye/hand), direction (left/right). If you are not using DataHigh toolbox (e.g., for scripting with GPFA's original [Matlab code](https://users.ece.cmu.edu/~byronyu/software/gpfa0203.tgz)), you can always generate your own data format from raw spike train data in RSG2prior/data/neurons/raw_spike_data.

Our paper's figure related to GPFA is in RSG2prior/figures/plot_scripts/figure6.mlx and associated data (i.e., GPFA's output) are in the followings:

RSG2prior/data/neurons/processed_data/[trajKS_161218_pr_H_condSpec_GPFA_CV_cellFullOnly.mat](https://github.com/hansem/RSG2prior/blob/main/data/neurons/processed_data/trajKS_161218_pr_H_condSpec_GPFA_CV_cellFullOnly.mat)

RSG2prior/data/neurons/processed_data/[trajKS_170511_pr_G_condSpec_GPFA_CV_cellFullOnly.mat](https://github.com/hansem/RSG2prior/blob/main/data/neurons/processed_data/trajKS_170511_pr_G_condSpec_GPFA_CV_cellFullOnly.mat)




