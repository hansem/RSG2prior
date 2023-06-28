function analRNN0(dataFileN,filen)
% analyzing RNN
% not matching behavior yet
%
% 1. PCA to look at population trajectory
% 2. look at single neurons

% analRNN0('ForHansem2Priors1.mat','RNN_20170302.mat');
% analRNN0('ForHansem2Priors_Train0TestInt001.mat','RNN_20170306_0_001.mat');
% analRNN0('ForHansem2Priors_Train0TestInt005.mat','RNN_20170306_0_005.mat');
% analRNN0('ForHansem2Priors_Train001TestInt001.mat','RNN_20170306_001_001.mat');
% analRNN0('ForHansem2Priors_Train001TestInt005.mat','RNN_20170306_001_005.mat'); %%%
% analRNN0('ForHansem2Priors_Train001TestInt008.mat','RNN_20170306_001_008.mat');

% analRNN0('ForHansem2Priors_Train001TestInt012.mat','RNN_20170308_001_012.mat');
% analRNN0('ForHansem2Priors_Train001TestInt016.mat','RNN_20170308_001_016.mat');
% analRNN0('ForHansem2Priors_Train001TestInt0012.mat','RNN_20170308_001_0012.mat');
% analRNN0('ForHansem2Priors_Train001TestInt001Ip01.mat','RNN_20170308_001_001_01.mat');
% analRNN0('ForHansem2Priors_Train001TestInt005Ip01.mat','RNN_20170308_001_005_01.mat');

% analRNN0('ForHansem2Priors_Noise.mat','RNN_20170428.mat');

idBehAnal=0; % 1; % 0;

% cc;
psthDir='/Users/hansem/Dropbox (MIT)/psthDataHigh/';

%% load data

%
%   Neurons           1x200            960022400  cell             {[200 x   3000]}   trials {neuron x time}
%   ReadyOnset      200x1                   1600  double              [100   200]
%   inputs            1x200              9622400  cell                {[2 x   3000]}
%   output            1x200              4822400  cell                
%   param             1x1                  32592  struct    
%
%                  N: 200
%                 dt: 1
%                tau: 10
%            ntrials: 200
%            ninputs: 1
%               time: 3000
%              times: [1x3000 double]
%             ntimes: 3000
%        burn_length: 50
%            noisesd: 0.0100
%         pulsewidth: 10
%                  g: 1
%              alpha: 1
%         readyonset: 100
%        readyoffset: 200
%         intervalno: 5
%            offbase: 0.4000
%            offsets: [0.3000 0.4000]
%             prior1: [480 560 640 720 800]
%             prior2: [800 900 1000 1100 1200]
%     input_duration: 2950
%              pulse: 20
%        pulseheight: 0.4000
%       ntrials_test: 200
%           noisearb: 1.0000e-04
%                vec: [1x400 double]
%         networkidx: 1
%                  A: 2.8500
%                eta: 3.3000
% load('/Users/hansem/Dropbox (MIT)/Jazlab/Mehrdad & Devika/Stuff for Hansem/ForHansem2Priors1.mat');
% filen='RNN_20170302.mat';

load(['/Users/hansem/Dropbox (MIT)/RNN/simDevika/' dataFileN]);
% load(['/Users/hansem/Dropbox (MIT)/Jazlab/Mehrdad & Devika/Stuff for Hansem/' dataFileN]);

v2struct(param);
th=0.99;

% simple transformation for intput output
y=cell2mat(output')'; % [time x trials]
x=zeros(size(inputs{1},2),ntrials,size(inputs{1},1));
for i=1:ntrials
    x(:,i,:)=inputs{i}'; % [time x trials x node]
end
r=zeros(size(Neurons{1},2),ntrials,N);
for i=1:N
    r(:,i,:)=Neurons{i}'; % [time x trials x neurons]
end
clear inputs outputs Neurons;

% finding ts
ts=zeros(ntrials,1);
dts=diff(x(:,:,2)); % [time x trials] for ts
[idT,idTr]=find(dts>pulseheight/2);
for i=1:ntrials
    ts(i)=diff(idT(idTr==i));
end
% ReadyOnset relative to input onset (50ms)
readyOnset2=zeros(ntrials,1); for i=1:ntrials,readyOnset2(i)=min(idT(idTr==i))+1;end;

% neural data from ReadyOnset (discard before ready)
nTpointMax=time-max(readyOnset2); % 3000-200
r2=zeros(nTpointMax,ntrials,N);
for i=1:ntrials
        r2(:,i,:)=r(readyOnset2(i):(readyOnset2(i)+nTpointMax-1),i,:); % [time x trials x neurons]
end

% finding prior
dPr=mean(x(:,:,1),1); % [time x trials] for prior cue
idPr=dPr>mean(dPr); idPr=idPr(:); % 0 for short

% finding tp
tp=zeros(ntrials,1);
for i=1:ntrials
    tSet=max(idT(idTr==i))+1;
    if ~isempty(find(y(tSet:end,i)>th,1,'first')) % no threshold crossing        
         tp(i)=find(y(tSet:end,i)>th,1,'first')-1; % finding after set
%         tp(i)=find(y(:,i)>th,1,'first')-max(idT(idTr==i));
    else
       tp(i)=-2*ts(i); 
    end
end

% check input, output
 figure; plot(y);ylabel('RNN output'); plotHorizon(gca,[0 1],[]);
xlabel('time from trial onset (ms)');

figure; subplot(2,1,1),plot(squeeze(x(:,:,1))); axis tight;
ylabel('RNN input (prior cue)'); 
subplot(2,1,2),plot(squeeze(x(:,:,2))); axis tight;
ylabel('RNN input (ready,set)'); 
xlabel('time from trial onset (ms)');

% %check neuron & trials
% neuronId=100; trialId=100; figure; plot(squeeze(r2(:,trialId,neuronId)));hold on; plot(mean(r2(:,tmpId,neuronId),2),'k-','linewidth',3);

% % all pupulation activity for single trials
% figure; 
% for i=1:round(ntrials/10):ntrials
%     trialId=i;
%     subplot(2,1,2), imagesc(squeeze(r(:,trialId,:))');colorbar;
%     xlabel('time from trial onset (ms)'); ylabel('neurons');
%     plotVertical(gca,readyOnset2(trialId)+[0 ts(trialId)],[]);
%     subplot(2,1,1),  plot(squeeze(x(:,trialId,:)));
%     ylabel('RNN input');
%     drawnow; 
%     waitforbuttonpress; clf;
% end

%% checking behavior
if idBehAnal
    if ~exist(['/Users/hansem/Dropbox (MIT)/fileFromNHPrig/dataMat/' filen(1:end-4) '.mat'],'file')
        t=tp; T=ts; idShortTrial=1-idPr;
        ballAlpha=zeros(size(t));
        radius=ones(size(t))*10;
        winF=0.2*ones(size(t));
        theta=180*(rand(size(t))>0.5);
        idHandEye=(rand(size(t))>0.5);
        bonusRewDur=zeros(size(t));
        rewardDur=70*ones(size(t));
        fixTimeDur=ReadyOnset;
        targetTimeDur=zeros(size(t));
        dtProdRew=zeros(size(t));
        FeedbackTimeDur=1000*ones(size(t));
        iti=zeros(size(t));
        timeout=zeros(size(t));
        save(['/Users/hansem/Dropbox (MIT)/fileFromNHPrig/dataMat/' filen(1:end-4) '.mat'],'t','T','ballAlpha','radius','idShortTrial','winF','theta','bonusRewDur','rewardDur','fixTimeDur','targetTimeDur','idHandEye',...
            'dtProdRew','FeedbackTimeDur','iti','timeout');
        mkdir('/Users/hansem/Dropbox (MIT)/fileFromNHPrig',filen(1:end-4));
        analHandEye2PriorRSG(filen);
    else
        analHandEye2PriorRSG(filen);
    end
end

%% single neurons
dirSingleRNN=[psthDir filen(1:end-4)];
idPlot=1;
for iN=1:N % 139:139 % 1:N
ws=v2struct;
plotSDF_RNN(ws); % x r r2[time x trial x neurons] y idPr ts tp readyOnset2
close
end
% 
%% population trajectory
% r [time x trials x node]
% ts [trials x1]
% idPr [trials x1] % 0 for short 1 for long
filen_dimred=['dimRed_' filen(1:end-4)];
if ~ exist([psthDir filen_dimred(1:end-4) '.mat']);
    if ~exist([psthDir 'psth_' filen(1:end-4) '.mat']);
        % define D data structure
        nPr=2;    prNm={'Short','Long'};
        nTspp=5;
        T=[];
        T{1}=linspace(500,820,nTspp); T{2}=linspace(820,1220,nTspp); Tall=[T{1} T{2}];
        
        % for now, avg across trial per condition
        
        % trial wise (session specific)
        % ex1_spikecounts: D(7*30).Data[61 400]
        % 	61 neurons, 400ms spike trains, 1ms resolution, 7 conditions, 30 trials/condition
        % D
        %     data
        %     condition
        %     epochStarts
        %     epochColors
        
        % condition wise (across days)
        % ex3_psths: D(2).Data[61 1018]
        % 	61 neurons, 1ms resolution, 2 conditions, so 122 PSTHs in total
        % D
        %     data
        %     traj
        %     epochStarts
        %     epochColors
        %     condition
        
        tmpCmap={flipud(autumn(nTspp));... % , flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
            winter(nTspp)}; % , tmp2(1:nTspp,:)}; % longEye, longHand
        iD=1;
        for i=1:nPr
            for j=1:length(T{i})
                D(iD).condition=[prNm{i} num2str(T{i}(j))];
                tmpId=T{i}(j)==ts & idPr==(i-1);
                D(iD).data=squeeze(mean(r2(:,tmpId,:),2))'; % [time x trials x neurons] > % [time x neuron]'
                % neuronId=100; figure; plot(squeeze(r2(:,tmpId,neuronId)));hold on; plot(mean(r2(:,tmpId,neuronId),2),'k-','linewidth',3);
                D(iD).traj='traj'; % not spiking data
                D(iD).epochColors=[tmpCmap{i}(j,:); tmpCmap{i}(j,:)]; % darker for 2nd phase of each epoch
                D(iD).epochStarts=[1 T{i}(j)];
                iD=iD+1;
            end % j
        end % i
        nCond=length(D); % 40
        save([psthDir 'psth_' filen(1:end-4) '.mat'],'D');
    else
        load([psthDir 'psth_' filen(1:end-4) '.mat']);
    end
    DataHigh(D,'DimReduce');
else
    load([psthDir filen_dimred(1:end-4) '.mat']);
    DataHigh(D);
end

%% fixed point analysis (dynamics)

