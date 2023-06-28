function analRNN(dataFileN) % ,filen)

% 2018/5/15
% - remove pulse (20ms)
% - save offset from PCA/PSTH2traj
% - also for tp
%
% remember ts is from ready offset to set onset
% ReadyOn vs readyOffset: estimated directly from input
% ReadyOn-readyOffset=-151(burn_length)
% readyOffset=burn_length(150)+readyOnset+pulse_length
%
% 2018/4/20
% plot state-space trajectory
% cc; analRNN('Weber_2Prior_Type2_000_v3.mat');
% cc; analRNN('Weber_2Prior_Type2_000_v4.mat');
% cc; analRNN('Weber_2Prior_Type2_005_v2.mat');
% cc; analRNN('Weber_2Prior_Type2_005_v5.mat');
%
% cc;analRNN('Weber_2Prior_Type2_005_v4_R2.mat');
% cc;analRNN('NetworkRun_Type2_w005_v4.mat')
% cc;analRNN('NetworkRun_Type2_w005_v3.mat')
% cc;analRNN('NetworkRun_Type2_w000_v1.mat')
%
% cc;analRNN('Weber_2Prior_Type2_005_v2.mat')
% cc;analRNN('Weber_2Prior_Type3_005_v1.mat')


%   ctx_idx       1x500                  4000  double              
%   d             1x1                   16704  struct              
%   inputs        1x500              28056000  cell                
%   ipidx         1x500                  4000  double              
%   net_name      1x30                     60  char                
%   neurons       1x500            2800056000  cell                
%   output        1x500              14056000  cell                
%   param         1x1                   36480  struct   

%% old
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

%%
idBehAnal=1; % 1; % 0;

% cc;
initRSG2prior;
% psthDir='/Users/hansem/Dropbox (MIT)/psthDataHigh/';

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

try
    if exist(['/Users/hansem/Dropbox (MIT)/DN & HS/WeberNetworks/' dataFileN])
        load(['/Users/hansem/Dropbox (MIT)/DN & HS/WeberNetworks/' dataFileN]);
    else
        load(['/Users/hansem/Dropbox (MIT)/DN & HS/WeberNetworks/OLD Networks/' dataFileN]);
    end
catch
    if exist(['/Users/seonminahn/Dropbox (MIT)/DN & HS/WeberNetworks/' dataFileN]);
        load(['/Users/seonminahn/Dropbox (MIT)/DN & HS/WeberNetworks/' dataFileN]);
    else
        load(['/Users/seonminahn/Dropbox (MIT)/DN & HS/WeberNetworks/OLD Networks/' dataFileN]);
    end
end
% load(['/Users/hansem/Dropbox (MIT)/Jazlab/Mehrdad & Devika/Stuff for Hansem/' dataFileN]);

v2struct(param); % N: # neurons
th=0.99; %8; % 1; % 0.99;

% simple transformation for intput output
y=cell2mat(output')'; % [time x trials]
clear output;
x=zeros(size(inputs{1},2),ntrials,size(inputs{1},1));
for i=1:ntrials
    x(:,i,:)=inputs{i}'; % [time x trials x node]
end
clear inputs
r=zeros(size(neurons{1},2),ntrials,N);
for i=1:ntrials
    r(:,i,:)=neurons{i}'; % [time x trials x neurons]
end
clear neurons;

% finding ts
dts=diff(x(:,:,2)); % [time x trials] for ts
[idT,idTr]=find(dts>pulseheight/2);
if exist('d')
    ts=d.ts(:);
else % ts from ready offset to set onset
    ts=zeros(ntrials,1);
    for i=1:ntrials
        ts(i)=diff(idT(idTr==i))-param.pulse;
    end
end
% ReadyOnset relative to input onset (50ms)
readyOffset=zeros(ntrials,1); 
if exist('ReadyOn')
    for i=1:ntrials,readyOffset(i)=param.burn_length+ReadyOn(i)+param.pulse;end;
else
    for i=1:ntrials,readyOffset(i)=min(idT(idTr==i))+1;end;
end

% % for new networks
% if min(ts)~=T{1}(1) % ts not taken into account pulse width
%     ts=ts-pulse;
%     
% %     neural data from ReadyOnset+pulse (discard before ready)
%     nTpointMax=max(Tall);
%     r2=nan(nTpointMax,ntrials,N);
%     for i=1:ntrials
%         r2(1:ts(i),i,:)=r((readyOffset(i)+param.pulse):(readyOffset(i)+ts(i)+param.pulse-1),i,:); % [time x trials x neurons]
%     end
% else
    % neural data from ReadyOnset+pulse (discard before ready)
    nTpointMax=max(Tall);
    r2=nan(nTpointMax,ntrials,N);
    for i=1:ntrials
        r2(1:ts(i),i,:)=r((readyOffset(i)):(readyOffset(i)+ts(i)-1),i,:); % [time x trials x neurons]
    end
% end

% finding prior
idPr=ctx_idx(:)-1; %  0 for short 1 for long
% dPr=mean(x(:,:,1),1); % [time x trials] for prior cue
% idPr=dPr>mean(dPr); idPr=idPr(:); % 0 for short

% check definition of tp: from set offset to t(R>th)
% for i=1:ntrials,tpTmp=[tpTmp;min(find(output{i}((500+ts(i)):end)>1)+500+ts(i)-1)-(readyOffset(i)+ts(i)+pulse)];end;
% figure; plot(tp,tpTmp,'o');plotIdentity(gca);figure; hist(tp-tpTmp,50);

% figure; plot(inputs{1}')
% ha; plot(output{1}')
% plotVertical(gca,readyOffset(1),[]);
% plotVertical(gca,readyOffset(1)+ts(1),[]);
% ha;plot(min(find(output{1}(1000:end)>th)+999),1,'kx');
% plotHorizon(gca,1,[]);

% finding tp
if exist('d')
    tp=d.tp(:);
else
    tp=zeros(ntrials,1);
    for i=1:ntrials
        tSet=max(idT(idTr==i))+1;
        if ~isempty(find(y(tSet:end,i)>=th,1,'first')) % no threshold crossing
            tp(i)=find(y(tSet:end,i)>=th,1,'first')-1; % finding after set
            %         tp(i)=find(y(:,i)>th,1,'first')-max(idT(idTr==i));
        else
            tp(i)=-2*ts(i);
        end
    end
end

% neural data tp: from set onset
nTpointMax=max(tp);
r3=nan(nTpointMax,ntrials,N);
for i=1:ntrials
    setOffset=readyOffset(i)+ts(i)+param.pulse;
    if ~isnan(tp(i)) % not reaching to threshold in Weber_2Prior_Type2_005_v2
        r3(1:tp(i),i,:)=r(setOffset:(setOffset+tp(i)-1),i,:); % [time x trials x neurons]
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
%     plotVertical(gca,readyOffset(trialId)+[0 ts(trialId)],[]);
%     subplot(2,1,1),  plot(squeeze(x(:,trialId,:)));
%     ylabel('RNN input');
%     drawnow; 
%     waitforbuttonpress; clf;
% end



%% single neurons
% dirSingleRNN=[psthDir filen(1:end-4)];
% idPlot=1;
% for iN=1:N % 139:139 % 1:N
% ws=v2struct;
% plotSDF_RNN(ws); % x r r2[time x trial x neurons] y idPr ts tp readyOffset
% close
% end
% % 
%% population trajectory

epNm={'ts','tp','prior'};
tBuffer=80; % for shortest ts

for iEp=1:length(epNm)

    % r2,r3 [time x trials x node]
    % ts [trials x1]
    % idPr [trials x1] % 0 for short 1 for long
    nPSTH=['PSTH_' epNm{iEp} '_' dataFileN];
    nTraj=['traj_' epNm{iEp} '_' dataFileN];
    
    % (1) PSTH
    if ~exist(nPSTH)
        D=[];
        iD=1;
        for i=1:nPr
            for j=1:length(T{i})
                D(iD).condition=[prNm{i} num2str(T{i}(j))];
                tmpId=T{i}(j)==ts & idPr==(i-1);
                switch iEp
                    case 1 % ts
                        D(iD).data=squeeze(mean(r2(1:T{i}(j),tmpId,:),2))'; % [time x trials x neurons] > % [time x neuron]'
                    case 2 % tp: given smoothing activity in RNN, include up to longest tp
                        D(iD).data=squeeze(nanmean(r3(:,tmpId,:),2))'; % [time x trials x neurons] > % [time x neuron]'
                        D(iD).data=D(iD).data(:,~isnan(D(iD).data(1,:))); % truncate nan at end
                    case 3 % 
                        D(iD).data=squeeze(mean(r2((T{i}(1)-tBuffer+1):T{i}(j),tmpId,:),2))'; % [time x trials x neurons] > % [time x neuron]'
                end
                % neuronId=100; figure; plot(squeeze(r2(:,tmpId,neuronId)));hold on; plot(mean(r2(:,tmpId,neuronId),2),'k-','linewidth',3);
                D(iD).traj='traj'; % not spiking data
                D(iD).epochColors=[tmpCmap{i,1}(j,:); tmpCmap{i,1}(j,:)]; % darker for 2nd phase of each epoch
                D(iD).epochStarts=1; % [1 T{i}(j)]; % for now only for ts
                iD=iD+1;
            end % j
        end % i
        save([psthDir nPSTH],'D');
        
        idShortTrial=1-idPr;
        save([psthDir nPSTH],'ts','tp','idShortTrial','param','-append');
        
        if exist('fps') % {20x1}(4).FP
            save([psthDir nPSTH],'fps','-append');
        end
        
        if exist('readyOffset')
            save([psthDir nPSTH],'readyOffset','-append');
        end
    else
        load([psthDir nPSTH]);
    end
    
    
    % (2) traj
    idScree=1; % 0;
    idPlotTraj=1; % 0;
    if ~exist(nTraj)
        binSize=20;
        smthWidth=40;
        use_sqrt=0;
        mean_thresh=-inf; % keeping all units
        
        [D,proj_matrix,keep_neurons,eigenvalues,offset]=PSTH2traj(D,...
            [],[],binSize,smthWidth,[],use_sqrt,mean_thresh); %  proj_matrix,keep_neurons optimD
        optimD=size(D(1).data,1);
        save([psthDir nTraj],'D','binSize','smthWidth','use_sqrt','optimD','proj_matrix','keep_neurons','eigenvalues','mean_thresh','offset');
        
        if idPlotTraj
            % plot
            plotTraj_figure3(nTraj); % plotTraj_SFN17(fnTmp);
        end
        
        
        
        if idScree % eigenvalues
            msize=2;
            lw=1;
            
            nPCmax=10;
            thVar=.75;
            
            figure;ha;setFigPos(1,6);
            plot(1:nPCmax,cumsum(eigenvalues(1:nPCmax))./sum(eigenvalues)*100,'o','color','k','markerfacecolor','w','markeredgecolor','k',...
                'markersize',msize,'linewidth',lw);
            optD=find(cumsum(eigenvalues(1:20))./sum(eigenvalues)>thVar,1,'first');
            disp(optD);
            xlabel('Number of PCs'); ylabel('Cum. var. explained');
            set(gca,'ytick',0:25:100,'yticklabel',{0;[];50;[];100},'ylim',[0 100],'xlim',[0 nPCmax],...
                'xtick',[0:3:nPCmax],'tickdir','out','ticklength',[.02 .02]);
        end % if idScree

        idShortTrial=1-idPr;
        save([psthDir nTraj],'ts','tp','idShortTrial','param','-append');
        
        if exist('fps') % {20x1}(4).FP
            save([psthDir nTraj],'fps','-append');
        end
        
        if exist('readyOffset')
            save([psthDir nTraj],'readyOffset','-append');
        end
        
    else
        load([psthDir nTraj]);
        plotTraj_figure3(nTraj); % plotTraj_SFN17(fnTmp);
        
    end  % if ~exist(nTraj)
    
end % for iEp=1:length(epNm)

%% OLD
% if ~ exist([psthDir filen_dimred(1:end-4) '.mat']);
%     if ~exist([psthDir 'psth_' filen(1:end-4) '.mat']);
%         % define D data structure
%         nPr=2;    prNm={'Short','Long'};
%         nTspp=5;
%         T=[];
%         T{1}=linspace(500,820,nTspp); T{2}=linspace(820,1220,nTspp); Tall=[T{1} T{2}];
%         
%         % for now, avg across trial per condition
%         
%         % trial wise (session specific)
%         % ex1_spikecounts: D(7*30).Data[61 400]
%         % 	61 neurons, 400ms spike trains, 1ms resolution, 7 conditions, 30 trials/condition
%         % D
%         %     data
%         %     condition
%         %     epochStarts
%         %     epochColors
%         
%         % condition wise (across days)
%         % ex3_psths: D(2).Data[61 1018]
%         % 	61 neurons, 1ms resolution, 2 conditions, so 122 PSTHs in total
%         % D
%         %     data
%         %     traj
%         %     epochStarts
%         %     epochColors
%         %     condition
%         
%         tmpCmap={flipud(autumn(nTspp));... % , flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
%             winter(nTspp)}; % , tmp2(1:nTspp,:)}; % longEye, longHand
%         iD=1;
%         for i=1:nPr
%             for j=1:length(T{i})
%                 D(iD).condition=[prNm{i} num2str(T{i}(j))];
%                 tmpId=T{i}(j)==ts & idPr==(i-1);
%                 D(iD).data=squeeze(mean(r2(:,tmpId,:),2))'; % [time x trials x neurons] > % [time x neuron]'
%                 % neuronId=100; figure; plot(squeeze(r2(:,tmpId,neuronId)));hold on; plot(mean(r2(:,tmpId,neuronId),2),'k-','linewidth',3);
%                 D(iD).traj='traj'; % not spiking data
%                 D(iD).epochColors=[tmpCmap{i}(j,:); tmpCmap{i}(j,:)]; % darker for 2nd phase of each epoch
%                 D(iD).epochStarts=[1 T{i}(j)];
%                 iD=iD+1;
%             end % j
%         end % i
%         nCond=length(D); % 40
%         save([psthDir 'psth_' filen(1:end-4) '.mat'],'D');
%     else
%         load([psthDir 'psth_' filen(1:end-4) '.mat']);
%     end
%     DataHigh(D,'DimReduce');
% else
%     load([psthDir filen_dimred(1:end-4) '.mat']);
%     DataHigh(D);
% end

%% fixed point analysis (dynamics)


%% checking behavior
if idBehAnal
%     if ~exist(['/Users/hansem/Dropbox (MIT)/fileFromNHPrig/dataMat/' filen([1:(end-4)]) '.mat'],'file')
%         t=tp; T=ts; idShortTrial=1-idPr;
%         ballAlpha=zeros(size(t));
%         radius=ones(size(t))*10;
%         winF=0.2*ones(size(t));
%         theta=180*(rand(size(t))>0.5);
%         idHandEye=(rand(size(t))>0.5);
%         bonusRewDur=zeros(size(t));
%         rewardDur=70*ones(size(t));
%         fixTimeDur=ReadyOnset;
%         targetTimeDur=zeros(size(t));
%         dtProdRew=zeros(size(t));
%         FeedbackTimeDur=1000*ones(size(t));
%         iti=zeros(size(t));
%         timeout=zeros(size(t));
%         save(['/Users/hansem/Dropbox (MIT)/fileFromNHPrig/dataMat/' filen([1:(end-4)]) '.mat'],'t','T','ballAlpha','radius','idShortTrial','winF','theta','bonusRewDur','rewardDur','fixTimeDur','targetTimeDur','idHandEye',...
%             'dtProdRew','FeedbackTimeDur','iti','timeout');
%         mkdir('/Users/hansem/Dropbox (MIT)/fileFromNHPrig',filen([1:(end-4)]));
%         analHandEye2PriorRSG(filen);
%     else
%         analHandEye2PriorRSG(filen);
%     end
    
% check FP
% cc; load('traj_ts_Weber_2Prior_Type2_000_v3.mat');disp('traj_ts_Weber_2Prior_Type2_000_v3'); % all non-nan FP
% for i=1:20,for j=1:4,disp(nnz(~isnan(fps{i}(j).FP)));end;end
% cc; load('traj_ts_Weber_2Prior_Type2_000_v4.mat');disp('traj_ts_Weber_2Prior_Type2_000_v4'); % no FP
% for i=1:20,for j=1:4,disp(nnz(~isnan(fps{i}(j).FP)));end;end
% cc; load('traj_ts_Weber_2Prior_Type2_005_v2.mat');disp('traj_ts_Weber_2Prior_Type2_005_v2'); % all non-nan FP
% for i=1:20,for j=1:4,disp(nnz(~isnan(fps{i}(j).FP)));end;end
% cc; load('traj_ts_Weber_2Prior_Type2_005_v5.mat');disp('traj_ts_Weber_2Prior_Type2_005_v5'); % had some
% for i=1:20,for j=1:4,disp(nnz(~isnan(fps{i}(j).FP)));end;end

cc;
initRSG2prior;
Tmat=T;
% load('traj_ts_Weber_2Prior_Type2_005_v5.mat');
% load('traj_ts_Weber_2Prior_Type2_005_v2.mat');
% load('traj_ts_Weber_2Prior_Type2_000_v4.mat');
load('traj_ts_Weber_2Prior_Type2_000_v3.mat');

ctidx=2-idShortTrial; % 1 for short, 2 for long
ipidx=nan(size(ctidx));
d.ts=T;
d.tp=t;

% data formatting: [# interval x maxRepetition]
% find max # trials 1st
for i=1:nPr
    for k=1:nTspp
        id=idShortTrial==(2-i) &... % idShortTrial
            T==Tmat{i}(k); % ts
        if exist('nMaxTrial')
            nMaxTrial=max([nMaxTrial; nnz(id)]);
        else
            nMaxTrial=nnz(id);
        end
        ipidx(id)=k; % for each trial in d.ts specifies which time bin or ts interval it belongs to. Values: 1-5
    end
end
% filling matrix
Ts=nan(length(Tall),nMaxTrial); % [K x R]
Tp=nan(length(Tall),nMaxTrial);
for i=1:nPr
    for k=1:nTspp
        id=idShortTrial==(2-i) &... % idShortTrial
            T==Tmat{i}(k); % ts
        Ts((i-1)*nTspp+k,1:nnz(id))=T(id);
        Tp((i-1)*nTspp+k,1:nnz(id))=t(id);
    end
end
d.Ts=Ts;
d.Tp=Tp;
% theory

    wm=.05; % wFit.w_m;
    wp=.05; % wFit.w_p;
    offset=[0;0]; % [wFit.offset2; wFit.offset1];

for i=1:nPr
    theory{i}.tm=[min(Tmat{i}):max(Tmat{i})]'; % theory{x}.tm - Tx1 - sample corresponding to above
    theory{i}.te=BLS(theory{i}.tm,wm,[min(Tmat{i}) max(Tmat{i})],'uniform')+offset(i);theory{i}.te=theory{i}.te(:);
    theory{i}.range=1:length(theory{i}.tm);
    theory{i}.wm=wm;
end
DnPlotResponses2P(d,theory,ctidx,ipidx);

% paper
% set(gca,'xtick',400:400:1200,'ytick',400:400:1200,'ylim',[300 1450]); 
xlim([450 1250]);
set(gca,'xtick',[480 640 800 1000 1200],'ytick',[480 640 800 1000 1200],'xticklabel',[],'yticklabel',[]);
ylim([350 1400]);
xlabel([]);ylabel([]);
title([])
    
    
    
end
