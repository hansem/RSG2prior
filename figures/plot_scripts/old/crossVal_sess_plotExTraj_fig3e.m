function crossVal_sess_plotExTraj_fig3e(varargin)

% 2018/7/30
% save traj with leave-1trial-out
% t# 102 (idTrial)
% fig3d1: rotating trajectory+readout vector in state space
% fig3d2: IC+readout vector in state space
% additional smoothing for visualization: idSmoothAfterAvgAtt
%   modifed from crossVal_sess.m: no parfor, expand plotting script

% crossVal_sess(1,12); crossVal_sess(2,1); crossVal_sess(1,16);             crossVal_sess(2,8); crossVal_sess(1,8);

% - for IC, PCA use only data around tIC(200ms) with tBuffer (only applied to tp, not pr): idUseTIC
% - estimating readout vector from PC of all sessions: idAllSessReadout
% debug for iNotFullRecord=cellId(:,3)-cellId(:,2)<nTrialIncNeuron; % >=nTrialIncNeuron;

% 2018/7/19
% after avgAttB4PCA (e.g. trajKS_prior_EL_H_supportOnly_avgAttB4PCA)
% other sessions: (criteria:At least, 40 cells with at least 2000 trials)
%         H: 161218 161222 (161211 >55 cells with >1700 trials)
%         G: 170817 170822 170818 170823 170511 170821 170507 170506 (more trials<>more neurons)
%               H/G: sess#12 16 (8) / 8 11 9 12 5 10 2 1 for varargin
%          if idNotFullRecord==0, include neurons with at least 2000 trials (nTrialIncNeuron=2000)
% clean up codes

% 2018/5/15
% remove cells without full recording all trials '_cellFullOnly.mat'  & ~idNotFullRecord
% H again
% parallerize with parfor iTrial=1:nTrial
%
% 2018/5/11
% G 5/6
%
% 2018/5/8
% usually PCs higher than 3rd are noisy (variable across ts during measurement)
% projection weights b/t sessPCA vs allSessPCA are different
% rotation visible with idPCpoolSess=1 for ex. iTrial=102
%
% 2018/5/3
% improving code using PSTHsess (not from raw data): idUseRawData  (N.B. if avg. from trialwise PSTHsess, inc. nonrecord)
% use all sessions to constrain projection weights: idPCpoolSess
% use scores for readout vector: idUseScoreReadout (N.B. nonrecord removed for PSTH to estimate PC scores)
% (TBD) if idUseScoreReadout==0, readout vector estimated from simple avg. (inc. nonrecord)
%
% trajKS_161218_pr_H_condSpec_sessPCA_CV_avgAtt_newCenter: idIncEarly=0(no ts in filename)? idUseScoreReadout=1?
% trajKS_161218_t_p_H_condSpec_sessPCA_CV_avgAtt_newCenter:
%
% 2018/5/2
% correcting centering for PSTH2traj.m (XX_newCenter.mat)
%
% 2018/4/20
% complete CV: estimating readout vector also should be done only from PC(allTrial-1)
% projection should be done independently across trials
% saving sSet (1D projection), sIC (1D projection), xro(readout vector), projection weight W for all trials
% (TBD) check corr b/t sSet & sIC, how much xro/W changes across trials
% depends on estReadout.m
%
% 2018/4/17
% debugging: 1) didn't remove cells<nTrailExc 2) PCA was also specific to prior
% PCA for [RS] rather than priorSupport (in avgTraj, no rotation): idIncEarly
% (TBD) use tBuffer for start & end (as done in runPSTH2traj_prior_conditionSpecific)

% 2018/3/20
% PCA with leave-one-trial-out cross validation for sessions
% period: priorSupport, tp
% choosing dimension up to PVAF>75%
% condition-specific
% combining runExportDataHigh_conditionSpecific.m & runPSTH2traj_ts_tp_conditionSpecific.m
% for now, only H 161218

%% aux. note
% data used
% - raw session data: e.g. 161218.mat
% - behavioral data: H_DMFC_RSG2prior.mat
% - other session PSTH: e.g. PSTH_t_s_H_full_SUMU.mat (dOther; used if idPCpoolSess; also for ID)
% - session PSTH: e.g. PSTH_161218_periSet_H_avgDir_SUMU.mat (psthSess used if ~idUseScoreReadout & if ~idUseRawData)
% 
% dMatK=D0; % population PSTH w/o test trial
% tmpD=[]; tmpD=dMatK; % {40x1}.data [53 x 384]
% tmpDold=[];tmpDold=tmpD;
% 
% tmpD1=[]; % a cell's PSTH w/o test trial; 40 cells [1 x time] 
% tmpD2=[]; % tmpD2=dOther.D; % clear dOther;  % parfor
% tmpD=tmpD2; % if idPCpoolSess
% 
% tmpDtraj=[]; % score from PSTH2traj(tmpD)
% Dtraj=[];

%% init
idTrial=175; % 152; % 175; % short hand - best for IC
% 102;
% 16; % theta0
%     36;
%     58;
%    986;
idTheta=180; % 0; % 180; % used for example trials
   %520; 
%    256 %%
%          989
%         1150
   % 1169;
   
   % short800
%             102 %%
%          137
%          385
%          759
%          967
%         1342
%         1501
%         1538
idCheckTraj=  1;

idPCpoolSess=1; % 0; % 1; % session-specific PCA if 0
idIncEarly= 0; % 1; % 0; %  if 1, include early measurement
idUseScoreReadout=1; % 0; use score from leave-one-trial-out PSTH> if 0, use psthSess w/o cross validation to estimate readout
idAllSessReadout= 0; % 1; % estimating readout vector from PC of all sessions: idAllSessReadout
idUseAvgAtt=1; % avgAtt for estimating readout vector
idAvgAttB4PCA=1;
idUseTIC=0; % 1;
% (keep for ref) idUseRawData=1; % 0; % PSTH.mat can't be used as it didn't take care of non-record period (0 in PSTH means no-recording or no-activity in record?)
idNotFullRecord=0; % 1; % 0; % 1; % if 1, including neurons with not fully recorded
if idNotFullRecord==1,
    nTrialIncNeuron_buffer=15; % threshold for including neurons with (total # trials - 15 trials)
    fnEnd=[]; 
    iNotFullRecord=[]; % parfor
else fnEnd='_cellFullOnly'; 
nTrialIncNeuron_buffer=15; % 0; % 15; % threshold for including neurons with (total # trials - 15 trials)
 iNotFullRecord=[]; % parfor
end
if idAllSessReadout, fnEnd=[fnEnd '_allSessReadout']; end
if idUseTIC, fnEnd=[fnEnd '_idUseTIC']; end
% neurons may not fire some early/late trials

idSmoothAfterAvgAtt=1;

param.idPCpoolSess=idPCpoolSess; % 0; % 1; % session-specific PCA if 0
param.idIncEarly= idIncEarly; % 1; %  if 1, include early measurement
param.idUseScoreReadout=idUseScoreReadout; % 0;
param.idUseAvgAtt=idUseAvgAtt;
% param.idUseRawData=idUseRawData;
param.idNotFullRecord=idNotFullRecord;
param.nTrialIncNeuron_buffer=nTrialIncNeuron_buffer;
param.idAvgAttB4PCA=idAvgAttB4PCA;
param.idAllSessReadout=idAllSessReadout;
param.idUseTIC=idUseTIC;

% deal with varargin
if nargin==2
    iAnimal=varargin{1};
    iSess=varargin{2};
else % default: H 161218
    iAnimal=1;
    iSess=12;
end

initRSG2prior; % psthDir neuDir2 animalNm prNm T ehNm targNm
% for parfor, variables need to be defined
T=T;
Tall=Tall;
animalNm=animalNm;
behDir=behDir;
bootDir=bootDir;
% cname=cname;
% cnameG=cnameG;
% cnameH=cnameH;
condNm=condNm;
dPCAdir=dPCAdir;
ehNm=ehNm;
epochNm=epochNm;
fNeuronDir=fNeuronDir;
figDir=figDir;
fn=fn;
fnTmp=fnTmp;
fnameG=fnameG;
fnameH=fnameH;
iAnimalTmp=iAnimalTmp;
iCond=iCond;
iEH=iEH;
iPr=iPr;
iTarg=iTarg;
iTs=iTs;
nAnimal=nAnimal;
nCell=nCell;
nEH=nEH;
nPr=nPr;
nSplit=nSplit;
nTarg=nTarg;
nTspp=nTspp;
neuDir=neuDir;
neuDir2=neuDir2;
optsAxes=optsAxes;
optsExpFig=optsExpFig;
optsFont=optsFont;
optsMarker=optsMarker;
pplot=pplot;
prNm=prNm;
psthDir=psthDir;
recordSiteDir=recordSiteDir;
sessDate=sessDate;
sessDir=sessDir;
singleNeuronDir=singleNeuronDir;
singleTDir=singleTDir;
sortQ=sortQ;
splitNm=splitNm;
stimDir=stimDir;
targNm=targNm;
tmpCmap=tmpCmap;
tmpCmap2=tmpCmap2;
tmpCmapOld=tmpCmapOld;

psthSess=[];

% PCA
binSize=20;
smthWidth=40;
use_sqrt=0;

tIC=200; iIC=round(tIC/binSize);
param.tIC=tIC;

%% runExportDataHigh_conditionSpecific.m 
% 2017/12/8
% condition-specific (effector x direction)
% for both H/G
% w/o outlier (noSet in G included in outlier)
idSUonly=0; % 1; % 0; % 1;

% exclusion crteria in exportDataHigh
nTrialExc=5; % if # trials per condition<nTrialExc, cell is removed from analysis
tBuffer=80; % include 80ms before shortest ts of each prior
    
%%
% initRSG2prior; % animalNm tmpCmap from pplot.mat
cd([neuDir2]); %     neuDir2='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/neuralData/kilosort';

nAnimal=length(animalNm);

% for iAnimal=1:1 %nAnimal:nAnimal %1:1 % nAnimal %nAnimal:nAnimal % 1:nAnimal % length(animalNm):(-1):1 
    % making D empty
    iD=1;
    for i=1:length(prNm) % now 3
        for j=1:length(T{i})
            for k=1:nEH
                for m=1:nTarg
                    D0(iD).condition=[prNm{i} num2str(T{i}(j)) '/' ehNm{k}  '/' targNm{m}]; %  '/' targNm{m}];
                    D0(iD).data=[];
                    D0(iD).id=[]; % sessionId unitId idSU
                    D0(iD).epochColors=[tmpCmap{i,1}(j,:)];
                    iD=iD+1;
                end % m
            end % k
        end % j
    end % i
%     D0=D;
    
    % epoch
    % 1) [fix-100 fix+500]
    % 2) [tOn-500 tOn+250]
    % 3) [ready-250 ready+ts]
    % 4) [set-ts set+tp]
    % 5) [production-tp production+100]
    % 6) measurement (ts)
    % 7) production (tp)
    % 8) prior support (pr)
    epochNm={'pr','t_p'}; % {'fixation','targetOn','periReady','periSet','production','t_s','t_p','pr'};
    
    epochId=3:4; % 7:8; 
    nEpoch=length(epochId);
    
    %% main
    % session
%     for iSess=12:12 %1:1 % 12:12 % 1:length(fn{iAnimal}) 161218 170506
        dMat=cell(nEpoch,1);
        for iEp=1:nEpoch
            dMat{iEp}=D0;
        end
%         tmpD=D0;
        
        fnm=fn{iAnimal}{iSess}; % H_RSGprior_20161203
        disp(['===== ' fnm ' =====']);
        fid=fnm(end-5:end); % '161203'
        fidData=load(fid); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo
        sp=fidData.sp;
        idClust=fidData.idClust;
        idSU=fidData.idSU;
        tInfo=fidData.tInfo;
        NEV=kilosort2NEV(sp,idClust,idSU);  % parfor
        NEV0=NEV;
        
        nTrialIncNeuron=max(tInfo(:,1))-nTrialIncNeuron_buffer; % 2097-15 for 161218
        
        % check start/end trials of each neurons: except 12 (2123,1174,1172,1170,1167,1164,1162,1159,1157,1147,1146,1138), almost all
        if ~idNotFullRecord
            % neurons are recorded most of the time
            %         for iidSU=1:size(idSU,1)
            %             disp([idSU(iidSU,1) ...
            %                 tInfo(find(sp(find(idClust==idSU(iidSU,1),1,'first'))/30000*1000<tInfo(:,3)*1000,1,'first'),1) ...
            %                 tInfo(find(sp(find(idClust==idSU(iidSU,1),1,'last'))/30000*1000>tInfo(:,3)*1000,1,'last'),1)]);
            %         end
            cellId=[]; % [cellId(XXXX) startTrial endTrial]
            for iidSU=1:size(idSU,1)
                cellId=[cellId;...
                    idSU(iidSU,1) ...
                    tInfo(find(sp(find(idClust==idSU(iidSU,1),1,'first'))/30000*1000<tInfo(:,3)*1000,1,'first'),1) ...
                    tInfo(find(sp(find(idClust==idSU(iidSU,1),1,'last'))/30000*1000>tInfo(:,3)*1000,1,'last'),1)];
            end
%             iNotFullRecord=(cellId(:,2)>12 | cellId(:,3)<2095); % 13 among 53 for H161218
            iNotFullRecord=cellId(:,3)-cellId(:,2)<nTrialIncNeuron; % >=nTrialIncNeuron;
        end
        
        % extract behavior
        beh=load([behDir animalNm{iAnimal} '_RSGprior_DMFC.mat']);
        
        % remove outlier
        idOut2=false(size(tInfo,1),1);  % [17913 x 1]
        movingBallTrials=unique(tInfo(:,1)); % size same as idNoise
        idSess=str2num(fid)==beh.sessId;
        if strcmp(fid,'161207'),idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
        if nnz(idSess)~=length(movingBallTrials),disp('smth wrong for trial align'); end;
        iOutSess=find(beh.idOut(idSess)); %disp(['outlier trials (incl. abort/wait): ' num2str(length(iOutSess)) '/' num2str(length(movingBallTrials)) ',' ...
        for iOut=1:length(iOutSess)
            idOut2(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
        end
        iNoOut2=unique(tInfo(~idOut2,1)); % [1610 x 1] 2 3 7 8 ...
        nTrial=length(iNoOut2);
        clear beh;
        
        for k=1:nEpoch %% 1:1 % nEpoch % nEpoch:(-1):1 % nEpoch:nEpoch % 1:nEpoch %%%%%
            
            %% leave-one-trial-out cross validation
            D(nTrial).data=[];
            D(nTrial).id=[];
            D(nTrial).condition=[];
            D(nTrial).epochColors=[];
            D(nTrial).epochStarts=[];
            D(nTrial).tp=[];
            D(nTrial).idShortTrial=[];
            D(nTrial).ts=[];
            D(nTrial).idHandEye=[];
            D(nTrial).theta=[];
            
            % saving sSet (1D projection), sIC (1D projection), xro(readout vector), projection weight W for all trials
            D(nTrial).W=[]; % projection weight [#neurons x # PC]
            D(nTrial).xro=[]; % [# PC x 1] normalized
            if k==1
                D(nTrial).sSet=[];
            else
                D(nTrial).sIC=[];
            end
            D(nTrial).offset=[];
            D(nTrial).keep_neurons=[];
            D(nTrial).Dcv=[];
            
            % load PSTH for no need for estimating PSTH of each trials (N.B. if avg. from trialwise PSTHsess, inc. nonrecord)
            if ~idUseScoreReadout
                psthSess=load(['PSTH_' fid '_periSet_' animalNm{iAnimal} '_avgDir_SUMU.mat']); % D(1:1610).data[53 x ts+tp]
                % truncate from periSet
                for iTrial=1:nTrial
                    if epochId(k)==3 % prior
                        if idIncEarly
                            psthSess.D(iTrial).data=psthSess.D(iTrial).data(:,1:(psthSess.D(iTrial).ts));
                        else
                            psthSess.D(iTrial).data=psthSess.D(iTrial).data(:,(T{2-psthSess.D(iTrial).idShortTrial}(1)-tBuffer+1):(psthSess.D(iTrial).ts));
                        end
                    else
                        psthSess.D(iTrial).data=psthSess.D(iTrial).data(:,(psthSess.D(iTrial).ts):end);
                    end
                end
            end
            % idPCpoolSess
            if k==1 % ts
                dOther=load(['PSTH_t_s_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D
            else
                dOther=load(['PSTH_t_p_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D
            end
            
            for iTrial=idTrial:idTrial % 1:nTrial % round(rand*(nTrial-1))+1 % idTrial:idTrial
%             delete(gcp('nocreate'));
%             parpool; % (1);
%             parfor iTrial=1:nTrial % length(iNoOut2) 
% % parfor iTrial=148:168 % 1:nTrial % length(iNoOut2) 
% % iTrial=102; % short800eyeLeft
% % disp('1');
dMatK=D0; % parfor
% % NEV=kilosort2NEV(sp,idClust,idSU);  % parfor
% disp('2');
tmpDold=[];
tmpDtraj=[];
Dtraj=[];
tmpD=[];
tmpD1=[];
tmpD2=[];


                disp(['----- Trial' num2str(iTrial) '/Trial' num2str(nTrial) '-----']);
                tmpIdOut2=idOut2; % [17913 x 1]
                trialId=tInfo(:,1)==iNoOut2(iTrial); % [17913 x 1] 1 for ~10 entries
                tmpIdOut2(trialId)=true; % setting a trial outlier
                
                % tInfo: trial#(movingBall_trials; not complete), stageCode, t(blackRock),
                %           idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)]

                % setting up other fields in D(iTrial)
                iTrialId=find(trialId,1,'first'); % picking only one time stampt of a trial
                tmpTs=tInfo(iTrialId,7); D(iTrial).ts=tmpTs;
                tmpIdPr=2-tInfo(iTrialId,4); D(iTrial).idShortTrial=2-tmpIdPr;
                D(iTrial).idHandEye=tInfo(iTrialId,5); 
                D(iTrial).theta=tInfo(iTrialId,6); 
                D(iTrial).tp=tInfo(iTrialId,8);
                D(iTrial).condition=[prNm{tmpIdPr} ...
                    num2str(tmpTs) ...
                    ehNm{2-D(iTrial).idHandEye} ...
                    targNm{(D(iTrial).theta==180)+1}];disp(D(iTrial).condition);
                D(iTrial).epochColors=tmpCmap{tmpIdPr,1}(tmpTs==T{tmpIdPr},:);
                D(iTrial).epochStarts=1;

                %% only checking traj if eye left
%                 if D(iTrial).idHandEye==0 & tmpIdPr==1 % & D(iTrial).theta==idTheta % 180 %
%                     idCheckTraj=1;
%                 else
%                     idCheckTraj=0;
%                 end
%                 if D(iTrial).idHandEye==0 & tmpIdPr==1 % D(iTrial).theta==idTheta % 180 % 
%                      disp(['      ts: ' num2str(D(iTrial).ts) ', tp: ' num2str(D(iTrial).tp)]);
% %                     idCheckTraj=1;
% %                 else
% %                     idCheckTraj=0;
% %                 end
                
                % loop through each unit
                nUnit=size(idSU,1); % both SU & MU
                for j=1:nUnit
                    if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
                    else
                        % consistent with plotSDFexport
                        
                        tmpId=[idSU(j,1) idSU(j,1) 1 length(movingBallTrials)]; % id(electrode,Unit,start trial,end trial)
%                         tic
                        tmpD1=exportDataHigh_conditionspecific(NEV0,tmpId,tInfo(~tmpIdOut2,:),epochId(k)); % 40 cells [1 x time] idOut2 NEV
%                         toc
                        if isempty(tmpD1) % &
%                             disp(['cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust)) '     # trial/condition < ' num2str(nTrialExc)]);
                            % still adding data to match matrix size
                            for m=1:(nPr*nEH*nTspp*nTarg) % length(tmpD1) % cond
                                tmpTs2=Tall(ceil(m/4));
                                if epochId(k)==4 %tp
                                    if isempty(dMatK(m).data)
                                        dMatK(m).data=nan(1,tmpTs2*2); % zeros
                                    else
                                        dMatK(m).data=[dMatK(m).data; nan(1,size(dMatK(m).data,2))]; % zeros
                                    end
                                elseif epochId(k)==3 % prior
                                    if strfind(dMatK(m).condition,'Short')
                                        if idIncEarly
                                            tInit=251;
                                        else
                                            tInit=251+T{1}(1)-tBuffer;
                                        end
                                        tEnd=251+tmpTs2;
                                    elseif strfind(dMatK(m).condition,'Long')
                                        if idIncEarly
                                            tInit=251;
                                        else
                                            tInit=251+T{2}(1)-tBuffer;
                                        end
                                        tEnd=251+tmpTs2;                                        
                                    end % if strcmp(dMatK(m).condition,'Short')
                                    dMatK(m).data=[dMatK(m).data; nan(1,tEnd-tInit)]; % zeros
                                end % k if
                                dMatK(m).epochStarts=[1];
                                dMatK(m).epochColors=dMatK(m).epochColors(1,:);
                                if isempty(dMatK(m).id), dMatK(m).id=[str2num(fid) idSU(j,:)]; else dMatK(m).id=[dMatK(m).id; str2num(fid) idSU(j,:)]; end
                            end % m
                        else % isempty(tmpD1)
%                             disp(['cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust))]);
                            for m=1:(nPr*nEH*nTspp*nTarg) % length(tmpD1) % cond
                                if epochId(k)==4 %tp
                                    tmpTs2=Tall(ceil(m/4)); % 480 480 480 480 540 ...
                                    if isempty(dMatK(m).data)
                                        dMatK(m).data=tmpD1{m}(tmpTs2:end);
                                    else
                                        minT=min([size(dMatK(m).data,2) length(tmpD1{m}(tmpTs2:end))]); % truncate with shorter tp
                                        dMatK(m).data=[dMatK(m).data(:,1:minT); tmpD1{m}(tmpTs2-1+[1:minT])];
                                    end
                                elseif epochId(k)==3 % prior
                                    if strfind(dMatK(m).condition,'Short')
                                        if idIncEarly
                                            tInit=251;
                                        else
                                            tInit=251+T{1}(1)-tBuffer; end
                                        tEnd=length(tmpD1{m});
                                    elseif strfind(dMatK(m).condition,'Long')
                                        if idIncEarly
                                            tInit=251;
                                        else
                                            tInit=251+T{2}(1)-tBuffer;end
                                        tEnd=length(tmpD1{m});
                                    end % if strcmp(dMatK(m).condition,'Short')
                                    dMatK(m).data=[dMatK(m).data; tmpD1{m}(tInit:tEnd)];
                                end % k if
                                dMatK(m).epochStarts=[1];
                                dMatK(m).epochColors=dMatK(m).epochColors(1,:);
                                if isempty(dMatK(m).id), dMatK(m).id=[str2num(fid) idSU(j,:)]; else dMatK(m).id=[dMatK(m).id; str2num(fid) idSU(j,:)]; end
                            end % m
                        end % isempty
                    end %                     if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
                end % j nUnit
%                 toc

                %% initial PCA
                
                tmpD=dMatK; % {40x1}.data [53 x 384]

                % resolving cell list difference as cells excluded in whole-session trajKS_XX
                % [538 neurons x(sessId,cellId,idSU)] as nTrialExc=5; % if # trials per condition<nTrialExc, cell is removed from analysis (exportDataHigh_conditionSpecific.m)
                % assuming cid is the same across conditions & epoches
    % 2018/7/20: remove use of dPSTH>use dOther
                cidInc=dOther.D(1).id(str2num(fid)==dOther.D(1).id(:,1),2); % 48 cells, w/o cells excluded in whole-session trajKS_XX
%                 clear dPSTH; % parfor
                nCellSess=size(tmpD(1).id,1);  %53
                idCSess=zeros(nCellSess,1);
                for iC=1:length(cidInc)
                    idCSess(cidInc(iC)==tmpD(1).id(:,2))=1;
                end
                idCSess=logical(idCSess); % 53-48=5 cells; 1147 1157 1162(SU) 1174 2123(SU) but SDF not that impressive
                for iD=1:length(tmpD)
                    tmpD(iD).data=tmpD(iD).data(idCSess,:);
                    tmpD(iD).id=tmpD(iD).id(idCSess,:);
                end

                % making it condition-specific; debug prior-specific
%                 if strcmp(epochNm{k},'pr')
                    % across ts
                    condId=[(0:(nTspp-1))*(nEH*nTarg)+... % ts
                    (2-D(iTrial).idHandEye-1)*(nTarg)+... % EH
                    (D(iTrial).theta==180)+1 ... % short
                    (nEH*nTspp*nTarg)+... % long prior
                    (0:(nTspp-1))*(nEH*nTarg)+... % ts
                    (2-D(iTrial).idHandEye-1)*(nTarg)+... % EH
                    (D(iTrial).theta==180)+1]; % e.g. 3     7    11    15    19    23    27    31    35    39
%                 else
%                     % finding this trial's condition
%                     condId=[(find(D(iTrial).ts==T{1})-1)*(nEH*nTarg)+... % ts
%                         (2-D(iTrial).idHandEye-1)*(nTarg)+... % EH
%                         (D(iTrial).theta==180)+1 ... % short
%                         (nEH*nTspp*nTarg)+... % long prior
%                         (find(D(iTrial).ts==T{2})-1)*(nEH*nTarg)+... % ts
%                         (2-D(iTrial).idHandEye-1)*(nTarg)+... % EH
%                         (D(iTrial).theta==180)+1];
%                 end
                %
                if idPCpoolSess

                    tmpD2=dOther.D; % clear dOther;  % parfor
                    sid=tmpD2(1).id(:,1)==str2num(fid); % [555 x 1]
                    for iOther=1:length(tmpD2) % across conditions
                        % assuming dOther has same size as tmpD
                        strTmp=tmpD2(iOther).condition(regexp(tmpD2(iOther).condition,'[^/]')); % to remove /
                        % truncate PSTH of other sessions
                        tmpSize=size(tmpD(iOther).data,2);
                        if epochId(k)==3 % prior % prior support
                            tmpD2(iOther).data=tmpD2(iOther).data(:,(end-tmpSize+1):end);
                            % putting leave-one-trial-out PSTH
                            if strcmp(D(iTrial).condition,strTmp)
                                tmpD2(iOther).data(sid,:)=tmpD(iOther).data;
%                                 disp('matched');
                            end
                        else % tp
                            if idUseTIC
                                % assuming min(tp)>tIC+tBuffer=200+80 >
                                % choose min among tIC+tBuffer vs #columns
                                try
                                    nColTmp=size(tmpD2(iOther).data,2);
                                    tmpD2(iOther).data=tmpD2(iOther).data(:,[tIC-(tBuffer-1):min(tIC+tBuffer,nColTmp)]); % [200-79 200+80]
                                catch
                                    disp('');
                                end
                                % putting leave-one-trial-out PSTH
                                if strcmp(D(iTrial).condition,strTmp)
                                    tmpD2(iOther).data(sid,:)=tmpD(iOther).data(:,[tIC-(tBuffer-1):min(tIC+tBuffer,nColTmp)]);
                                    %                                 disp('matched');
                                end
                            else
                                tmpSize=min([tmpSize;size(tmpD2(iOther).data,2)]);
                                tmpD2(iOther).data=tmpD2(iOther).data(:,1:tmpSize);
                                % putting leave-one-trial-out PSTH
                                if strcmp(D(iTrial).condition,strTmp)
                                    tmpD2(iOther).data(sid,1:tmpSize)=tmpD(iOther).data(:,1:tmpSize);
                                    %                                 disp('matched');
                                end
                            end % idUseTIC
                        end % epochId(k)==3 % prior % prior support
%                         if idCheckTraj
%                             figure; setFigPos(1,1);imagesc(tmpD2(iOther).data(sid,:));title(tmpD2(iOther).condition);
%                             figure; setFigPos(1,2);imagesc(tmpD(iOther).data);title(tmpD(iOther).condition);
%                             figure; setFigPos(2,2);imagesc(tmpD2(iOther).data(sid,:)-tmpD(iOther).data); colorbar;
%                             waitforbuttonpress;close all;
%                         end% should be almost same for condition of leave-out-trial; should be the same for other conditoins

                    end % for iOther=1:length(tmpD2)
                    tmpDold=tmpD; % this session
                    tmpD=tmpD2; % all session
%                     clear tmpD2;  % parfor
                end % idPCpoolSess

                % 2018/7/20 avgAttB4PCA
                if idAvgAttB4PCA & epochId(k)==3 % prior % prior support
                    for iPr=1:nPr
                        idTmpPr=(iPr-1)*nTspp+[1:nTspp]; % 1:5 for short, 6:10 for long
                        dTmp=struct2mat(tmpD(condId(idTmpPr)),'data'); % [neuron x time x 5ts]
                        dTmp=nanmean(dTmp,3); % [neuron x time]
                        
                        tmpD((iPr-1)*nTspp+1).data=dTmp;
                    end
                    d4PSTH2traj=tmpD([1 nTspp+1]);
                else
                    d4PSTH2traj=tmpD(condId);
                end % if idAvgAttB4PCA
                
                try
                    [tmpDtraj,proj_matrix,keep_neurons,eigenvalues,offsetCenter]=PSTH2traj(d4PSTH2traj,...
                        [],[],binSize,smthWidth,[],use_sqrt);
                    if isempty(tmpDtraj)
                        disp('data for pca may have nan');
                    end
                    
                    % for meas, % variance explained for its own PCA: 87.8208
                    % for prod, 
                    
%                     if idCheckTraj,plotPC(tmpDtraj,1);title('PC<PC(all sessions)');
%                     end
                    if idPCpoolSess
                        % removing other sessions
                        proj_matrix=proj_matrix(sid(keep_neurons),:); % [513 x 6]> [39 x 6]
                        offsetCenter=offsetCenter(sid (keep_neurons)); % [1 x 513]> [39 x 1]
                        keep_neurons=keep_neurons(sid,:); % [555 x 1]> [48 x 1]
                    end % if idPCpoolSess

                    D(iTrial).W=proj_matrix; % projection weight [#neurons x # PC]
                    D(iTrial).offset=offsetCenter(:); % offset for centering [# neurons x 1]
                    D(iTrial).keep_neurons=keep_neurons;
                    
                    proj_matrixTmp=cat(3,proj_matrix,[offsetCenter(:) nan(size(proj_matrix,1),size(proj_matrix,2)-1)]);
                    optimD=size(tmpDtraj(1).data,1);
                    
                    %% readout vector: xro
                    if idAllSessReadout & idPCpoolSess % estimating readout vector from PC of all sessions: idAllSessReadout
                         if k==1 % prior support
                             [tmpDtraj,proj_matrix,keep_neurons,eigenvalues,offsetCenter]=PSTH2traj(tmpD2(condId),...
                                 proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt);
                             D(iTrial).xro=estReadout(tmpDtraj([1:nTspp]+(tmpIdPr-1)*nTspp)); % now implemented with idUseAvgAtt=1
                         else % k==1 % prior support
                             tmpDtrajMat=struct2mat(tmpDtraj,'data'); % [#PC x time x 10prTs]
                             D(iTrial).xro=squeeze(mean(diff(tmpDtrajMat(:,iIC,[1:nTspp]+(tmpIdPr-1)*nTspp),1,3),3)); % mean difference vector
                             D(iTrial).xro=normalize(D(iTrial).xro(:),1);
                         end % k==1 % prior support
                    else % idAllSessReadout % estimating readout vector from PC of all sessions: idAllSessReadou
                        if k==1 % prior support
                            % tmpIdPr==1 % short; 2 for long
                            if idUseScoreReadout
                                if idPCpoolSess
                                    % N.B. just using tmpDtraj would be the same if idAvgAttB4PCA
                                    % use PC score of this sessions's PSTH
                                    [tmpDtraj,proj_matrix,keep_neurons,eigenvalues,offsetCenter]=PSTH2traj(tmpDold(condId),...
                                        proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt);
%                                     if idCheckTraj,plotPC(tmpDtraj,1);title('PC(this session with CV)<PC(all sessions)'); % #23(2007) with max W
%                                         [tmpDtraj2,P,K,~,O]=PSTH2traj(tmpDold(condId),... % if ~idPCpoolSess
%                                             [],[],binSize,smthWidth,[],use_sqrt);
%                                         plotPC(tmpDtraj2,1);title('PC(this session with CV)<PC(this session with CV)'); % #24(2022) with max W
%                                         figure; plot(D(iTrial).W(:,1),'b-');ha; plot(P(:,1),'r-'); xlabel('cell'); ylabel('loadings');
%                                     end
                                end
                                D(iTrial).xro=estReadout(tmpDtraj([1:nTspp]+(tmpIdPr-1)*nTspp)); % now implemented with idUseAvgAtt=1
                                
                            else % mean trajectory across trials; condition&prior specific
                                % N.B. no cross validation here
                                iTrialCond=([psthSess.D.idHandEye]==D(iTrial).idHandEye) &...
                                    ([psthSess.D.idShortTrial]==D(iTrial).idShortTrial) &...
                                    ([psthSess.D.theta]==D(iTrial).theta); % ~1600/4/2
                                [tmpDtraj2,tmpProj_matrix,tmpKeep_neurons,eigenvalues]=PSTH2traj(psthSess.D(iTrialCond),...
                                    proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt); %  proj_matrix,keep_neurons optimD
                                
                                tmpDtraj3=struct2mat(tmpDtraj2,'data'); % [dim x Maxtime x ntrial/cond/pr]
                                mD=squeeze(nanmean(tmpDtraj3,3)); % [dim Maxtime]
                                % get t(set)
                                nTime5=size(tmpDtraj(5+(tmpIdPr-1)*nTspp).data,2); % ts(5)-ts(1)
                                nTime1=size(tmpDtraj(1+(tmpIdPr-1)*nTspp).data,2);
                                nTime4=size(tmpDtraj(4+(tmpIdPr-1)*nTspp).data,2); %  ts(4)-ts(2)
                                nTime2=size(tmpDtraj(2+(tmpIdPr-1)*nTspp).data,2);
%                                 if idCheckTraj,figure;setFigPos(2,2); plot3(mD(1,:),mD(2,:),mD(3,:),'k-'); ha;grid on;plot3(mD(1,[nTime1 nTime2 nTime4 nTime5]),mD(2,[nTime1 nTime2 nTime4 nTime5]),mD(3,[nTime1 nTime2 nTime4 nTime5]),'kx');end;
                                
                                xro=mD(:,nTime5)-mD(:,nTime1);
                                xro2=mD(:,nTime4)-mD(:,nTime2);
                                D(iTrial).xro=normalize((xro+xro2)/2,1); % [# PC x 1] normalized
                            end
                            % % - not use all data by averaging across ts (e.g. earlier data from longest ts)
                            %                             nTime2=size(tmpDtraj(5+(tmpIdPr-1)*nTspp).data,2); % ts(5)-ts(1)
                            %                             nTime1=size(tmpDtraj(1+(tmpIdPr-1)*nTspp).data,2);
                            %                             xro=tmpDtraj(5+(tmpIdPr-1)*nTspp).data(:,nTime2)-tmpDtraj(1+(tmpIdPr-1)*nTspp).data(:,nTime1);
                            %
                            %                             nTime2=size(tmpDtraj(4+(tmpIdPr-1)*nTspp).data,2); %  ts(4)-ts(2)
                            %                             nTime1=size(tmpDtraj(2+(tmpIdPr-1)*nTspp).data,2);
                            %                             xro2=tmpDtraj(4+(tmpIdPr-1)*nTspp).data(:,nTime2)-tmpDtraj(2+(tmpIdPr-1)*nTspp).data(:,nTime1);
                            %
                            %                             D(iTrial).xro=normalize((xro+xro2)/2,1); % [# PC x 1] normalized
                            
                        else % @IC
                            if idUseScoreReadout
                                if idPCpoolSess
                                    % use PC score of this sessions's PSTH
                                    [tmpDtraj,proj_matrix,keep_neurons,eigenvalues,offsetCenter]=PSTH2traj(tmpDold(condId),...
                                        proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt);
%                                     if idCheckTraj,plotPC(tmpDtraj,1);title('this session with CV');end
                                end
                                tmpDtrajMat=struct2mat(tmpDtraj,'data'); % [#PC x time x 10prTs]
                                D(iTrial).xro=squeeze(mean(diff(tmpDtrajMat(:,iIC,[1:nTspp]+(tmpIdPr-1)*nTspp),1,3),3)); % mean difference vector
                                D(iTrial).xro=normalize(D(iTrial).xro(:),1);
                            else
                                % N.B. no cross validation here
                                iTrialCond=([psthSess.D.idHandEye]==D(iTrial).idHandEye) &...
                                    ([psthSess.D.idShortTrial]==D(iTrial).idShortTrial) &...
                                    ([psthSess.D.theta]==D(iTrial).theta); % ~1600/4/2
                                [tmpDtraj2,tmpProj_matrix,tmpKeep_neurons,eigenvalues]=PSTH2traj(psthSess.D(iTrialCond),...
                                    proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt); %  proj_matrix,keep_neurons optimD
                                
                                % averaging for each ts
                                tmpDtrajMat=[];
                                for iTs=1:nTspp
                                    tmpTmpDtrajMat=struct2mat(tmpDtraj2([psthSess.D(iTrialCond).ts]==T{tmpIdPr}(iTs)),'data'); %% [#PC x time x trials/ts/pr/cond]
                                    mtmpTmpDtrajMat=squeeze(nanmean(tmpTmpDtrajMat,3));% [#PC x Maxtime
                                    if iTs==1
                                        tmpDtrajMat=mtmpTmpDtrajMat;
                                    else
                                        tmpMin=min([size(tmpDtrajMat,2) size(mtmpTmpDtrajMat,2)]);
                                        tmpDtrajMat=cat(3,tmpDtrajMat(:,1:tmpMin,:),mtmpTmpDtrajMat(:,1:tmpMin));
                                    end
                                end
                                D(iTrial).xro=squeeze(mean(diff(tmpDtrajMat(:,iIC,[1:nTspp]),1,3),3)); % mean difference vector
                                D(iTrial).xro=normalize(D(iTrial).xro(:),1);
                            end
                            D(iTrial).Dcv=tmpDtraj; %%%%%
                            
                        end % if k==1
                        %                     clear tmpDold; %  % parfor
                        
                        %                     Dtmp=struct2mat(tmpDtraj,'data'); % [3dim 24time 10conditions]
                        %                     Dtmp=permute(Dtmp,[3 2 1]); % [10cond x 25time x 3]
                        %                     figure(999); ha
                        %                     DtmpMean=cat(1,nanmean(Dtmp(1:nTspp,:,:),1),nanmean(Dtmp(nTspp+[1:nTspp],:,:),1)); % [2pr x 25time x 3]
                        %                     plot3Dmultiple(DtmpMean,[tmpCmap{1,1}(1,:);tmpCmap{2,1}(1,:)],[]); grid on;
                        %                     plot3(DtmpMean(1,1,1),DtmpMean(1,1,2),DtmpMean(1,1,3),'^','color',tmpCmap{1,1}(1,:));
                        %                     plot3(DtmpMean(2,1,1),DtmpMean(2,1,2),DtmpMean(2,1,3),'^','color',tmpCmap{2,1}(1,:));
                        %                     if idIncEarly
                        %                         plot3(DtmpMean(1,24,1),DtmpMean(1,24,2),DtmpMean(1,24,3),'x','color',tmpCmap{1,1}(1,:)); % 480
                        %                         plot3(DtmpMean(2,40,1),DtmpMean(2,40,2),DtmpMean(2,40,3),'x','color',tmpCmap{2,1}(1,:));
                        %                     end
                    
                    end %                     if idAllSessReadout % estimating readout vector from PC of all sessions: idAllSessReadout

                catch
                    disp('error');
                end
                
                
                %% apply PCA to leaved-out trial      
                
                % loop through each unit
                nUnit=size(idSU,1); % both SU & MU
                for j=1:nUnit
                    if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
                    else
                        % consistent with plotSDFexport
                        NEV=kilosort2NEV(sp,idClust,idSU);
                        tmpId=[idSU(j,1) idSU(j,1) 1 length(movingBallTrials)]; % id(electrode,Unit,start trial,end trial)
                        
                        tmpD=exportDataHigh_conditionspecific(NEV,tmpId,tInfo(~idOut2,:),epochId(k),tInfo(iTrialId,1)); % 40 cells [1 x time]
                        if isempty(tmpD) % &
%                             disp(['cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust)) '     # trial/condition < ' num2str(nTrialExc)]);
                        else
%                             disp(['cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust))]);
                            
                            if epochId(k)==4 %tp
                                if isempty(D(iTrial).data)
                                    D(iTrial).data=tmpD(D(iTrial).ts:end);
                                else
                                    minT=min([size(D(iTrial).data,2) length(tmpD(D(iTrial).ts:end))]);
                                    D(iTrial).data=[D(iTrial).data(:,1:minT); tmpD(D(iTrial).ts-1+[1:minT])];
                                end
                                if isempty(D(iTrial).id), D(iTrial).id=[str2num(fid) idSU(j,:)]; else D(iTrial).id=[D(iTrial).id; str2num(fid) idSU(j,:)]; end
                            elseif epochId(k)==3 % prior
                                if strfind(D(iTrial).condition,'Short')
                                    if idIncEarly
                                        tInit=251;
                                    else
                                        tInit=251+T{1}(1)-tBuffer;end
                                    tEnd=length(tmpD);
                                    
                                    D(iTrial).data=[D(iTrial).data; tmpD(tInit:tEnd)];
                                elseif strfind(D(iTrial).condition,'Long')
                                    if idIncEarly
                                        tInit=251;
                                    else
                                        tInit=251+T{2}(1)-tBuffer;end
                                    tEnd=length(tmpD);
                                    
                                    D(iTrial).data=[D(iTrial).data; tmpD(tInit:tEnd)];
                                end % if strcmp(dMat{k}(m).condition,'Short')
                                if isempty(D(iTrial).id), D(iTrial).id=[str2num(fid) idSU(j,:)]; else D(iTrial).id=[D(iTrial).id; str2num(fid) idSU(j,:)]; end
                            end % k if

                        end % isempty
                    end %                     if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
                end % j nUnit
%                 toc % ~1s

                
                % resolving cell list difference as cells excluded in whole-session trajKS_XX
                D(iTrial).data=D(iTrial).data(idCSess,:);
                D(iTrial).id=D(iTrial).id(idCSess,:);
%                 disp('3');

                %  & ~idNotFullRecord
                if ~idNotFullRecord
                    iNotFullRecord2=iNotFullRecord(idCSess); % [48 x 1]
                    proj_matrixTmp=proj_matrixTmp(~iNotFullRecord2(keep_neurons),:,:); % [39 x 7 x 2] > [32 x 7 x 2]
                    keep_neurons(iNotFullRecord2)=false; % [48 x 1] 32 nnz
                    D(iTrial).keep_neurons=keep_neurons;
                    D(iTrial).offset=D(iTrial).offset(~iNotFullRecord2(keep_neurons)); % [39 x 1 > [32 x 1]
                    D(iTrial).W=D(iTrial).W(~iNotFullRecord2(keep_neurons),:);% [39 x 7 > [32 x 7]
                end
%                 disp('4');

%                 Dtraj=[];
                try
                [Dtraj,tmpProj_matrix,tmpKeep_neurons,eigenvalues]=PSTH2traj(D(iTrial),...
                    proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt); %  proj_matrix,keep_neurons optimD
%                 disp(iTrial);
                catch
                    disp('error');
%                     disp(['----- Trial' num2str(iTrial) '/Trial' num2str(nTrial) '-----']);
                end
                D(iTrial).data=Dtraj(1).data;
%                 disp('5');

%%                % check traj %%%%%
% fig3c1: rotating trajectory+readout vector in state space
% fig3c2: IC+readout vector in state space
                if idCheckTraj
                    Dtmp=struct2mat(tmpDtraj,'data'); % [3dim 24time 10conditions]

                    Dtmp=permute(Dtmp,[3 2 1]); % [10cond x 25time x 3]
                    figure; setFigPos(1,k); ha; grid on;
                    
                    if k==1 % avg-ts traj
                        DtmpMean=cat(1,nanmean(Dtmp(1:nTspp,:,:),1),nanmean(Dtmp(nTspp+[1:nTspp],:,:),1)); % [2pr x 25time x 3]
                        
                        if idSmoothAfterAvgAtt
                            binWidth=1; % 20 ms bin in PC
                            kernSD=2; % 40ms smoothing was used in PC
                            for iDtmpMean=1:size(DtmpMean,1)
                                tmpTmpTmp=squeeze(DtmpMean(iDtmpMean,:,:))'; % [3pc 24time]
                                isNanTmp=~isnan(tmpTmpTmp(1,:));
                                DtmpMean(iDtmpMean,isNanTmp,:)=smoother(tmpTmpTmp(:,isNanTmp),kernSD,binWidth)';
                            end
                        end
                                
                        plot3Dmultiple(DtmpMean,[tmpCmap{1,1}(1,:);tmpCmap{2,1}(1,:)],[],lw); % avg-ts traj
%                         plot3(DtmpMean(1,1,1),DtmpMean(1,1,2),DtmpMean(1,1,3),'^','color',tmpCmap{1,1}(1,:)); % start
%                         plot3(DtmpMean(2,1,1),DtmpMean(2,1,2),DtmpMean(2,1,3),'^','color',tmpCmap{2,1}(1,:));
                    end
                    
                    % states for each ts
                    cmapAll=cell2mat(tmpCmap(:,1)); % [10ts x 3]
                    for iTsTmp=1:size(Dtmp,1)
                        idPrTmp=floor((iTsTmp-1)./nTspp)+1; % 1111122222
                        if k==1 % rotation
                            tTmp=nnz(~isnan(Dtmp(iTsTmp,:,1)));
                            plot3(DtmpMean(idPrTmp,tTmp,1),DtmpMean(idPrTmp,tTmp,2),DtmpMean(idPrTmp,tTmp,3),...
                                'o','markerfacecolor',cmapAll(iTsTmp,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2); % states for each ts
                        else
                            tTmp=iIC;
%                             plotTraj_figure3('trajKS_tp_EL_H_bin20_smth40.mat',tmpDtraj);
                            if rem(iTsTmp,nTspp)~=0 % connecting lines
                                plot3(Dtmp(iTsTmp+[0:1],tTmp,1),Dtmp(iTsTmp+[0:1],tTmp,2),Dtmp(iTsTmp+[0:1],tTmp,3),...
                                    '-','color',tmpCmap{idPrTmp}(1,:),'linewidth',lw); % states for each ts
                            end
                            plot3(Dtmp(iTsTmp,tTmp,1),Dtmp(iTsTmp,tTmp,2),Dtmp(iTsTmp,tTmp,3),...
                                'o','markerfacecolor',cmapAll(iTsTmp,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2); % states for each ts
                        end
                    end
                    
                    % readout
                    if k==1
                        tmpX=nanmean(DtmpMean(tmpIdPr,[1 nnz(~isnan(DtmpMean(tmpIdPr,:,1)))],1)); % plot readout around center of rotation
                        tmpY=nanmean(DtmpMean(tmpIdPr,[1 nnz(~isnan(DtmpMean(tmpIdPr,:,2)))],2));
                        tmpZ=nanmean(DtmpMean(tmpIdPr,[1 nnz(~isnan(DtmpMean(tmpIdPr,:,3)))],3));
                    else
                        tmpX=nanmean(Dtmp([1:nTspp]+(tmpIdPr-1)*nTspp,iIC,1)); % plot readout around center of rotation
                        tmpY=nanmean(Dtmp([1:nTspp]+(tmpIdPr-1)*nTspp,iIC,2));
                        tmpZ=nanmean(Dtmp([1:nTspp]+(tmpIdPr-1)*nTspp,iIC,3));
                    end
                    plot3(tmpX+[-1;1]*D(iTrial).xro(1),...
                        tmpY+[-1;1]*D(iTrial).xro(2),...
                        tmpZ+[-1;1]*D(iTrial).xro(3),'--','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % readout
%                     if idIncEarly
%                         plot3(DtmpMean(1,24,1),DtmpMean(1,24,2),DtmpMean(1,24,3),'x','color',tmpCmap{1,1}(1,:)); % 480
%                         plot3(DtmpMean(2,40,1),DtmpMean(2,40,2),DtmpMean(2,40,3),'x','color',tmpCmap{2,1}(1,:));
%                     end
%                     plot3Dmultiple(Dtmp,cell2mat(tmpCmap(:,1)),[]); grid on; % individual ts traj
                    
%                     % this trial's state        
%                     cmap3=tmpCmap{tmpIdPr,1}(D(iTrial).ts==T{tmpIdPr},:);
% %                     plot3(Dtraj.data(1,:),Dtraj.data(2,:),Dtraj.data(3,:),'-','linewidth',2,'color',cmap3); % leaved trial's trajectory
%                     if k==1
%                         plot3(Dtraj.data(1,end),Dtraj.data(2,end),Dtraj.data(3,end),'o','linewidth',lw2,'markeredgecolor',cmap3,'markerfacecolor',cmap3,'markersize',msize); % start leaved trial's trajectory
%                     else
%                         plot3(Dtraj.data(1,iIC),Dtraj.data(2,iIC),Dtraj.data(3,iIC),'o','linewidth',lw2,'markeredgecolor',cmap3,'markerfacecolor',cmap3,'markersize',msize); % start leaved trial's trajectory
%                     end
%                     
%                     % trial state to readout
% %                     % xb=x(x'x)^(-1)x'y
% %                     x=D(iTrial).xro(1:3);
% %                     if k==1,y=Dtraj.data(1:3,end); else y=Dtraj.data(1:3,iIC); end
% %                     xb=x*((x'*x)^(-1))*x'*y;
%                     % (a-p)-((a-p)'*n)*n
%                     a=[tmpX;tmpY;tmpZ]; n=D(iTrial).xro(1:3);
%                     if k==1,p=Dtraj.data(1:3,end); else p=Dtraj.data(1:3,iIC); end
%                     xb=(a-p)-((a-p)'*n)*n;
%                     if k==1
%                         %                         plot3([Dtraj.data(1,end);xb(1)],...
%                         %                             [Dtraj.data(2,end);xb(2)],...
%                         %                             [Dtraj.data(3,end);xb(3)],':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
%                         plot3(Dtraj.data(1,end)+[0;1]*xb(1),...
%                             Dtraj.data(2,end)+[0;1]*xb(2),...
%                             Dtraj.data(3,end)+[0;1]*xb(3),':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
%                     else
%                         %                         plot3([Dtraj.data(1,iIC);xb(1)],...
%                         %                             [Dtraj.data(2,iIC);xb(2)],...
%                         %                             [Dtraj.data(3,iIC);xb(3)],':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
%                         plot3(Dtraj.data(1,iIC)+[0;1]*xb(1),...
%                             Dtraj.data(2,iIC)+[0;1]*xb(2),...
%                             Dtraj.data(3,iIC)+[0;1]*xb(3),':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
%                     end
                    
                    
                    axis tight;
                    if k==1
                        set(gca,'view',[-21 83],... % 40 74],... % 66 48],...
                            'xtick',-1:1,'ytick',-1:1,'ztick',-3:1:3);
                        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
                    else
%                         axis tight;
                        set(gca,'view',[-18.8000  -56.4000],...  % [147.6000   71.6000] % [-108 34],...
                            'xtick',-3:1:3,'ytick',-2:1:3,'ztick',-2:1:1);
                        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
                    end
%                     waitforbuttonpress; close;
                end
%                 clear Dtraj; % 4M %  % parfor
                if k==1 % Sset
                    D(iTrial).sSet=D(iTrial).xro(:)'*D(iTrial).data(:,end);
                else
                    D(iTrial).sIC=D(iTrial).xro(:)'*D(iTrial).data(:,iIC);
                end
%                  disp('6');
                % initializing
                proj_matrix=[];
                keep_neurons=[];
%                  disp('7');
%                 dMat=cell(nEpoch,1);
%                 for iEp=1:nEpoch
%                     dMat{iEp}=D0;
%                 end

title(['trial#' num2str(iTrial)]);drawnow;
                
%                 end % if D(iTrial).idHandEye==1 & D(iTrial).theta==180 %

            end % for iTrial=1:length(iNoOut2) % leave-one-trial-out cross validation
%             disp('8');
            %% save
%             if idPCpoolSess %%%%%
%                 if idIncEarly
%                     if idUseAvgAtt
%                         if idUseScoreReadout
%                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_ts_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat']); % psthDir
%                         else
%                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_ts_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAtt_newCenter_noUseScoreReadout' fnEnd '.mat']); % psthDir
%                         end
%                     else
%                         fnTmp=fullfile(singleTDir,['trajKS_' fid '_ts_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_newCenter' fnEnd '.mat']); % psthDir
%                     end
%                 else %%%%%
%                     if idUseAvgAtt %%%%%
% %                         fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat']); % psthDir
%                         fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAttB4PCA' fnEnd '.mat']); % psthDir
%                     else
%                         fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_newCenter' fnEnd '.mat']); % psthDir
%                     end
%                 end
%             else
%                 if idIncEarly
%                     if idUseAvgAtt
%                         if idUseScoreReadout
%                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_ts_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_sessPCA_CV_avgAtt_newCenter' fnEnd '.mat']); % psthDir
%                         else
%                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_ts_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_sessPCA_CV_avgAtt_newCenter_noUseScoreReadout' fnEnd '.mat']); % psthDir
%                         end
%                     else
%                         fnTmp=fullfile(singleTDir,['trajKS_' fid '_ts_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_sessPCA_CV_newCenter' fnEnd '.mat']); % psthDir
%                     end
%                 else
%                     if idUseAvgAtt
%                         fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_sessPCA_CV_avgAtt_newCenter' fnEnd '.mat']); % psthDir
%                     else
%                         fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_sessPCA_CV_newCenter' fnEnd '.mat']); % psthDir
%                     end
%                 end
%             end
%                         
% % for whole ts
% %             fnTmp=fullfile(singleTDir,['trajKS_' fid '_ts_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_sessPCA_CV.mat']); % psthDir
%             tp=[D.tp]';ts=[D.ts]';idPr=2-[D.idShortTrial]';idHandEye=[D.idHandEye]';theta=[D.theta]';        % 0hand1eye, 0right180left
%              save(fnTmp,...
%                 'D','binSize','smthWidth','use_sqrt','tp','ts','idPr','idHandEye','theta','param');
%             
%             % sSet sIC
%             if k==1
%                 [S{1:length(D)}]=deal(D.data);
%                 S=cellfun(@(x)x(:,end),S,'UniformOutput',0); % @set
%                 
%                 sSet=fillNanCell(S); % [nTrial x max#PC]      
%                 save(fnTmp,'sSet','-append');
%             else % k==2
%                 
%                 [S{1:length(D)}]=deal(D.data);
%                 S=cellfun(@(x)x(:,iIC),S,'UniformOutput',0); % 
%                 
%                 sIC=fillNanCell(S); % [nTrial x max#PC]          
%                 save(fnTmp,'sIC','tIC','-append');
%             end
            clear D;
        end % k epoch
%     end % for iSess=1:length(fn{iAnimal})
%     
% 
% 
% end % iAnimal



%% runPSTH2traj_ts_tp_conditionSpecific.m


return;

%% find out best trial with straight readout vector for production
initRSG2prior;
load('trajKS_161218_t_p_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly.mat'); % D(trial).Dcv ([1:10])
nPC=3;
binSize=20;tIC=200/20;

x=nan(nPC,nPr,nTspp);
distP=squeeze(nan(size(D)));

for iT=1:length(D)
    tmp=struct2mat(D(iT).Dcv,'data'); % [7dim x 46time x 10prTs]
    
    x(:,1,:)=tmp(1:nPC,tIC,1:nTspp); % short
    x(:,2,:)=tmp(1:nPC,tIC,nTspp+[1:nTspp]); % long
    
    tmpV=normalize(D(iT).xro(1:3),1); % readout vector % [3pc x 1]
    
    if D(iT).idShortTrial==1 % short
        dTmp=squeeze(x(:,1,:)); % [3pc x 5ts]
        pX=tmpV'*dTmp; % projection [1 x 5ts]
        distP(iT)=sum(sum((dTmp-tmpV*pX).^2,1),2); % distance to projection [3pc x 5ts]
        
    else % long
        dTmp=squeeze(x(:,2,:)); % [3pc x 5ts]
        pX=tmpV'*dTmp; % projection [1 x 5ts]
        distP(iT)=sum(sum((dTmp-tmpV*pX).^2,1),2); % distance to projection [3pc x 5ts]
    end % if D(iT).idShortTrial==1 % short
    
end % for iT=1:length(D)
figure; hist(distP,100);

idHandEye=[D.idHandEye];
idPr=2-[D.idShortTrial];
theta=[D.theta];
ts=[D.ts];
tp=[D.tp];
disp([distP distP<7 idPr idHandEye theta ts]); % short hand is best
