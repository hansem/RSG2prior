function estR2BLSrotation_fig3c(varargin)

% 2019/3/7
% final fig 3c: tp/Xu =f(ts), short in long included
% inset: RMSE difference (BLS-Linear) histogram 
% use original data set (before removing buffer after set for speed as this can
% generate artificial curvature-PCA edge truncation effect)
% check runPSTH2traj_prior_ts_conditionSpecific.m
% result hold idSmoothAfterAvgAtt=1 or 0 (but use 1 for better visual)
% G handRight

% 2019/1/29: speed with big buffer (idBigBuffer)
% smthWidth 20 %%%%%% % no buffer after Set before PCA

% 2019/1/24:
% test across-prior projection

% 2019/1/14
% - test sigmoid: simplest (dXu(extreme) vs dXu(middle)
% - analyze short in long: idSiL update
%   for now, use PCA(S,L,SiL): trajKS_prior_HR_G_avgAttB4PCA
% - measure speed: idSiL=0

% 2018/10/8
% check how uts' (60 deg from uts) can generate high R2
% check a and b from a*XUts+b
% optimal decoding: idDecode
% checking angle of uts/vtp across priors

% 2018/9/19
% (TBD) checking angle between uts across conditions
% plotyy
% use animal-, condition-specific wM > supplementary figures
% (TBD) decoding
%
% representative: estR2BLSrotation_fig3c('HL_H');
% estR2BLSrotation_fig3c('ER_H'); estR2BLSrotation_fig3c('EL_H'); estR2BLSrotation_fig3c('HR_H'); estR2BLSrotation_fig3c('HL_H'); estR2BLSrotation_fig3c('ER_G'); estR2BLSrotation_fig3c('EL_G'); estR2BLSrotation_fig3c('HR_G'); estR2BLSrotation_fig3c('HL_G'); 
% idSaveFig

% 2018/7/17
% new PCA: "supportOnly_avgAttB4PCA"
% idSiL: apply to [short in long]
% bootstrap with avgAttB4PCA
% (TBD) decode readout vector: idDecode
% (TBD) Srdjan's vector [R+200 ts(5)]: need to use PCA for ts

% 2018/6/23
% improve neural data fit to BLS
% issue: undershoot before prior mean (overshoot after) when plotted te-ts
% 1) (done) use all mid-points for readout vector, uts; update estReadout.m
% 2) (done) no avgAtt: weigted fitting for each ts (more weight to earlier) using fitlm
% 3) (done) find best readout vector (max corr with BLS: regression)
% -conclude: use 1),2) but plot te vs ts

% 2018/6/19
% fig3c2: going back to cartesian but dot product (cos)
%           subplot for short/long
% fig3d1: error bar from bootstrapping (not for worst eg for now)
%           tp-ts format (not tp) - idTpMTs
%           even fitting is done teBLS-ts & detrending - idDetrend
% smoothing after averaging across ts with attrition

% 2018/6/13
% polar plot for R2: left for short, right for long

% 2018/6/11
% fig3c1: aV'X+b vs ts (V:readout vector), a & b fitted to BLS)
% fig3c2: R2(BLS) vs angle(V,V')
% example data set: H EL

% estimating R^2 b/t rotation states projected to readout vector and BLS
% repeat with randomizing readout vector (nRand)
% across animals and conditions
% use trajKS_prior_(CondNm)_(H/G)_supportOnly.mat
% modified analSet_randProj.m
% depends on estReadout.m, randVect.m, BLS.m
% check estR2BLSIC.m too
% 2018/6/3

%% initial
initRSG2prior;
cd(psthDir);
load pplot.mat; % pplot.cmap for separate figure
nWin=50; % window for moving average
boxWin=ones(nWin,1)/nWin;

idSiL=1; % 0; % for speed, 1; % 0; % also analyzing SiL data
idBigBuffer=0; % 1; % for speed
idTpMTs=0; % 1; % plot te-ts vs ts
idDetrend=0; % 1;
idAvgAtt=1; % 0; % if 1, averaging across ts with attrition and use WLS with fitlm; if 0, do multiple linear regression
idDecode=0;% 1; % if 1, find best readout vector (max corr with BLS: regression)
idSmoothAfterAvgAtt=1; % 1; % smoothing after avg. across ts w/ attrition %%%%%
idBoot=0; % 1; % 0; % 1; % avgAttB4PCA
idBLSts=1; % use BLSts after marginzaling over tm
idSaveFig=0; % 1;

idPlotSep=true; % plot separately across animals and conditions
nRand=1000;

idDim=0; % 3
% H/G, ER/EL/HR/HL: 8 7 10 9 7 7 9 6 for priorIncSL, 8 7 9 8 10 9 11 7

% % selecting one data sets
% if isempty(varargin)
%     dsName='ER_H'; %'ER_G'; %'ER_H'; % 'HL_G'; % EL_H'; % eye left H
% else
%     dsName=varargin{1};
% end

% BLS
% wm=0.05;

% bootstrap
nBoot=1000;

% rotation specific?

% plot
cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
cmapRot=[rgb('DeepPink'); rgb('Aqua')];
msize=4; % 10; % 4; % 6; % 2;
msize2=1;
lw=.5; % 1.2; % 1; % 3; % 1.5;
lw2=2;

idDebug=0; % 1;

%% main
if idBigBuffer
    fnEnd='_bigBuffer';
else
    fnEnd='';
end
if idSiL
%     d=dir('trajKS_prior*_avgAttB4PCA_1half.mat'); %%%%%
    d=dir(['trajKS_prior*_avgAttB4PCA' fnEnd '.mat']); % _smth20 smthWidth 20 %%%%%% no buffer after Set
    tmpIdD=[];
    for iDS=1:length(d)
        if contains(d(iDS).name,'_supportOnly_')
            tmpIdD=[tmpIdD; false];
        else
            tmpIdD=[tmpIdD; true];
        end
    end
    d=d(logical(tmpIdD));
    nPr=nPr+1;
    T{end+1}=T{1};
else
    d=dir(['trajKS_prior*_supportOnly_avgAttB4PCA' fnEnd '.mat']); % _smth20 smthWidth 20 %%%%%% no buffer after Set
end
% d=dir('trajKS_prior*_supportOnly_avgAttB4PCA.mat');
nDS=length(d);

% summary variables

% if exist('R2BLSrotation.mat')
%     load('R2BLSrotation.mat'); % dAngle r2
% else
% output
dAngle=nan(nDS,nRand,nPr);
r2=nan(nDS,1+nRand,nPr); % last for actual data
a=nan(nDS,1+nRand,nPr); % last for actual data
b=nan(nDS,1+nRand,nPr); % last for actual data
te=cell(2,nPr); % actualData/worstExample
teBoot=cell(2,nPr); % te from bootstrap
pmDBoot=cell(2,nPr);

% 2019/1/14: test sigmoid: simplest (dXu(extreme) vs dXu(middle)
xu=cell(nDS,nPr); % [17/21 x 1]
nxu=cell(nPr,1); % # observation
RMSE=nan(nDS,nPr,2); % last for BLS vs linear
teMat=cell(nDS,nPr); % te from BLS
fxu=cell(nDS,nPr); % fitted neural projections [17/21 x 1]

% speed measurement
sp=cell(nPr,1);

aU=nan(nAnimal,nEH,nTarg);
idBootAngle=1;
nBootAngle=1000; % 000;
aUBoot=nan(nAnimal,nEH,nTarg,nBootAngle); % angle between priors
aTmp=nan(nAnimal,nEH,nTarg,nPr,nBootAngle); % angle between original vs bootstrapped data

for iDS=1:nDS
    %% checking angle of uts/vtp across priors % H/G, Hand/Eye, Right/Left
% %     disp(['===== ' d(iDS).name ' =====']);
%     % figure out which data set
%     if contains(d(iDS).name,'_ER_H_')
%         iAnimal=1; iEH=1; iTarg=1;iAnimalNm='H';
%         tmpD=load(['trajKS_prior_ER_H_supportOnly_avgAttB4PCA.mat']); % D
%     elseif contains(d(iDS).name,'_EL_H_')
%         iAnimal=1; iEH=1; iTarg=2;iAnimalNm='H';
%         tmpD=load(['trajKS_prior_EL_H_supportOnly_avgAttB4PCA.mat']); % D
%     elseif contains(d(iDS).name,'_HR_H_')
%         iAnimal=1; iEH=2; iTarg=1;iAnimalNm='H';
%         tmpD=load(['trajKS_prior_HR_H_supportOnly_avgAttB4PCA.mat']); % D
%     elseif contains(d(iDS).name,'_HL_H_')
%         iAnimal=1; iEH=2; iTarg=2;iAnimalNm='H';
%         tmpD=load(['trajKS_prior_HL_H_supportOnly_avgAttB4PCA.mat']); % D
%     elseif contains(d(iDS).name,'_ER_G_')
%         iAnimal=2; iEH=1; iTarg=1;iAnimalNm='G';
%         tmpD=load(['trajKS_prior_ER_G_supportOnly_avgAttB4PCA.mat']); % D
%     elseif contains(d(iDS).name,'_EL_G_')
%         iAnimal=2; iEH=1; iTarg=2;iAnimalNm='G';
%         tmpD=load(['trajKS_prior_EL_G_supportOnly_avgAttB4PCA.mat']); % D
%     elseif contains(d(iDS).name,'_HR_G_')
%         iAnimal=2; iEH=2; iTarg=1;iAnimalNm='G';
%         tmpD=load(['trajKS_prior_HR_G_supportOnly_avgAttB4PCA.mat']); % D
%     elseif contains(d(iDS).name,'_HL_G_')
%         iAnimal=2; iEH=2; iTarg=2;iAnimalNm='G';
%         tmpD=load(['trajKS_prior_HL_G_supportOnly_avgAttB4PCA.mat']); % D
%     end
%     
%     for iPr=1:nPr  % prior-specific
%         idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
%         uts{iAnimal,iEH,iTarg,iPr}=estReadout(tmpD.D(idTs)); % [dim x 1]
%         
%         % bootstrap
%         bootDir=['/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'...
%             'trajKS_prior_supportOnly_avgAttB4PCA/']; % trajKS_prior_supportOnly/']; % '/Volumes/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'; % macmini
%         
%         epSelect='prior'; % 'ts'; % 'tp';
%         for iBoot=1:nBootAngle
%             
%             tsNm=fullfile(bootDir,...
%                 ['trajKS_' epSelect '_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_supportOnly_avgAttB4PCA_' num2str(iBoot) '.mat']); % '_bin20_smth40.mat']; % ['trajKS_ts_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40.mat'];
%             tmpBoot=load(tsNm); % D
%             utsBoot{iAnimal,iEH,iTarg,iPr,iBoot}=estReadout(tmpBoot.D(idTs)); % [dim x 1]
%             
%             % check angle betwen uts and uts(bootstrap)
%             %             disp(tsNm);
%             %             disp(angleVectors(uts{iAnimal,iEH,iTarg,iPr},utsBoot{iAnimal,iEH,iTarg,iPr,iBoot}));
%             aTmp(iAnimal,iEH,iTarg,iPr,iBoot)=angleVectors(uts{iAnimal,iEH,iTarg,iPr},utsBoot{iAnimal,iEH,iTarg,iPr,iBoot});
%             %         % check trajectory
%             %         [diffV,sSet]=plotTraj_figure3(tsNm);
%             %         hSS=gcf;setFigPos(2,2);
%             %         if iPr==nPr
%             %             ha; arrow(sSet(1:3,1+nTspp),sSet(1:3,1+nTspp)+uts{iAnimal,iEH,iTarg,iPr}(1:3),'Type','line','Color','b','linewidth',2); % original u
%             %         else
%             %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+uts{iAnimal,iEH,iTarg,iPr}(1:3),'Type','line','Color','r','linewidth',2);  % original u
%             %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+utsBoot{iAnimal,iEH,iTarg,iPr,iBoot},'Type','line','Color',[0.5 0 0],'linewidth',1);
%             %         end
%             % %         waitforbuttonpress;
%             % %         close;
%             
%             if iPr==nPr
%                 aUBoot(iAnimal,iEH,iTarg,iBoot)=angleVectors(utsBoot{iAnimal,iEH,iTarg,1,iBoot}(:),utsBoot{iAnimal,iEH,iTarg,2,iBoot}(:));
%             end
%         end
%         
%     end % for iPr=1:nPr
%     
%     %     % check trajectory
%     %         [diffV,sSet]=plotTraj_figure3(d(iDS).name);
%     %         hSS=gcf;setFigPos(2,1);
%     % %         if iPr==nPr % long
%     %             ha; arrow(sSet(1:3,1+nTspp),sSet(1:3,1+nTspp)+uts{iAnimal,iEH,iTarg,nPr}(1:3),'Type','line','Color','b','linewidth',2);  % original u
%     %             ha; arrow(sSet(1:3,1+nTspp),sSet(1:3,1+nTspp)+ utsBoot{iAnimal,iEH,iTarg,nPr,iBoot}(1:3),'Type','line','Color',[0 0 0.5],'linewidth',1);
%     % %             else % short
%     %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+uts{iAnimal,iEH,iTarg,1}(1:3),'Type','line','Color','r','linewidth',2);  % original u
%     %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+utsBoot{iAnimal,iEH,iTarg,1,iBoot},'Type','line','Color',[0.5 0 0],'linewidth',1);
%     % %         end
%     
%     aU(iAnimal,iEH,iTarg)=angleVectors(uts{iAnimal,iEH,iTarg,1}(:),uts{iAnimal,iEH,iTarg,2}(:));
    
    %% main
%     if contains(d(iDS).name,['_' dsName '_']) % strcmp(d(iDS).name,['trajKS_prior_' dsName '_supportOnly.mat'])
        
        disp(['===== ' d(iDS).name ' =====']);
        load(d(iDS).name); % binSize smthWidth optimD use_sqrt proj_matrix keep_neurons D eigenvalues meanPSTH
        
        % animal-specific wM
        if strfind(d(iDS).name,'_H_')
            load('H_RSGprior_DMFC.mat','wFitCond');
            if contains(d(iDS).name,'_ER_') % strcmp(dsName(1:2),'ER')
                wm=wFitCond(1,1).w_m; % 0.0508
                dsName='ER_H'; disp(wm);
            elseif contains(d(iDS).name,'_EL_') % strcmp(dsName(1:2),'EL')
                wm=wFitCond(1,2).w_m; % 0.0508
                dsName='EL_H'; disp(wm);
            elseif contains(d(iDS).name,'_HR_') % strcmp(dsName(1:2),'HR')
                wm=wFitCond(2,1).w_m; % 0.0508
                dsName='HR_H'; disp(wm);
            elseif contains(d(iDS).name,'_HL_') % strcmp(dsName(1:2),'HL')
                wm=wFitCond(2,2).w_m; % 0.0508
                dsName='HL_H'; disp(wm);
            end
            %         wm=wFit.w_m; % 0.0508
            iAnimalNm='H';
        else
            load('G_RSGprior_DMFC.mat','wFitCond');
            if contains(d(iDS).name,'_ER_') % strcmp(dsName(1:2),'ER')
                wm=wFitCond(1,1).w_m; % 0.0508
                dsName='ER_G'; disp(wm);
            elseif contains(d(iDS).name,'_EL_')% strcmp(dsName(1:2),'EL')
                wm=wFitCond(1,2).w_m; % 0.0508
                dsName='EL_G'; disp(wm);
            elseif contains(d(iDS).name,'_HR_') % strcmp(dsName(1:2),'HR')
                wm=wFitCond(2,1).w_m; % 0.0508
                dsName='HR_G'; disp(wm);
            elseif contains(d(iDS).name,'_HL_') % strcmp(dsName(1:2),'HL')
                wm=wFitCond(2,2).w_m; % 0.0508
                dsName='HL_G'; disp(wm);
            end
            %         wm=wFit.w_m; % 0.0481
            iAnimalNm='G';
        end
        
        % estimate Vreadout
        if idSiL
%             tmpD=load(['trajKS_prior_' dsName '_avgAttB4PCA.mat']); % D
%             D((nPr-1)*nTs+[1:nTs])=tmpD.D((nPr-1)*nTs+[1:nTs]); % SiL
            %     dVro=load(d(iDS).name); % ,'Vro');
            %     if ~isfield(dVro,'Vro')
            for iPr=1:nPr  % prior-specific
                idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
                Vro{iPr}=estReadout(D(idTs)); % [dim x 1]
            end % for iPr=1:nPr
            %         save(d(iDS).name,'Vro','-append');
            %     else
            %         Vro=dVro.Vro; % {pr}[dim x 1]
            %     end
        else %%%%% used now
            for iPr=1:nPr  % prior-specific
                idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
                Vro{iPr}=estReadout(D(idTs)); % [dim x 1]
            end % for iPr=1:nPr
        end % if idSiL
        
        % debug
        if idDebug
            [diffV,sSet]=plotTraj_figure3(d(iDS).name);
            hSS=gcf;setFigPos(1,3);
            ha; arrow(sSet(1:3,1),sSet(1:3,1)+Vro{1}(1:3),'Type','line','Color','r','linewidth',2);
            ha; arrow(sSet(1:3,1+nTspp),sSet(1:3,1+nTspp)+Vro{2}(1:3),'Type','line','Color','b','linewidth',2);
            %         ha; plot3([0;Vro{1}(1)],[0;Vro{1}(2)],[0;Vro{1}(3)],'r-','linewidth',2);
            %         ha; plot3([0;Vro{2}(1)],[0;Vro{2}(2)],[0;Vro{2}(3)],'b-','linewidth',2);
        end
        
        % modify dimension if needed
        if idDim==0 % use optimD (PVAF>75%)
            dim=optimD;
        else % trim data
            %         try
            proj_matrix=proj_matrix(:,1:dim);
            %         catch
            %             disp('');
            %         end
            for iD=1:length(D),D(iD).data=D(iD).data(1:dim,:);end
        end
        
        %% speed
        tmpSp{1}=cell(1,1); % short
        for iD=1:nTspp % length(D)
            tmpSp{1}=[tmpSp{1}; sqrt(sum(diff(squeeze(D(iD).data),1,2).^2,1))];
        end
        tmpSp{2}=cell(1,1); % long
        for iD=(nTspp+1):(2*nTspp) % length(D)
            tmpSp{2}=[tmpSp{2}; sqrt(sum(diff(squeeze(D(iD).data),1,2).^2,1))];
        end
        if idSiL
            tmpSp{3}=cell(1,1); % short in long
            for iD=(2*nTspp+1):(3*nTspp) % length(D)
                tmpSp{3}=[tmpSp{3}; sqrt(sum(diff(squeeze(D(iD).data),1,2).^2,1))];
            end
        end
        for iPr=1:nPr
            sp{iPr}=[sp{iPr}; fillNanCell(tmpSp{iPr})];
        end
        
        %% projection: after averaging across ts
        for iPr=1:nPr  % prior-specific
            idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
            dTmp=struct2mat(D(idTs),'data'); % [PC x time x ts]
            if idAvgAtt %%%%% actually used
                mD{iPr}=nanmean(dTmp,3); % [PC x time x ts] > [PC x time]
                
%                 % speed measurement before additional smoothing
%                 sp{iPr}=[sp{iPr}; sqrt(sum(diff(squeeze(mD{iPr}),1,2).^2,1))];
                
                if ~(idSiL & iPr==nPr) % no need for this for SiL
                    % smoothing after averaging across ts with attrition
                    if idSmoothAfterAvgAtt
                        binWidth=1; % 20 ms bin in PC
                        kernSD=2; % 40ms smoothing was used in PC
                        mD{iPr} = smoother(mD{iPr}, kernSD, binWidth);
                    end
                    
                    % removing points out of prior support (shortest ts has two data points: 4 20-ms bin before & after > remove 3 for both)
                    nTShortest=size(D(idTs(1)).data,2); % 2 for both short/long
                    mD{iPr}=mD{iPr}(:,nTShortest:end);
                    
%                     % speed measurement
%                     sp{iPr}=[sp{iPr}; sqrt(sum(diff(squeeze(mD{iPr}),1,2).^2,1))];
                
                    %             pmD{iPr}=pmD{iPr}(nTShortest:end);
                    dTmp=dTmp(:,nTShortest:end,:);
                else % still removing points out of prior support for SiL
                    nTShortest=size(D(1).data,2); % 2 for both short/long
                    mD{iPr}=mD{iPr}(:,nTShortest:end);
                    dTmp=dTmp(:,nTShortest:end,:);
                    
%                     % speed measurement
%                     sp{iPr}=[sp{iPr}; sqrt(sum(diff(squeeze(mD{iPr}),1,2).^2,1))];
                end % if idSiL & iPr~=nPr
                
                % do actual projection
                pmD{iPr}=Vro{iPr}(:)'*mD{iPr}; % [1 x time]
                
                %%%%%%
                % 2023/6/25: override if Xu exist (from 2019/3/9)
                %%%%%
                if exist('Xu')
                    pmD=Xu;
                end
                
                % 2019/1/14: test sigmoid: simplest (dXu(extreme) vs dXu(middle)
                xu{iDS,iPr}=pmD{iPr}(:); % [17/21 x 1]
                
                % R2 with BLS
                tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:); % 17/21 for short/long
                if idBLSts
                    teBLS{iPr}=BLSts(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                else
                    teBLS{iPr}=BLS(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                end
                teMat{iDS,iPr}=teBLS{iPr};
                if idDetrend
                    stat{iPr}=regstats(teBLS{iPr}-tm(:),detrend(pmD{iPr}),'linear'); % rsquare
                    te{1,iPr}=stat{iPr}.yhat(:)+tm(:);
                else %%%%% actually used
                    nObs{iPr}=sum(~isnan(squeeze(dTmp(1,:,:))),2); % [time x 1]
                    nxu{iPr}=nObs{iPr}(:);
                    % lsocv(pmD{iPr}(:),teBLS{iPr}(:),nObs) will do the job too
                    LM=fitlm(pmD{iPr}(:),teBLS{iPr}(:),'linear',...
                        'Weights',nObs{iPr});
                    LMM{iPr}=LM;
                    te{1,iPr}=LM.predict;
                    fxu{iDS,iPr}=te{1,iPr};
                    %                 stat{iPr}=regstats(teBLS{iPr},pmD{iPr},'linear'); % rsquare
                    %                 te{1,iPr}=stat{iPr}.yhat(:);
                    
                    % 2019/1/14: test sigmoid: simplest (dXu(extreme) vs dXu(middle)
                    % switch xy
                    LMBLS=fitlm(teBLS{iPr}(:),pmD{iPr}(:),'linear',...
                        'Weights',nObs{iPr});
                    RMSE(iDS,iPr,1)=LMBLS.RMSE;
                    LMlinear=fitlm(tm,pmD{iPr}(:),'linear',...
                        'Weights',nObs{iPr});
                    RMSE(iDS,iPr,2)=LMlinear.RMSE;
                    
                    %%
                    if idDecode % decode best readout vector with max. corr. with BLS
                        % do L2/1-regularization lasso ridge
                        [B,sLasso]=lasso(squeeze(mD{iPr})',teBLS{iPr}(:),...
                            'Alpha',0.1,... % if 1, lasso; near 0, ridge
                            'Weights',nObs{iPr}); %                       'Lambda',
                        lassoPlot(B,'PredictorNames',{'PC1','PC2','PC3'});
                        lassoPlot(B,sLasso,'PlotType','Lambda','XScale','log');
                        %                     lassoPlot(B,sLasso,'PlotType','CV');
                        
                        if iPr==1, iBest=15; B2=B(:,iBest);
                        else
                            iBest=13; B2=B(:,iBest);
                        end
                        X=squeeze(mD{iPr})';
                        teFit=[ones(size(X,1),1) X] * [sLasso.Intercept(iBest); B2];
                        figure; plot(teFit,teBLS{iPr}(:),'o'); xlabel('t_e from regularized regression');  ylabel('t_e'); plotIdentity(gca);
                        disp(['angle: ' num2str(angleVectors(B2,Vro{iPr}),3)]);
                        decVro=normalize(B2,1);
                        
                        % linear regression
                        %                     decLM=fitlm(squeeze(mD{iPr})',teBLS{iPr}(:),'linear',... % [time x PC]
                        %                         'Weights',nObs{iPr});
                        %                     betaTmp=table2array(decLM.Coefficients(:,1));
                        %                     disp(['beta: ' num2str(betaTmp(:)',3)]);
                        %                     disp(['angle: ' num2str(angleVectors(betaTmp(2:end),Vro{iPr}),3)]);
                        %                     decVro=normalize(betaTmp(2:4),1);
                        % plot
                        [diffV,sSet]=plotTraj_figure3(d(iDS).name);
                        hSS=gcf;setFigPos(1,3);
                        ha; arrow(sSet(1:3,1+nTspp*(iPr-1)),sSet(1:3,1+nTspp*(iPr-1))+Vro{iPr}(1:3),'Type','line','Color',pplot.cmap{2*(iPr-1)+1},'linewidth',2); % uts
                        ha; arrow(sSet(1:3,1+nTspp*(iPr-1)),sSet(1:3,1+nTspp*(iPr-1))+decVro(1:3),'Type','line','Color',pplot.cmap{4-(iPr-1)*2},'linewidth',2);
                    end % idDecode
                end
                if exist('stat','var')
                    r2(iDS,end,iPr)=stat{iPr}.rsquare;
                    disp(['R^2: ' num2str(stat{iPr}.rsquare,3)]);
                else % use fitlm
                    r2(iDS,end,iPr)=LM.Rsquared.Ordinary;
                    disp(['R^2: ' num2str(LM.Rsquared.Ordinary,3)]);
                    
                    b(iDS,end,iPr)=table2array(LM.Coefficients(1,1)); disp(['intercept: ' num2str(b(iDS,end,iPr),3)]);
                    a(iDS,end,iPr)=table2array(LM.Coefficients(2,1)); disp(['slope: ' num2str(a(iDS,end,iPr),3)]);
                end
                
                %         % debug
                if idDebug
                    if idTpMTs
                        %                     figure; setFigPos(iPr,1);plotyy(tm,teBLS{iPr}-tm,tm,pmD{iPr}-tm);
                        if exist('stat','var')
                            figure; setFigPos(iPr,2);plot(tm,teBLS{iPr}-tm,'-',tm,stat{iPr}.yhat-tm,'x','color',tmpCmap{iPr,1}(1,:));
                        else
                            figure; setFigPos(iPr,2);plot(tm,teBLS{iPr}-tm,'-',tm,LM.predict-tm,'x','color',tmpCmap{iPr,1}(1,:));
                        end
                        axis tight; plotHorizon(gca,0,[]); %plotIdentity(gca);
                        figure; plotResiduals(LM,'fitted');setFigPos(iPr,3); % probability
                        %                     plotDiagnostics(LM,'cookd'); setFigPos(iPr,3);
                    else
                        figure; setFigPos(iPr,1);plotyy(tm,teBLS{iPr},tm,pmD{iPr});
                        figure; setFigPos(iPr,2);plot(tm,teBLS{iPr},'-',tm,te{1,iPr},'x','color',tmpCmap{iPr,1}(1,:));
                        axis tight; plotIdentity(gca);
                    end
                end % idDebug
                
            else % w/o avgAtt > ts-specific; not used
                %              %%
                %             % goal: formal treat of different time length across ts
                %             pmD{iPr}=[]; % [ts x time]
                %
                %             for iTs=1:nTspp
                %                 pmD{iPr}=[pmD{iPr};...
                %                     Vro{iPr}(:)'*squeeze(dTmp(:,:,iTs))]; % [1 x PC] x [PC x time] > [1 x time]
                %             end
                %             pmD{iPr}=pmD{iPr}'; % [time x ts]
                %             pmD{iPr}(isnan(pmD{iPr}))=0;
                %
                %             % removing points out of prior support (shortest ts has two data points: 4 20-ms bin before & after > remove 3 for both)
                %             nTShortest=size(D(idTs(1)).data,2); % 2 for both short/long
                %             pmD{iPr}=pmD{iPr}(nTShortest:end,:);
                %
                %             % R2 with BLS
                %             tm=linspace(T{iPr}(1),T{iPr}(end),size(pmD{iPr},1)); tm=tm(:); % 17/21 for short/long
                %             teBLS{iPr}=BLS(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                %                 stat{iPr}=regstats(teBLS{iPr},pmD{iPr},'linear'); % rsquare, yhat, beta
                %                 te{1,iPr}=stat{iPr}.yhat(:);
                %                 disp(['beta: ' num2str(stat{iPr}.beta(:)',3)]);
                %             r2(iDS,end,iPr)=stat{iPr}.rsquare;
                %             disp(['R^2: ' num2str(stat{iPr}.rsquare,3)]);
                %
                %             %         % debug
                %             if idDebug
                %                 if idTpMTs
                %                     figure; setFigPos(iPr,1);plotyy(tm,teBLS{iPr}-tm,tm,pmD{iPr}-tm);
                %                     figure; setFigPos(iPr,2);plot(tm,teBLS{iPr}-tm,'-',tm,stat{iPr}.yhat-tm,'x','color',tmpCmap{iPr,1}(1,:));
                %                     axis tight; plotHorizon(gca,0,[]); %plotIdentity(gca);
                %                 else
                %                     figure; setFigPos(iPr,1);plotyy(tm,teBLS{iPr},tm,pmD{iPr});
                %                     figure; setFigPos(iPr,2);plot(tm,teBLS{iPr},'-',tm,stat{iPr}.yhat,'x','color',tmpCmap{iPr,1}(1,:));
                %                     axis tight; plotIdentity(gca);
                %                 end
                %             end % idDebug
            end % idAvgAtt
        end % iPr
        
        %% bootstrap
        if idBoot
            bootDir=['/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'...
                'trajKS_prior_supportOnly_avgAttB4PCA/']; % trajKS_prior_supportOnly/']; % '/Volumes/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'; % macmini
            
            epSelect='prior'; % 'ts'; % 'tp';
            for iBoot=1:nBoot
                tsNm=fullfile(bootDir,...
                    [d(iDS).name(1:(end-4)) '_' num2str(iBoot) '.mat']);
%                     ['trajKS_' epSelect '_' dsName(1:2) '_' iAnimalNm '_supportOnly_avgAttB4PCA_' num2str(iBoot) '.mat']); % '_bin20_smth40.mat']; % ['trajKS_ts_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40.mat'];
                [teBoot,pmDBoot]=estTeRv(teBoot,tsNm,nPr,nTs,idDim,T,idDetrend,1,idBLSts,pmDBoot); % (teBoot,tsNm,nPr,T); last idRand 1 for actual data
                
                % short in long
                if idSiL
                    %                 bootDir=['/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'...
                    %                     'trajKS_prior_avgAttB4PCA/'];
                    %                 tsNm=fullfile(bootDir,...
                    %                     ['trajKS_' epSelect '_' dsName(1:2) '_' iAnimalNm '_avgAttB4PCA_' num2str(iBoot) '.mat']); % '_bin20_smth40.mat']; % ['trajKS_ts_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40.mat'];
                    %                 teBoot=estTeRv(teBoot,tsNm,nPr,nTs,idDim,T,idDetrend,1,idBLSts,idSiL); % (teBoot,tsNm,nPr,T); last idRand 1 for actual data,
                    
                end % if idSiL
            end % for iBoot=1:nBoot
        end % if idBoot
        
        %%    % random readout vectors: uniformly distributed [0 90] deg
%         if idSiL, nPrTmp=nPr-1; else nPrTmp=nPr; end % not for SiL
%         for iPr=1:nPrTmp  % prior-specific
%             aRand=rand(nRand,1)*90; % (0 90)
%             dAngle(iDS,:,iPr)=aRand;
%             
%             tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr}));  tm=tm(:);% 17/21 for short/long
%             for iRand=1:nRand
%                 rVro=randVect(Vro{iPr},aRand(iRand),1);
%                 
%                 % projection
%                 pmDr{iPr}=rVro(:)'*mD{iPr}; % [1 x time]
%                 %             pmDr{iPr}=pmDr{iPr}(nTShortest:end); % removing points out of prior support
%                 if idDetrend
%                     statTmp=regstats(teBLS{iPr}-tm(:),detrend(pmDr{iPr}),'linear'); % rsquare
%                     r2(iDS,iRand,iPr)=statTmp.rsquare;
%                 else
%                     %                 statTmp=regstats(teBLS{iPr},pmDr{iPr},'linear'); % rsquare
%                     %                 r2(iDS,iRand,iPr)=statTmp.rsquare;
%                     LM2=fitlm(pmDr{iPr}(:),teBLS{iPr}(:),'linear',...
%                         'Weights',nObs{iPr});
%                     LMM2{iPr}=LM2;
%                     r2(iDS,iRand,iPr)=LM2.Rsquared.Ordinary;
%                     
%                     b(iDS,iRand,iPr)=table2array(LM2.Coefficients(1,1));
%                     a(iDS,iRand,iPr)=table2array(LM2.Coefficients(2,1));
%                 end
%                 %             r2(iDS,iRand,iPr)=statTmp.rsquare;
%                 %             disp(['dAngle: ' num2str(aRand(iRand),3) ', R^2: ' num2str(r2(iDS,iRand,iPr),3)]);
%                 
%                 % check how uts' (60 deg from uts) can generate high R2
%                 %             if r2(iDS,iRand,iPr)<0.9 & aRand(iRand)>75
%                 %                 [diffV,sSet]=plotTraj_figure3(d(iDS).name);setFigPos(2,3);
%                 %                 hSS=gcf;
%                 %                 ha; arrow(sSet(1:3,1+nTspp*(iPr-1)),sSet(1:3,1+nTspp*(iPr-1))+rVro(1:3),'Type','line','Color',cmapPr(iPr,:),'linewidth',2);
%                 %
%                 %                 tm=linspace(T{iPr}(1),T{iPr}(end),length(pmDr{iPr})); % 17/21 for short/long
%                 %                 figure; setFigPos(iPr,5); drawnow;plotyy(tm,teBLS{iPr},tm,pmDr{iPr}); title(aRand(iRand));
%                 %                 figure; setFigPos(iPr,6);drawnow;plot(tm,teBLS{iPr},'-',tm(:),LM2.predict,'x','color',tmpCmap{iPr,1}(1,:));  title(r2(iDS,iRand,iPr));
%                 %                 axis tight; plotIdentity(gca);
%                 %
%                 %                  disp(['intercept: ' num2str(b(iDS,iRand,iPr),3)]);
%                 %                   disp(['slope: ' num2str(a(iDS,iRand,iPr),3)]);
%                 %                 waitforbuttonpress; close;
%                 %             end
%                 
%                 % pick up worst example
%                 if r2(iDS,iRand,iPr)<0.4 & aRand(iRand)>70 & isempty(te{2,iPr})  & aRand(iRand)<85
%                     if idDetrend
%                         te{2,iPr}=LM2.predict+tm(:); % statTmp.yhat(:)+tm(:);
%                         iRandEx{iPr}=iRand;
%                     else
%                         te{2,iPr}=LM2.predict; % statTmp.yhat(:);
%                         iRandEx{iPr}=iRand;
%                     end
%                     if idBoot
%                         % bootstrap
%                         bootDir=['/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'...
%                             'trajKS_prior_supportOnly_avgAttB4PCA/'];
%                         for iBoot=1:nBoot
%                             tsNm=fullfile(bootDir,...
%                                 ['trajKS_' epSelect '_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_supportOnly_avgAttB4PCA_' num2str(iBoot) '.mat']); % '_bin20_smth40.mat']; % ['trajKS_ts_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40.mat'];
%                             [teBoot,pmDBoot]=estTeRv(teBoot,tsNm,nPrTmp,nTs,idDim,T,idDetrend,2,idBLSts,pmDBoot); % (teBoot,tsNm,nPr,T); last 1 for actual data
%                         end
%                     end
%                     disp(['dAngle: ' num2str(aRand(iRand),3) ', R^2: ' num2str( r2(iDS,iRand,iPr),3)]);
%                     
%                     % debug
%                     if idDebug
%                         [diffV,sSet]=plotTraj_figure3(d(iDS).name);setFigPos(2,3);
%                         hSS=gcf;
%                         ha; arrow(sSet(1:3,1+nTspp*(iPr-1)),sSet(1:3,1+nTspp*(iPr-1))+rVro(1:3),'Type','line','Color',cmapPr(iPr,:),'linewidth',2);
%                         
%                         tm=linspace(T{iPr}(1),T{iPr}(end),length(pmDr{iPr})); % 17/21 for short/long
%                         figure; setFigPos(iPr,5); drawnow;plotyy(tm,teBLS{iPr},tm,pmDr{iPr}); title(aRand(iRand));
%                         figure; setFigPos(iPr,6);drawnow;plot(tm,teBLS{iPr},'-',tm(:),LM2.predict,'x','color',tmpCmap{iPr,1}(1,:));  title(r2(iDS,iRand,iPr));
%                         axis tight; plotIdentity(gca);
%                     end
%                     
%                 end % worst case
%             end % for iRand=1:nRand
%         end % for iPr=1:nPr  % prior-specific
        
        
%     end % if strcmp(d(iDS).name,'trajKS_prior_EL_H_supportOnly.mat')
end % for iDS=1:nDS
% save('R2BLSrotation.mat','dAngle','r2');
% end % if exist('R2BLSrotation.mat')
%% speed measurement
% figure; ha; setFigPos(2,1);
% for iPr=1:nPr
%      tm=linspace(T{iPr}(1),T{iPr}(end),size(sp{iPr},2));  tm=tm(:);% 17/21 for short/long
%      for iDS=1:(nDS*nTspp)
%         plot(tm(:)',sp{iPr}(iDS,:),'-','color',cmapPr(iPr,:),'linewidth',0.5);
%      end
%      shadedErrorBar(tm(:)',nanmean(sp{iPr},1),sem(sp{iPr},1,'omitnan'),{'o-','color',cmapPr(iPr,:),'linewidth',2},1); drawnow; % 2 ,'markersize',msize2
% end
% xlabel('time from ready');
% ylabel('speed');

%% 2019/1/14: test sigmoid: simplest (dXu(extreme) vs dXu(middle)
% model-based approach: RMSE(BLS) vs RMSE(linear) both have two free params
figure; setFigPos(1,1); ha;
for iPr=1:nPr
    plot(RMSE(:,iPr,1),RMSE(:,iPr,2),'o','color',cmapPr(iPr,:));
    p=signrank(RMSE(:,iPr,1),RMSE(:,iPr,2));
    disp(['p(signrank): ' num2str(p)]);
end
% plot(RMSE(:,:,1),RMSE(:,:,2),'ko'); 
axisEqual;
plotIdentity(gca);
xlabel('RMSE(BLS)');
ylabel('RMSE(linear)');
% stat: friedman for nonparametric of 2 way anova2
% reps=nDS;
% [pFried,tFried,statFried]=friedman(reshape(RMSE,nDS*nPr,2),reps); % p for p value of H0(different BLS/linear), t{4,6} for interaction p value
RMSE4tbl=reshape(RMSE,nDS,nPr*2); %[8 x 6]
Tbl=array2table(RMSE4tbl,'VariableNames',{'y1','y2','y3','y4','y5','y6'}); % {'shortBLS';'longBLS';'SiLBLS';'shortlinear';'longlinear';'SiLlinear'});
% Tbl=table(RMSE(:),repmat(1:nDS,1,nPr*2)',repmat([ones(nDS,1);2*ones(nDS,1);3*ones(nDS,1)],2,1),[ones(nDS*nPr,1);ones(nDS*nPr,1)*2],...
%     'VariableNames',{'RMSE','dataId','prior','model'});
rm=fitrm(Tbl,'y1-y6~1','WithinDesign',array2table([1 1;2 1;3 1;1 2;2 2;3 2],'VariableNames',{'prior','model'}));
rtbl=ranova(rm,'WithinModel','prior:model') % prior+model prior*model
% 2019/3/7, inset: RMSE difference (BLS-Linear) histogram 
figure; setFigPos(1,1); ha;
nBin=8;
barWidth=0.3;
markersize=10;
% getting bin locations
RMSE2=squeeze(RMSE(:,:,1)-RMSE(:,:,2)); % BLS-linear
RMSE3=squeeze(mean(RMSE2,1)); % [3pr x 1]
x=linspace(min(RMSE2(:)),max(RMSE2(:)),nBin);
% plot
for iPr=1:nPr
    h=histogram(RMSE2(:,iPr),x); % nBin); % ,'DisplayStyle',DisplayStyle);
     Xval    = h.BinEdges + h.BinWidth*0.5 +h.BinWidth*(iPr-2)*0.3; % (i-(nPr+1)/2)/2; % binWidth*0.25 for short, 0.75 for long > binWidth*0.2 for short, 0.5 for long, 0.8 for SiL
     Xval    = Xval(1:end-1);
     Yval    = h.Values;
     delete(h);
     hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',cmapPr(iPr,:),'EdgeColor','none');
end
axis tight; 
maxYlim=max(ylim);
% plot mean
for iPr=1:nPr
    plot(mean(RMSE3(iPr)),1.1*maxYlim,'color',cmapPr(iPr,:),'markerfacecolor',cmapPr(iPr,:),'marker','v','markersize',markersize);
end % for i=1:nPr
xlabel('RMSE(Bayesian)-RMSE(linear)'); ylabel('# data sets'); set(gca,'tickdir','out','ticklength',[0.02 0.02],'xtick',-0.04:0.04:0.04,'ytick',1:5); % 'xtick','','ytick',''
xlim([-0.05 0.05]);
plotVertical(gca,0,[]);

% model-free
% dTp=nan(nDS,nPr,2); % % last for middle vs extreme
% slope=nan(nDS,nPr);
% for iDS=1:nDS
%     for iPr=1:nPr
%         xu2=xu{iDS,iPr}(:); % [17/21 x 1]
%         tm=linspace(T{iPr}(1),T{iPr}(end),length(xu2));  tm=tm(:);% 17/21 for short/long
%         
%         tmpLM=fitlm(tm,xu2);
%         slope(iDS,iPr)=table2array(tmpLM.Coefficients(2,1));
%         
%        % model-free test b/t extreme vs middle (only for points having corresponding ts)
%        nBin=unique(diff(T{iPr}))/binSize; % 4 for short, 5 for long
%        muTp=xu2(1:nBin:end);
%        dTp(iDS,iPr,1)=mean(diff(muTp(2:4)));
%        dTp(iDS,iPr,2)=mean([diff(muTp(1:2)) diff(muTp(4:5))]);
%         
%     end % for iPr=1:nPr
% end % for iDS=1:nDS
% 
% figure; ha; setFigPos(1,2);
%     tmpX=dTp(:,:,1); tmpX=tmpX(:);
%     tmpY=dTp(:,:,2); tmpY=tmpY(:);
% %     plot(tmpX,tmpY,'ko');
% for iPr=1:nPr
%     plot(dTp(:,iPr,1),dTp(:,iPr,2),'o','color',cmapPr(iPr,:));
% end
% axis tight; plotIdentity(gca);
% ylabel('difference in mean t_p for extreme t_s');
% xlabel('difference in mean t_p for middle t_s');
% 
% % histogram
% figure; ha; setFigPos(1,3);
% % hHist=histogram(tmpY-tmpX,25);
% % set(hHist,'Facecolor','w');
% for iPr=1:nPr
%     hHist=histogram(-dTp(:,iPr,1)+dTp(:,iPr,2),10);
%     set(hHist,'Edgecolor',cmapPr(iPr,:),'Facecolor','w');
% end
% plotVertical(gca,0,[]);
% ylabel('# data sets');
% xlabel('degree of being sigmoidal');
% 
% % check slope
% if idSiL
%     figure; ha; setFigPos(1,4);
%     plot(slope(:,1),slope(:,3),'ko');
%     plotIdentity(gca);
%     ylabel('slope(short in long)');
%     xlabel('slope(short)');
% end



%% plot
% dAngle(2an x 2EH x 2targ, nRand, 2pr)
% r2(2an x 2EH x 2targ, nRand+actual data, 2pr)

if ~idPlotSep
%     % use pplot.cmap: rgbomcpy
%     % two figures for two priors but all data sets in one figure 
%     for iPr=1:nPr
%         hFig=figure; ha;
%         for iDS=1:nDS
%             tmpX=dAngle(iDS,:,iPr);
%             tmpY=r2(iDS,1:(end-1),iPr);
%             plot(tmpX,tmpY,'.','color',pplot.cmap{iDS},'markersize',8);ha; % individual rand
%             s=plotReg(tmpX,tmpY,hFig,pplot.cmap{iDS}); % regression line
%                 disp([d(iDS).name ': slope=' num2str(s.beta(2),3) ',  p=' num2str(s.tstat.pval(2))]);
%             plot(0,r2(iDS,end,iPr),'o','color',pplot.cmap{iDS},'markersize',10,'markerfacecolor','w');ha; % acutal data
%         end % for iDS=1:nDS
%         axis tight; xlabel('angle from actual V(readout)'); ylabel('R^2 with BLS');
%     end % for iPr=1:nPr
else % if ~idPlotSep
    
    
    % separate figure for each data set but two priors together
    for iDS=1:nDS
%         if contains(d(iDS).name,['_' dsName '_']) % strcmp(d(iDS).name,['trajKS_prior_' dsName '_supportOnly.mat'])
            
            %% fig3c1: aV'X+b vs ts (V:readout vector), a & b fitted to BLS)
            figure; ha; setFigPos(1,iDS); % 1);
            dsName=d(iDS).name(14:17);
%             title(d(iDS).name(14:17));
             for iPr=1:nPr  % prior-specific
                 tm=linspace(T{iPr}(1),T{iPr}(end),length(xu{iDS,iPr}(:)));  tm=tm(:);% 17/21 for short/long
                 % bootstrap
                 if idBoot
                     sTe=squeeze(std(teBoot{1,iPr},0,1)); % [nBoot x timePoints] > [timePoints x 1]%  pmDBoot{1,iPr}(:)
                     if isempty(sTe)
                         sTe=zeros(size(tm));
                     end
                 else
                     sTe=zeros(size(tm));
                 end

                 if idTpMTs
                     % BLS
                     plot(tm,teBLS{iPr}-tm(:),'-','linewidth',lw2,'color',cmapPr(iPr,:));
                     % neural te
                     h=shadedErrorBar(tm(:),te{1,iPr}(:)-tm(:),sTe(:),{'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2},1); drawnow; % 2
%                      plot(tm(:),te{1,iPr}(:)-tm(:),'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % acutal
                     %                 plot(tm,te{2,iPr},'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % rotated
                 else%%%%% used
                     % BLS
                     yyaxis left;
                     if ~(idSiL & iPr==nPr)
                         plot(tm,teMat{iDS,iPr},'-','linewidth',lw2,'color',cmapPr(iPr,:));
                     end
                     % neural te
                     yyaxis right;
                     h=shadedErrorBar(tm(:),fxu{iDS,iPr}(:),sTe(:),{'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2},1); drawnow; % 2 pmD{iPr}(:)
%                      plot(tm,te{1,iPr},'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % acutal
                     %                 plot(tm,te{2,iPr},'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % rotated
                 end
             end
             if idTpMTs
                 axis tight; % 
                 xlim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                 plotHorizon(gca,0,[]); % plotIdentity(gca);
                 set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                     'ytick',-200:25:200); % [T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
             else%%%%% used
                 %                  axis tight;
                 yyaxis left;
                 ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                 plotIdentity(gca);
                 set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                     'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)});
                 
                 yyaxis right;
                 ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                 %                  plotIdentity(gca);
                 c=table2array(LMM{1}.Coefficients); cs1=c(1,1); cs2=c(2,1);
                 c=table2array(LMM{2}.Coefficients); cl1=c(1,1); cl2=c(2,1);
                 y1=(T{1}(1)-cs1)/cs2;
                 y2=(T{1}(3)-cs1)/cs2;
                 y3=(T{1}(end)-cs1)/cs2;
                 y4=(T{2}(3)-cl1)/cl2;
                 y5=(T{2}(end)-cl1)/cl2;
                 y1=round(y1*100)/100;y2=round(y2*100)/100;y3=round(y3*100)/100;
                 y4=round(y4*100)/100;y5=round(y5*100)/100;
                 set(gca,'ticklength',[0.01 0.01],'tickdir','out',...
                     'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{y1;y2;y3;y4;y5});
                 drawnow;
             end
             applytofig(gcf,optsExpFig);
% applytofig4keynote;
if idSaveFig
        savefig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/revision/figure3/'...
        dsName '.fig']);
    remLabel;yyaxis left;remLabel;
    exportfig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/revision/figure3/'...
        dsName '.eps'],optsExpFig);
%     savefig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/3/new/c/'...
%         dsName '.fig']);
end
             
             %%  fig3c12
%              figure; ha; setFigPos(2,1); % u'_ts
%              for iPr=1:nPr  % prior-specific
%                  tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr}));  tm=tm(:);% 17/21 for short/long
%                  % bootstrap
%                  if idBoot
%                      sTe=squeeze(std(teBoot{2,iPr},0,1)); % [nBoot x timePoints] > [timePoints x 1] %  pmDBoot{2,iPr}(:)
%                      if isempty(sTe)
%                          sTe=zeros(size(tm));
%                      end
%                  else
%                      sTe=zeros(size(tm));
%                  end
%                  
%                  if idTpMTs
%                      % BLS
%                      plot(tm,teBLS{iPr}-tm(:),'-','linewidth',lw2,'color',cmapPr(iPr,:));
%                      % neural te
%                      %                         plot(tm,te{1,iPr},'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % acutal
%                      plot(tm(:),te{2,iPr}-tm(:),'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % rotated
%                  else%%%%% used
%                      % BLS
%                      yyaxis left;
%                      if ~(idSiL & iPr==nPr)
%                          plot(tm,teBLS{iPr},'-','linewidth',lw2,'color',cmapPr(iPr,:));
%                      end
%                      yyaxis right;
%                      h=shadedErrorBar(tm(:),te{2,iPr}(:),sTe(:),{'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2},1); drawnow; % 2 pmDr{iPr}(:)
%                      % neural te
%                      %                         plot(tm,te{1,iPr},'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % acutal
% %                      plot(tm,te{2,iPr},'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % rotated
%                  end
%              end
%              if idTpMTs
%                  axis tight; % 
%                  xlim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
%                  plotHorizon(gca,0,[]); % plotIdentity(gca);
%                  set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
%                      'ytick',-500:50:500); % 'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
%              else%%%%% used
%                  yyaxis left;
%                  %                  axis tight; 
%                  ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
%                  plotIdentity(gca);
%                  set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
%                      'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)});
%                  
%                   yyaxis right;
% %                  axis tight; 
%                  ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
%                  c=table2array(LMM2{1}.Coefficients); cs1=c(1,1); cs2=c(2,1);
%                  c=table2array(LMM2{2}.Coefficients); cl1=c(1,1); cl2=c(2,1);
%                  y1=(T{1}(1)-cs1)/cs2; 
%                  y2=(T{1}(3)-cs1)/cs2;
%                  y3=(T{1}(end)-cs1)/cs2;
%                  y4=(T{2}(3)-cl1)/cl2;
%                  y5=(T{2}(end)-cl1)/cl2;
%                  y1=round(y1*100)/100;y2=round(y2*100)/100;y3=round(y3*100)/100;
%                  y4=round(y4*100)/100;y5=round(y5*100)/100;
%                  set(gca,'ticklength',[0.01 0.01],'tickdir','out',...
%                      'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{y1;y2;y3;y4;y5});
%                  drawnow;
%              end
% %             applytofig(gcf,optsExpFig);
% if idSaveFig
%     savefig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/3/new/c/'...
%         dsName '2.fig']);
% end

            
            %% fig3c2: R2(BLS) vs angle(V,V')
%             dAngleDot=cosd(dAngle);
% %             deg2polar=@(x,iPr) deg2rad(((iPr==1)+(iPr==nPr)*(-1))*x+90); % left for short, right for long
%             hFig=figure; setFigPos(floor((iDS-1)/4)+1,rem(iDS-1,4)+1); % eye for 1st row; LeftG, LeftH, RightG, RightH
%             for iPr=1:nPr
%                 subplot(1,nPr,iPr);
%                 
%                 tmpX=squeeze(dAngleDot(iDS,:,iPr)); % *2*(iPr-1.5);  % short in left
%                 tmpY=squeeze(r2(iDS,1:(end-1),iPr));
%                 
%                 % individual rand
% %                 polarplot(deg2polar(tmpX,iPr),tmpY,'.','color',tmpCmap{iPr,1}(1,:),'markersize',msize2);ha; 
%                 plot(tmpX,tmpY,'.','color',tmpCmap{iPr,1}(1,:),'markersize',msize2);ha; % individual rand
% 
%                 % moving average
%                 [tmpXsort,idSort]=sort(tmpX(:));
%                 x2=smooth(tmpXsort(:),nWin);y2=smooth(tmpY(idSort),nWin); 
% %                 x2=conv(tmpXsort(:),boxWin(:),'full');y2=conv(tmpY(idSort),boxWin(:),'full'); % valid
% %                 polarplot(deg2polar(x2,iPr),y2,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2);
%                 plot(x2,y2,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2);
%                 %             s=plotReg(tmpX,tmpY,hFig,tmpCmap{iPr,1}(1,:)); % regression line
%                 %             disp([d(iDS).name ': slope=' num2str(s.beta(2),3) ',  p=' num2str(s.tstat.pval(2))]);
%                 plot(1,r2(iDS,end,iPr),'o','color',cmapPr(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % acutal data 0 2*(iPr-1.5)
%                 plot(dAngleDot(iDS,iRandEx{iPr},iPr),r2(iDS,iRandEx{iPr},iPr),'d','color',cmapRot(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % *2*(iPr-1.5)
% 
% %                 % acutal data & example
% %                 polarplot(deg2polar(0,iPr),r2(iDS,end,iPr),'o','color',cmapPr(iPr,:),'markersize',msize2,'markerfacecolor','w','linewidth',lw);ha; % acutal data
% %                 polarplot(deg2polar(dAngle(iDS,iRandEx{iPr},iPr),iPr),r2(iDS,iRandEx{iPr},iPr),'d','color',cmapRot(iPr,:),'markersize',msize2,'markerfacecolor','w','linewidth',lw);ha;
%                 
% %                 set(gca,'thetalim',[0 180],'rlim',[0 1],'Rtick',0:.25:1,'ThetaTick',[0 30 60 90 120 150 180],'tickdir','out','TickLength',[0.02 0.02],...
% %                     'ThetaTickLabel',{'90','60','30','0','30','60','90'});
% %                 axis([0 90 0 1]); %                 axis tight; 
% axis([0 1 0 1]); % axis([-1 1 0 1]);
%                 xlabel('u\prime_t_s projected to u_t_s'); if iPr==1, ylabel('Similarity to Bayesian model'); end
%                 set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[-1:0.5:1],'ytick',0:0.5:1,...
%                     'xticklabel',[1 0.5 0 0.5 1]); % 0:30:90
%                 if iPr==2, set(gca,'yticklabel',[]); end
%                 
%                 box off;
%                 remAxLabel; remTickLabel;
%                 
%             end % for iPr=1:nPr
% %             applytofig(gcf,optsExpFig);
% if idSaveFig
%     savefig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/3/new/c/'...
%         dsName '3.fig']);
% end

            %% check a and b from te=a*XUts+b
%             figure; setFigPos(floor((iDS-1)/4)+1,rem(iDS-1,4)+1); % eye for 1st row; LeftG, LeftH, RightG, RightH
%             for iPr=1:nPr
%                 subplot(1,nPr,iPr);
%                 
%                 tmpX=squeeze(dAngleDot(iDS,:,iPr)); % *2*(iPr-1.5);  % short in left
%                 tmpY=squeeze(a(iDS,1:(end-1),iPr));
%                 
%                 % individual rand
%                 plot(tmpX,tmpY,'.','color',tmpCmap{iPr,1}(1,:),'markersize',msize2);ha; % individual rand
%                 
%                 % moving average
%                 [tmpXsort,idSort]=sort(tmpX(:));
%                 x2=smooth(tmpXsort(:),nWin);y2=smooth(tmpY(idSort),nWin);
%                 plot(x2,y2,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2);
%                 plot(1,a(iDS,end,iPr),'o','color',cmapPr(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % acutal data 0 2*(iPr-1.5)
% %                 plot(dAngleDot(iDS,iRandEx{iPr},iPr),a(iDS,iRandEx{iPr},iPr),'d','color',cmapRot(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % *2*(iPr-1.5)
% %                 axis([0 1 0 1]); % axis([-1 1 0 1]);
%                 xlabel('u\prime_t_s projected to u_t_s'); if iPr==1, ylabel('slope to Bayesian model'); end
% %                 set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[-1:0.5:1],'ytick',0:0.5:1,...
% %                     'xticklabel',[1 0.5 0 0.5 1]); % 0:30:90
% %                 if iPr==2, set(gca,'yticklabel',[]); end
%                 box off;
%                 
%             end % for iPr=1:nPr
%             
%             figure; setFigPos(floor((iDS-1)/4)+1,rem(iDS-1,4)+1); % eye for 1st row; LeftG, LeftH, RightG, RightH
%             for iPr=1:nPr
%                 subplot(1,nPr,iPr);
%                 
%                 tmpX=squeeze(dAngleDot(iDS,:,iPr)); % *2*(iPr-1.5);  % short in left
%                 tmpY=squeeze(b(iDS,1:(end-1),iPr));
%                 
%                 % individual rand
%                 plot(tmpX,tmpY,'.','color',tmpCmap{iPr,1}(1,:),'markersize',msize2);ha; % individual rand
%                 
%                 % moving average
%                 [tmpXsort,idSort]=sort(tmpX(:));
%                 x2=smooth(tmpXsort(:),nWin);y2=smooth(tmpY(idSort),nWin);
%                 plot(x2,y2,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2);
%                 plot(1,b(iDS,end,iPr),'o','color',cmapPr(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % acutal data 0 2*(iPr-1.5)
%                 %                 plot(dAngleDot(iDS,iRandEx{iPr},iPr),b(iDS,iRandEx{iPr},iPr),'d','color',cmapRot(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % *2*(iPr-1.5)
%                 %                 axis([0 1 0 1]); % axis([-1 1 0 1]);
%                 xlabel('u\prime_t_s projected to u_t_s'); if iPr==1, ylabel('intercept to Bayesian model'); end
%                 %                 set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[-1:0.5:1],'ytick',0:0.5:1,...
%                 %                     'xticklabel',[1 0.5 0 0.5 1]); % 0:30:90
%                 %                 if iPr==2, set(gca,'yticklabel',[]); end
%                 box off;
%                 
%             end % for iPr=1:nPr

%         end % if strcmp(d(iDS).name,'trajKS_prior_EL_H_supportOnly.mat')
    end % for iDS=1:nDS
    
    
    %% checking angle of uts/vtp across priors % H/G, Hand/Eye, Right/Left
%     figure; ha; setFigPos(1,5);
%     aUBoot2=reshape(aUBoot,nAnimal*nEH*nTarg,nBootAngle)'; % [1000 x 8]
%     violinPlot(aUBoot2,'showMM',2); plotHorizon(gca,90,[]);
%     plot(aU(:),'o');
%     set(gca,'xticklabel',[]);
%     xlabel('data sets'); ylabel('angle(u(short),u(long))');
%     
%     % check bootstrap
%     figure;setFigPos(2,5);violinPlot(reshape(aTmp,2*2*2*2,nBootAngle)','showMM',2);
%     set(gca,'xticklabel',[]); plotHorizon(gca,90,[]);
%     xlabel('data sets (animal,effector,direction,prior)'); ylabel('angle(u(original),u(bootstrap))');
%     
%     % compare to vtp
%     aV=cat(3,[48.1756   39.9877;...
%    45.6345   45.0398],[...
%    31.5721   41.4173;...
%    55.5135   48.5545]);
%     
% figure; plot(aU(:),aV(:),'ko');
% plotIdentity(gca); xlabel('angle(u(short),u(long))');
%     ylabel('angle(v(short),v(long))');
    
end % if ~idPlotSep

%%
function [teBoot,pmDBoot]=estTeRv(teBoot,nm,nPr,nTs,idDim,T,idDetrend,idRand,idBLSts,pmDBoot,varargin)
% only for bootstraping
if ~isempty(varargin)
    idSiL=varargin{1}; % nPr=3
else
    idSiL=0; % nPr=2
end

load(nm); % D optimD
% nm=d(iDS).name;

% animal-specific wM
if strcmp(nm,'_H_')
    load('H_RSGprior_DMFC.mat','wFit');
    wm=wFit.w_m; % 0.0508
else
    load('G_RSGprior_DMFC.mat','wFit');
    wm=wFit.w_m; % 0.0481
end

% estimate Vreadout
if idSiL
    for iPr=1:nPr  % prior-specific
        idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
        Vro{iPr}=estReadout(D(idTs)); % [dim x 1]
    end % for iPr=1:nPr
    %         save(d(iDS).name,'Vro','-append');
    %     else
    %         Vro=dVro.Vro; % {pr}[dim x 1]
    %     end
else
    for iPr=1:nPr  % prior-specific
        idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
        Vro{iPr}=estReadout(D(idTs)); % [dim x 1]
    end % for iPr=1:nPr
end % if idSiL

% modify dimension if needed
if idDim==0 % use optimD (PVAF>75%)
    dim=optimD;
else % trim data
    %         try
    proj_matrix=proj_matrix(:,1:dim);
    %         catch
    %             disp('');
    %         end
    for iD=1:length(D),D(iD).data=D(iD).data(1:dim,:);end
end

     %% projection: after averaging across ts
    for iPr=1:nPr  % prior-specific
        idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
        dTmp=struct2mat(D(idTs),'data'); % [PC x time x ts]
%         if idAvgAtt %%%%% actually used
            mD{iPr}=nanmean(dTmp,3); % [PC x time x ts] > [PC x time]
            
            if ~(idSiL & iPr==nPr) % no need for this for SiL
                % smoothing after averaging across ts with attrition
%                 if idSmoothAfterAvgAtt
                    binWidth=1; % 20 ms bin in PC
                    kernSD=2; % 40ms smoothing was used in PC
                    mD{iPr} = smoother(mD{iPr}, kernSD, binWidth);
%                 end
                
                % removing points out of prior support (shortest ts has two data points: 4 20-ms bin before & after > remove 3 for both)
                nTShortest=size(D(idTs(1)).data,2); % 2 for both short/long
                mD{iPr}=mD{iPr}(:,nTShortest:end);
                %             pmD{iPr}=pmD{iPr}(nTShortest:end);
                dTmp=dTmp(:,nTShortest:end,:);
            else % still removing points out of prior support for SiL
                nTShortest=size(D(1).data,2); % 2 for both short/long
                mD{iPr}=mD{iPr}(:,nTShortest:end);
                dTmp=dTmp(:,nTShortest:end,:);
            end % if idSiL & iPr~=nPr
            
            % do actual projection
            pmD{iPr}=Vro{iPr}(:)'*mD{iPr}; % [1 x time]
            pmDBoot{idRand,iPr}=cat(1,pmDBoot{idRand,iPr},pmD{iPr}(:)');
            
            % R2 with BLS
            if idSiL
                if iPr==nPr
                    % R2 with BLS
                    tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:); % 17/21 for short/long
                    if idBLSts
                        teBLS{iPr}=BLSts(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                    else
                        teBLS{iPr}=BLS(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                    end
                    if idDetrend
                        stat{iPr}=regstats(teBLS{iPr}-tm(:),detrend(pmD{iPr}),'linear'); % rsquare
                        teBoot{idRand,iPr}=cat(1,teBoot{idRand,iPr},stat{iPr}.yhat(:)'+tm(:)'); % [nBoot x timePoints] to Avg
                    else %%%%% actually used
                        nObs=sum(~isnan(squeeze(dTmp(1,:,:))),2); % [time x 1]
                        % lsocv(pmD{iPr}(:),teBLS{iPr}(:),nObs) will do the job too
                        LM=fitlm(pmD{iPr}(:),teBLS{iPr}(:),'linear',...
                            'Weights',nObs);
                        teBoot{idRand,iPr}=cat(1,teBoot{idRand,iPr},LM.predict');
                        %                 stat{iPr}=regstats(teBLS{iPr},pmD{iPr},'linear'); % rsquare
                        %                 te{1,iPr}=stat{iPr}.yhat(:);
                    end % idDetrend
                    %     r2(iDS,end,iPr)=stat{iPr}.rsquare;
                    %     disp(['R^2: ' num2str(stat{iPr}.rsquare,3)]);
                end % iPr
            else % idSiL
                % R2 with BLS
                tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:); % 17/21 for short/long
                if idBLSts
                    teBLS{iPr}=BLSts(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                else
                    teBLS{iPr}=BLS(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                end
                if idDetrend
                    stat{iPr}=regstats(teBLS{iPr}-tm(:),detrend(pmD{iPr}),'linear'); % rsquare
                    teBoot{idRand,iPr}=cat(1,teBoot{idRand,iPr},stat{iPr}.yhat(:)'+tm(:)'); % [nBoot x timePoints] to Avg
                else %%%%% actually used
                    nObs=sum(~isnan(squeeze(dTmp(1,:,:))),2); % [time x 1]
                    % lsocv(pmD{iPr}(:),teBLS{iPr}(:),nObs) will do the job too
                    LM=fitlm(pmD{iPr}(:),teBLS{iPr}(:),'linear',...
                        'Weights',nObs);
                    teBoot{idRand,iPr}=cat(1,teBoot{idRand,iPr},LM.predict');
                    %                 stat{iPr}=regstats(teBLS{iPr},pmD{iPr},'linear'); % rsquare
                    %                 te{1,iPr}=stat{iPr}.yhat(:);
                end % idDetrend
                %     r2(iDS,end,iPr)=stat{iPr}.rsquare;
                %     disp(['R^2: ' num2str(stat{iPr}.rsquare,3)]);
            end % if idSiL
            %         end % idAvgAtt
    end % iPr
%%
% 1. RNN must reflect the same task analysis pipeline as implemented on data
% 2. Same RNN - different wm - should show increased bias - we reverse engineer measurement mechanism
% 3. Different RNNs either producing or not producing bias - what is different in measurement epoch.
% 
% Other analysis:
% 
% Random projections around Set + 200 ms IC space (v1) and at the end of measurement phase (v2). Find the cone of vectors plus spurious vector sets between v1 and v2 that have a high R2. Use cross validation to determine meaningful subset or at least get closer to the cone.