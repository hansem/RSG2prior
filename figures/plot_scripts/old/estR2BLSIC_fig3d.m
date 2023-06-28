function estR2BLSIC_fig3d(varargin)

% 2019/1/24:
% test across-prior projection

% 2019/1/20:
% - test coding at IC: ts(linear), te(BLS), 1/ts, 1/te
% - test sigmoid: simplest > dXu(extreme) vs dXu(middle)
% - measure set transient [Set IC]

% 2018/10/8
% checking angle of uts/vtp across priors

% 2018/9/19
% plotyy
% use animal-, condition-specific wM > supplementary figures
%
% estR2BLSIC_fig3d('ER_H'); estR2BLSIC_fig3d('EL_H'); estR2BLSIC_fig3d('HR_H'); estR2BLSIC_fig3d('HL_H'); estR2BLSIC_fig3d('ER_G'); estR2BLSIC_fig3d('EL_G'); estR2BLSIC_fig3d('HR_G'); estR2BLSIC_fig3d('HL_G'); 
% idSaveFig=1;

% 2018/6/19
% fig3d2: going back to cartesian but dot product (cos)
%           subplot for short/long
% fig3d1: error bar from bootstrapping (not for worst eg for now)
%            tp-ts format (not tp) - idTpMTs
%           BLS more fine-grained

% 2018/6/13
% polar plot for R2: left for short, right for long

% 2018/6/11
% fig3e1: aW'X+b vs ts (W:readout vector), a & b fitted to BLS)
% fig3e2: R2(BLS) vs angle(W,W')
% example data set: H EL

% difference from estR2BLSrotation.m
% 1) readout vector: average difference vector (c.f. mean(ts5-ts1 ts4-ts2))
% 2) IC specific: PCA for tp, tIC=200ms
%
% estimating R^2 b/t IC states projected to readout vector and BLS
% repeat with randomizing readout vector (nRand)
% across animals and conditions
% use trajKS_tp_(CondNm)_(H/G)_bin20_smth40.mat
% modified analSet_randProj.m
% depends on estReadout.m, randVect.m, BLS.m
% 2018/6/3


%% initial
initRSG2prior;
cd(psthDir);
load pplot.mat; % pplot.cmap for separate figure
nWin=50; % window for moving average
boxWin=ones(nWin,1)/nWin;

idTpMTs=0; % 1; % plot te-ts vs ts
idBLSts=0; % 1; % use BLSts after marginzaling over tm
idSaveFig=0; % 1;

idPlotSep=true; % plot separately across animals and conditions
nRand=1000;

idDim=0; % 3
% H/G, ER/EL/HR/HL: 8 7 10 9 7 7 9 6 for priorIncSL, 8 7 9 8 10 9 11 7

% if isempty(varargin)
%     dsName='ER_H'; %'ER_G'; %'ER_H'; % 'HL_G'; % EL_H'; % eye left H
% else
%     dsName=varargin{1};
% end

% BLS
% wm=0.05;

% bootstrap
idBoot=0; % 
nBoot=1000;

% IC specific
tIC=200;

idDebug=0;

% plot
cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
cmapRot=[rgb('DeepPink'); rgb('Aqua')];
msize=4; % 10; % 4; % 6; % 2;
msize2=1;
lw=.5;% 1.2; % 1; % 3; % 1.5;
lw2=2;

%% main
d=dir('trajKS_tp_*.mat'); % prior*_supportOnly.mat'); % IC specific % _bin20_smth40
nDS=length(d);

% if exist('R2BLS_IC.mat') % rotation.mat') % IC specific
%     load('R2BLS_IC.mat'); % rotation.mat'); % dAngle r2
% else
% output
dAngle=nan(nDS,nRand,nPr);
r2=nan(nDS,1+nRand,nPr); % last for actual data
neuralTe=cell(2,nPr); % actualData-neural prediction/worstExample
teBoot=cell(2,nPr); % te from bootstrap

% 2019/1/20: test coding at IC: ts(linear), te(BLS), 1/ts, 1/te
nModel=4;
xu=nan(nDS,nPr,nTspp);
RMSE=nan(nDS,nPr,nModel); 

% 2019/1/24: test across-prior projection
tePr=nan(nDS,nPr,nTspp); % across-prior prediction
te=nan(nDS,nPr,nTspp); % behavioral
v=cell(nDS,nPr); % decoding vector
c1=nan(nDS,nPr); % slope
c2=nan(nDS,nPr); % intercept
hFig=figure; setFigPos(1,2);

aU=nan(nAnimal,nEH,nTarg);
idBootAngle=1;
nBootAngle=1000; % 000;
aUBoot=nan(nAnimal,nEH,nTarg,nBootAngle); % angle between priors
aTmp=nan(nAnimal,nEH,nTarg,nPr,nBootAngle); % angle between original vs bootstrapped data

for iDS=1:nDS
    %% checking angle of uts/vtp across priors % H/G, Hand/Eye, Right/Left
%     binSize=20;
%         kIC=round(tIC/binSize); % IC specific
%         
%     disp(['===== ' d(iDS).name ' =====']);
%     % figure out which data set
%     if contains(d(iDS).name,'_ER_H_')
%         iAnimal=1; iEH=1; iTarg=1;iAnimalNm='H';
%         tmpD=load(['trajKS_tp_ER_H_bin20_smth40.mat']); % D
%         dsName='ER_H';
%     elseif contains(d(iDS).name,'_EL_H_')
%         iAnimal=1; iEH=1; iTarg=2;iAnimalNm='H';
%         tmpD=load(['trajKS_tp_EL_H_bin20_smth40.mat']); % D
%         dsName='EL_H';
%     elseif contains(d(iDS).name,'_HR_H_')
%         iAnimal=1; iEH=2; iTarg=1;iAnimalNm='H';
%         tmpD=load(['trajKS_tp_HR_H_bin20_smth40.mat']); % D
%         dsName='HR_H';
%     elseif contains(d(iDS).name,'_HL_H_')
%         iAnimal=1; iEH=2; iTarg=2;iAnimalNm='H';
%         tmpD=load(['trajKS_tp_HL_H_bin20_smth40.mat']); % D
%         dsName='HL_H';
%     elseif contains(d(iDS).name,'_ER_G_')
%         iAnimal=2; iEH=1; iTarg=1;iAnimalNm='G';
%         dsName='ER_G';
%         tmpD=load(['trajKS_tp_ER_G_bin20_smth40.mat']); % D
%     elseif contains(d(iDS).name,'_EL_G_')
%         iAnimal=2; iEH=1; iTarg=2;iAnimalNm='G';
%         tmpD=load(['trajKS_tp_EL_G_bin20_smth40.mat']); % D
%         dsName='EL_G';
%     elseif contains(d(iDS).name,'_HR_G_')
%         iAnimal=2; iEH=2; iTarg=1;iAnimalNm='G';
%         tmpD=load(['trajKS_tp_HR_G_bin20_smth40.mat']); % D
%         dsName='HR_G';
%     elseif contains(d(iDS).name,'_HL_G_')
%         iAnimal=2; iEH=2; iTarg=2;iAnimalNm='G';
%         tmpD=load(['trajKS_tp_HL_G_bin20_smth40.mat']); % D
%         dsName='HL_G';
%     end
%         
%     for iPr=1:nPr  % prior-specific
%         idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
%         
%         Dtmp=struct2mat(tmpD.D(idTs),'data'); % [dim Maxtime 5conditions]
%             mDtmp=squeeze(nanmean(diff(Dtmp(:,kIC,:),1,3),3)); % [dim x 1] mean-difference
%             vtp{iAnimal,iEH,iTarg,iPr}=normalize(mDtmp,1); % IC specific
%         
%         % bootstrap
%         bootDir=['/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'...
%             'trajKS_tp/']; % trajKS_prior_supportOnly/']; % '/Volumes/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'; % macmini
%         
%         epSelect='tp'; % prior'; % 'ts'; % 'tp';
%         for iBoot=1:nBootAngle
%             
%             tsNm=fullfile(bootDir,...
%                 ['trajKS_' epSelect '_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40_' num2str(iBoot) '.mat']); % '_bin20_smth40.mat']; % ['trajKS_ts_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40.mat'];
%             tmpBoot=load(tsNm); % D
%             Dtmp=struct2mat(tmpBoot.D(idTs),'data'); % [dim Maxtime 5conditions]
%             mDtmp=squeeze(nanmean(diff(Dtmp(:,kIC,:),1,3),3)); % [dim x 1] mean-difference
%             vtpBoot{iAnimal,iEH,iTarg,iPr,iBoot}=normalize(mDtmp,1); % IC specific
%             
%             % check angle betwen vtp and vtp(bootstrap)
% %             disp(tsNm);
% %             disp(angleVectors(vtp{iAnimal,iEH,iTarg,iPr},vtpBoot{iAnimal,iEH,iTarg,iPr,iBoot}));
%             aTmp(iAnimal,iEH,iTarg,iPr,iBoot)=angleVectors(vtp{iAnimal,iEH,iTarg,iPr},vtpBoot{iAnimal,iEH,iTarg,iPr,iBoot});
% %         % check trajectory
% %         [diffV,sSet]=plotTraj_figure3(tsNm);
% %         hSS=gcf;setFigPos(2,2);
% %         if iPr==nPr
% %             ha; arrow(sSet(1:3,1+nTspp),sSet(1:3,1+nTspp)+vtp{iAnimal,iEH,iTarg,iPr}(1:3),'Type','line','Color','b','linewidth',2); % original u
% %         else
% %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+vtp{iAnimal,iEH,iTarg,iPr}(1:3),'Type','line','Color','r','linewidth',2);  % original u
% %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+vtpBoot{iAnimal,iEH,iTarg,iPr,iBoot},'Type','line','Color',[0.5 0 0],'linewidth',1);
% %         end
% % %         waitforbuttonpress;
% % %         close;
%         
%             if iPr==nPr
%                 aUBoot(iAnimal,iEH,iTarg,iBoot)=angleVectors(vtpBoot{iAnimal,iEH,iTarg,1,iBoot}(:),vtpBoot{iAnimal,iEH,iTarg,2,iBoot}(:));
%             end   
%         end
%     
%     end % for iPr=1:nPr
%     
% %     % check trajectory
% %         [diffV,sSet]=plotTraj_figure3(d(iDS).name);
% %         hSS=gcf;setFigPos(2,1);
% % %         if iPr==nPr % long
% %             ha; arrow(sSet(1:3,1+nTspp),sSet(1:3,1+nTspp)+vtp{iAnimal,iEH,iTarg,nPr}(1:3),'Type','line','Color','b','linewidth',2);  % original u
% %             ha; arrow(sSet(1:3,1+nTspp),sSet(1:3,1+nTspp)+ vtpBoot{iAnimal,iEH,iTarg,nPr,iBoot}(1:3),'Type','line','Color',[0 0 0.5],'linewidth',1);
% % %             else % short
% %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+vtp{iAnimal,iEH,iTarg,1}(1:3),'Type','line','Color','r','linewidth',2);  % original u
% %             ha; arrow(sSet(1:3,1),sSet(1:3,1)+vtpBoot{iAnimal,iEH,iTarg,1,iBoot},'Type','line','Color',[0.5 0 0],'linewidth',1);
% % %         end
%     
%     aU(iAnimal,iEH,iTarg)=angleVectors(vtp{iAnimal,iEH,iTarg,1}(:),vtp{iAnimal,iEH,iTarg,2}(:));
    
    
    %% main
%     if strcmp(d(iDS).name,['trajKS_tp_' dsName '_bin20_smth40.mat'])
        
    disp(['===== ' d(iDS).name ' =====']);
    load(d(iDS).name); % binSize smthWidth optimD use_sqrt proj_matrix keep_neurons D eigenvalues meanPSTH
    kIC=round(tIC/binSize); % IC specific
    
    % animal-specific wM
    if strfind(d(iDS).name,'_H.mat')
        load('H_RSGprior_DMFC.mat','wFitCond');
        if contains(d(iDS).name,'_ER_') %strcmp(dsName(1:2),'ER')
            wm=wFitCond(1,1).w_m; % 0.0508
            dsName='ER_H';
        elseif contains(d(iDS).name,'_EL_') %strcmp(dsName(1:2),'EL')
            wm=wFitCond(1,2).w_m; % 0.0508
            dsName='EL_H';
        elseif contains(d(iDS).name,'_HR_') %strcmp(dsName(1:2),'HR')
            wm=wFitCond(2,1).w_m; % 0.0508
            dsName='HR_H';
        elseif contains(d(iDS).name,'_HL_') %strcmp(dsName(1:2),'HL')
            wm=wFitCond(2,2).w_m; % 0.0508
            dsName='HL_H';
        end
        %         wm=wFit.w_m; % 0.0508
        iAnimalNm='H';
    else
        load('G_RSGprior_DMFC.mat','wFitCond');
        if contains(d(iDS).name,'_ER_') %strcmp(dsName(1:2),'ER')
            wm=wFitCond(1,1).w_m; % 0.0508
            dsName='ER_G';
        elseif contains(d(iDS).name,'_EL_')%strcmp(dsName(1:2),'EL')
            wm=wFitCond(1,2).w_m; % 0.0508
            dsName='EL_G';
        elseif contains(d(iDS).name,'_HR_') %strcmp(dsName(1:2),'HR')
            wm=wFitCond(2,1).w_m; % 0.0508
            dsName='HR_G';
        elseif contains(d(iDS).name,'_HL_') %strcmp(dsName(1:2),'HL')
            wm=wFitCond(2,2).w_m; % 0.0508
            dsName='HL_G';
        end
%         wm=wFit.w_m; % 0.0481
        iAnimalNm='G';
    end
    
    % estimate Vreadout
%     dVro=load(d(iDS).name); % ,'Vro');
%     if ~isfield(dVro,'Vro')
        for iPr=1:nPr  % prior-specific
            idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
            Dtmp=struct2mat(D(idTs),'data'); % [dim Maxtime 5conditions]
            mDtmp=squeeze(nanmean(diff(Dtmp(:,kIC,:),1,3),3)); % [dim x 1] mean-difference
            Vro{iPr}=normalize(mDtmp,1); % IC specific
            v{iDS,iPr}=Vro{iPr}(:);
%             Vro{iPr}=estReadout(D(idTs)); % [dim x 1] 
        end % for iPr=1:nPr
%         save(d(iDS).name,'Vro','-append');
%     else
%         Vro=dVro.Vro; % {pr}[dim x 1]
%     end
    
    % debug
    if idDebug
        plotTraj_figure3(d(iDS).name); setFigPos(2,3);
        ha; plot3([0;Vro{1}(1)],[0;Vro{1}(2)],[0;Vro{1}(3)],'r--')
        ha; plot3([0;Vro{2}(1)],[0;Vro{2}(2)],[0;Vro{2}(3)],'b--')
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
    
        %% bootstrap
        if idBoot
            bootDir=['/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'...
                'trajKS_tp/']; % '/Volumes/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'; % macmini
            
            epSelect='tp'; % 'ts'; % 'prior';
            for iBoot=1:nBoot
                tsNm=fullfile(bootDir,...
                    ['trajKS_' epSelect '_' dsName(1:2) '_' iAnimalNm '_bin20_smth40_' num2str(iBoot) '.mat']); % trajKS_tp_EL_G_bin20_smth40_20
                teBoot=estTeRv(teBoot,tsNm,nPr,nTs,idDim,T,kIC,1,idBLSts); % (teBoot,tsNm,nPr,T); last 1 for actual data
            end
        end
    
    %% projection: after averaging across ts
    for iPr=1:nPr  % prior-specific
        idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
        Dtmp=struct2mat(D(idTs),'data'); % [PC x time x ts] 
        mD{iPr}=squeeze(Dtmp(:,kIC,:)); % [PC x ts]
        pmD{iPr}=Vro{iPr}(:)'*mD{iPr}; % [1 x ts]
        
        % R2 with BLS
        tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:); % 5ts
        if idBLSts
            teBLS{iPr}=BLSts(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
        else
            teBLS{iPr}=BLS(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
        end
        te(iDS,iPr,:)=teBLS{iPr};
        stat{iPr}=regstats(teBLS{iPr},pmD{iPr},'linear'); % rsquare
        neuralTe{1,iPr}=stat{iPr}.yhat(:);
        r2(iDS,end,iPr)=stat{iPr}.rsquare;
        disp(['R^2: ' num2str(stat{iPr}.rsquare,3)]);
        
        c2(iDS,iPr)=stat{iPr}.beta(1); % intercept
        c1(iDS,iPr)=stat{iPr}.beta(2); % slope
        
        % 2019/1/20: test coding at IC: ts(linear), te(BLS), 1/ts, 1/te
        xu(iDS,iPr,:)=pmD{iPr}; % RMSE=nan(nDS,nPr,nModel);
        sTmp=regstats(pmD{iPr},T{iPr},'linear'); % ts(linear) mse
        RMSE(iDS,iPr,1)=sqrt(sTmp.mse);
        sTmp=regstats(pmD{iPr},teBLS{iPr},'linear'); %te(BLS) mse
        RMSE(iDS,iPr,2)=sqrt(sTmp.mse);
        sTmp=regstats(pmD{iPr},1./(T{iPr}),'linear'); % 1/ts mse
        RMSE(iDS,iPr,3)=sqrt(sTmp.mse);
        sTmp=regstats(pmD{iPr},1./(teBLS{iPr}),'linear'); % 1/te mse
        RMSE(iDS,iPr,4)=sqrt(sTmp.mse);
        
        % debug
        if idDebug
            figure; setFigPos(2,1);plotyy(tm,teBLS{iPr},tm,pmD{iPr});
            figure; setFigPos(2,2);plot(tm,teBLS{iPr},'-',tm,stat{iPr}.yhat,'x','color',tmpCmap{iPr,1}(1,:));
            axis tight; plotIdentity(gca);
        end
    end % iPr
    
    for iPr=1:nPr  % prior-specific
        % 2019/1/24: test across-prior projection
%         tePr(iDS,iPr,:)=c2(iDS,3-iPr)+c1(iDS,3-iPr)*(v{iDS,3-iPr}'*mD{iPr});% 
        tePr(iDS,iPr,:)=c2(iDS,iPr)+c1(iDS,iPr)*(v{iDS,3-iPr}'*mD{iPr}); % use scaling/offset of original
        
        figure(hFig); 
        for iTs=1:nTspp
            plot(te(iDS,iPr,iTs),tePr(iDS,iPr,iTs),'color',tmpCmap{iPr,1}(iTs,:),'marker',pplot.marker8{iDS}); ha;
        end % for iTs=1:nTspp
        
        % plot example
        if contains(d(iDS).name,'_EL_') & contains(d(iDS).name,'_H.mat')
            if iPr==1, hFig99=figure; setFigPos(2,3); else figure(hFig99); end
            tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:); % 5ts
            plot(tm,squeeze(tePr(iDS,iPr,:)),'color',cmapPr(iPr,:),'marker',pplot.marker8{iDS}); ha;
            tm=T{iPr}(1):T{iPr}(end); tm=tm(:); % 5ts
            plot(tm,BLSts(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform'),'-','color',cmapPr(iPr,:)); ha;
        end % tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:); % 5ts
    end
    
    % random readout vectors: uniformly distributed [0 90] deg
    for iPr=1:nPr  % prior-specific
        aRand=rand(nRand,1)*90; % (0 90)
        dAngle(iDS,:,iPr)=aRand;
        for iRand=1:nRand
            rVro=randVect(Vro{iPr},aRand(iRand),1);
            
            % projection
            pmDr{iPr}=rVro(:)'*mD{iPr}; % [1 x time]
%             pmDr{iPr}=pmDr{iPr}(nTShortest:end); % removing points out of prior support
            statTmp{iPr}=regstats(teBLS{iPr},pmDr{iPr},'linear'); % rsquare
            r2(iDS,iRand,iPr)=statTmp{iPr}.rsquare;
%             disp(['dAngle: ' num2str(aRand(iRand),3) ', R^2: ' num2str(statTmp.rsquare,3)]);
            
            % pick up worst example 
            if statTmp{iPr}.rsquare<0.4 & aRand(iRand)>70 & isempty(neuralTe{2,iPr}) & aRand(iRand)<85
                neuralTe{2,iPr}=statTmp{iPr}.yhat(:); iRandEx{iPr}=iRand;
                disp(['dAngle: ' num2str(aRand(iRand),3) ', R^2: ' num2str( statTmp{iPr}.rsquare,3)]);
                
                if idBoot
                    % bootstrap
                     bootDir=['/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'...
                         'trajKS_tp/'];
                    for iBoot=1:nBoot
                        tsNm=fullfile(bootDir,...
                            ['trajKS_' epSelect '_' dsName(1:2) '_' iAnimalNm '_bin20_smth40_' num2str(iBoot) '.mat']); % '_bin20_smth40.mat']; % ['trajKS_ts_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40.mat'];
                        teBoot=estTeRv(teBoot,tsNm,nPr,nTs,idDim,T,kIC,2,idBLSts); % (teBoot,tsNm,nPr,T); last 1 for actual data
                    end
                end
                
            end
            
            % debug
            if idDebug
                tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); 
                figure; setFigPos(1,1);plotyy(tm,teBLS{iPr},tm,pmDr{iPr}); title(aRand(iRand));
                figure; setFigPos(1,2);plot(tm,teBLS{iPr},'-',tm,statTmp{iPr}.yhat,'x','color',tmpCmap{iPr,1}(1,:));  title(statTmp{iPr}.rsquare);
                axis tight; plotIdentity(gca);
            end
        end % for iRand=1:nRand
    end % for iPr=1:nPr  % prior-specific
%     end % . if strcmp(d(iDS).name,'trajKS_prior_EL_H_supportOnly.mat')
end % for iDS=1:nDS

figure(hFig); 
axisEqual;
plotIdentity(gca);
xlabel('behavioral t_e'); ylabel('neural t_e predicted across prior');

% save('R2BLS_IC.mat','dAngle','r2','tIC'); % rotation.mat','dAngle','r2');
% end % if exist('R2BLSrotation.mat')

 % 2019/1/20: test coding at IC: ts(linear), te(BLS), 1/ts, 1/te
% model-based approach: RMSE(BLS) vs RMSE(linear) both have two free params
RMSE2=reshape(RMSE,nDS*nPr,nModel); % RMSE=nan(nDS,nPr,nModel);
figure; setFigPos(1,1); ha;
plot(RMSE2','k'); 
shadedErrorBar(1:nModel,nanmean(RMSE2,1),sem(RMSE2,1,'omitnan'),{'k-','linewidth',2});
xlabel('model'); ylabel('RMSE');
set(gca,'xtick',1:4,'xticklabel',{'t_s','t_e','1/t_s','1/t_e'});

% compare te vs 1/te
figure; setFigPos(2,1); ha;
for iPr=1:nPr
    plot(RMSE(:,iPr,2),RMSE(:,iPr,4),'o','color',cmapPr(iPr,:));
    p=signrank(RMSE(:,iPr,2),RMSE(:,iPr,4));
    disp(['p(signrank): ' num2str(p)]);
end
axisEqual;
plotIdentity(gca);
xlabel('RMSE(t_e)');
ylabel('RMSE(1/t_e)');

% compare ts vs te
figure; setFigPos(2,2); ha;
for iPr=1:nPr
    plot(RMSE(:,iPr,1),RMSE(:,iPr,2),'o','color',cmapPr(iPr,:));
    p=signrank(RMSE(:,iPr,1),RMSE(:,iPr,2));
    disp(['p(signrank): ' num2str(p)]);
end
axisEqual;
plotIdentity(gca);
xlabel('RMSE(t_s)');
ylabel('RMSE(t_e)');

% compare 1/ts vs 1/te
figure; setFigPos(2,3); ha;
for iPr=1:nPr
    plot(RMSE(:,iPr,3),RMSE(:,iPr,4),'o','color',cmapPr(iPr,:));
    p=signrank(RMSE(:,iPr,3),RMSE(:,iPr,4));
    disp(['p(signrank): ' num2str(p)]);
end
axisEqual;
plotIdentity(gca);
xlabel('RMSE(1/t_s)');
ylabel('RMSE(1/t_e)');

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
        if strcmp(d(iDS).name,['trajKS_tp_' dsName '.mat']) % _bin20_smth40
        
            %% fig3e1: aV'X+b vs ts (V:readout vector), a & b fitted to BLS)
            figure; ha; setFigPos(1,1);
            for iPr=1:nPr  % prior-specific
                tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr}));  tm=tm(:); % 17/21 for short/long
                % bootstrap
                if idBoot
                    sTe=squeeze(std(teBoot{1,iPr},0,1)); % [nBoot x timePoints] > [timePoints x 1]
                    if isempty(sTe)
                        sTe=zeros(size(tm));
                    end
                else
                    sTe=zeros(size(tm));
                end
                 
                % BLS fine grained
                tmBLS=linspace(T{iPr}(1),T{iPr}(end),1000); tmBLS=tmBLS(:); % 5ts
                if idBLSts
                    teBLS2{iPr}=BLSts(tmBLS,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                else
                teBLS2{iPr}=BLS(tmBLS,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                end
                
                if idTpMTs
                    % BLS
                    plot(tmBLS,teBLS2{iPr}-tmBLS,'-','linewidth',lw2,'color',cmapPr(iPr,:));
                    % neural te
                    h=shadedErrorBar(tm(:),neuralTe{1,iPr}(:)-tm(:),sTe(:),{'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2},1); drawnow; % 2
%                     plot(tm,neuralTe{1,iPr}-tm,'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % acutal
                    %                 plot(tm,neuralTe{2,iPr},'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize,'markerfacecolor','w'); drawnow; % rotated
                else %%%%% used
                    % BLS
                    yyaxis left;
                    plot(tmBLS,teBLS2{iPr},'-','linewidth',lw2,'color',cmapPr(iPr,:));
                    % neural te
                    yyaxis right;
                    h=shadedErrorBar(tm(:),neuralTe{1,iPr}(:),sTe(:),{'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2},1); drawnow; % 2
%                     plot(tm,neuralTe{1,iPr},'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % acutal
                    %                 plot(tm,neuralTe{2,iPr},'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize,'markerfacecolor','w'); drawnow; % rotated
                end
            end
            if idTpMTs
                axis tight; %
                xlim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                plotHorizon(gca,0,[]); % plotIdentity(gca);
                set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                    'ytick',-200:25:200);
            else
                yyaxis left;
%                 axis tight; 
                ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                plotIdentity(gca);
                set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                     'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)});
                 
                 yyaxis right;
                 ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                 xlim([480 1200]);
%                  plotIdentity(gca);
                 cs1=stat{1}.beta(2); cs2=stat{1}.beta(1);
                cl1=stat{2}.beta(2); cl2=stat{2}.beta(1);
                 y1=(T{1}(1)-cs1)/cs2; 
                 y2=(T{1}(3)-cs1)/cs2;
                 y3=(T{1}(end)-cs1)/cs2;
                 y4=(T{2}(3)-cl1)/cl2;
                 y5=(T{2}(end)-cl1)/cl2;
                 y1=round(y1*100)/100;y2=round(y2*100)/100;y3=round(y3*100)/100;
                 y4=round(y4*100)/100;y5=round(y5*100)/100;
                 set(gca,'ticklength',[0.01 0.01],'tickdir','out',...
                     'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{y1;y2;y3;y4;y5});
            end
%             applytofig(gcf,optsExpFig);
            if idSaveFig
    savefig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/3/new/d/'...
        dsName '.fig']);
end
            %% fig3e12
            figure; ha; setFigPos(2,1); % v'_tp
            for iPr=1:nPr  % prior-specific
                tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:);% 17/21 for short/long
                
                 % bootstrap
                 if idBoot
                     sTe=squeeze(std(teBoot{2,iPr},0,1)); % [nBoot x timePoints] > [timePoints x 1]
                     if isempty(sTe)
                         sTe=zeros(size(tm));
                     end
                 else
                     sTe=zeros(size(tm));
                 end
                 
                tmBLS=linspace(T{iPr}(1),T{iPr}(end),1000); tmBLS=tmBLS(:); % 5ts
                if idTpMTs
                    % BLS
                    plot(tmBLS,teBLS2{iPr}-tmBLS,'-','linewidth',lw2,'color',cmapPr(iPr,:));
                    % neural neuralTe
                    %                         plot(tm,neuralTe{1,iPr},'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize,'markerfacecolor','w'); drawnow; % acutal
                    plot(tm,neuralTe{2,iPr}-tm,'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % rotated
                else %%%%% used
                    % BLS
                    plot(tmBLS,teBLS2{iPr},'-','linewidth',lw2,'color',cmapPr(iPr,:));
                     h=shadedErrorBar(tm(:),neuralTe{2,iPr}(:),sTe(:),{'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2},1); drawnow; % 2
                    % neural te
                    %                         plot(tm,neuralTe{1,iPr},'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize,'markerfacecolor','w'); drawnow; % acutal
%                     plot(tm,neuralTe{2,iPr},'d','color',cmapRot(iPr,:),'linewidth',lw,'markersize',msize2,'markerfacecolor','w'); drawnow; % rotated
                end
            end
            if idTpMTs
                axis tight; % 
                 xlim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                 plotHorizon(gca,0,[]); % plotIdentity(gca);
                 set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                     'ytick',-500:50:500); % 'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
            else
                yyaxis left;
%                 axis tight; 
                ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                plotIdentity(gca);
                set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                     'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)});
                 
                 yyaxis right;
                 ylim([T{1}(1)-10 T{2}(end)+10]); %              axis tight;
                 xlim([480 1200]);
%                  plotIdentity(gca);
                 cs1=statTmp{1}.beta(2); cs2=statTmp{1}.beta(1);
                cl1=statTmp{2}.beta(2); cl2=statTmp{2}.beta(1);
                 y1=(T{1}(1)-cs1)/cs2; 
                 y2=(T{1}(3)-cs1)/cs2;
                 y3=(T{1}(end)-cs1)/cs2;
                 y4=(T{2}(3)-cl1)/cl2;
                 y5=(T{2}(end)-cl1)/cl2;
                 y1=round(y1*10)/10;y2=round(y2*10)/10;y3=round(y3*10)/10;
                 y4=round(y4*10)/10;y5=round(y5*10)/10;
                 set(gca,'ticklength',[0.01 0.01],'tickdir','out',...
                     'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{y1;y2;y3;y4;y5});

            end
%             applytofig(gcf,optsExpFig);
if idSaveFig
    savefig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/3/new/d/'...
        dsName '2.fig']);
end
            
            %% fig3e2: R2(BLS) vs angle(V,V')
            dAngleDot=cosd(dAngle);
%             deg2polar=@(x,iPr) deg2rad(((iPr==1)+(iPr==nPr)*(-1))*x+90); % left for short, right for long
        hFig=figure; setFigPos(floor((iDS-1)/4)+1,rem(iDS-1,4)+1); % eye for 1st row; LeftG, LeftH, RightG, RightH
        for iPr=1:nPr
            subplot(1,nPr,iPr);
            
            tmpX=squeeze(dAngleDot(iDS,:,iPr)); % 2*(iPr-1.5)* % short in left
            tmpY=squeeze(r2(iDS,1:(end-1),iPr));
            
            % individual rand
%             polarplot(deg2polar(tmpX,iPr),tmpY,'.','color',tmpCmap{iPr,1}(1,:),'markersize',msize2);ha; 
            plot(tmpX,tmpY,'.','color',tmpCmap{iPr,1}(1,:),'markersize',msize2);ha; % individual rand
            
            % moving average
            [tmpXsort,idSort]=sort(tmpX(:));
            x2=smooth(tmpXsort(:),nWin);y2=smooth(tmpY(idSort),nWin); 
%             x2=conv(tmpXsort(:),boxWin(:),'full');y2=conv(tmpY(idSort),boxWin(:),'full'); % valid
%             polarplot(deg2polar(x2,iPr),y2,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2);
            plot(x2,y2,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2);
            %             s=plotReg(tmpX,tmpY,hFig,tmpCmap{iPr,1}(1,:)); % regression line
            %             disp([d(iDS).name ': slope=' num2str(s.beta(2),3) ',  p=' num2str(s.tstat.pval(2))]);
            plot(1,r2(iDS,end,iPr),'o','color',cmapPr(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % acutal data 0 2*(iPr-1.5)
            plot(dAngleDot(iDS,iRandEx{iPr},iPr),r2(iDS,iRandEx{iPr},iPr),'d','color',cmapRot(iPr,:),'markersize',msize,'markerfacecolor','w','linewidth',lw);ha; % *2*(iPr-1.5)
                
%             % acutal data & example
%             polarplot(deg2polar(0,iPr),r2(iDS,end,iPr),'o','color',cmapPr(iPr,:),'markersize',msize2,'markerfacecolor','w','linewidth',lw);ha; % acutal data
%             polarplot(deg2polar(dAngle(iDS,iRandEx{iPr},iPr),iPr),r2(iDS,iRandEx{iPr},iPr),'d','color',cmapRot(iPr,:),'markersize',msize2,'markerfacecolor','w','linewidth',lw);ha;
            
%             set(gca,'thetalim',[0 180],'rlim',[0 1],'Rtick',0:.25:1,'ThetaTick',[0 30 60 90 120 150 180],'tickdir','out','TickLength',[0.02 0.02],...
%                     'ThetaTickLabel',{'90','60','30','0','30','60','90'});
%             axis([0 90 0 1]); %                 axis tight;
axis([0 1 0 1]); %axis([-1 1 0 1]);
            xlabel('v\prime_t_p projected to v_t_p'); if iPr==1, ylabel('Similarity to Bayesian model'); end
            set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[-1:0.5:1],'ytick',0:0.5:1,...
                'xticklabel',[1 0.5 0 0.5 1]); % 0:30:90
            if iPr==2, set(gca,'yticklabel',[]); end
            box off;
            remAxLabel; remTickLabel;
            
        end % for iPr=1:nPr
%         applytofig(gcf,optsExpFig);
if idSaveFig
    savefig(gcf,['/Users/hansem/Dropbox (MIT)/figuresRSG2prior/3/new/d/'...
        dsName '3.fig']);
end
        
        end % if strcmp(d(iDS).name,'trajKS_prior_EL_H_supportOnly.mat')
    end % for iDS=1:nDS
    
     %% checking angle of uts/vtp across priors % H/G, Hand/Eye, Right/Left
    figure; ha; setFigPos(1,5);
    aUBoot2=reshape(aUBoot,nAnimal*nEH*nTarg,nBootAngle)'; % [1000 x 8]
    violinPlot(aUBoot2,'showMM',2); plotHorizon(gca,90,[]);
    plot(aU(:),'o');
    set(gca,'xticklabel',[]);
    xlabel('data sets'); ylabel('angle(v(short),v(long))');
    
    % check bootstrap
    figure;setFigPos(2,5);violinPlot(reshape(aTmp,2*2*2*2,nBootAngle)','showMM',2);
    set(gca,'xticklabel',[]); plotHorizon(gca,90,[]);
    xlabel('data sets (animal,effector,direction,prior)'); ylabel('angle(v(original),v(bootstrap))');
    
    
end % if ~idPlotSep


%%
function teBoot=estTeRv(teBoot,nm,nPr,nTs,idDim,T,kIC,idRand,idBLSts)

load(nm); % D optimD
% nm=d(iDS).name;

% animal-specific wM
if strcmp(nm,'_H.mat')
    load('H_RSGprior_DMFC.mat','wFit');
    wm=wFit.w_m; % 0.0508
else
    load('G_RSGprior_DMFC.mat','wFit');
    wm=wFit.w_m; % 0.0481
end

% estimate Vreadout
% dVro=load(nm); % ,'Vro');
% if ~isfield(dVro,'Vro')
    for iPr=1:nPr  % prior-specific
        idTs=(iPr-1)*nTs+[1:nTs]; % 1:5 for short 6:10 for long
        Dtmp=struct2mat(D(idTs),'data'); % [dim Maxtime 5conditions]
        mDtmp=squeeze(nanmean(diff(Dtmp(:,kIC,:),1,3),3)); % [dim x 1] mean-difference
        Vro{iPr}=normalize(mDtmp,1); % IC specific
    end % for iPr=1:nPr
%     save(nm,'Vro','-append');
% else
%     Vro=dVro.Vro; % {pr}[dim x 1]
% end

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
    Dtmp=struct2mat(D(idTs),'data'); % [PC x time x ts]
    mD{iPr}=squeeze(Dtmp(:,kIC,:)); % [PC x ts]
    pmD{iPr}=Vro{iPr}(:)'*mD{iPr}; % [1 x ts]
    
    % R2 with BLS
    tm=linspace(T{iPr}(1),T{iPr}(end),length(pmD{iPr})); tm=tm(:); % 5ts
    if idBLSts
    teBLS{iPr}=BLSts(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
    else
        teBLS{iPr}=BLS(tm,wm,[min(T{iPr}) max(T{iPr})],'uniform');
    end
    stat{iPr}=regstats(teBLS{iPr},pmD{iPr},'linear'); % rsquare
    teBoot{idRand,iPr}=cat(1,teBoot{idRand,iPr},stat{iPr}.yhat(:)');
%     r2(iDS,end,iPr)=stat{iPr}.rsquare;
%     disp(['R^2: ' num2str(stat{iPr}.rsquare,3)]);
    
%     r2(iDS,end,iPr)=stat{iPr}.rsquare;
%     disp(['R^2: ' num2str(stat{iPr}.rsquare,3)]);
    
end

%%
% 1. RNN must reflect the same task analysis pipeline as implemented on data
% 2. Same RNN - different wm - should show increased bias - we reverse engineer measurement mechanism
% 3. Different RNNs either producing or not producing bias - what is different in measurement epoch.
% 
% Other analysis:
% 
% Random projections around Set + 200 ms IC space (v1) and at the end of measurement phase (v2). Find the cone of vectors plus spurious vector sets between v1 and v2 that have a high R2. Use cross validation to determine meaningful subset or at least get closer to the cone.