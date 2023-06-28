function fig4adef

% 2019/3/5
% check d' result with sqrt transformation (idSqrt): done
% use GPFA
% no IC/production > TBT speed from GPFA
% check corr. b/t tp and sSet with color map plot + corr=f(ts); partial corr

% 2018/10/5
% debug: swithc X & Y for surf
% trying stat on dprime

% 2018/9/19
% repeat analysis for G
% corr=f(t(IC)) for supplementary figure
% variance decrease by changing uts connecting only ts1-ts2, ts4-ts5? idUtsNew
%    trajKS_161218_pr_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly_newUts
%    D.sSet1, D.sSet2 (D.uts1, D.uts2)

% 2018/7/19
% after avgAttB4PCA (e.g. trajKS_prior_EL_H_supportOnly_avgAttB4PCA)
% other sessions: (criteria:At least, 40 cells with at least 2000 trials)
%         H: 161218 161222 (161211 >55 cells with >1700 trials)
%         G: 170817 170822 170818 170823 170511 170821 170507 170506 (more trials<>more neurons)
%               H/G: sess#12 16 (8) / 8 11 9 12 5 10 2 1 for varargin
% crossVal_sess(1,12); crossVal_sess(2,1); crossVal_sess(1,16); crossVal_sess(2,8); crossVal_sess(1,8);
%          if idNotFullRecord==0, include neurons with at least 2000 trials (nTrialIncNeuron=2000)
%
% 2018/6/26
% plot histogram of projected states
%
% 2018/6/12
%
% fig3c3: correlation b/t V'X(rot) & W'X(IC)
% fig3c4: correlation map with angle(V,V') and angle(W,W')
%
% (TBD) check G, tIC
%
% 2018/6/4
% modified from crossVal_sess_anal.m
% randomly perturb readout vectors: xro
% z scoring projected states to combine across conditions (1600 trials/40
% conditions=40 trials/conditions): idCondSpec=3
% parallelized
%
% 2018/5/16
% remove cells without full recording all trials '_cellFullOnly.mat'  & ~idNotFullRecord
%
% 2018/5/14
% varying t(IC) to argue against autocorr
%
% 2018/5/3
% correcting centering for PSTH2traj.m (XX_newCenter.mat): idNewCenter
%
% 2018/4/28
% analyzing results from crossVal_sess.m
%
% 2018/4/20
% complete CV: estimating readout vector also should be done only from PC(allTrial-1)
% projection should be done independently across trials
% saving sSet (1D projection), sIC (1D projection), xro(readout vector), projection weight W for all trials
% (TBD) check corr b/t sSet & sIC, how much xro/W changes across trials
%
% 2018/4/17
% debugging: 1) didn't remove cells<nTrailExc 2) PCA was also specific to prior
% (TBD) PCA for [RS] rather than priorSupport (in avgTraj, no rotation): idIncEarly
% (TBD) use tBuffer for start & end (as done in runPSTH2traj_prior_conditionSpecific)

initRSG2prior; % singleNeuronDir

% specifying session
anNm='H'; % 'H'; % 'G'; %  H: 161218 161222 161211
sessId='161218'; % '161218'; % '170506'; % G: 170817 170506

idGPFA=1;
idSqrt=1; % 0; % 1;
if idSqrt
    fnEnd2='_newUts_sqrtSingleT';
else
    fnEnd2='';
end

%% init
% perturbation parameters
 idRandPerturb=0; % 1;
nRand=10000;
nCmap=256;
cmap=flipud(copper(nCmap)); % parula(nCmap);
lw=2;
msize=10;

% plot
cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
cmapRot=[rgb('DeepPink'); rgb('Aqua')];
msize=4; % 10; % 4; % 6; % 2;
msize2=2;
lw=1.2; % 1; % 3; % 1.5;
optsExpFig.Width=6; optsExpFig.Height=6;

% params in correlation
idCondSpec=  3; % 1; % if 2, pr/ts-specific; if 0, ignore conditions; if 3, z score projected states to combine across conditions
idRobustReg=false; % false; % true;
idRegModel=1; % for plotReg - idRegModel: 0 for no error in variable x, 1 otherwise

% params in crossVal_sess
idPCpoolSess=1; % 0; % 1; % session-specific PCA if 0
idIncEarly= 0; % idIncEarly in crossVal_sess.m
idUseScoreReadout= 1; % use score from leave-one-trial-out PSTH> if 0, use psthSess w/o cross validation to estimate readout
idAllSessReadout=0; % 1; % estimating readout vector from PC of all sessions: idAllSessReadout
idUseAvgAtt=1;

idAvgAttB4PCA=1; % 0; % 1; %%%%%%%%%%
idUseTIC=0; % 1;
idNotFullRecord= 0; % 1; % 0; % 1; %    if 1, including neurons with not fully recorded
if idNotFullRecord==1,
    nTrialIncNeuron_buffer=15; % threshold for including neurons with (total # trials - 15 trials)
    fnEnd=[];
else fnEnd='_cellFullOnly';
    nTrialIncNeuron_buffer=15; % 0; % 15; % threshold for including neurons with (total # trials - 15 trials)
end
if idAllSessReadout, fnEnd=[fnEnd '_allSessReadout']; end
nPC= 0; % 3; % if 0, use nPC(PVAF>75%)
idNewCenter=1;
idSave=1; % crossVal_sess_anal_randPerturb

if idPCpoolSess
    if idIncEarly
        if idUseScoreReadout
            d1=load(fullfile(singleTDir,'trajKS_161218_ts_pr_H_condSpec_poolSessNew_CV_avgAtt_newCenter.mat')); % D tp ts idPr idHandEye theta tIC sSet
        else
            d1=load(fullfile(singleTDir,'trajKS_161218_ts_pr_H_condSpec_poolSessNew_CV_avgAtt_newCenter_noUseScoreReadout.mat')); % D tp ts idPr idHandEye theta tIC sSet
        end
    else % idIncEarly
        if idUseAvgAtt
            if idNewCenter==1
                if idAvgAttB4PCA
                    if idNotFullRecord 
                        d1=load(fullfile(singleTDir,['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA.mat'])); % D tp ts idPr idHandEye theta tIC sSet _newCenter
                        disp(['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA.mat']);disp(size(d1.D(1).W));
                        %                                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAttB4PCA_' fnEnd '.mat']); % psthDir
                    else  %%%%% currently used
                        if idGPFA
                            d1=load(fullfile(singleTDir,['trajKS_' sessId '_pr_' anNm '_condSpec_GPFA_CV' fnEnd '.mat'])); % D tp ts idPr idHandEye theta tIC sSet _newCenter _cellFullOnly
                            disp(['trajKS_' sessId '_pr_' anNm '_condSpec_GPFA_CV' fnEnd '.mat']);
                        else
                            d1=load(fullfile(singleTDir,['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA' fnEnd fnEnd2 '.mat'])); % D tp ts idPr idHandEye theta tIC sSet _newCenter _cellFullOnly
                            disp(['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA' fnEnd fnEnd2 '.mat']);disp(size(d1.D(1).W));
                        end
                        %                                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAttB4PCA_' fnEnd '.mat']); % psthDir
                    end
                else % idAvgAttB4PCA
                    if idNotFullRecord % OLD
                        d1=load(fullfile(singleTDir,['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter.mat'])); % D tp ts idPr idHandEye theta tIC sSet _newCenter
                        disp(['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter.mat']);disp(size(d1.D(1).W));
                        %                                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAttB4PCA_' fnEnd '.mat']); % psthDir
                    else  
                        d1=load(fullfile(singleTDir,['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat'])); % D tp ts idPr idHandEye theta tIC sSet _newCenter _cellFullOnly
                        disp(['trajKS_' sessId '_pr_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat']);disp(size(d1.D(1).W));
                        %                                             fnTmp=fullfile(singleTDir,['trajKS_' fid '_' epochNm{k} '_' animalNm{iAnimal} '_condSpec_poolSessNew_CV_avgAttB4PCA_' fnEnd '.mat']); % psthDir
                    end
                end % idAvgAttB4PCA
            end
        else
        end
    end
    % tp
    if idUseTIC, fnEnd=[fnEnd '_idUseTIC']; end
    if idNewCenter==1
        if idUseScoreReadout
            if idAvgAttB4PCA
                if idNotFullRecord
                    d2=load(fullfile(singleTDir,['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA' fnEnd '.mat'])); % D tp ts idPr idHandEye theta tIC sIC
                    disp(['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA' fnEnd '.mat']);disp(size(d2.D(1).W));
                else  %%%%% currently used
                    if idGPFA
                        d2=load(fullfile(singleTDir,['trajKS_' sessId '_tp_' anNm '_condSpec_GPFA_CV' fnEnd '.mat'])); % D tp ts idPr idHandEye theta tIC sSet _newCenter _cellFullOnly
                        disp(['trajKS_' sessId '_tp_' anNm '_condSpec_GPFA_CV' fnEnd '.mat']);
                    else
                        d2=load(fullfile(singleTDir,['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA' fnEnd '.mat'])); % D tp ts idPr idHandEye theta tIC sIC _newCenter_cellFullOnly
                        disp(['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA' fnEnd '.mat']);disp(size(d2.D(1).W));
                    end
                end
            else % if idAvgAttB4PCA
                if idNotFullRecord % OLD
                    d2=load(fullfile(singleTDir,['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat'])); % D tp ts idPr idHandEye theta tIC sIC
                    disp(['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat']);disp(size(d2.D(1).W));
                else  
                    d2=load(fullfile(singleTDir,['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat'])); % D tp ts idPr idHandEye theta tIC sIC _newCenter_cellFullOnly
                    disp(['trajKS_' sessId '_t_p_' anNm '_condSpec_poolSessNew_CV_avgAtt_newCenter' fnEnd '.mat']);disp(size(d2.D(1).W));
                end
            end
        else
            d2=load(fullfile(singleTDir,'trajKS_161218_ts_t_p_H_condSpec_poolSessNew_CV_avgAtt_newCenter_noUseScoreReadout.mat')); % D tp ts idPr idHandEye theta tIC sIC
        end
    else
    end
else % idPCpoolSess
    if idIncEarly
        d1=load(fullfile(singleTDir,'trajKS_161218_ts_pr_H_condSpec_sessPCA_CV.mat')); % D tp ts idPr idHandEye theta tIC sSet
    else
        if idUseAvgAtt
            if idNewCenter==1
                d1=load(fullfile(singleTDir,'trajKS_161218_pr_H_condSpec_sessPCA_CV_avgAtt_newCenter.mat')); % D tp ts idPr idHandEye theta tIC sSet
            else
                d1=load(fullfile(singleTDir,'trajKS_161218_pr_H_condSpec_sessPCA_CV_avgAtt.mat')); % D tp ts idPr idHandEye theta tIC sSet
            end
        else
            d1=load(fullfile(singleTDir,'trajKS_161218_pr_H_condSpec_sessPCA_CV.mat')); % D tp ts idPr idHandEye theta tIC sSet
        end
    end
%     if idNewCenter==1
%         d2=load(fullfile(singleTDir,'trajKS_161218_t_p_H_condSpec_sessPCA_CV_avgAtt_newCenter.mat')); % D tp ts idPr idHandEye theta tIC sIC
%     else
%         d2=load(fullfile(singleTDir,'trajKS_161218_t_p_H_condSpec_sessPCA_CV.mat')); % D tp ts idPr idHandEye theta tIC sIC
%     end
end

if idGPFA
    % readout vectors CV
    xro1=d1.xro'; % squeeze(struct2mat(d1.D,'xro')); % [6(maxDim) x 1610 trials]
%     xro2=squeeze(struct2mat(d2.D,'xro'));
    
    idPr=d1.idPr; % 1 short 2 long
    tp=d1.tp;
    ts=d1.ts; tsId=nan(size(ts)); for iPr=1:nPr,for j=1:nTspp,tsId(T{iPr}(j)==ts & idPr==iPr)=j;end;end
    idHE=2-d1.idHandEye; % 0 hand 1 eye
    theta=d1.theta-1; % d1.theta==180; % 0 right 1 left    
else
    % readout vectors CV
    xro1=squeeze(struct2mat(d1.D,'xro')); % [6(maxDim) x 1610 trials]
    xro2=squeeze(struct2mat(d2.D,'xro'));
    
    idPr=d1.idPr; % 1 short 2 long
    tp=d1.tp;
    ts=d1.ts; tsId=nan(size(ts)); for iPr=1:nPr,for j=1:nTspp,tsId(T{iPr}(j)==ts & idPr==iPr)=j;end;end
    idHE=d1.idHandEye; % 0 hand 1 eye
    theta=d1.theta==180; % 0 right 1 left
end
idCond=(nEH*nTspp*nTarg)*(idPr-1)+(nEH*nTarg)*(tsId-1)+(1-idHE)*nTarg+theta+1; % short 480 eye right ...
% idCond=(nPr*nTspp*nTarg)*idHE+(nPr*nTspp)*theta+(idPr-1)*nTspp+(tsId-1)+1; % hand right short 480 ...

%% speed estimation
% for now, measure instantanious distance (c.f. KiNeT)
nTrial=length(d2.D);
v=nan(nTrial,1);
for iT=1:nTrial
    v(iT)=mean(sqrt(sum(diff(d2.D(iT).data,1,2).^2,1)));
end % iT=1:nTrial

%% main figure: scatterplot(X(Set)*u,X(IC)*v), histogram of states (1,3; 1,4; 2,4)
% % choose example trials: t 102
% id=(d1.ts==T{1}(end)) & (d1.idHandEye==1) & (d1.theta==180) & (d1.idPr==1) ... % eye left long 800
%     & (([d1.D.sSet]')>mean([d1.D.sSet])) & (([d2.D.sIC]')>mean([d2.D.sIC])) & d1.tp<T{1}(end); %
% find(id)
% % go back to crossVal_sess.m and find traj with leave-one-trial-out

% % fig3c3: correlation b/t V'X(rot) & W'X(IC)
%     % estimating sSet
% if nPC==0 % if 0, use nPC(PVAF>75%)
%     sSet=[d1.D.sSet];
%     sIC=[d2.D.sIC];
% else
%     sSet=sum(d1.sSet(:,1:nPC).*xro1(1:nPC,:)',2); % [trials x dim].* [trials x dim] > [trials x 1]
%     sIC=sum(d2.sIC(:,1:nPC).*xro2(1:nPC,:)',2);    
% end
sSet=[d1.sSet];
sSet=sSet(:); 
% sIC=sIC(:);

calcCorr(sSet,tp,ts,v,idCondSpec,idCond,idPr,tsId,idRobustReg,idRegModel,1); % V); last for idPlot
applytofig(gcf,optsExpFig);
% calcCorr(sSet,sIC,1,idCond,idPr,tsId,1,0,1);

return;

%% varying t(IC) to argue against autocorr (2,2)
% idCondSpec=3; % 1;
tICmax=360;
if strcmp(sessId,'170506'), tICmax=340; end;
iIClist=1:(tICmax/d2.binSize);
nDimIC=size(d2.D(1).data,1);
nTrial=length(d2.D);
sICv=nan(nTrial,length(iIClist));
rMatv=[];

for iiIC=1:length(iIClist)
    iIC=iIClist(iiIC);
    for iTrial=1:nTrial
        tmpIdPr=2-d2.D(iTrial).idShortTrial;
        
        % new readout
        tmpDtraj=d2.D(iTrial).Dcv; %%%%% from tmpDtraj: score from leave-one-trial-out
        tmpDtrajMat=struct2mat(tmpDtraj,'data'); % [#PC x time x 10prTs]
        xroTmp=squeeze(mean(diff(tmpDtrajMat(:,iIC,[1:nTspp]+(tmpIdPr-1)*nTspp),1,3),3)); % mean difference vector
        xroTmp=normalize(xroTmp,1);
        
        % state(IC)
        try
        sICv(iTrial,iiIC)=xroTmp(:)'*d2.D(iTrial).data(:,iIC);
        catch
            disp('');
        end
    end % for iTrial=1:length(D)

    sIC=sICv(:,iiIC);
%     V=v2struct;
    [rMatTmp,stdR]=calcCorr(sSet,sIC,idCondSpec,idCond,idPr,tsId,idRobustReg,1,0); % V); lasts for idRegModel, idPlot 
    rMatv=[rMatv; rMatTmp stdR]; % [length(iIClist) x #conditions]
    
end % for iiIC=1:length(iIClist)
figure; setFigPos(2,2); 
shadedErrorBar(iIClist,rMatv(:,1),rMatv(:,2),{'-','color','k','linewidth',2},1);
xlabel('time after Set (ms)'); ylabel('correlation');
set(gca,'xtick',iIClist,'xticklabel',(d2.binSize/2):d2.binSize:(tICmax+d2.binSize/2),'tickdir','out','ytick',0:0.04:0.12);
axis tight; plotHorizon(gca,0,[]);plotVertical(gca,d2.tIC/d2.binSize,[]);
applytofig(gcf,optsExpFig); box off;
% set(gca,'xtick',1+[0:3:15]);

% errorbar(iIClist,mean(rMatv,2),sem(rMatv,2));
% now errorbar : SD across bootstrap across trials but across conditions
% %     cBoot=bootstrp(1000,@corr,sSet(:),sIC(:));

% return;

%% angle between readout vectors u & v in original state space: not single neurons?
% % check projection weight (W [39neurons x nDim]) and readout vector (xro [nDim x 1])
% w1=struct2mat(d1.D,'W'); % [42 x 6(maxDim) x 1610trials]
% w2=struct2mat(d2.D,'W');
% nPCtmp=3;
% nBin=15;
% 
% % angles b/t Set vs IC readout vectors
% IC=unique(idCond); % 1:40
% angle=cell(length(IC),1); % {40 x 1}
% 
% h=figure; setFigPos(1,5);  % PC coeff
% iPlot=1;
% for iPC=1:nPCtmp
%     for jPC=1:nPCtmp
%         subplot(nPCtmp,nPCtmp,iPlot);
%         for iC=1:length(IC)
%             tmpId=IC(iC)==idCond; iTmpId=find(tmpId);
%             sw1=normalize(squeeze(w1(:,iPC,tmpId)),1); % [neurons x 1PC x trials/cond]
%             sw2=normalize(squeeze(w2(:,jPC,tmpId)),1); % [neurons x 1PC x trials/cond]
%             histogram(squeeze(nansum(sw1.*sw2,1)),nBin,'edgecolor','k','displaystyle','stairs'); ha; % trials/cond
%         end % for iC=1:length(IC)
%         iPlot=iPlot+1;
%     end % for jPC=1:nPCtmp
% end % for iPC=1:nPCtmp
% 
% h2=figure; setFigPos(2,5); % original state space
% for iC=1:length(IC)
%     tmpId=IC(iC)==idCond; iTmpId=find(tmpId);
%     tmpU=xro1(:,tmpId); % [3dim x ~40 trials]
%     tmpW1=w1(:,:,tmpId); % [42 x 6(maxDim) x ~40]
%     tmpV=xro2(:,tmpId); % [3dim x ~40 trials]
%     tmpW2=w2(:,:,tmpId); % [42 x 6(maxDim) x ~40]
%     
%     for iX=1:nnz(tmpId)
%         idTmp=~isnan(tmpU(:,iX));
%         idTmp2=~isnan(tmpV(:,iX));
%         cu=squeeze(tmpW1(:,idTmp,iX))*tmpU(idTmp,iX);
%         cv=squeeze(tmpW2(:,idTmp2,iX))*tmpV(idTmp2,iX);
%         angle{iC}=[angle{iC}; ...
%             angleVectors(cu,cv)];
% %         if angleVectors(tmpXro(:,iX),tmpXro(:,jX))>90
% %             disp([num2str(iTmpId(iX)) ' vs ' num2str(iTmpId(jX))]);
% %         end
%     end
%     histogram(angle{iC},nBin,'edgecolor','k','displaystyle','stairs'); ha; % trials/cond
% end % for iC=1:length(unique(idCond))

%% random perturbation
% if idRandPerturb
%     idCondSpec=3;
%     aRand=rand(nRand,2)*90; % (0 90)
%     rMat=nan(nRand,1);
%     
%     if exist([anNm '_' sessId '_crossVal_sess_anal_randPerturb.mat']) % 'crossVal_sess_anal_randPerturb.mat')
%         load([anNm '_' sessId '_crossVal_sess_anal_randPerturb.mat']); % 'crossVal_sess_anal_randPerturb.mat');
%     else
%         
%         delete(gcp('nocreate'));
%         parpool;
%         parfor iRand=1:nRand
%             %     tstart=tic;
%             xro1r=nan(size(xro1));
%             xro2r=nan(size(xro2));
%             for iTrial=1:size(xro1,2)
%                 xro1r(~isnan(xro1(:,iTrial)),iTrial)=randVect(xro1(~isnan(xro1(:,iTrial)),iTrial),aRand(iRand,1),1);
%                 xro2r(~isnan(xro2(:,iTrial)),iTrial)=randVect(xro2(~isnan(xro2(:,iTrial)),iTrial),aRand(iRand,2),1);
%                 %         if rem(iTrial,100)==0
%                 %             disp(iTrial);
%                 %         end
%             end % for iTrial=1:size(xro1,2)
%             %     toc(tstart) ~ 10sec
%             % estimating sSet
%             if nPC==0 % if 0, use nPC(PVAF>75%)
%                 sSet=nansum(d1.sSet.*xro1r',2); % [trials x dim].* [trials x dim] > [trials x 1]
%                 sIC=nansum(d2.sIC.*xro2r',2);
%             else
%                 sSet=sum(d1.sSet(:,1:nPC).*xro1r(1:nPC,:)',2); % [trials x dim].* [trials x dim] > [trials x 1]
%                 sIC=sum(d2.sIC(:,1:nPC).*xro2r(1:nPC,:)',2);
%             end
%             sSet=sSet(:); sIC=sIC(:);
%             
%             %     V=v2struct;
%             rMat(iRand)=calcCorr(sSet,sIC,idCondSpec,idCond,idPr,tsId,idRobustReg,1,0); % V); lasts for idRegModel, idPlot
%             % input: sSet,sIC,idCondSpec,idCond,idPr,tsId,tmpCmap,idRobustReg,condNm,xro1
%             %     toc(tstart)
%         end % for iRand=1:nRand
%         if idSave
%             save([anNm '_' sessId '_crossVal_sess_anal_randPerturb.mat'],'aRand','nRand','rMat','cmap','nCmap','msize','lw','d1','d2'); % key var: aRand, rMat
%         end
%     end
%     
%     %% plot corr=f(angle@rotation,angle@set) color-coded (1,1)
%     % load crossVal_sess_anal_randPerturb_old.mat aRand nRand rMat cmap nCmap msize lw;
%     
%     % convert into cosd (inner product)
%     aRand=cosd(aRand); % [0 90 deg]>[1 0]
%     
%     % getting color: low corr>blue
%     nRmat=(rMat-min(rMat))./(max(rMat)-min(rMat)); % [0 1]
%     cnRmat=cmap(round(nRmat*(nCmap-1)+1),:); % round([0 255]+1) > [trials x 3]
%     
%     % % 2D plot
%     % hFig=figure;ha; setFigPos(1,1);
%     % for iRand=1:nRand
%     %     plot(aRand(iRand,1),aRand(iRand,2),'.','markersize',msize,'linewidth',lw,'color',cnRmat(iRand,:)); ha; % drawnow;
%     % end % for iRand=1:nRand
%     % colormap(cmap);
%     % hBar=colorbar; hBar.Label.String='corr. coef.'; set(hBar,'Ticks',linspace(0,1,5),'TickLabels',round(linspace(min(rMat),max(rMat),5)*100)/100);
%     % set(gca,'tickdir','out','ticklength',[0.02 0.02]);
%     % axis tight; xlabel('angle from actual V(readout)@Set'); ylabel('angle from actual V(readout)@IC');
%     
%     % fig3c4: correlation map with angle(V,V') and angle(W,W')
%     % contour plot
%     hFig=figure;ha; setFigPos(1,1); colormap(cmap);
%     nBin=15; nContour=15;
%     corrMat=nan(nBin,nBin);
%     [N,xedges,yedges]=histcounts2(aRand(:,1),aRand(:,2),nBin);
%     for iX=1:(length(xedges)-1)
%         for iY=1:(length(yedges)-1)
%             tmpId=xedges(iX)<=aRand(:,1) & aRand(:,1)<xedges(iX+1) & ...
%                 yedges(iY)<=aRand(:,2) & aRand(:,2)<yedges(iY+1);
%             corrMat(iX,iY)=mean(rMat(tmpId));
%         end % iY
%     end % iX
%     tmpX=xedges(1:(end-1))+diff(xedges)/2;
%     tmpY=yedges(1:(end-1))+diff(yedges)/2;
%     % contour(tmpX,tmpY,corrMat,nContour,'linewidth',lw); % ,'ShowText','on')
%     surf(tmpY,tmpX,corrMat,'linewidth',lw); % ,'ShowText','on')
%     axis tight;
%     colormap(cmap);
%     hBar2=colorbar; hBar2.Label.String='corr. coef.';
%     set(hBar2,'Ticks',round(linspace(min(corrMat(:)),max(corrMat(:)),5)*100)/100); % ,'TickLabels',round(linspace(min(rMat),max(rMat),5)*100)/100);
%     xlabel('u_t_s^Tu_t_s`'); ylabel('v_t_p^Tv_t_p`');
%     % xlabel('angle from actual V(readout)@Set'); ylabel('angle from actual V(readout)@IC');
%     set(gca,'tickdir','out','ticklength',[0.02 0.02],'xtick',0:0.2:1,'ytick',0:0.2:1,'xticklabel',0:0.2:1,'yticklabel',0:0.2:1); % 0:20:80,'ytick',0:20:80);
%     applytofig(gcf,optsExpFig);
% end % if idRandPerturb

return;

%% inspect projection weight W & readout vector xro1
% % conclude: readout vectors are well aligned
% % check projection weight (W [39neurons x nDim]) and readout vector (xro [nDim x 1])
% % matrix similarity?
% w1=struct2mat(d1.D,'W'); % [42 x 6(maxDim) x 1610trials]
% w2=struct2mat(d2.D,'W');
% 
% % finding minimum dim
% minD1=min(sum(squeeze(~isnan(w1(1,:,:))),1));
% w1=w1(:,1:minD1,:);
% minD2=min(sum(squeeze(~isnan(w2(1,:,:))),1));
% w2=w2(:,1:minD2,:);
% 
% % angles b/t trials' readout vectors
% IC=unique(idCond); % 1:40
% angle=cell(length(IC),1); % {40 x 1}
% h=figure; setFigPos(1,5); h2=figure; setFigPos(2,5);
% for iC=1:length(IC)
%     tmpId=IC(iC)==idCond; iTmpId=find(tmpId);
%     tmpXro=xro1(1:minD1,tmpId); % [3dim x ~40 trials]
%     for iX=1:(size(tmpXro,2)-1)
%         for jX=2:(size(tmpXro,2))
%             angle{iC}=[angle{iC}; angleVectors(tmpXro(:,iX),tmpXro(:,jX))];
%             if angleVectors(tmpXro(:,iX),tmpXro(:,jX))>90
%                 disp([num2str(iTmpId(iX)) ' vs ' num2str(iTmpId(jX))]);
%             end
%         end
%     end
%     histStairs(angle{iC},25,0,h);
%     figure(h2); plot(iC,angle{iC},'.'); ha;
% end % for iC=1:length(unique(idCond))
% 
% return;

%% main: calculating correlation
function [rMat,varargout]=calcCorr(sSet,tp,ts,v,idCondSpec,idCond,idPr,tsId,idRobustReg,idRegModel,idPlot) % V
% input: sSet,sIC,idCondSpec,idCond,idPr,tsId,tmpCmap,idRobustReg,condNm,xro1

% v2struct(V);
initRSG2prior; % condNm

if idCondSpec==3 % z scoring projected values to comibne across conditions
    IC=unique(idCond); % [1:40]
    
    if idPlot
        load('pplot.mat','tmpCmap');
        
        %% plot histogram (eye right)
        idStart=1; % eye right
        idHist=1; % 0; % 1;
        hFig=figure(123);setFigPos(1,4);ha; % rotation
        
        % get bin locations (common b/t short and long)
        nBin=20; % 10;
        barWidth=0.4;
        x=linspace(min(sSet),max(sSet),nBin);
        
        % set
        tmpIC=idStart:4:length(IC); % 1 5 9 ... 37
        xmeshDense=linspace(-4,4,500);
        for iC=1:length(tmpIC) % 1 2 3 ... 10
            tmpId=idCond==tmpIC(iC);
            iPr=1+(length(tmpIC)/2<iC); % 1 for short, 2 for long
            iTs=iC-length(tmpIC)/2*(iPr-1);
            subplot(2,1,iPr); ha;
            if iTs==1 | iTs==nTspp % rem(iTs,2)~=0 % only plot ts1, ts3, ts5
                if idHist
                    hH=histogram(sSet(tmpId),x);
                    Xval    = hH.BinEdges + hH.BinWidth*0.5 +hH.BinWidth*(iTs-(nTspp+1)/2)/8; % binWidth*0.25 for short, 0.75 for long
                    Xval    = Xval(1:end-1);
                    Yval    = hH.Values;
                    delete(hH);
                    hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',tmpCmap{iPr,1}(iTs,:),'EdgeColor','none');     ha;
%                     histogram(sSet(tmpId),nBin,'DisplayStyle','stairs','edgecolor',tmpCmap{iPr,1}(iTs,:)); ha;
                else
                    [bandwidth,density,xmesh,cdf]=kde(sSet(tmpId),2^5); % ,tmpMin,tmpMax);
                    densityDense=ksdensity(sSet(tmpId),xmeshDense,'bandwidth',bandwidth/2);
                    plot(xmeshDense,densityDense,'-','linewidth',0.5,'color',tmpCmap{iPr,1}(iTs,:)); ha;
                end
            end
        end % for iC=1:length(tmpIC)
        %%%%%
        xtickTmp=-.3:.3:.3;xlimTmp=[-0.45 0.4];
        subplot(2,1,1); set(gca,'xtick',xtickTmp,'xticklabel',[],'tickdir','out','ticklength',[0.01 0.01],'xlim',xlimTmp); box off; % ,'ytick',0:0.2:0.4
        subplot(2,1,2); set(gca,'xtick',xtickTmp,'tickdir','out','ticklength',[0.01 0.01],'xlim',xlimTmp); box off; % ,'ytick',0:0.2:0.4
        
        % IC
%         x=linspace(min(sIC),max(sIC),nBin);
%         hFig=figure(124);setFigPos(2,4);ha; % IC
%         xmeshDense=linspace(-6,4,500);
%         for iC=1:length(tmpIC) % 1 2 3 ... 10
%             tmpId=idCond==tmpIC(iC);
%             iPr=1+(length(tmpIC)/2<iC); % 1 for short, 2 for long
%             iTs=iC-length(tmpIC)/2*(iPr-1);
%             subplot(2,1,iPr); ha;
%             if iTs==1 | iTs==nTspp % rem(iTs,2)~=0 % only plot ts1, ts3, ts5
%                 if idHist
%                     hH=histogram(sIC(tmpId),x);
%                     Xval    = hH.BinEdges + hH.BinWidth*0.5 +hH.BinWidth*(iTs-(nTspp+1)/2)/8; % binWidth*0.25 for short, 0.75 for long
%                     Xval    = Xval(1:end-1);
%                     Yval    = hH.Values;
%                     delete(hH);
%                     hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',tmpCmap{iPr,1}(iTs,:),'EdgeColor','none');     ha;
% % histogram(sIC(tmpId),nBin,'DisplayStyle','stairs','edgecolor',tmpCmap{iPr,1}(iTs,:)); ha;
%                 else
%                     [bandwidth,density,xmesh,cdf]=kde(sIC(tmpId),2^5); % ,tmpMin,tmpMax);
%                     densityDense=ksdensity(sIC(tmpId),xmeshDense,'bandwidth',bandwidth/2);
%                     plot(xmeshDense,densityDense,'-','linewidth',0.5,'color',tmpCmap{iPr,1}(iTs,:)); ha;
%                 end
%             end
%         end % for iC=1:length(tmpIC)
%         subplot(2,1,1); set(gca,'xtick',-8:4:4,'xticklabel',[],'tickdir','out','ticklength',[0.01 0.01],'xlim',[-8 4]); box off; % ,'ytick',0:0.2:0.4
%         subplot(2,1,2); set(gca,'xtick',-8:4:4,'tickdir','out','ticklength',[0.01 0.01],'xlim',[-8 4]); box off; % ,'ytick',0:0.2:0.4
        
    
    %% dprime vs ts: set
    % 18/10/5: do 
    hF=figure; setFigPos(1,1); ha; % combined across conditions
    mDp=zeros(nEH*nTarg,nPr,nTspp);
    sdDp=zeros(nEH*nTarg,nPr,nTspp);
    
    pDp=zeros(nEH*nTarg,nPr,nTspp-1); % paired d prime
    
    bDp=zeros(nEH*nTarg,nPr);
    bDp2=zeros(nEH*nTarg,nPr);
%     hF2=figure; setFigPos(1,2); ha; % residual combined across conditions
    
    for idStart=1:(nEH*nTarg) % eye right
        tmpIC=idStart:(nEH*nTarg):length(IC); % 1 5 9 ... 37
        
        for iC=1:length(tmpIC) % 1 2 3 ... 10
            iPr=1+(nTspp<iC); % 1 for short, 2 for long
            iTsPrM=(iPr-1)*nTspp+round(nTspp/2); % 3 or 8
            iTs=iC-nTspp*(iPr-1); % 12345
            if iC~=iTsPrM % no comparison for priorMean itself
                tmpId=idCond==tmpIC(iC);
                
                sPrM=sSet(idCond==tmpIC(iTsPrM)); % state for prior mean
%                 if mean(sSet(tmpId))>mean(sPrM) % choose which is signal/noise distribution
                    [mDp(idStart,iPr,iTs),sdDp(idStart,iPr,iTs)]=dprime(sPrM,sSet(tmpId));
%                 else
%                     [mDp(iPr,iTs),sdDp(iPr,iTs)]=dprime(sSet(tmpId),sPrM);
%                 end
            end % if iC~=3 & iC~=7
             % paired neighbors
            if iC~=length(tmpIC) % not for the last
                tmpId=idCond==tmpIC(iC);
                [pDp(idStart,iPr,iTs),sdDp(idStart,iPr,iTs)]=dprime(sSet(tmpId),sSet(idCond==tmpIC(iC+1))); % 1,4 vs 2,3
            end
        end % for iC=1:length(tmpIC) % 1 2 3 ... 10
        
%         figure; setFigPos(1,idStart); ha; % dprime vs ts for each         condition
%         for iPr=1:nPr
%             shadedErrorBar(T{iPr},squeeze(mDp(idStart,iPr,:)),squeeze(sdDp(idStart,iPr,:)),{'-','color',pplot.cmap{(iPr-1)*2+1}},1);
%             for iTs=1:nTspp
%                 plot(T{iPr}(iTs),mDp(idStart,iPr,iTs),'o','markerfacecolor','w','color',tmpCmap{iPr,1}(iTs,:));
%             end % for iTs=1:nTspp
%         end %for iPr=1:nPr
%         axis tight; 
%         plotHorizon(gca,0,[]); 
%         plotVertical(gca,[median(T{1}) median(T{2})],[]);
%         xlabel('t_s (ms)');
%         ylabel('d^\prime');
%         set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ytick',-3:0.5:3,'ticklength',[0.01 0.01],'tickdir','out');
%         applytofig(gcf,optsExpFig);
        
        for iPr=1:nPr
            figure(hF);
            plot(T{iPr},squeeze(mDp(idStart,iPr,:)),'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',0.5);
            
            % for stat
             regstat=regstats(squeeze(mDp(idStart,iPr,:)),T{iPr},'linear');
             bDp(idStart,iPr)=regstat.beta(2); % second for slope
%              figure(hF2); plot(T{iPr}(:),regstat.r(:),'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',0.5);
             regstat2=regstats(regstat.r,T{iPr},'linear');
             bDp2(idStart,iPr)=regstat2.beta(2); % second for slope
        end %for iPr=1:nPr
        
    end % for idStart=1; % eye right
    
%     % stat for paired neighbors
%     figure; setFigPos(1,2);
%     for iPr=1:nPr
%         a=squeeze(pDp(:,iPr,[1 4])); % [4 conditions x 2 pairs]
%         b=squeeze(pDp(:,iPr,[2 3])); % [4 conditions x 2 pairs]
% %         [p,h,stat]=signrank(a(:),b(:))
% %         stats=testt(a(:),b(:))
%         plot(a(:),b(:),'o','color',pplot.cmap{(iPr-1)*nPr+1});ha;
%     end
%     plotIdentity(gca);
%             a=squeeze(pDp(:,:,[1 4])); % [4 conditions x 2 pairs]
%         b=squeeze(pDp(:,:,[2 3])); % [4 conditions x 2 pairs]
%         [p,h,stat]=signrank(a(:),b(:))
%         stats=testt(a(:),b(:))
    
    % stat for slopes
%     [p,h,stat]=signrank(bDp2(:)) % H0: negative slope?
    [p,h,stat]=signrank(bDp(:,1),bDp(:,2))
    stats=testt(bDp(:,1),bDp(:,2))
    
%     dmDp=[reshape(mDp(:,:,5)-mDp(:,:,4),(nTspp-1)*nPr,1); reshape(mDp(:,:,2)-mDp(:,:,1),(nTspp-1)*nPr,1)];
%     mDpMid=[reshape(mDp(:,:,4),(nTspp-1)*nPr,1); -reshape(mDp(:,:,2),(nTspp-1)*nPr,1)];
%     [p,h,stat]=signrank(bDp(:,1),bDp(:,2))
%     stats=testt(bDp(:,1),bDp(:,2))
    
    figure(hF); % mean across conditions
    for iPr=1:nPr
        shadedErrorBar(T{iPr},squeeze(mean(mDp(:,iPr,:),1)),squeeze(sem(mDp(:,iPr,:),1)),{'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',2},1);
        for iTs=1:nTspp
            plot(T{iPr}(iTs),squeeze(mean(mDp(:,iPr,iTs),1)),'o','markerfacecolor','w','color',tmpCmap{iPr,1}(iTs,:));
        end % for iTs=1:nTspp
        
%         regstat=regstats(squeeze(mean(mDp(:,iPr,:),1)),T{iPr},'linear'); 
%         res=squeeze(mDp(:,iPr,:))-(regstat.beta(1)+regstat.beta(2)*repmat(T{iPr}(:)',size(mDp,1),1)); % [4 cond x 5 ts]
%         figure; plot(res');
%         s=regstats(reshape(res',numel(res'),1),repmat(T{iPr}(:),4,1),'linear'); disp(s.beta(2));
        
        bDp(idStart,iPr)=regstat.beta(2); % second for slope
%         figure(hF2); plot(T{iPr}(:),regstat.r(:),'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',0.5);
        
    end %for iPr=1:nPr
    axis tight;
    plotHorizon(gca,0,[]);
    plotVertical(gca,[median(T{1}) median(T{2})],[]);
    xlabel('t_s (ms)');
    ylabel('d^\prime');
    %%%%%
    set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ytick',-1.5:0.5:1.5,'ticklength',[0.01 0.01],'tickdir','out');
    applytofig(gcf,optsExpFig);
    
    
    
    end % idPlot
    
    %% plot scatter sSet vs speed, speed vs tp    
    % outlier %%%%%
%     nSD=2;
%     [~,idNoOutSet,pOutSet]=removeOutlier(sSet,nSD); disp(['p(out|set): ' num2str(pOutSet)]);
%     [~,idNoOutIC,pOutIC]=removeOutlier(tp,nSD); disp(['p(out|tp): ' num2str(pOutIC)]);
%     disp(['p(out|both): ' num2str(mean(~idNoOutSet | ~idNoOutIC))]);
%     sSet=sSet(idNoOutSet & idNoOutIC);
%     tp=tp(idNoOutSet & idNoOutIC);
%     ts=ts(idNoOutSet & idNoOutIC);
%     idCond=idCond(idNoOutSet & idNoOutIC);
    
    % partial correlation
    [pc,ppc]=partialcorr(sSet(:),v(:),ts(:));
    disp(['partial correlation(Xu,Speed): ' num2str(pc) ', ' num2str(ppc)]);
    [pc,ppc]=partialcorr(tp(:),v(:),ts(:));
    disp(['partial correlation(Speed,tp): ' num2str(pc) ', ' num2str(ppc)]);
    
    % cond-spec correlation
    cCond1=nan(length(IC),1);
    cCond2=nan(length(IC),1);
    
    % condition-specific z-scoring
    for iC=1:length(IC)
        tmpId=IC(iC)==idCond;
        
        sSet(tmpId)=zscore(sSet(tmpId));
        tp(tmpId)=zscore(tp(tmpId));
        v(tmpId)=zscore(v(tmpId));
        
        cCond1(iC)=corr(sSet(tmpId),v(tmpId));
        cCond2(iC)=corr(v(tmpId),tp(tmpId));
        
        % ROC analysis
        
    end
    figure; setFigPos(2,3); ha; % cond-spec correlation
    for iStartCond=1:4
    plot(iStartCond:4:40,cCond1(iStartCond:4:40),'o-','color',pplot.cmap{iStartCond},'markerfacecolor','w'); 
    end
    plotVertical(gca,4.5:4:40,[]); applytofig4keynote; %  set(gca,'xtick',2.5:4:38.5);
    figure; setFigPos(2,3); ha; % cond-spec correlation
    for iStartCond=1:4
    plot(iStartCond:4:40,cCond2(iStartCond:4:40),'o-','color',pplot.cmap{iStartCond},'markerfacecolor','w'); 
    end
    plotVertical(gca,4.5:4:40,[]); applytofig4keynote; %  set(gca,'xtick',2.5:4:38.5);
    
    rMat1=[];
    [r1,p1]=corr(sSet(:),v(:));
    rMat1=[rMat1;r1];
    rMat2=[];
    [r2,p2]=corr(v(:),tp(:));
    rMat2=[rMat2;r2];
%     if idRobustReg
%          [B,s]=robustfit(sSet(:),tp(:));
%          s.beta=B;
%         rMat=[rMat;s.beta(2)]; % r];
%     else % used
%         s=regstats(tp(:),sSet(:),'linear',{'beta','yhat','r','rsquare','tstat'});
%         rMat=[rMat;r];
%     end
    % plot
            % for colormap
        nBin=25; nCmap=256;
cmap=parula(nCmap); % flipud(copper(nCmap));
    if idPlot
        hFig=figure;setFigPos(1,4);ha;
        disp(['Xu vs Speed, pearsonC: ' num2str(r1,3) ', p:' num2str(p1,3)]);
        [N,xedges,yedges]=histcounts2(sSet(:),v(:),nBin);
        colormap(cmap);
        tmpX=xedges(1:(end-1)); %+diff(xedges)/2;
        tmpY=yedges(1:(end-1)); % +diff(yedges)/2;
        surf(tmpY,tmpX,N,'linewidth',lw,'EdgeAlpha',0); % ,'ShowText','on')
        axis tight;
        s=plotReg(sSet(:),v(:),hFig,zeros(1,3),idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
        colormap(cmap);
        hBar2=colorbar; hBar2.Label.String='# Trials';
        set(hBar2,'Ticks',round(linspace(min(N(:)),max(N(:)),5)*100)/100); % ,'TickLabels',round(linspace(min(rMat),max(rMat),5)*100)/100);
        xlabel('z-scored projection (zXu)'); ylabel('z-scored speed (zV)'); axis tight; plotVertical(gca,0,[]); plotHorizon(gca,0,[]);
        set(gca,'xtick',-2:2:2,'ytick',-2:2:4,'tickdir','out'); %%%%%
        applytofig(gcf,optsExpFig);
        
        hFig2=figure;ha; setFigPos(1,5);
        disp(['Speed vs tp, pearsonC: ' num2str(r2,3) ', p:' num2str(p2,3)]);
        [N,xedges,yedges]=histcounts2(v(:),tp(:),nBin);
        colormap(cmap);
        tmpX=xedges(1:(end-1)); %+diff(xedges)/2;
        tmpY=yedges(1:(end-1)); % +diff(yedges)/2;
        surf(tmpY,tmpX,N,'linewidth',lw,'EdgeAlpha',0); % ,'ShowText','on')
        axis tight;
        s=plotReg(v(:),tp(:),hFig2,zeros(1,3),idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
        colormap(cmap);
        hBar2=colorbar; hBar2.Label.String='# Trials';
        set(hBar2,'Ticks',round(linspace(min(N(:)),max(N(:)),5)*100)/100); % ,'TickLabels',round(linspace(min(rMat),max(rMat),5)*100)/100);
        xlabel('z-scored speed (zV)'); ylabel('z-scored tp'); axis tight; plotVertical(gca,0,[]); plotHorizon(gca,0,[]);
        set(gca,'xtick',-2:2:2,'ytick',-2:2:4,'tickdir','out'); %%%%%
        applytofig(gcf,optsExpFig);
        
    end
%     

% % for main figure
% % xlim([-3.5 3.1]); ylim([-3.5 4.1]); title([]); remTickLabel; remAxLabel;
% 
% %     % bootstrap
%     cBoot=bootstrp(1000,@corr,sSet(:),tp(:));
%     varargout{1}=std(cBoot);
% %     figure; setFigPos(1,2);histogram(cBoot,35,'DisplayStyle','stairs');plotVertical(gca,0,[]);
% %     plot(rMat,0,'ko','markersize',10,'markerfacecolor','w');
% %     xlabel('corr b/t projS(Set) & projS(IC)'); ylabel('# bootstrapped corr. coef.');
    
    




%% plot scatter sSet & tp
%     %(old) % plot scatter X(set)*u vs X(IC)*v
%     
%     % outlier %%%%%
%     nSD=2;
%     [~,idNoOutSet,pOutSet]=removeOutlier(sSet,nSD); disp(['p(out|set): ' num2str(pOutSet)]);
%     [~,idNoOutIC,pOutIC]=removeOutlier(tp,nSD); disp(['p(out|tp): ' num2str(pOutIC)]);
%     disp(['p(out|both): ' num2str(mean(~idNoOutSet | ~idNoOutIC))]);
%     sSet=sSet(idNoOutSet & idNoOutIC);
%     tp=tp(idNoOutSet & idNoOutIC);
%     ts=ts(idNoOutSet & idNoOutIC);
%     idCond=idCond(idNoOutSet & idNoOutIC);
%     
%     % partial correlation
%     [pc,ppc]=partialcorr(sSet(:),tp(:),ts(:));
%     disp(['partial correlation: ' num2str(pc) ', ' num2str(ppc)]);
%     
%     % cond-spec correlation
%     cCond=nan(length(IC),1);
%     
%     % condition-specific z-scoring
%     for iC=1:length(IC)
%         tmpId=IC(iC)==idCond;
%         
%         sSet(tmpId)=zscore(sSet(tmpId));
%         tp(tmpId)=zscore(tp(tmpId));
%         
%         cCond(iC)=corr(sSet(tmpId),tp(tmpId));
%         
%         % ROC analysis
%         
%     end
%     figure; setFigPos(2,3); ha; % cond-spec correlation
%     for iStartCond=1:4
%     plot(iStartCond:4:40,cCond(iStartCond:4:40),'o-','color',pplot.cmap{iStartCond},'markerfacecolor','w'); 
%     end
%     plotVertical(gca,4.5:4:40,[]); applytofig4keynote; %  set(gca,'xtick',2.5:4:38.5);
%     
%     rMat=[];
%     [r,p]=corr(sSet(:),tp(:));
%     if idRobustReg
%          [B,s]=robustfit(sSet(:),tp(:));
%          s.beta=B;
%         rMat=[rMat;s.beta(2)]; % r];
%     else % used
%         s=regstats(tp(:),sSet(:),'linear',{'beta','yhat','r','rsquare','tstat'});
%         rMat=[rMat;r];
%     end
%     % plot
%     if idPlot
%         hFig=figure;setFigPos(1,3);ha;
%         plot(sSet(:),tp(:),'k.'); axis tight;
%         s=plotReg(sSet(:),tp(:),hFig,zeros(1,3),idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
%         title([' pearsonC: ' num2str(r,3) ', p:' num2str(p,3) ', slope: ' num2str(s.beta(2),3)]);
%         xlabel('z-scored projection (Set)'); ylabel('z-scored projection (IC)'); axis tight; plotVertical(gca,0,[]); plotHorizon(gca,0,[]);
%         disp(signrank(rMat));
%         set(gca,'xtick',-2:2:2,'ytick',-2:2:4,'tickdir','out');
%         
%         % for colormap
%         nBin=25; nCmap=256;
%         [N,xedges,yedges]=histcounts2(sSet(:),tp(:),nBin);
% %         hFig=figure;ha; setFigPos(1,2);
%         cmap=parula(nCmap); % flipud(copper(nCmap));
%         colormap(cmap);
% %         imagesc(fliplr(N'));
% %         s=plotReg(sSet(:),tp(:),hFig,zeros(1,3),idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
%         
%         hFig=figure;ha; setFigPos(1,4);
%         tmpX=xedges(1:(end-1)); %+diff(xedges)/2;
%         tmpY=yedges(1:(end-1)); % +diff(yedges)/2;
%         surf(tmpY,tmpX,N,'linewidth',lw,'EdgeAlpha',0); % ,'ShowText','on')
%         axis tight;
%         s=plotReg(sSet(:),tp(:),hFig,zeros(1,3),idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
%         colormap(cmap);
%         hBar2=colorbar; hBar2.Label.String='# Trials';
%         set(hBar2,'Ticks',round(linspace(min(N(:)),max(N(:)),5)*100)/100); % ,'TickLabels',round(linspace(min(rMat),max(rMat),5)*100)/100);
%         xlabel('z-scored projection (Set)'); ylabel('z-scored projection (IC)'); axis tight; plotVertical(gca,0,[]); plotHorizon(gca,0,[]);
%         set(gca,'xtick',-2:2:2,'ytick',-2:2:4,'tickdir','out'); %%%%%
% applytofig(gcf,optsExpFig);
%         
%         
%     end
% %     
% 
% % for main figure
% % xlim([-3.5 3.1]); ylim([-3.5 4.1]); title([]); remTickLabel; remAxLabel;
% 
% %     % bootstrap
%     cBoot=bootstrp(1000,@corr,sSet(:),tp(:));
%     varargout{1}=std(cBoot);
% %     figure; setFigPos(1,2);histogram(cBoot,35,'DisplayStyle','stairs');plotVertical(gca,0,[]);
% %     plot(rMat,0,'ko','markersize',10,'markerfacecolor','w');
% %     xlabel('corr b/t projS(Set) & projS(IC)'); ylabel('# bootstrapped corr. coef.');
    
    
%% skewness
%         %% set
%     hF=figure; setFigPos(1,5); ha; % combined across conditions
%     sk=zeros(nEH*nTarg,nPr,nTspp);
%     
%     zSSet=cell(nPr*nTspp,1); % pooling across conditions
%     
%     for idStart=1:(nEH*nTarg) % eye right
%         tmpIC=idStart:(nEH*nTarg):length(IC); % 1 5 9 ... 37
%         
%         for iC=1:length(tmpIC) % 1 2 3 ... 10
%             iPr=1+(nTspp<iC); % 1 for short, 2 for long
%             iTs=iC-nTspp*(iPr-1); % 12345
%             tmpId=idCond==tmpIC(iC);
%             
%             zSSet{iC}=[zSSet{iC}; zscore(sSet(tmpId))];
%             
%             sk(idStart,iPr,iTs)=skewness(sSet(tmpId));
%             figure(hF); plot(T{iPr}(iTs),sk(idStart,iPr,iTs),'o','color',tmpCmap{iPr,1}(iTs,:));
%         end % for iC=1:length(tmpIC) % 1 2 3 ... 10
% 
%     end % for idStart=1; % eye right
%     for iPr=1:nPr
%         shadedErrorBar(T{iPr},squeeze(nanmean(sk(:,iPr,:),1)),squeeze(nanstd(sk(:,iPr,:),0,1)),{'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',2},1);
%         %             for iTs=1:nTspp
%         %                 plot(T{iPr}(iTs),mDp(idStart,iPr,iTs),'o','markerfacecolor','w','color',tmpCmap{iPr,1}(iTs,:));
%         %             end % for iTs=1:nTspp
%     end %for iPr=1:nPr
%     
%     figure; setFigPos(2,5); ha; % pooled across conditions
%     mzSSet=cellfun(@skewness,zSSet);
%     for iPr=1:nPr
%         sdzSSet=[]; % bootstrap
%         for iTs=1:nTspp
%             iTmp=(iPr-1)*nTspp+iTs;
%             sdzSSet=[sdzSSet; quantile(bootstrp(1000,@skewness,zSSet{iTmp}),[.025 .975])];
%         end
%         iTmp=(iPr-1)*nTspp+(1:nTspp);
%         shadedErrorBar(T{iPr}',mzSSet(iTmp),[mzSSet(iTmp)-sdzSSet(:,1) sdzSSet(:,2)-mzSSet(iTmp)],{'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',2},1);
%     end
%     axis tight;
%     plotHorizon(gca,0,[]);
%     xlabel('t_s (ms)');
%     ylabel('skewness pooled across conditions');
%     set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out');
%     
%     figure(hF); % not pooling
%     axis tight;
%     plotHorizon(gca,0,[]);
%     xlabel('t_s (ms)');
%     ylabel('skewness');
%     set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out');
%     applytofig(gcf,optsExpFig);
%         
%         %% IC
%         hF=figure; setFigPos(1,6); ha; % combined across conditions
%     sk=zeros(nEH*nTarg,nPr,nTspp);
%         zSSet=cell(nPr*nTspp,1); % pooling across conditions
% 
%     
%     for idStart=1:(nEH*nTarg) % eye right
%         tmpIC=idStart:(nEH*nTarg):length(IC); % 1 5 9 ... 37
%         
%         for iC=1:length(tmpIC) % 1 2 3 ... 10
%             iPr=1+(nTspp<iC); % 1 for short, 2 for long
%             iTs=iC-nTspp*(iPr-1); % 12345
%             tmpId=idCond==tmpIC(iC);
%             
%              zSSet{iC}=[zSSet{iC}; zscore(sIC(tmpId))];
%             
%             sk(idStart,iPr,iTs)=skewness(sIC(tmpId));
%             figure(hF);plot(T{iPr}(iTs),sk(idStart,iPr,iTs),'o','color',tmpCmap{iPr,1}(iTs,:));
%         end % for iC=1:length(tmpIC) % 1 2 3 ... 10
%     end % for idStart=1; % eye right
%     for iPr=1:nPr
%         shadedErrorBar(T{iPr},squeeze(nanmean(sk(:,iPr,:),1)),squeeze(nanstd(sk(:,iPr,:),0,1)),{'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',2},1);
%         %             for iTs=1:nTspp
%         %                 plot(T{iPr}(iTs),mDp(idStart,iPr,iTs),'o','markerfacecolor','w','color',tmpCmap{iPr,1}(iTs,:));
%         %             end % for iTs=1:nTspp
%     end %for iPr=1:nPr
%     
%     figure; setFigPos(2,6); ha; % pooled across conditions
%     mzSSet=cellfun(@skewness,zSSet);
%     for iPr=1:nPr
%         sdzSSet=[]; % bootstrap
%         for iTs=1:nTspp
%             iTmp=(iPr-1)*nTspp+iTs;
%             sdzSSet=[sdzSSet; quantile(bootstrp(1000,@skewness,zSSet{iTmp}),[.025 .975])];
%         end
%         iTmp=(iPr-1)*nTspp+(1:nTspp);
%         shadedErrorBar(T{iPr}',mzSSet(iTmp),[mzSSet(iTmp)-sdzSSet(:,1) sdzSSet(:,2)-mzSSet(iTmp)],{'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',2},1);
%     end
%     axis tight;
%     plotHorizon(gca,0,[]);
%     xlabel('t_s (ms)');
%     ylabel('skewness pooled across conditions');
%     set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out');
%     
%     figure(hF); % not pooling
%     axis tight;
%     plotHorizon(gca,0,[]);
%     xlabel('t_s (ms)');
%     ylabel('skewness');
%     set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out');
%     applytofig(gcf,optsExpFig);

    
    %% dprime vs ts: IC
%     hF=figure; setFigPos(2,1); ha; % combined across conditions
%     mDp=zeros(nEH*nTarg,nPr,nTspp);
%     sdDp=zeros(nEH*nTarg,nPr,nTspp);
%     
%         pDp=zeros(nEH*nTarg,nPr,nTspp-1); % paired d prime
% 
%     bDp=zeros(nEH*nTarg,nPr);
%     bDp2=zeros(nEH*nTarg,nPr);
% %     hF2=figure; setFigPos(1,2); ha; % residual combined across conditions
%         
%     for idStart=1:(nEH*nTarg) % eye right
%         tmpIC=idStart:(nEH*nTarg):length(IC); % 1 5 9 ... 37
%         
%         for iC=1:length(tmpIC) % 1 2 3 ... 10
%             iPr=1+(nTspp<iC); % 1 for short, 2 for long
%             iTsPrM=(iPr-1)*nTspp+round(nTspp/2); % 3 or 8
%             iTs=iC-nTspp*(iPr-1); % 12345
%             if iC~=iTsPrM % no comparison for priorMean itself
%                 tmpId=idCond==tmpIC(iC);
%                 
%                 sPrM=sIC(idCond==tmpIC(iTsPrM)); % state for prior mean
% %                 if mean(sIC(tmpId))>mean(sPrM) % choose which is signal/noise distribution
%                     [mDp(idStart,iPr,iTs),sdDp(idStart,iPr,iTs)]=dprime(sPrM,sIC(tmpId));
% %                 else
% %                     [mDp(iPr,iTs),sdDp(iPr,iTs)]=dprime(sIC(tmpId),sPrM);
% %                 end
%             end % if iC~=3 & iC~=7
%                          % paired neighbors
%             if iC~=length(tmpIC) % not for the last
%                 tmpId=idCond==tmpIC(iC);
%                 [pDp(idStart,iPr,iTs),sdDp(idStart,iPr,iTs)]=dprime(sIC(tmpId),sIC(idCond==tmpIC(iC+1))); % 1,4 vs 2,3
%             end
%         end % for iC=1:length(tmpIC) % 1 2 3 ... 10
%         
% %         figure; setFigPos(2,idStart); ha; % dprime vs ts for each condition
% %         for iPr=1:nPr
% %             shadedErrorBar(T{iPr},squeeze(mDp(idStart,iPr,:)),squeeze(sdDp(idStart,iPr,:)),{'-','color',pplot.cmap{(iPr-1)*2+1}},1);
% %             for iTs=1:nTspp
% %                 plot(T{iPr}(iTs),mDp(idStart,iPr,iTs),'o','markerfacecolor','w','color',tmpCmap{iPr,1}(iTs,:));
% %             end % for iTs=1:nTspp
% %         end %for iPr=1:nPr
% %         axis tight; 
% %         plotHorizon(gca,0,[]); 
% %         plotVertical(gca,[median(T{1}) median(T{2})],[]);
% %         xlabel('t_s (ms)');
% %         ylabel('d^\prime');
% %         set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ytick',-3:0.5:3,'ticklength',[0.01 0.01],'tickdir','out');
% %         applytofig(gcf,optsExpFig);
%         
%         for iPr=1:nPr
%             figure(hF);
%             plot(T{iPr},squeeze(mDp(idStart,iPr,:)),'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',0.5);
%             
%             % for stat
%              regstat=regstats(squeeze(mDp(idStart,iPr,:)),T{iPr},'linear');
%              bDp(idStart,iPr)=regstat.beta(2); % second for slope
% %              figure(hF2); plot(T{iPr}(:),regstat.r(:),'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',0.5);
%              regstat2=regstats(regstat.r,T{iPr},'linear');
%              bDp2(idStart,iPr)=regstat2.beta(2); % second for slope
% 
%         end %for iPr=1:nPr
%         
%     end % for idStart=1; % eye right
%     
%         % stat for paired neighbors
%     figure; setFigPos(2,2);
%     for iPr=1:nPr
%         a=squeeze(pDp(:,iPr,[1 4])); % [4 conditions x 2 pairs]
%         b=squeeze(pDp(:,iPr,[2 3])); % [4 conditions x 2 pairs]
% %         [p,h,stat]=signrank(a(:),b(:))
% %         stats=testt(a(:),b(:))
%         plot(a(:),b(:),'o','color',pplot.cmap{(iPr-1)*nPr+1});ha;
%     end
%     plotIdentity(gca);
%             a=squeeze(pDp(:,:,[1 4])); % [4 conditions x 2 pairs]
%         b=squeeze(pDp(:,:,[2 3])); % [4 conditions x 2 pairs]
%         [p,h,stat]=signrank(a(:),b(:))
%         stats=testt(a(:),b(:))
%         
%     % stat for slopes
% %     [p,h,stat]=signrank(bDp2(:)) % H0: negative slope?
%     [p,h,stat]=signrank(bDp(:,1),bDp(:,2))
%     stats=testt(bDp(:,1),bDp(:,2))
%     
% %     dmDp=[reshape(mDp(:,:,5)-mDp(:,:,4),(nTspp-1)*nPr,1); reshape(mDp(:,:,2)-mDp(:,:,1),(nTspp-1)*nPr,1)];
% %     mDpMid=[reshape(mDp(:,:,4),(nTspp-1)*nPr,1); -reshape(mDp(:,:,2),(nTspp-1)*nPr,1)];
% %     [p,h,stat]=signrank(dmDp,mDpMid)
% %     stats=testt(dmDp,mDpMid)
%     
%     figure(hF); % mean across conditions
%     for iPr=1:nPr
%         shadedErrorBar(T{iPr},squeeze(mean(mDp(:,iPr,:),1)),squeeze(sem(mDp(:,iPr,:),1)),{'-','color',pplot.cmap{(iPr-1)*2+1},'linewidth',2},1);
%         for iTs=1:nTspp
%             plot(T{iPr}(iTs),squeeze(mean(mDp(:,iPr,iTs),1)),'o','markerfacecolor','w','color',tmpCmap{iPr,1}(iTs,:));
%         end % for iTs=1:nTspp
%     end %for iPr=1:nPr
%     axis tight;
%     plotHorizon(gca,0,[]);
%     plotVertical(gca,[median(T{1}) median(T{2})],[]);
%     xlabel('t_s (ms)');
%     ylabel('d^\prime');
%     set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ytick',-3:0.5:3,'ticklength',[0.01 0.01],'tickdir','out');
%     applytofig(gcf,optsExpFig);



elseif idCondSpec==1
    % condition-specific
    IC=unique(idCond); % [1:40]
    angle=cell(length(IC),1); % {40 x 1}

    rMat=[];
    if idPlot
        load('pplot.mat','tmpCmap');
%     hFig=figure;setFigPos(1,1);ha;
    else
        hFig=[];
    end
    for iC=1:length(IC)
              hFig=figure;ha;
        tmpId=IC(iC)==idCond; iTmpId=find(tmpId);
        if idPlot
            cmap3=tmpCmap{unique(idPr(tmpId)),1}(unique(tsId(tmpId)),:);
            plot(sSet(tmpId),sIC(tmpId),'.','color',cmap3);
                
            s=plotReg(sSet(tmpId),sIC(tmpId),hFig,cmap3,idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
                
        end
        [r,p]=corr(sSet(tmpId),sIC(tmpId));
        if idRobustReg
            rMat=[rMat;s.beta(2)]; % r];
        else
            rMat=[rMat;r];
        end
        title([condNm{iC} ' pearsonC: ' num2str(r,3) ', p:' num2str(p,3) ', slope: ' num2str(s.beta(2),3)]);waitforbuttonpress; close;
%         if r<-0.1
%             % check outlier: conclude no outlier
% %             robustdemo(sSet(tmpId),sIC(tmpId));
%             % check readout vector
%             nDimTmp=sum(~isnan(xro1(:,find(tmpId,1,'first'))));
%             tmpXro=xro1(1:nDimTmp,tmpId); % [3dim x ~40 trials]
%             for iX=1:(size(tmpXro,2)-1)
%                 for jX=2:(size(tmpXro,2))
%                     angle{iC}=[angle{iC}; angleVectors(tmpXro(:,iX),tmpXro(:,jX))];
%                     if angleVectors(tmpXro(:,iX),tmpXro(:,jX))>90
%                         disp([num2str(iTmpId(iX)) ' vs ' num2str(iTmpId(jX))]);
%                     end
%                 end
%             end
%             h=figure; histStairs(angle{iC},25,0,h);
%         end % if r<-0.1
    end % for iC=1:length(unique(idCond))
    
    if idPlot
    xlabel('projected state (Set)'); ylabel('projected state (IC)'); axis tight; plotVertical(gca,0,[]); plotHorizon(gca,0,[]);
    figure; setFigPos(1,2);histogram(rMat,25,'DisplayStyle','stairs'); plotVertical(gca,0,[]);
    plot(mean(rMat),0,'ko','markerfacecolor','w');
%     plotVertical(gca,mean(rMat),'r');
    xlabel('corr b/t projS(Set) & projS(IC)'); ylabel('# conditions');
    end
    disp(signrank(rMat));
    
    
elseif idCondSpec==2
    % pr/ts-specific
    IC=unique(idCond);
    rMat=[];
    hFig=figure;setFigPos(2,1);ha;
    for iC=1:(nPr*nTspp) % length(IC)
        %     hFig=figure;ha;
%         tmpId=(iC==(rem(idCond-1,nPr*nTspp)+1)); % [1:10 1:10...
        tmpId=(iC==(floor((idCond-1)./(nEH*nTarg))+1)); % 1111 2222 3333 ...
        cmap3=tmpCmap{unique(idPr(tmpId)),1}(unique(tsId(tmpId)),:);
        plot(sSet(tmpId),sIC(tmpId),'.','color',cmap3);
        s=plotReg(sSet(tmpId),sIC(tmpId),hFig,cmap3,idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
        [r,p]=corr(sSet(tmpId),sIC(tmpId));
        if idRobustReg
            rMat=[rMat;s.beta(2)]; % r];
        else
            rMat=[rMat;r];
        end
        %     title([num2str(r,3) ', p:' num2str(p,3)]);
        %     waitforbuttonpress; clf;
    end % for iC=1:length(unique(idCond))
    xlabel('projected state (Set)'); ylabel('projected state (IC)'); axis tight; plotVertical(gca,0,[]); plotHorizon(gca,0,[]);
    figure; setFigPos(1,2);histogram(rMat,25,'DisplayStyle','stairs');plotVertical(gca,0,[]);
    xlabel('corr b/t projS(Set) & projS(IC)'); ylabel('# conditions');
    disp(signrank(rMat));
    
%%





else % idCondSpec
    % 2pr x 5ts x 2HE x 2RL
    [r,p]=corr(sSet(:),sIC(:)); % r=~0.09
    hFig=figure;setFigPos(2,1);ha; plot(sSet(:),sIC(:),'k.');    
    s=plotReg(sSet(:),sIC(:),hFig,zeros(1,3),idRobustReg,idRegModel); % 'beta','yhat','r','rsquare','tstat'
    title(['pearsonC: ' num2str(r,3) ', p:' num2str(p,3) ', slope: ' num2str(s.beta(2),3)]);
end

return;

%% script back up
d20=load('trajKS_161218_t_p_H_condSpec_poolSessNew_CV_avgAtt_newCenter.mat');
d1=load('trajKS_161218_pr_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly.mat');
d2=load('trajKS_161218_t_p_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly.mat');
d10=load('trajKS_161218_pr_H_condSpec_poolSessNew_CV_avgAtt_newCenter_cellFullOnly.mat');
initRSG2prior
for iPr=1:2
for iTs=1:5
tmpId=d1.ts==T{iPr}(iTs)&d1.idPr==iPr&d1.idHandEye==1&d1.theta==0;
ha; plot([d1.D(tmpId).sSet],[d10.D(tmpId).sSet],'o','color',tmpCmap{iPr}(iTs,:));
end
end
xlabel('old projection to u');ylabel('projection to u');
axis tight;
plotHorizon(gca,0,[]); plotVertical(gca,0,[]);
