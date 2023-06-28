function runReg

% new regression % 2018/8/28
% issue: prior is intrinsically confounded by ts 
% 1) use only overlap+SiL for prior (previously short800 and long800 was separate regressors, wrong)
% 2) do ts-regression for short/long separately (not using SiL)
% 3) still regressed out effector, direction
% link in glmfit: linear (previously log, default for poisson)
%
% save as regFinal.mat
% change subfunction reg.m > regFinal.m

% single-neuron peak/curvature analysis @ 2018/7/11

% re-run after correcting tInfo (12/7, 8/22, 8/23) @2018/3/27
%     log(r|t(set))=sum_i [a(i)timeAfterReady(i)] + c*idPrior + d*idEH + e*theta
%       idAtt: use attrition (use data in prior span for longer ts): if idAtt>0, tp is tricky to include (setting tp empty in reg.m)
%       use short in long (idAtt=2)
% use idOut2.m to remove behOutlier
% advantage over PCA: use all TBT info, truncate trials customized for each cells, no smoothing

% additional analysis @2018/1/23
% getting single neuorns' tuning parameters (amplitude, width, baseline) to compare b/t priors

% runReg: update @ 2017/12/13

% goal 
% 1) report % cells modulated by ts, tp, prior (by model comparison)
% 2) reject static representation scheme showing @ cells with sig. interaction b/t overlap ts x prior (or use shortInLong?)
% 3) peak of ts-mod cells: distributed for prior support
% 4) population mean: curvature during prior support?
% 5) is ramping neurons' set activity correlated with tp?
% (old) running regression on set activity to find out 1) relation to tp 2) targetedDimRedu 3) decoding PPC

%  log(r|t(set))=sum_i [a(i)ts(i)] + sum_i [b(i)tp(i)] + c*idPrior + d*idEH + e*theta

%
% spec of regressors
%   log(r|t(set)): 80ms bin # spikes before set for all prior support (inc. short in long>enabling idPrior)
%
%   - ts(i): dummy variable (1 for each ts i) = mean for each ts
%   - tp(i): zscored tp for each ts i
%   - idPrior: 0 for short, 1 for long
%   - idEH: 0 for eye, 1 for hand
%   - theta: 0 for right, 1 for left
%
% idAtt: use attrition (use data in prior span for longer ts)
% idPlot: plot raster, PSTH, fitted PSTH with GLM
%
% depends on reg.m
%
% note
% - if idAtt>0, tp is tricky to include (setting tp empty in reg.m)
% - now using redo-kilosort results: for both H/G
% - poisson regression
% - ignoring higher interactions for now

% H, 17 sessions, 34 proble penetrations, 196 single 421 multi > 617 total > 555 in PSTH mat (nTrialExc=5) > 536 in regFinal.mat (thNT(trial/cond):5, thFR(prSup):1)
% G, 12 sessions, 35 proble penetrations, 260 single 481 multi > 741 total > 717 in PSTH mat (nTrialExc=5) > 636 in regFinal.mat (thNT:5, thFR:1)

%% init
initRSG2prior;
try
cd(psthDir);
catch
    psthDir='/Users/seonminahn/Dropbox (MIT)/psthDataHigh/';
end
load pplot.mat;

if ~exist(fullfile(psthDir,'regFinal.mat')) % exist('regFinal.mat')

    % time window t(set)+[-winSize dAfterSet]
    regP.winSize=80;
    regP.stepSize=[80 100]; % short/long
    regP.dAfterSet=0; % 70;
    
    regP.idAtt=2; % 0; % if 1, use data during prior support; if 2, use even short in long
    regP.idInteract=1; % 0;
    regP.idNormalize=0;
    
    % thresholding with # trials for the longest ts (attrition)
    regP.thNT=5; % 10; % 25; % 0; % 15;% minimum 15 trials
    regP.thFR=1; % 2; % 0; % 2; % 1; % [sp/s] if trial-averaged FR <thFR, discard
    
    nAnimal= length(animalNm); % 1; % length(animalNm)
    
    idOldKilosort=0;
    if idOldKilosort
        neuDir2=[neuDir2 '_old'];
    end
    
    %% model setting - # params
    % full_prior: data(short with attrition+SiL),
    %      beta(5ts)+beta(prior)+beta(HE)+beta(RL) - 8
    % reduced_prior: data(short with attrition+SiL),
    %      beta(5ts)+beta(HE)+beta(RL) - 7
    % full_tsShort: data(short with attrition)
    %      beta(5ts)+beta(HE)+beta(RL) - 7
    % reduced_tsShort: data(short with attrition)
    %      beta0+beta(HE)+beta(RL) - 3
    % full_tsLong: data(long with attrition)
    %      beta(5ts)+beta(HE)+beta(RL) - 7 
    % reduced_tsLong: data(long with attrition)
    %      beta0+beta(HE)+beta(RL) - 3
    
    regP.nModel=6; % 'full_prior','reduced_prior','full_tsShort','reduced_tsShort','full_tsLong','reduced_tsLong'
    regP.modelDescript={'full_prior: data(short with attrition+SiL) beta(5ts)+beta(prior)+beta(HE)+beta(RL)';....
        'reduced_prior: data(short with attrition+SiL) beta(5ts)+beta(HE)+beta(RL)';...
        'full_tsShort: data(short with attrition) beta(5ts)+beta(HE)+beta(RL)';...
        'reduced_tsShort: data(short with attrition) beta0+beta(HE)+beta(RL)';...
        'full_tsLong: data(long with attrition) beta(5ts)+beta(HE)+beta(RL)';...
        'reduced_tsLong: data(long with attrition) beta0+beta(HE)+beta(RL)'};
    regP.nP=[8 7 7 3 7 3];
    
    b=cell(regP.nModel,1); % beta regression coefficient {nModel x 1} [# neurons x nReg]
    p=cell(regP.nModel,1); % p value of regression {nModel x 1} [# neurons x nReg]
    c=[]; % cell ID [session cell idSU]
    d=[]; % deviance: -2*(LL(model)-LL(saturated model)): [full_prior reduced_prior full_tsShort reduced_tsShort full_tsLong reduced_tsLong]
    nT=[]; % number of data points/trials
    nP=[]; % [full noTs noTp noPrior wInteract noEH noTarg]
    idA=[]; % animal ID
    
    regP.bInfo='beta regression coefficient {nModel x 1}[# neurons x nParam]';
    regP.pInfo='p value of regression{nModel x 1} [# neurons x nParam]';
    regP.cInfo='cell ID [session cell idSU]';
    regP.dInfo='deviance: -2*(LL(model)-LL(saturated model)): [full_prior reduced_prior full_tsShort reduced_tsShort full_tsLong reduced_tsLong]';
    regP.nTInfo='number of data points/trials';
    regP.nPInfo='[full_prior reduced_prior full_tsShort reduced_tsShort full_tsLong reduced_tsLong]';
    regP.idAInfo='animal ID';
    regP.idTsModInfo='idTsMod=dBICts<0; disp([% ts-modulation:  num2str(mean(idTsMod)*100)]);'; % % ts-modulation: 39.7611
    regP.idPrModInfo='idPrMod=dBIC(:,3)<0; disp([% prior-modulation:  num2str(mean(idPrMod)*100)])';% prior-modulation: 68.3447
%     save regFinal.mat regInfo -append;

    % load cell ID from reg.mat to remove neurons with # trial/cond<5, meanFR<1
    tmp=load('reg.mat','c'); cid=tmp.c; % [neuron x (session,id)]
    
    % misc setting
    idSUonly=0;
    eventId=4; % set
    
    idPlot=0; %1;
    iPlot=1;
    
    load pplot.mat; % tmpCmap
    
    for iAnimal=1:nAnimal % length(animalNm) % length(animalNm):(-1):1
        % session
        for i=1:length(fn{iAnimal})
            fnm=fn{iAnimal}{i}; % H_RSGprior_20161203
            disp(['===== ' fnm ' =====']);
            fid=fnm(end-5:end); % '161203'
            load(fullfile(neuDir2,fid)); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo
            
            % extract behavior
            beh=load([behDir animalNm{iAnimal} '_RSGprior_DMFC.mat']);
            
            % remove outlier
            idOut2=makeIdOut(tInfo,beh,fid);
            %         idOut2=false(size(tInfo,1),1);
            %         movingBallTrials=unique(tInfo(:,1)); % size same as idNoise
            %         idSess=str2num(fid)==beh.sessId;
            %         if strcmp(fid,'161207'),idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
            %         if nnz(idSess)~=length(movingBallTrials),disp('smth wrong for trial align'); end;
            %         iOutSess=find(beh.idOut(idSess)); disp(['outlier trials (incl. abort/wait): ' num2str(length(iOutSess)) '/' num2str(length(movingBallTrials)) ',' ...
            %             num2str(length(iOutSess)/length(movingBallTrials)*100) '%']);
            %         idValid=beh.t(idSess)>0&beh.t(idSess)<3*beh.T(idSess);
            %         disp(['outlier trials (w/o abort/wait): ' num2str(nnz(beh.idOut(idSess)&idValid)) '/' num2str(length(movingBallTrials)) ',' ...
            %             num2str(nnz(beh.idOut(idSess)&idValid)/length(movingBallTrials)*100) '%']);
            %         for iOut=1:length(iOutSess)
            %             idOut2(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
            %         end
            
            % loop through each unit
            nUnit=size(idSU,1); % both SU & MU
            for j=1:nUnit
                
                if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
                    % SU only
                else % SU & MU
                    
                    if sum(str2num(fid)==cid(:,1) & idSU(j,1)==cid(:,2))>0 % if included in the original regression (meanFR(Short/Long with attrition, Short in Long)>1)
                        
                        %                 iCn=find(str2num(fid)==cn(:,1) & idSU(j,1)==cn(:,2));
                        %                 if  iCn % only for intersting neurons
                        
                        tmpCellId=idSU(j,1); % XXXX
                        disp(['cluster#' num2str(tmpCellId) ': #sp=' num2str(nnz(tmpCellId==idClust))]);
                        % consistent with plotSDFexport
                        NEV=kilosort2NEV(sp,idClust,idSU); % sp [#spikeTiming x 1], idClust [cellId(XXXX) of # spikes x 1], idSU[cellID idSU]
                        tmpId=[idSU(j,1) idSU(j,1) 1 length(unique(tInfo(:,1)))]; % id(electrode,Unit,start trial,end trial)
                        NEV.MetaTags.DateTime=fid;
                        
                        % run actual regression
                        %                 try % # trials too small
                        statTmp=regFinal(NEV,tmpId,tInfo(~idOut2,:),eventId,idPlot,regP); % {[2(b,p) x nReg]; meanFR[10 x1]}
                        if ~isempty(statTmp)
                            
                            for iModel=1:length(statTmp.b)
                                b{iModel}=[b{iModel}; statTmp.b{iModel}(:)']; % cell
                                p{iModel}=[p{iModel}; statTmp.p{iModel}(:)']; % cell
                            end
                            
                            c=[c; str2num(fid) tmpCellId idSU(j,2)]; % last for idSU
                            d=[d; statTmp.dev(:)'];
                            nT=[nT; statTmp.nT(:)'];
                            nP=[nP; statTmp.nP(:)'];
                            idA=[idA; iAnimal];
                            %                         r=cat(3,r,statTmp{2}); % [2priors x 5 ts x #cell]
                        end
                        
                    end % if str2num(fid)==cid(:,1) & idSU(j,1)==cid(:,2)
                    
                end %                     if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
            end % j nUnit
        end % for i=1:length(fn{iAnimal})
    end % animal
    clear NEV beh sp tInfo idClust
    save('regFinal.mat'); % b p c d nT nP idA regP

else
    load('regFinal.mat');
end
%% determine neurons modulated by ts: model comparison
% check distribution of peaks (cmap), also distribution of peak times
% if peak is at extreme ts>ramping?
% if peak is in the middle>curved

disp(['H # neurons: ' num2str(nnz(c(:,1)<170000))]);
disp(['G # neurons: ' num2str(nnz(c(:,1)>170000))]);
%  H # neurons: 536
% G # neurons: 636

% % check distribution of b
% bts=b(:,1:10); % ts
% % btp=b(:,11:20); % tp
% 
% % figure; set(gcf,'position',pplot.rect1_1);imagesc(bts);colorbar;xlabel('ts'); ylabel('cell');
% % figure; set(gcf,'position',pplot.rect2_1);imagesc(btp);colorbar;xlabel('tp'); ylabel('cell');
% % figure; set(gcf,'position',pplot.rect1_2);hist(bts(:),50);xlabel('reg. coeff. for ts'); ylabel('cell*10ts');
% 
% [maxTmp,iMax]=max(bts,[],2); % sorting by peak times
% [siMax,iiMax]=sort(iMax);

% figure; set(gcf,'position',pplot.rect1_3);imagesc(bts(iiMax,:));colorbar;xlabel('ts'); ylabel('sorted cell');
% 
% figure; set(gcf,'position',pplot.rect1_4);imagesc(p<0.05);xlabel('regressors ([ts tp idEH])'); ylabel('sorted cell');

    % full_prior: data(short with attrition+SiL),
    %      beta(5ts)+beta(prior)+beta(HE)+beta(RL) - 8
    % reduced_prior: data(short with attrition+SiL),
    %      beta(5ts)+beta(HE)+beta(RL) - 7
    % full_tsShort: data(short with attrition)
    %      beta(5ts)+beta(HE)+beta(RL) - 7
    % reduced_tsShort: data(short with attrition)
    %      beta0+beta(HE)+beta(RL) - 3
    % full_tsLong: data(long with attrition)
    %      beta(5ts)+beta(HE)+beta(RL) - 7 
    % reduced_tsLong: data(long with attrition)
    %      beta0+beta(HE)+beta(RL) - 3
    
% determine neurons modulated by tp: model comparison
modelN={'full_prior','reduced_prior','full_tsShort','reduced_tsShort','full_tsLong','reduced_tsLong'};
nModel=length(modelN); % 7; % [full noTs noTp noPrior wInteract noEH noTarg]
% bic difference: bic = -2*logL + numParam*log(numObs) > % cell(ts modulated)
dBICpr=d(:,1)-d(:,2)+(nP(:,1).*log(nT(:,1))-nP(:,2).*log(nT(:,2))); % -2*LLR(full/ts-reduced)
dBICtsS=d(:,3)-d(:,4)+(nP(:,3).*log(nT(:,3))-nP(:,4).*log(nT(:,4))); % -2*LLR(full/ts-reduced)
dBICtsL=d(:,5)-d(:,6)+(nP(:,5).*log(nT(:,5))-nP(:,6).*log(nT(:,6))); % -2*LLR(full/ts-reduced)
% dBICtp=d(:,1)-d(:,3)+(nP(:,1).*log(nT)-nP(:,3).*log(nT)); % -2*LLR(full/tp-reduced)
% dBIC=bsxfun(@minus,d(:,1),d(:,2:end))+... % repmat(d(:,1),1,nModel-1)-d(:,2:end)+...
%     bsxfun(@minus,nP(:,1).*log(nT),bsxfun(@times,nP(:,2:end),log(nT))); %     (repmat(nP(:,1),1,nModel-1).*log(nT)-nP(:,2:end).*log(nT)); % -2*LLR(full/tp-reduced)
idTsShortMod=dBICtsS<0; disp(['% ts-modulation(short): ' num2str(mean(idTsShortMod)*100)]); % % ts-modulation: 39.7611
idTsLongMod=dBICtsL<0; disp(['% ts-modulation(long): ' num2str(mean(idTsLongMod)*100)]); % % ts-modulation: 39.7611
idPrMod=dBICpr<0; disp(['% prior-modulation: ' num2str(mean(idPrMod)*100)]); % % prior-modulation: 68.3447

% idEffMod=dBIC(:,5)<0; disp(['% effector-modulation: ' num2str(mean(idEffMod)*100)]);
% idDirectMod=dBIC(:,6)<0; disp(['% direction-modulation: ' num2str(mean(idDirectMod)*100)]);
% % disp(c(find(min(dBICts)==dBICts),:)); % pick up ts-modulated (curved) and tp-modulated cells > PSTH
% % idTpMod=dBICtp<0; disp(['% tp-modulation: ' num2str(mean(idTpMod)*100)]);
% % disp(c(find(min(dBICtp)==dBICtp),:));

% save regFinal.mat idTsMod idPrMod -append;

% % find representative cells
% for iM=1:(nModel-1)
%     if iM~=2 % no tp-mod        
%         if iM==4 % wInteract b/t ts & prior
%             idMod=dBIC(:,iM)>0; disp(['% better ' modelN{iM+1} ' (H/G) : ' num2str(mean(idMod)*100) '(' num2str(mean(idMod(idA==1))*100) '/' num2str(mean(idMod(idA==2))*100) ')']);
% %             disp(c(find(max(dBIC(:,iM))==dBIC(:,iM)),:)); % pick up ts-modulated (curved) and tp-modulated cells > PSTH
%             [YdBIC,iDBIC]=sort(dBIC(:,iM),1,'descend');
%             disp(c(iDBIC(1:15),:));
%         else
%             idMod=dBIC(:,iM)<0; disp(['% worse w/ ' modelN{iM+1} ' (H/G): ' num2str(mean(idMod)*100) '(' num2str(mean(idMod(idA==1))*100) '/' num2str(mean(idMod(idA==2))*100) ')']);
%             disp(c(find(min(dBIC(:,iM))==dBIC(:,iM)),:)); % pick up ts-modulated (curved) and tp-modulated cells > PSTH
%         end
%     end
% end % for iM=1:nModel
% % % worse w/ noTs (H/G): 39.7611(40.2985/39.3082)
% %       170818        3088
% % 
% % % worse w/ noPrior (H/G): 68.3447(69.5896/67.2956)
% %       161218        2022
% % 
% % % better wInteract (H/G) : 20.2218(19.0299/21.2264)
% %       170823        2111
% %       170823        1104
% %       161221        1116
% %       170818        2153
% %       170822        3137
% %       161211        1109
% %       170818        2017
% %       170823        1114
% %       170823        3167
% %       161207        2055
% %       161211        1044
% %       161221        2155
% %       161211        1048
% %       161211        1104
% %       170821        3054
% % 
% % % worse w/ noEH (H/G): 71.3311(70.5224/72.0126)
% %       161211        1104
% % 
% % % worse w/ noTarg (H/G): 57.2526(53.1716/60.6918)
% %       170822        3053

% % cells with ts-modulation is more likely to be cells with context dependence (prior)
% figure; setFigPos(2,4);
% scatterhist(dBICts,dBIC(:,3),'Style','stairs','color','k','marker','.'); % ,'k.','markersize',1,'linewidth',2,'markerfacecolor','w');
% axis tight; plotVertical(gca,0,[]); plotHorizon(gca,0,[]);
% xlabel('BIC(t_s-reduced)-BIC(full model)');
% ylabel('BIC(prior-reduced)-BIC(full model)');

% % test independence b/t ts-modulated and prior-modulated neurons
% nMat=[nnz(~idTsMod&~idPrMod) nnz(~idTsMod&idPrMod);...
%     nnz(idTsMod&~idPrMod) nnz(idTsMod&idPrMod)];
% 
%  [chi2,pval]=chi2testIndep(nMat)
% % chi2 =143.9912    pval =0

 %% pie chart showing % cells (separately for each animal)
% test independence b/t ts-modulated and prior-modulated neurons
 
 for i=1:nAnimal
     % report # cells & %
     tmp=[nnz(~idTsShortMod&~idTsLongMod&idPrMod&idA==i);... % P
            nnz(idTsShortMod&~idTsLongMod&~idPrMod&idA==i);... % TS
            nnz(~idTsShortMod&idTsLongMod&~idPrMod&idA==i);... % TL
            nnz(idTsShortMod&~idTsLongMod&idPrMod&idA==i);... % P*TS
            nnz(~idTsShortMod&idTsLongMod&idPrMod&idA==i);... % P*TL
            nnz(idTsShortMod&idTsLongMod&~idPrMod&idA==i);... % TS*TL
            nnz(idTsShortMod&idTsLongMod&idPrMod&idA==i);... % P*TS*TL
         ];    
%      disp(tmp); % # cells
%      disp(tmp./nnz(idA==i)*100); % percent
     
% %      % H
%    228
%     18
%      6
%     60
%     17
%      5
%     41
% 
%    42.5373
%     3.3582
%     1.1194
%    11.1940
%     3.1716
%     0.9328
%     7.6493
% chi2 =
%    20.1340
% pval =
%    7.2201e-06
% chi2 =
%    13.1680
% pval =
%    2.8476e-04
% chi2 =
%    21.8741
% pval =
% 2.9114e-06
%    
% %      % G
% 
%    257
%     24
%     19
%     58
%     32
%     21
%     45
% 
%    40.4088
%     3.7736
%     2.9874
%     9.1195
%     5.0314
%     3.3019
%     7.0755
% chi2 =
%     5.1677
% pval =
%     0.0230
% chi2 =
%     1.0578
% pval =
%     0.3037
% chi2 =
%     4.7142
% pval =
%     0.0299
    
% % interaction b/t ts-short & prior
%      nMatAn{i}=[nnz(~idTsShortMod&~idPrMod&idA==i) nnz(~idTsShortMod&idPrMod&idA==i);... % 2nd col: prior-dep
%          nnz(idTsShortMod&~idPrMod&idA==i) nnz(idTsShortMod&idPrMod&idA==i)];     % 2nd row: ts-dep
%      [chi2,pval]=chi2testIndep(nMatAn{i})
%      
% % interaction b/t ts-long & prior
%      nMatAn{i}=[nnz(~idTsLongMod&~idPrMod&idA==i) nnz(~idTsLongMod&idPrMod&idA==i);... % 2nd col: prior-dep
%          nnz(idTsLongMod&~idPrMod&idA==i) nnz(idTsLongMod&idPrMod&idA==i)];     % 2nd row: ts-dep
%      [chi2,pval]=chi2testIndep(nMatAn{i})
     
     %% interaction b/t (ts-short|ts-long) & prior

     idTsMod=idTsShortMod | idTsLongMod;
     nMatAn{i}=[nnz(~idTsMod&~idPrMod&idA==i) nnz(~idTsMod&idPrMod&idA==i);... % 2nd col: prior-dep
         nnz(idTsMod&~idPrMod&idA==i) nnz(idTsMod&idPrMod&idA==i)];     % 2nd row: ts-dep
     [chi2,pval]=chi2testIndep(nMatAn{i})
     
     % H
%         161   228
%     29   118
% 
%    30.0373   42.5373
%     5.4104   22.0149
     
     % chi2 =
     %    21.8741
     % pval =
     % 2.9114e-06
     
     % G
%         180   257
%     64   135
% 
%    28.3019   40.4088
%    10.0629   21.2264
     
     % chi2 =
     %     4.7142
     % pval =
     %     0.0299

          disp(nMatAn{i}); % # cells
     disp(nMatAn{i}./nnz(idA==i)*100); % percent
     
     % reshape so that both-dependent in the middle; both>ts-only>N/A>prior-only
     tmp=[nMatAn{i}(2,2);nMatAn{i}(2,1);nMatAn{i}(1,1);nMatAn{i}(1,2)];
     
     figure(i+1000); setFigPos(1,i);
     p=pie(tmp);
     
     % color & text setting
     for j=2:2:length(tmp)*2        
        p(j).String=''; % [num2str(tmp(round((j-1)/2))) '(' num2str(tmp(round((j-1)/2))/sum(tmp)*100,2) '%)']; % '';
     end
     p(1).FaceColor=[.7 .7 .7];
     p(3).FaceColor=[.55 .55 .55];
     p(5).FaceColor=[1 1 1];
     p(7).FaceColor=[.85 .85 .85];
%      
%      
% %      explode
%      
 end
 
 close all;

%% check distribution of peaks (cmap), also distribution of peak times c(idTsMod,:)
%% 3) peak of ts-mod cells: distributed for prior support

% single-neuron peak/curvature analysis @ 2018/7/11
%   condition-specific
%   measure: 1) curvature (estCurvature)
%                   2) t(peak): indirect (findStatP)
%                   3) polynomial order:  dynamics complexity (polyorder)

% TBD:
%   both H & G
% use PSTH_periSet to avoid edge effect for prior support
%   how to deal with effector & directions>avg with PSTH? # trial may differ across conditions
%   how to deal with negative peak?
%   map: only prior span (choose ts-mod neurons that peak during prior support) whole ts? but our focus not visually responsive neurons
%   corr b/t peak(short) & peak(long)

%   D=combineD('PSTH_periSet_H_full_SUMU.mat','PSTH_periSet_G_full_SUMU.mat');
%   save('PSTH_periSet_HG_full_SUMU.mat','D');
load('PSTH_periSet_HG_full_SUMU.mat'); % D(1:40)
nCell=size(D(1).data,1);

% parameters
idMed=0; % 1; % summary stat for curvature: if 1, median, 0 for mean
tBuffer=80; % set+80ms included
tmpId=idTsMod; % true(size(c,1),1); % idTsMod; % idTpMod; % idTsMod;
smthWidth=40; % 20; %%%%%
binW=10;
curv=cell(nPr+1,nEH,nTarg);
mCurv=nan(nPr+1,nEH,nTarg,nCell);
maxFR=nan(nPr+1,nEH,nTarg,nCell);

% cell selection
% full data from periReady
iPSTH=findNeuronId(D(1).id,c(tmpId,:)); % index in D.data
for i=1:length(D) % selecting only cells modulated by ts
    D(i).data=D(i).data(iPSTH,:);
    D(i).id=D(i).id(iPSTH,:);
end

% preprocessing: binning>truncate@set > meanAcrConds > smoothing
% binning 
[D2(1:length(D)).spikes] = deal(D.data);
for itrial = 1:length(D)
    D2(itrial).trialId = itrial;
%     D(itrial).epochStarts = ceil(D(itrial).epochStarts./ binW);    
end
s = getSeq(D2, binW, 'useSqrt', 0);
[D.data] = deal(s.y);

% truncate [shortestTs set+buffer] incl. SiL
for i = 1:length(D) 
    if i<=nEH*nTarg*nTspp % short
        d0=floor((T{1}(1)-tBuffer-1)/binW);
        d1=floor((D(i).epochStarts(2)+tBuffer)/binW);
        D(i).data=D(i).data(:,d0:d1);  % [480-80-1 set+80]
    else % long: no buffer in between L & SiL
        d0=floor((T{1}(1)-tBuffer-1)/binW);
        d1=floor(T{1}(end)/binW);
        D(i+nEH*nTarg*nTspp).data=D(i).data(:,d0:d1); % SiL [480-80-1 800] 
        d0=floor((T{2}(1))/binW);
        d1=floor((D(i).epochStarts(2)+tBuffer)/binW);
        D(i).data=D(i).data(:,d0:d1); % long [800-1 set+80]
    end
end

% time vectors in sec
tV{1}=((T{1}(1)-tBuffer-binW):binW:(T{1}(end)+tBuffer))./1000;
tV{2}=((T{2}(1)):binW:(T{2}(end)+tBuffer))./1000; % L
tV{3}=((T{1}(1)-tBuffer-binW):binW:T{1}(end))./1000; % SiL

% mean across ts w/ attirtion (condition-specific) > smoothing
mD=cell(nPr+1,nEH,nTarg);
D2=reshape(D,[nTarg nEH nTspp nPr+1]);
for iEH=1:nEH
    for iTarg=1:nTarg
        for iPr=1:(nPr+1)
            mD{iPr,iEH,iTarg}=...
                nanmean(struct2mat(D2(iTarg,iEH,:,iPr),'data'),3); % [neuron x time]
            % smoothing done later for t(peak)
%             mD{iPr,iEH,iTarg}=smoother(mD{iPr,iEH,iTarg},smthWidth, binW)*1000/binW;
        end
        % debug
        %             figure; plot(tV{iPr},mD{iPr,iEH,iTarg}(1:25:end,:)'); axis tight; ylabel('spike rate (/s)'); xlabel('time after ready (s)');
        
        % remove buffer
        
        for iPr=1:(nPr+1)
            %   measure: 1) curvature 2) t(peak): indirect 3) polynomial order: dynamics complexity
            for iCell=1:size(mD{iPr,iEH,iTarg},1)
                dataCurv=[tV{iPr}(:)'; mD{iPr,iEH,iTarg}(iCell,:)];
                tmpCurv=estCurvature(dataCurv);
                curv{iPr,iEH,iTarg}=[curv{iPr,iEH,iTarg}; tmpCurv];
                
                % 2) t(peak): indirect (findStatP)
                if iPr~=(nPr+1) % finding t(peak) for both SiL & L for long prior
                    if iPr==1
                        tmpData=smoother(dataCurv(2,:),smthWidth, binW)*1000/binW;
                        tmpT=findStatP(tmpData,3); % 3-point stencil
                        tmpTv=tV{iPr};
                    else % long
                        tmpData=[mD{iPr+1,iEH,iTarg}(iCell,:) mD{iPr,iEH,iTarg}(iCell,:)];
                        tmpData=smoother(tmpData,smthWidth, binW)*1000/binW;
                        tmpT=findStatP(tmpData,3); % 3-point stencil
                        tmpTv=[tV{iPr+1} tV{iPr}];
                    end
                    if ~isempty(tmpT)
                        try
                            mCurv(iPr,iEH,iTarg,iCell)=tmpTv(tmpT);
                            maxFR(iPr,iEH,iTarg,iCell)=tmpData(tmpT);
                        catch
                            disp('');
                        end
                    end
                end % if iPr~=(nPr+1)
                
                % 3) polynomial order:  dynamics complexity (polyorder)
%                 try
%                     tmpOrder=polyorder(dataCurv(1,:),dataCurv(2,:)); % maximum degree: 8 in polyorder.m
%                     if ~isempty(tmpOrder)
%                         mCurv(iPr,iEH,iTarg,iCell)=tmpOrder;
%                     end
%                 catch
%                     disp('');
%                 end
                %  1) curvature (estCurvature)
%                 if idMed
%                     mCurv(iPr,iEH,iTarg,iCell)=nanmedian(tmpCurv);
%                 else
%                     mCurv(iPr,iEH,iTarg,iCell)=nanmean(tmpCurv);
%                 end
            end
            
            % debug for curv
%             idCell=round(rand*(size(c,1)-1)+1);
%             tmpCurv=curv{iPr,iEH,iTarg}(idCell,:);
%             figure;
%             yyaxis left;
%             plot(tV{iPr}(:)',mD{iPr,iEH,iTarg}(idCell,:));ylabel('spike rate (/s)');
%             yyaxis right;
%             plot(dataCurv(1,round(length(tV{iPr})-length(tmpCurv))/2+1+[1:length(tmpCurv)]),tmpCurv(:)');ylabel('curvature');
%             title(['mean curv: ' num2str(nanmean(tmpCurv(:)),2)]);
%             xlabel('time after ready (s)');
%             axis tight;

            % debug for t(peak)
%             if iPr==2
%                 idCell=round(rand*(nnz(tmpId)-1)+1);
%                 tmpTv=[tV{iPr+1} tV{iPr}];
%                 tmpCurv=mCurv(iPr,iEH,iTarg,idCell);
%                 tmpData=smoother([mD{iPr+1,iEH,iTarg}(idCell,:) mD{iPr,iEH,iTarg}(idCell,:)],smthWidth, binW)*1000/binW;
%                 figure;ha;
%                 plot(tmpTv,tmpData);ylabel('spike rate (/s)');
%                 plot(tmpCurv,tmpData(tmpTv==tmpCurv),'ro');
%                 xlabel('time after ready (s)');
%                 axis tight;
%             end
            
        end % iPr
    end % iTarg
end % iEH

% plot
optsExpFig.Width=3.5; % for keynote
optsExpFig.Height=3.5;
idH=D(1).id(:,1)<170000;
for iAnimal=1:nAnimal
    for iEH=1:nEH
        for iTarg=1:nTarg
            figure; setFigPos(iAnimal,(iEH-1)*nTarg+iTarg); ha;
            % scatter plot
%             tmpX=squeeze(mCurv(2,iEH,iTarg,idH==(2-iAnimal)));
%             tmpY=squeeze(mCurv(1,iEH,iTarg,idH==(2-iAnimal)));
%             plot(tmpX,... % long
%                 tmpY,... % short
%                 'ro','markersize',msize,'linewidth',lw,'markerfacecolor','w');
%             tmpX2=squeeze(mCurv(2,iEH,iTarg,idH==(2-iAnimal)));
%             tmpY2=squeeze(mCurv(3,iEH,iTarg,idH==(2-iAnimal)));
%             plot(tmpX2,... % long
%                 tmpY2,... % short in long
%                 'o','markersize',msize,'linewidth',lw,'markerfacecolor','w','color',[0 0.6 0]);
%             plotIdentity(gca);
% if idMed, xlabel('median curvature (long)'); ylabel('median curvature (short/SiL)'); 
%         else    xlabel('mean curvature (long)'); ylabel('mean curvature (short/SiL)'); end
            
            % difference histogram
%             nH=75;
%             tmpX=squeeze(mCurv(2,iEH,iTarg,idH==(2-iAnimal)));
%             tmpY=squeeze(mCurv(1,iEH,iTarg,idH==(2-iAnimal)));
%             histogram(tmpY-tmpX,nH,'DisplayStyle','stairs','EdgeColor','r');ha;
%             tmpX2=squeeze(mCurv(2,iEH,iTarg,idH==(2-iAnimal)));
%             tmpY2=squeeze(mCurv(3,iEH,iTarg,idH==(2-iAnimal)));
%             histogram(tmpY2-tmpX2,nH,'DisplayStyle','stairs','EdgeColor',[0 0.6 0]);
%             plotVertical(gca,0,[]);
%             if idMed,  xlabel('median curvature (short/SiL - long)')
%             else xlabel('mean curvature (short/SiL - long)'); end
%             ylabel('# cells');
%             applytofig(gcf,optsExpFig);
            
            % for t(peak), scatterhist
            nH=50;
            tmpX=squeeze(mCurv(1,iEH,iTarg,idH==(2-iAnimal))); % short
            tmpY=squeeze(mCurv(2,iEH,iTarg,idH==(2-iAnimal)));% long
            scatterhist(tmpX,tmpY,'NBins',nH,'Style','stairs','marker','.','color','k');
            xlabel('t(peak) (short)');
            ylabel('t(peak) (long)');
           applytofig(gcf,'Width',3.5,'Height',3.5);
            figure; setFigPos(iAnimal,(iEH-1)*nTarg+iTarg); ha;
            tmpX=squeeze(mCurv(1,iEH,iTarg,idH==(2-iAnimal))); % short
            tmpY=squeeze(maxFR(1,iEH,iTarg,idH==(2-iAnimal)));
            scatterhist(tmpX,tmpY,'NBins',nH,'Style','stairs','marker','.','color','k');
            xlabel('t(peak) (short)');
            ylabel('spike rate(peak) (short)');applytofig(gcf,'Width',3.5,'Height',3.5);
            figure; setFigPos(iAnimal,(iEH-1)*nTarg+iTarg); ha;
            tmpX=squeeze(mCurv(2,iEH,iTarg,idH==(2-iAnimal))); % short
            tmpY=squeeze(maxFR(2,iEH,iTarg,idH==(2-iAnimal)));
            scatterhist(tmpX,tmpY,'NBins',nH,'Style','stairs','marker','.','color','k');
            xlabel('t(peak) (long)');
            ylabel('spike rate(peak) (long)');applytofig(gcf,'Width',3.5,'Height',3.5);
        end
    end
end

% (old code)

% mean across conditions(effector/direction) w/ attrition
% assuming # trials similar across conditions
% DS=fillNan(D(1:(nTspp*nEH*nTarg))); % 1:20
% DL=fillNan(D((nTspp*nEH*nTarg+1):end)); % 21:40
% rs=cat(3,DS.data); % [cells x time x 5ts*EH*RL)
% mrs=nanmean(rs,3); % [cells x time]
% zBtsMod=zscore(mrs,0,2);
% rl=cat(3,DL.data); % [cells x time x 5ts*EH*RL)
% mrl=nanmean(rl,3); % [cells x time]
% zBtlMod=zscore(mrl,0,2);
% 
% % smoothing (must be last step always)
% smthWidth=40; % 20; %%%%%
% mrs=smoother(mrs, smthWidth, binW)*1000/binW; % spikes/s
% mrl=smoother(mrl, smthWidth, binW)*1000/binW; % spikes/s
% % for i = 1:length(D) % old smoothing codes
% %     D(i).data = smoother(D(i).data, smthWidth, binW)*1000/binW; % spikes/s
% %     D(i).epochStarts = ceil(D(i).epochStarts./ binW);    
% %     D(i).data=D(i).data(:,1:D(i).epochStarts(2)); % truncate at set
% % end
% 
% % % sorting including whole ts>remove cells with peak<min(ts)
% % [rMaxS,iMaxS]=max(mrs,[],2);[rMinS,iMinS]=min(mrs,[],2); % sorting by peak times
% % [siMaxS,iiMaxS]=sort(iMaxS);
% % % figure;setFigPos(1,4); hist(iMaxS,50); xlabel('peak times (short)'); ylabel('# cells');
% % [rMaxL,iMaxL]=max(mrl,[],2); [rMinL,iMinL]=min(mrl,[],2);% sorting by peak times
% % [siMaxL,iiMaxL]=sort(iMaxL);
% % % figure; setFigPos(2,4);hist(iMaxL,50); xlabel('peak times (long)'); ylabel('# cells');
% % % % cmap for whole ts
% % % figure;  set(gcf,'position',pplot.rect1_1); % short data sorted by short
% % % imagesc(zBtsMod(iiMaxS,:)); colorbar; 
% % % set(gca,'xtick',[0 T{1}/binW],'xticklabel',[0 T{1}]);
% % % figure;  set(gcf,'position',pplot.rect1_2); % short data sorted by long
% % % imagesc(zBtsMod(iiMaxL,:)); colorbar; 
% % % set(gca,'xtick',[0 T{1}/binW],'xticklabel',[0 T{1}]);
% % % figure;  set(gcf,'position',pplot.rect2_2); % long data sorted by long
% % % imagesc(zBtlMod(iiMaxL,:)); colorbar; 
% % % set(gca,'xtick',[0 T{1}(1)/binW T{2}/binW],'xticklabel',[0 T{1}(1) T{2}]);
% % % figure;  set(gcf,'position',pplot.rect2_1); % long data sorted by short
% % % imagesc(zBtlMod(iiMaxS,:)); colorbar; 
% % % set(gca,'xtick',[0 T{1}(1)/binW T{2}/binW],'xticklabel',[0 T{1}(1) T{2}]);
% 
% % grab only for prior support: truncate before min(ts)
% t0s=T{1}(1)-(T{1}(2)-T{1}(1)); t1s=T{1}(end); %[ 480-80 800]
% t0L=T{2}(1)-(T{2}(2)-T{2}(1)); t1L=T{2}(end); % [800-100 1200]
% % remove cells with peak<min(ts)
% [rMaxS,iMaxS]=max(mrs,[],2);[rMinS,iMinS]=min(mrs,[],2); % sorting by peak times
% [rMaxL,iMaxL]=max(mrl,[],2); [rMinL,iMinL]=min(mrl,[],2);% sorting by peak times
% tmpIdPeak=iMaxS>t0s &iMaxL>t0s;
% mrs=mrs(tmpIdPeak,:);mrl=mrl(tmpIdPeak,:);
% % tmpIdPeakS=abs(iMaxS-median(T{1}))<(T{1}(2)-T{1}(1)); % only peak around prior mean
% % tmpIdPeakL=abs(iMaxL-median(T{2}))<(T{2}(2)-T{2}(1)); % only peak around prior mean
% % mrs=mrs(tmpIdPeakS,:);mrl=mrl(tmpIdPeakL,:);
% [rMaxS,iMaxS]=max(mrs,[],2);[rMinS,iMinS]=min(mrs,[],2); % sorting by peak times
% [rMaxL,iMaxL]=max(mrl,[],2); [rMinL,iMinL]=min(mrl,[],2);% sorting by peak times
% % sorting
% [siMaxS,iiMaxS]=sort(iMaxS);[siMaxL,iiMaxL]=sort(iMaxL);
% dS=mrs(iiMaxS,t0s:t1s); zdS=zscore(dS,0,2);msdS=minSubtract(dS,2);% [cell x 400:800]
% dSL=mrl(iiMaxS,t0s:t1s);  zdSL=zscore(dSL,0,2);msdSL=minSubtract(dSL,2);% [cell x 400:800] % for SL, use sort by short
% dL=mrl(iiMaxL,t0L:t1L); zdL=zscore(dL,0,2);msdL=minSubtract(dL,2);% [cell x 700:1200]
% 
% % figure; for i=1:size(zdS,1),plot(mrs(iiMaxS(i),:));drawnow;pause(0.1);end;
% 
%  % single neurons' tuning parameters (amplitude, width, baseline) to compare b/t priors
% % amplitude=|max FR-min FR|
% % baseline=minFR
% % width=FWHM=t(max FR)-t((max FR+min FR)/2)
% 
% bS=min(dS,[],2); aS=max(dS,[],2)-bS;
% bL=min(dL,[],2); aL=max(dL,[],2)-bL;
% bSL=min(dSL,[],2); aSL=max(dSL,[],2)-bSL;
% % baseline
% figure; setFigPos(1,1); plotScatter(bS,bL); xlabel('baseline FR (short)'); ylabel('baseline FR (long)'); set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% figure; setFigPos(1,2); plotScatter(bS,bSL); xlabel('baseline FR (short)'); ylabel('baseline FR (short in long)');set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% figure; setFigPos(1,3); plotScatter(bSL,bL); xlabel('baseline FR (short in long)'); ylabel('baseline FR (long)');set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% % amplitude
% figure; setFigPos(2,1); plotScatter(aS,aL); xlabel('amplitude FR (short)'); ylabel('amplitude FR (long)'); set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% figure; setFigPos(2,2); plotScatter(aS,aSL); xlabel('amplitude FR (short)'); ylabel('amplitude FR (short in long)');set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% figure; setFigPos(2,3); plotScatter(aSL,aL); xlabel('amplitude FR (short in long)'); ylabel('amplitude FR (long)');set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% 
% figure; histogram(aS-aL,'DisplayStyle','stairs','edgecolor','k'); ha; ytmp=get(gca,'ylim'); %errorbar(mean(aS-aL),mean(ytmp),sem(aS-aL),'color','k');
% plot(mean(aS-aL),max(ytmp)-2,'kv','markersize',10); 
% signrank(aS,aL)
% 
% % max FR
% figure; setFigPos(1,1); plotScatter(max(dS,[],2),max(dL,[],2)); xlabel('max FR (short)'); ylabel('max FR (long)'); set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% figure; setFigPos(1,2); plotScatter(max(dS,[],2),max(dSL,[],2)); xlabel('max FR (short)'); ylabel('max FR (short in long)');set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% figure; setFigPos(1,3); plotScatter(max(dSL,[],2),max(dL,[],2)); xlabel('max FR (short in long)'); ylabel('max FR (long)');set(gca,'xscale','log','yscale','log');plotIdentity(gca);
% 
% % % cmap for prior support
% idRaw=2; % 2; % 1; % 0;
% figure;  setFigPos(1,1); % short data sorted by short
% if idRaw==1,imagesc(dS); elseif idRaw==2, imagesc(msdS); else imagesc(zdS); end; colorbar; 
% set(gca,'xtick',[T{1}/binW-t0s/binW],'xticklabel',[T{1}/binW]);
% figure;   setFigPos(1,2);  % shortInLong data sorted by short
% if idRaw==1,imagesc(dSL); elseif idRaw==2, imagesc(msdSL);  else imagesc(zdSL); end; colorbar; 
% set(gca,'xtick',[T{1}/binW-t0s/binW],'xticklabel',[T{1}/binW]);
% figure;   setFigPos(1,3);  % long data sorted by short
% if idRaw==1,imagesc(dL);elseif idRaw==2, imagesc(msdL); else imagesc(zdL); end; colorbar; % hBar=colorbar; set(hBar.Label,'String','firing rate (sp/s)');
% set(gca,'xtick',[T{2}/binW-t0L/binW],'xticklabel',[T{2}/binW]);
% 
% close all;
% 
% [X{1:3}]=deal(dS,dSL,dL);
% [zX{1:3}]=deal(zdS,zdSL,zdL);

%% 4) population mean: curvature during prior support?
% % curved activity shown even in population mean but how?
% % curvature exists in z-scored: 1) ruling out higher mean FR for peak-at-priorMean neurons 2) supporting higher gain/activation relative to baseline in peak-at-priorMean neurons
% % histogram of peak times: ruling out denser distribution around priorMean
% % corr. b/t peak times & max FR:
% % broader tuning
% %
% % PSTH_periReady, periSet, periGo
% % (TBD) separately for each cond (effector/direciton)
%     % rs=cat(3,DS.data); % [cells x time x 5ts*EH*RL)
%     % mrs=nanmean(rs,3); % [cells x time]
%     % rl=cat(3,DL.data); % [cells x time x 5ts*EH*RL)
%     % mrl=nanmean(rl,3); % [cells x time]
% 
% % histogram of peak times: ruling out denser distribution around priorMean
% figure; histogram(iMaxS,50,'Normalization','count','DisplayStyle','stairs','edgecolor','r'); % ,'BinMethod','auto'
% set(gca,'xtick',[0 T{1}/binW],'xticklabel',[0 T{1}]); xlabel('Time from Ready'); ylabel('# cells');
% figure; histogram(iMaxL,50,'Normalization','count','DisplayStyle','stairs','edgecolor','b'); % 
% set(gca,'xtick',[0 T{2}/binW],'xticklabel',[0 T{2}]);    xlabel('Time from Ready'); ylabel('# cells');
% 
% % higher max/min FR for peak-at-priorMean neurons?
% tmpY=rMaxS(iiMaxS);
% figure; setFigPos(1,1);plot(siMaxS,tmpY,'ro-','markerfacecolor','w'); axis tight;set(gca,'xtick',[0 T{1}/binW],'xticklabel',[0 T{1}]); xlabel('Time from Ready');ylabel('max spike rate @ peak');
% tmpY=rMinS(iiMaxS);figure;setFigPos(1,2); plot(siMaxS,tmpY,'ro-','markerfacecolor','w'); axis tight;set(gca,'xtick',[0 T{1}/binW],'xticklabel',[0 T{1}]); xlabel('Time from Ready');ylabel('min spike rate @ peak');
% tmpY=rMaxL(iiMaxL);
% figure;setFigPos(2,1); plot(siMaxL,tmpY,'bo-','markerfacecolor','w'); axis tight;set(gca,'xtick',[0 T{2}/binW],'xticklabel',[0 T{2}]);    xlabel('Time from Ready');ylabel('max spike rate @ peak');
% tmpY=rMinL(iiMaxL);figure;setFigPos(2,2); plot(siMaxL,tmpY,'bo-','markerfacecolor','w'); axis tight;set(gca,'xtick',[0 T{2}/binW],'xticklabel',[0 T{2}]);    xlabel('Time from Ready');ylabel('min spike rate @ peak');
% % movWin?
% 
% % broader tuning around priorMean?
%     
% % figure setting
% load pplot.mat;optsExpFig.Format='png';
% optsExpFig.Height=1.2; % '2'; % 7;
% optsExpFig.Width=2.4; % '4';
% optsExpFig.FontSize='12';%         rmfield(optsExpFig,'Width');%         rmfield(optsExpFig,'Height');
% optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';%         optsExpFig.Color='cmyk';
% optsExpFig.LockAxes=1;
% optsExpFig.LineMode='scaled';
% lw=1; linestyle={'-','--'};
% 
% % raw_mean
% figure; setFigPos(1,1); % raw_mean,eye ready, 
% shadedErrorBar(t0s:t1s,mean(dS,1),sem(dS,1),{'-','color',pplot.cmap{1},'linewidth',2},1);ha;
% shadedErrorBar(t0s:t1s,mean(dSL,1),sem(dSL,1),{'-','color',pplot.cmap{2},'linewidth',2},1);ha;
% shadedErrorBar(t0L:t1L,mean(dL,1),sem(dL,1),{'-','color',pplot.cmap{3},'linewidth',2},1);ha;
% set(gca,'xtick',[t0s t1s t1L],'xticklabel',[t0s t1s t1L]);
% axis tight;plotVertical(gca,t1s,[]);
% 
% % figure; setFigPos(1,2); % raw_mean,eye set,
% % idEH=1;idEpoch=2;
% % ws=v2struct;plotRaw(ws)
% % 
% % figure; setFigPos(1,3); % raw_mean,eye go,
% % idEH=1;idEpoch=3;
% % ws=v2struct;plotRaw(ws)
% 
% % normalized_mean
% figure; setFigPos(2,1);% normalized_mean,eye ready, 
% shadedErrorBar(t0s:t1s,mean(zdS,1),sem(zdS,1),{'-','color',pplot.cmap{1},'linewidth',2},1);ha;
% shadedErrorBar(t0s:t1s,mean(zdSL,1),sem(zdSL,1),{'-','color',pplot.cmap{2},'linewidth',2},1);ha;
% shadedErrorBar(t0L:t1L,mean(zdL,1),sem(zdL,1),{'-','color',pplot.cmap{3},'linewidth',2},1);ha;
% set(gca,'xtick',[t0s t1s t1L],'xticklabel',[t0s t1s t1L]);
% axis tight;plotVertical(gca,t1s,[]);
% 
% % figure; setFigPos(2,2);% normalized_mean,eye set,
% % idEH=1;idEpoch=2;
% % ws=v2struct;plotNorm(ws)
% % 
% % figure; setFigPos(2,3);% normalized_mean,eye go,
% % idEH=1;idEpoch=3;
% % ws=v2struct;plotNorm(ws)
% 
% % % variance across neurons, not across trials
% % figure; set(gcf,'position',pplot.rect2_1);% var,eye ready, 
% % idEH=1;idEpoch=1;
% % ws=v2struct;plotVar(ws)
% % 
% % figure; set(gcf,'position',pplot.rect2_2);% var,eye set,
% % idEH=1;idEpoch=2;
% % ws=v2struct;plotVar(ws)
% % 
% % figure; set(gcf,'position',pplot.rect2_3);% var,eye go,
% % idEH=1;idEpoch=3;
% % ws=v2struct;plotVar(ws)

%% 5) is ramping neurons' set activity correlated with tp?
 
 
 
%% tDR


% % check stat b,p [nCell x 10 tp]
% load pplot.mat;
% tmpCmap=[flipud(autumn(nTspp)); winter(nTspp)];
% uniCn=unique(cn(:,3));
% for iCn=1:length(uniCn)
%     % hist of beta
%     figure; ha; set(gcf,'position',pplot.(['rect1_' num2str(iCn)]));
%     for iB=1:size(b,2)    
%         histStairs(b(cn(:,3)==uniCn(iCn),iB),15,0,gcf,tmpCmap(iB,:));
%     end % for iB=1:size(b,2)
%     
%     % p(sign.)
%     figure(56); ha; set(gcf,'position',pplot.rect2_1);
%     plot(mean(p(cn(:,3)==uniCn(iCn),:)<0.05),'color',pplot.cmap{iCn});
%     
% end % for iCn=1:length(unique(cn(:,3)))
% figure(56);
% legend('NA','ramp','curved');
% 
% idNA=(cn(:,3)==0);
% idRamp=(cn(:,3)==1);
% idCurved=(cn(:,3)==2);


% save('regFinal.mat','b','p','c','r','statTmp');
% 
%% checking proportion of prior, ramp, curved neurons
% load('regFinal.mat');
% 
% %
% %% check p(prior), p(ramp), p(curve)
% % # valid cells: 592
% % # modulation cells: 450
% % % modulation cells: 0.76014
% % i(ts)
% % # cells: 258
% % % cells: 0.43581
% % i(ts)^2
% % # cells: 124
% % % cells: 0.20946
% % prior
% % # cells: 358
% % % cells: 0.60473
% 
% % setting
% regNm={'i(ts)','i(ts)^2','prior'}; % ,'i(ts) x prior','i(ts)^2  x prior'};
% idCorrectMC=1; % controling for FDR
% 
% nReg=length(regNm);
% 
% % controling for FDR
% if idCorrectMC
%     H=[];
%     P=p;
%     thP=nan(1,nReg);
% %     P=nan(size(p));
%     for i=2:size(p,2) % across regressors        
%         [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p(:,i),0.05,'pdep','no'); % desired false discovery rate, method:assuming indep/pos dep, no report
%         thP(i-1)=crit_p;
%         H=[H h(:)];
%     end    
% else
%     P=p;
%     thP=repmat(0.05,1,nReg);
% end
% save('regFinal.mat','H','thP','-append');
% %     figure;plot(sort(p(:,2)));
% %     ha; plot([1:size(p,1)]*0.05/size(p,1),'r-');
% %     legend('sorted p value','FDR criteria'); xlabel('id for multiple comparison'); ylabel('p value');
%     
% % window
% winSize=80;
% stepSize=[80 100]; % short/long
% % attrition
% 
% % thresholding with # trials for the longest ts (attrition)
% thNT=15;% minimum 15 trials
% thFR=1; % [sp/s] if trial-averaged FR <thFR, discard
% 
% disp('check p(prior), p(ramp), p(curve)');
% 
% % total number of celss
% % 364 single 623 multi
% % -18 46;... % 12/7
% % -17 20;...
% %
% % 364+623-18-46-17-20=886
% 
% % 162 cells too low FR or too small # trials
% nValidCell=size(b,1); % 724
% disp(['# valid cells: ' num2str(nValidCell)]);
% 
% iPr=4;iRamp=2;iCurv=3;iRampPr=5;iCurvPr=6;
% 
% % task modulation
% idTaskMod=sum(P(:,2:end)<repmat(thP,size(P,1),1),2)>0;
% disp(['# modulation cells: ' num2str(nnz(idTaskMod))]); % 534
% disp(['% modulation cells: ' num2str(mean(idTaskMod))]); % 74
% 
% for iReg=2:(nReg+1) % const 1st columbn
%     disp(regNm{iReg-1});
%     
%     idSig=P(:,iReg)<thP(iReg-1); % TBD: ramp& curv separately
%     disp(['# cells: ' num2str(nnz(idSig))]);
%     disp(['% cells: ' num2str(mean(idSig))]);    
%     
% end % iReg
% 
% % check mean of curved neurons
% r3=permute(r,[3 2 1]);
% r3=r3(:,:); % [neurons x 10]
% idCurved=p(:,iCurv)<thP(iCurv-1);
% idConcave=b(:,iCurv)>0;
% figure; 
% errorbar(1:size(r3(idCurved&idConcave,1:5),2),mean(r3(idCurved&idConcave,1:5),1),sem(r3(idCurved&idConcave,1:5),1),'color','r');ha;
% errorbar(1:size(r3(idCurved&idConcave,6:10),2),mean(r3(idCurved&idConcave,6:10),1),sem(r3(idCurved&idConcave,6:10),1),'color','b');ha;
% set(gca,'xtick',1:5,'xticklabel',[]);xlabel('t(set)'); ylabel('mean FR (sp/s)');
% figure; 
% errorbar(1:size(r3(idCurved&~idConcave,1:5),2),mean(r3(idCurved&~idConcave,1:5),1),sem(r3(idCurved&~idConcave,1:5),1),'color','r');ha;
% errorbar(1:size(r3(idCurved&~idConcave,6:10),2),mean(r3(idCurved&~idConcave,6:10),1),sem(r3(idCurved&~idConcave,6:10),1),'color','b');
% set(gca,'xtick',1:5,'xticklabel',[]);xlabel('t(set)'); ylabel('mean FR (sp/s)');
% 
% 
% % % ramp: scatter for slope
% % % figure; hist(b(:,iRamp),100);
% % tmpX=normalize(([1:5]-3),2); % regressor: change into
% % tmpX(4); % 0.3162
% % 
% % figure; ha;
% % nBin=15;cmap=[.5 .5 .5];markersize=10;
% % idSigRamp=P(:,iRamp)<thP;
% % b2=b(idSigRamp,iRampPr)*tmpX(4)/(winSize/1000);
% % [n,x]=hist(b2,nBin);
% % bar(x,n,'FaceColor',cmap,'EdgeColor','none');
% % axis tight; maxYlim=max(ylim);
% % plot(mean(b2),1.1*maxYlim,'color',cmap,'markerfacecolor',cmap,'marker','v','markersize',markersize);
% % plotVertical(gca,0,[]);
% % xlabel('differential ramping slope (spike/s)'); ylabel('# neurons');
% % 
% %         set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
% %         'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
% %         'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
% %         'LineWidth', 1)   
% % 
% % % curv: scatter for curv
% % % figure; hist(b(:,iCurv),100);
% % tmpX=normalize(([1:5]-3).^2,2); % regressor: change into
% % 
% % figure; ha;
% % nBin=15;cmap=[.5 .5 .5];markersize=10;
% % idSigCurv=P(:,iCurv)<thP;
% % b2=b(idSigCurv,iCurvPr)*tmpX(4)/(winSize/1000);
% % [n,x]=hist(b2,nBin);
% % bar(x,n,'FaceColor',cmap,'EdgeColor','none');
% % axis tight; maxYlim=max(ylim);
% % plot(mean(b2),1.1*maxYlim,'color',cmap,'markerfacecolor',cmap,'marker','v','markersize',markersize);
% % plotVertical(gca,0,[]);
% % xlabel('differential curvature (spike/s)'); ylabel('# neurons');
% % 
% %         set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
% %         'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
% %         'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
% %         'LineWidth', 1)   
% 
% % colormap of mean FR: plot with log
% % figure; hist(log(r(:)),100);
% r=permute(r,[2 1 3]);
% r2=reshape(r,size(r,1)*size(r,2),size(r,3)); % [10(tsShort,tsLong) x #neuron]
% % r2=log(r2); % transform Poisson to more of normal
% %PCA
% [coeff,score]=pca(r2');
% figure; plot(coeff(:,1:3));
% figure; hist(score(:,1),50);
% % TBD: sorting & smoothing, 
% figure; imagesc(r2');
% colorbar;
% % mean: concluding pretty noisy
% figure; errorbar(1:10,mean(r2,2),std(r2,0,2));
% % check representative prior cell > code okay but issue is log transform????
% i=find(c(:,1)==161210 & c(:,2)==3011); % 104);
%  figure; plot(r2(:,i))

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
%%  backup
% idPlot=0; % 1;
%         
% nAnal=6; % 12; % 16; % 11; % 10; % treating sliding window seprately
% % nCoeff=67;
% Formula=cell(nAnal,1);
% coeff=[];coeffTmp=[];
% pval=[];pvalTmp=[];
% 
% statSlide=[];
% % iCell=1;
% % formatting regressors appropriate for regression (centering for continuous/0&1 for discrete dummy)
% 
% % FINAL
% % log(nSpikesWin100_fixOn[m200_500]) ~ 1 + idLong + idHand + idSwitch +idLong*idHand (5)
% % log(nSpikesWin100_tOn[m250_250]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
% % log(nSpikesWin100_ready[m250_480]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
% % log(nSpikesWin100_set[m480_300]) ~ 1 + catTs (8)+ idLong + tp + idHand + idLeft + catTs(800)*idLong (1) + catTs*tp (8)+ tp*idHand + tp*idHand*idLeft (24) 
% % log(nSpikesWin100_go[m300_50]) ~ 1 + t_p + idHand + idLeft + t_p*idHand + t_p*idHand*idLeft (6)
% % log(nSpikesWin100_reward[m100_100]) ~ 1 + rewarded + rewDur:rewarded (3)
% 
% % new stat: poisson regression
% %     '1)    ' prior x eyeHand 22GLM [fix fix+500]'                                  0.24361      0.06015
% %     '2)    ' prior x eyeHand 22GLM [tOn-250 tOn]'                                  0.16692     0.021053
% %     '3)    ' prior x eyeHand x target 222GLM [tOn tOn+250]'                        0.06015    0.0015038
% %     '4)    ' prior x eyeHand x target 222GLM [ready-250 ready]'                   0.046617            0
% %     '5) ' prior x eyeHand x target 222GLM   '  [ready ready+480]'                              0.26617     0.082707
% %     '6)   ' ts + ts^2 + idPrior regression ' [set-200 set]'                                  0.26917     0.096241
% %     '7) prior x eyeHand x target 222GLM [ready+480 ready+800] for overlap ts inc. all long prior?                           0.29624      0.14436       
% %     ?8) ' prior x eyeHand x target 222GLM   [set-200 set]'
% %     '9) intercept tp + eyeHand x directions regression [set set+min(tp)]'      0.8812      0.73835
% %     '10) intercept tp + eyeHand x directions regression [production-min(tp) production]';
% %     '11) intercept reward regression [rewardOnset onset+100]'                 0.99098      0.95489
% %     ' rewarded regression [rewardOnset onset+100]'                             0.1203     0.013534
% %     ' rewDur*rewarded regression [rewardOnset onset+100]'                    0.079699    0.0045113
% % stat=[];
% 
% % 12)    ' prior x eyeHand 22GLM [fix-100 fix]'  
% % 13) 100 sliding window  ' prior x eyeHand 22GLM [fix-100 fix+500]'
% % 14) 100 sliding window   prior x eyeHand x target 222GLM [tOn-250 tOn+250]' 
% % 15) 100 sliding window   prior x eyeHand x target 222GLM [ready-250 ready+480]'    
% % 16 ) 100 sliding window  [set-480 set+100]
% %         % log(nSpikes_100_set) ~ b1 + b2*t_s +  b3*idShort + b4*idEye + b5*idRight ...
% %         %                                         + b6*t_s*idShort + b7*t_s*idEye + b8*t_s*idRight ...
% %         %                                           + b9*t_s^2 + b10*t_s^2* *idShort + b11*t_s^2*idEye + b12*t_s^2*idRight
% 
% for iAnimal=1:length(animalNm)
%     
% for i=1:length(fn{iAnimal})
%     disp(['===== ' fn{iAnimal}{i} ' =====']);
%     
%     % MAT file: RSG prior
%     fid=fn{iAnimal}{i}(end-5:end);
%     if iAnimal==1
%         load([neuDir fid '_1.mat']);
%         beh=load([behDir 'H_RSGprior_DMFC.mat']);
%     else
%         load([neuDir fid '_1.mat']);
%         beh=load([behDir 'G_RSGprior_DMFC.mat']); % idOut sessId
%     end
%     % remove outlier
%     idOut2=false(size(tInfo,1),1);
%     movingBallTrials=unique(tInfo(:,1)); % size same as idNoise
%     idSess=str2num(fid)==beh.sessId;
%     if iAnimal==1 & i==5,idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
%     if nnz(idSess)~=length(movingBallTrials),disp('smth wrong for trial align'); end;
%     iOutSess=find(beh.idOut(idSess)); disp(['outlier trials (incl. abort/wait): ' num2str(length(iOutSess)) '/' num2str(length(movingBallTrials)) ',' ...
%         num2str(length(iOutSess)/length(movingBallTrials)*100) '%']);
%     idValid=beh.t(idSess)>0&beh.t(idSess)<3*beh.T(idSess);
%     disp(['outlier trials (w/o abort/wait): ' num2str(nnz(beh.idOut(idSess)&idValid)) '/' num2str(length(movingBallTrials)) ',' ...
%         num2str(nnz(beh.idOut(idSess)&idValid)/length(movingBallTrials)*100) '%']);
%     for iOut=1:length(iOutSess)
%         idOut2(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
%     end
%     CoefficientNames=[]; % =cell(nCoeff,1);
% 
%     for j=1:size(cname{iAnimal}{i},1)
%         if chkSDF{iAnimal}{i}(j)==1
% %             tic
%             disp(['ch unit t(start) t(end): ' num2str(cname{iAnimal}{i}(j,:))]);
%             
%             statTmp=plotSDFexport(NEV,cname{iAnimal}{i}(j,:),tInfo(~idOut2,:),0,idPlot); % p values cell
% %             toc
%                 % GLM: Coefficients (table: regressors x estimate/SE/tStat/pValue), ModelCriterion.AIC/BIC, Deviance, Rsquared.Ordinary/Adjusted/Deviance, ResponseName, PredictorNames, Variables, SSE/SST/SSR, plotResiduals
% 
% %             for k=1:nAnal
% %                 % for fitglm
% % %                 if j==1
% % %                 Formula{k}=statTmp{k}.Formula;
% % %                 CoefficientNames=[CoefficientNames(:); statTmp{k}.CoefficientNames(:)];
% % %                 end
% %                 coeffTmp=[coeffTmp table2array(statTmp{k}.Coefficients(:,1))']; % 1 x nCoeff
% %                 pvalTmp=[pvalTmp table2array(statTmp{k}.Coefficients(:,end))'];
% %             end
% %             coeff=[coeff;coeffTmp];coeffTmp=[];
% %             pval=[pval;pvalTmp ];pvalTmp=[];
% % %             stat=[stat cell2mat(statTmp)];
% 
% iCell=cellId{iAnimal}{i}(j);
% % % [statSlide{iCell,:}]=deal(statTmp{13:16});
% % % [statSlide{iCell,1}]=deal(statTmp{16});[statSlide{iCell,2}]=deal(statTmp{17});
% %             statSlide{iCell,1}=statTmp{13}; % 2(beta,pvalue) x 4([intercept idShort idEye idShort_Long:idEye_Hand]) x time
% %             statSlide{iCell,2}=statTmp{14}; % 2(beta,pvalue) x 7([intercept idShort idEye idRight_Left idShort_Long:idEye_Hand idEye_Hand:idRight_Left idShort_Long:idRight_Left]) x time
% %             statSlide{iCell,3}=statTmp{15};% 2(beta,pvalue) x 7([intercept idShort idEye idRight_Left idShort_Long:idEye_Hand idEye_Hand:idRight_Left idShort_Long:idRight_Left]) x time
% %             statSlide{iCell,4}=statTmp{16};% 2(beta,pvalue) x 12([intercept %  ts+t_s2+idShort+idEye+idRight+ idShort:ts+idEye:ts+idRight:ts+idShort:t_s2+idEye:t_s2+idRight:t_s2]) x time
% %             statSlide{iCell,5}=statTmp{17};
%             
%             % FINAL
%             stat{iCell,1}=statTmp{13}; 
%             stat{iCell,2}=statTmp{14}; 
%             stat{iCell,3}=statTmp{15}; 
%             stat{iCell,4}=statTmp{16}; 
%             stat{iCell,5}=statTmp{17}; 
%             stat{iCell,6}=statTmp{18}; 
% % log(nSpikesWin100_fixOn[m200_500]) ~ 1 + idLong + idHand + idSwitch +idLong*idHand (5)
% % log(nSpikesWin100_tOn[m250_250]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
% % log(nSpikesWin100_ready[m250_480]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
% % log(nSpikesWin100_set[m480_300]) ~ 1 + catTs (8)+ idLong + tp + idHand + idLeft + catTs*idLong (8) + catTs*tp (8)+ tp*idHand + tp*idHand*idLeft (31) 
% % log(nSpikesWin100_go[m300_50]) ~ 1 + t_p + idHand + idLeft + t_p*idHand + t_p*idHand*idLeft (6)
% % log(nSpikesWin100_reward[m100_100]) ~ 1 + rewarded + rewDur:rewarded (3)
% 
% % iCell=iCell+1;
%         end
%         if idPlot
%         set(gcf,'PaperPositionMode','auto');
%         saveas(gcf,[neuDir fid '_' num2str(cname{iAnimal}{i}(j,1)) '_' num2str(cname{iAnimal}{i}(j,2)) '.png']);
%         end
% %         keydown=waitforbuttonpress; %chkSDF{i}=[chkSDF{i}; keydown];
% %         close all;
% % exportfig(gcf,[neuDir fid '_' num2str(cname{i}(j,1)) '_' num2str(cname{i}(j,2)) '.png'],'color','rgb','Height',13,'Width',11,'format','png');
% 
% close all;
% %         save([neuDir 'chkSDF.mat'],'chkSDF');
%     end % for j=1:size(cname{iAnimal}{i},1)
%     
% %     statName={'1) prior 22ANOVA [fix fix+500]';
% %     ' eyeHand 22ANOVA [fix fix+500]';
% %     ' prior x eyeHand 22ANOVA [fix fix+500]';
% %     '2) prior 22ANOVA [tOn-250 tOn]';
% %     ' eyeHand 22ANOVA [tOn-250 tOn]';
% %     ' prior x eyeHand 22ANOVA [tOn-250 tOn]';
% %     '3) prior 222ANOVA [tOn tOn+250]';
% %     ' eyeHand 222ANOVA [tOn tOn+250]';
% %     ' target 222ANOVA [tOn tOn+250]';
% %     ' prior x eyeHand 222ANOVA [tOn tOn+250]';
% %     ' prior x target 222ANOVA [tOn tOn+250]';
% %     ' eyeHand x target 222ANOVA [tOn tOn+250]';
% %     ' prior x eyeHand x target 222ANOVA [tOn tOn+250]';
% %     '4) prior 222ANOVA [ready-250 ready]';
% %     ' eyeHand 222ANOVA [ready-250 ready]';
% %     ' target 222ANOVA [ready-250 ready]';
% %     ' prior x eyeHand 222ANOVA [ready-250 ready]';
% %     ' prior x target 222ANOVA [ready-250 ready]';
% %     ' eyeHand x target 222ANOVA [ready-250 ready]';
% %     ' prior x eyeHand x target 222ANOVA [ready-250 ready]';
% %     '5) prior 22ANOVA [ready ready+400]';
% %     ' eyeHand 22ANOVA [ready ready+400]';
% %     ' prior x eyeHand 22ANOVA [ready ready+400]';
% %     '6) prior 22ANOVA [set-400 set]';
% %     ' eyeHand 22ANOVA [set-400 set]';
% %     ' prior x eyeHand 22ANOVA [set-400 set]';
% %     '7) prior 22ANOVA [set-320 set] for overlap ts';
% %     ' eyeHand 22ANOVA [set-320 set] for overlap ts';
% %     ' prior x eyeHand 22ANOVA [set-320 set] for overlap ts';
% %     '8) intercept tp + eyeHand x directions regression [set set+min(tp)]';
% %     ' tp regression [set set+min(tp)]';
% %     ' eyeHand regression [set set+min(tp)]';
% %     ' directions regression [set set+min(tp)]';
% %     ' tp x eyeHand regression [set set+min(tp)]';
% %     ' tp x directions regression [set set+min(tp)]';
% %     ' eyeHand x directions regression [set set+min(tp)]';
% %     ' tp x eyeHand x directions regression [set set+min(tp)]';
% %     '9) intercept tp + eyeHand x directions regression [production-min(tp) production]';
% %     ' tp regression [production-min(tp) production]';
% %     ' eyeHand regression [production-min(tp) production]';
% %     ' directions regression [production-min(tp) production]';
% %     ' tp x eyeHand regression [production-min(tp) production]';
% %     ' tp x directions regression [production-min(tp) production]';
% %     ' eyeHand x directions regression [production-min(tp) production]';
% %     ' tp x eyeHand x directions regression [production-min(tp) production]';
% %     '10) intercept reward regression [rewardOnset onset+100]';
% %     ' rewarded regression [rewardOnset onset+100]';
% %     ' rewDur*rewarded regression [rewardOnset onset+100]';};
% %     save([neuDir 'stat.mat'],'stat','statName'); 
%     
% %     Generalized Linear regression model:
% %     log(numSpikes) ~ 1 + idShort*idEye
% %     Distribution = Poisson
% % 
% % Estimated Coefficients:
% %                                Estimate      SE        tStat       pValue  
% %                                ________    _______    _______    __________
% % 
% %     (Intercept)                 -1.2879      0.125    -10.303    6.8415e-25
% %     idShort_Long                 0.5025    0.16428     3.0588     0.0022225
% %     idEye_Hand                  0.37914     0.1583     2.3951      0.016617
% %     idShort_Long:idEye_Hand    -0.35537    0.21156    -1.6798      0.092998
% % 
% % 
% % 945 observations, 941 error degrees of freedom
% % Dispersion: 1
% % Chi^2-statistic vs. constant model: 14.2, p-value = 0.00262
%     
% end % for i=1:length(fname)
% end % for iAnimal=1:length(animalNm)
% 
%  save([neuDir 'singleNeuronStat_interaction_final.mat'],'stat','nAnal','statTmp');
%  
% % iAnalPerCoeff=[];n=[];for i=1:nAnal % length(statTmp),
% %     n=n+statTmp{i}.NumCoefficients;iAnalPerCoeff=[iAnalPerCoeff;repmat(i,statTmp{i}.NumCoefficients,1)];end
% % nCoeff=length(CoefficientNames);
% %  save([neuDir 'singleNeuronStat_interaction_final.mat'],'statSlide','CoefficientNames','Formula','coeff','pval','nCoeff','nAnal','iAnalPerCoeff','statTmp'); % singleNeuronStat_interaction_noSlide
 
 %%
%  chkRegSlide
 
 %% old singleNeuronStat_interaction_final.mat
% + neurons calculated from only among neurons with significant task modulation in any period? not from all
% - use glmfit (=fitglm confirmed) 
% - sliding window (windowSize=100 'winSize=100;' in plotSDFexport, small enough to separate different set in short prior)
% - categorical ts (catTs) to accomodate complex tuning
% - centering for continuous variables: ts, tp,rewDur
% - responseVar: # spikes for Poisson regression
% 
% log(nSpikesWin100_fixOn[m200_500]) ~ 1 + idLong + idHand + idSwitch +idLong*idHand (5)
% log(nSpikesWin100_tOn[m250_250]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
% log(nSpikesWin100_ready[m250_480]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
% log(nSpikesWin100_set[m480_300]) ~ 1 + catTs (8)+ idLong + tp + idHand + idLeft + catTs(800)*idLong (1) + catTs*tp (8)+ tp*idHand + tp*idHand*idLeft (24) 
% log(nSpikesWin100_go[m300_50]) ~ 1 + t_p + idHand + idLeft + t_p*idHand + t_p*idHand*idLeft (6)
% log(nSpikesWin100_reward[m100_100]) ~ 1 + rewarded + rewDur:rewarded (3)
% 
% Hypothesis testing
% before fixOn, effect of context?
% prior modulation (idShort) around ready,set?
% ts modulation (catTs) around set?
% interval tuning: catTs beta, linear or quadratic?
% static representation coding: catTs*idShort
% link for measure & produce: catTs & catTs*tp
% modulation by idSwitch?
% after targetOn, target?

%% old old data description

% singleNeuronStat.mat: 11Anal, 67 Coeffi, H, #spikes, fitglm
% singleNeuronStat_smallWindow.mat: reduced ts timeWindow, same as singleNeuronStat
% singleNeuronStat_interaction_noSlide.mat: 12 Anal (+ts,ts^2), 70 Coeff, H&G(1 neuron removed in G), coeff & pval
% singleNeuronStat_interaction_slide.mat: statSlide{969 cell x 4 epoch}(coeff/pval x #regreesor x time)
% singleNeuronStat_interaction_slide_perisetOnly.mat:coeff & pval
% singleNeuronStat_interaction_glmfit.mat

%% oldold
% run Poisson regression for peri-set activity
% log(nSpikes_100_set) ~ b1 + b2*t_s +  b3*idShort + b4*idEye + b5*idRight ...
%                                         + b6*t_s*idShort + b7*t_s*idEye + b8*t_s*idRight ...
%                                           + b9*t_s^2 + b10*t_s^2* *idShort + b11*t_s^2*idEye + b12*t_s^2*idRight

% also 100ms sliding window for those
% - peri-fix: 2pr x 2EH [-100 500]
% - peri-tOn: 2x2x2target [-250 250]
% - peri-ready: 2x2x2 [-250 480]

% 100ms window: to avoid set flash response, only before set (but perhaps ok with up to 100ms after set)

 
 
 
 %%
% thP=0.05;
% for i=1:length(iAnalPerCoeff)
%     if i==1
%         disp(['===== ']);
%         disp(Formula{iAnalPerCoeff(i)});
%     else
%         if iAnalPerCoeff(i)-iAnalPerCoeff(i-1)==1
%             disp(['===== ' num2str(iAnalPerCoeff(i))]);
%             disp(Formula{iAnalPerCoeff(i)});
%         end
%     end
%     pSig=mean(pval(:,i)<thP)*100;    
%     disp([num2str(i)  '. ' CoefficientNames{i} '         :           ' num2str(pSig,3)]);
% end
% % 
% % i=1;
% % h=figure; load pplot.mat;
% %  histStairs(coeff(:,i),40,0,h,'k');
% %  axis tight; ylabel('# neurons');plotVertical(gca);xlabel('regression coefficient');
% %  title([CoefficientNames{i} ' in ' Formula{iAnalPerCoeff(i)}]);
 
 %% checking statSlide
% statSlide{iCell,1}=statTmp{13}; % 2(beta,pvalue) x 4([intercept idShort idEye idShort_Long:idEye_Hand]) x [fix-100 fix+500]
%  statSlide{iCell,2}=statTmp{14}; % 2(beta,pvalue) x 7([intercept idShort  idEye idRight_Left idShort_Long:idEye_Hand idEye_Hand:idRight_Left  idShort_Long:idRight_Left]) x [tOn-250 tOn+250]
%  statSlide{iCell,3}=statTmp{15};% 2(beta,pvalue) x 7([intercept idShort  idEye idRight_Left idShort_Long:idEye_Hand idEye_Hand:idRight_Left  idShort_Long:idRight_Left]) x [ready-250 ready+480]
%  statSlide{iCell,4}=statTmp{16};% 2(beta,pvalue) x 12([intercept  ts+t_s2+idShort+idEye+idRight+ idShort:ts+idEye:ts+idRight:ts+idShort:t_s2+idEye:t_s2+idRight:t_s2]) x [set-480 set+100]

% plot figures of % cells
% 1) fixation: idShort idEye idShortxidEye
% 2) tOn:idShort idEye idRight
% 3) ready:idShort idEye idRight
% 4) set: ts ts2 idShort idShort:ts idShort:t_s2

%  thP=0.05;
% load pplot.mat;
% for iFig=1:size(statSlide,2)
%     figure; set(gcf,'position',pplot.(['rect1_' num2str(iFig)]));
%  tmpStat=cat(4,statSlide{:,iFig}); % 2(beta,pValue) x 12 (regressors) x time x 969 neurons
%  
%  if iFig==1
%      tmpStatNm={'(Intercept)';...
%          'idShort';'idEye';'idShort:idEye'};
%      idReg=[2 3 4]; xL='time from fixation (ms)';
%  elseif iFig==4
%      tmpStatNm={'(Intercept)';...
%          'ts';'t_s2';'idShort';'idEye';'idRight';... % 1st order
%          'idShort:ts';'idEye:ts';'idRight:ts';'idShort:t_s2';'idEye:t_s2';'idRight:t_s2'}; % interaction
%      idReg=[2 3 4 7 10]; xL='time from set (ms)';
%  else
%      tmpStatNm={'(Intercept)';...
%          'idShort';'idEye';'idRight';...
%          'idShort:idEye';'idShort:idRight';'idEye:idRight '};
%      idReg=[2 3 4]; 
%      if iFig==2,xL='time from target on (ms)';else xL='time from ready (ms)'; end;
%  end
%  
%  idCmap=1;hmat=[];
%  for iReg=1:size(tmpStat,2)
%      pSig=squeeze(meanWONan(tmpStat(2,iReg,:,:)<thP,4)*100); % [time x 1] ,'omitnan'
%      
%      xTmp=squeeze(sumWONan(tmpStat(2,iReg,:,:)<thP,4)); % [time x 1] % ,'omitnan'
%      nTmp=squeeze(sum(~isnan(tmpStat(2,iReg,:,:)),4)); % [time x 1]     
%      [pHat,pCI]=binofit(xTmp,nTmp); % pCI: 1st column for lower bounds
%      
%       if iFig==1 & sum(iReg==idReg)~=0 % (iReg==2|iReg==3|iReg==4)
%         tmpT=-100+[0:(size(tmpStat,3)-1)];
%         hTmp=plot(tmpT,pSig,'color',pplot.cmap{idCmap});drawnow;ha;
% %         hTmp=shadedErrorBar(tmpT(:),pSig,[pCI(:,2)'-pSig(:)'; pCI(:,1)'-pSig(:)'],{'-','color',pplot.cmap{idCmap},'linewidth',2},1); drawnow;ha;
% %         hTmp=hTmp.mainLine;
%         hmat=[hmat;hTmp];
%         idCmap=idCmap+1;
%         
%       elseif iFig==2& sum(iReg==idReg)~=0 %& (iReg==2|iReg==3|iReg==4)
%           tmpT=-250+[0:(size(tmpStat,3)-1)];          
%         hTmp=plot(tmpT,pSig,'color',pplot.cmap{idCmap});drawnow;ha;
% %         hTmp=shadedErrorBar(tmpT(:),pSig,[pCI(:,2)'-pSig(:)'; pCI(:,1)'-pSig(:)'],{'-','color',pplot.cmap{idCmap},'linewidth',2},1); drawnow;ha;
% %         hTmp=hTmp.mainLine;
%         hmat=[hmat;hTmp];
%         idCmap=idCmap+1;
%         
%       elseif iFig==3& sum(iReg==idReg)~=0 %& (iReg==2|iReg==3|iReg==4)
%           tmpT=-250+[0:(size(tmpStat,3)-1)];          
%         hTmp=plot(tmpT,pSig,'color',pplot.cmap{idCmap});drawnow;ha;
% %         hTmp=shadedErrorBar(tmpT(:),pSig,[pCI(:,2)'-pSig(:)'; pCI(:,1)'-pSig(:)'],{'-','color',pplot.cmap{idCmap},'linewidth',2},1); drawnow;ha;
% %         hTmp=hTmp.mainLine;
%         hmat=[hmat;hTmp];
%         idCmap=idCmap+1;
%       elseif iFig==4& sum(iReg==idReg)~=0 %& (iReg==2|iReg==3|iReg==4|iReg==7|iReg==10)
%           tmpT=-480+[0:(size(tmpStat,3)-1)];
%         hTmp=plot(tmpT,pSig,'color',pplot.cmap{idCmap});drawnow;ha;
% %         hTmp=shadedErrorBar(tmpT(:),pSig,[pCI(:,2)'-pSig(:)'; pCI(:,1)'-pSig(:)'],{'-','color',pplot.cmap{idCmap},'linewidth',2},1); drawnow;ha;
% %         hTmp=hTmp.mainLine;
%         hmat=[hmat;hTmp];
%         idCmap=idCmap+1;
%       end
%       
%  end %  for iReg=1:size(tmpStat,2)
% axis tight;
%  legend(hmat,tmpStatNm(idReg),'location','best'); legend boxoff;
%  xlabel(xL); ylabel('% significant neurons');
%  
%  
% end % iFig (refernce time)
 

%% noSliding resutls using glmfit & sum(# spikes)
% ===== 
% log(spikeRate_fix_500) ~ 1 + idShort*idEye
% 1. (Intercept)         :           94.8
% 2. idShort_Long         :           41.3
% 3. idEye_Hand         :           44.2
% 4. idShort_Long:idEye_Hand         :           37.4
% ===== 2
% log(spikeRate_250_tOn) ~ 1 + idShort*idEye
% 5. (Intercept)         :           92.5
% 6. idShort_Long         :           29.1
% 7. idEye_Hand         :           34.7
% 8. idShort_Long:idEye_Hand         :           25.4
% ===== 3
% log(spikeRate_tOn_250) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 9. (Intercept)         :           88.9
% 10. idShort_Long         :           25.8
% 11. idEye_Hand         :           30.7
% 12. idRight_Left         :           17.4
% 13. idShort_Long:idEye_Hand         :           24.6
% 14. idShort_Long:idRight_Left         :           13.8
% 15. idEye_Hand:idRight_Left         :           14.7
% ===== 4
% log(spikeRate_ready_m250) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 16. (Intercept)         :           89.6
% 17. idShort_Long         :           27.2
% 18. idEye_Hand         :           36.5
% 19. idRight_Right         :           20.7
% 20. idShort_Long:idEye_Hand         :           26
% 21. idShort_Long:idRight_Right         :           12.5
% 22. idEye_Hand:idRight_Right         :           21.8
% ===== 5
% log(spikeRate_ready_480) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 23. (Intercept)         :           93.1
% 24. idShort_Long         :           48.6
% 25. idEye_Hand         :           56.1
% 26. idRight_Right         :           33.4
% 27. idShort_Long:idEye_Hand         :           42.1
% 28. idShort_Long:idRight_Right         :           20.9
% 29. idEye_Hand:idRight_Right         :           42.9
% ===== 6
% log(spikeRate_200_set) ~ 1 + t_s*idShort + t_s2*idShort
% 30. (Intercept)         :           89.3
% 31. t_s         :           46.4
% 32. t_s2         :           23.2
% 33. idShort_Long         :           43.1
% 34. t_s:idShort_Long         :           30.1
% 35. t_s2:idShort_Long         :           17.2
% ===== 7
% log(spikeRate_shortPrior) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 36. (Intercept)         :           82.4
% 37. idShort_Long         :           37.4
% 38. idEye_Hand         :           37.7
% 39. idRight_Right         :           19.8
% 40. idShort_Long:idEye_Hand         :           28.1
% 41. idShort_Long:idRight_Right         :           19
% 42. idEye_Hand:idRight_Right         :           37.3
% ===== 8
% log(spikeRate_200_set) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 43. (Intercept)         :           88.5
% 44. idShort_Long         :           39.8
% 45. idEye_Hand         :           49.8
% 46. idRight_Right         :           22.4
% 47. idShort_Long:idEye_Hand         :           30.5
% 48. idShort_Long:idRight_Right         :           18.5
% 49. idEye_Hand:idRight_Right         :           37.7
% ===== 9
% log(spikeRate_set_mintp) ~ 1 + t_p*idEye + t_p*idRight + idEye*idRight
% 50. (Intercept)         :           83.3
% 51. t_p         :           46.6
% 52. idEye_Hand         :           52.1
% 53. idRight_Right         :           25
% 54. t_p:idEye_Hand         :           46.5
% 55. t_p:idRight_Right         :           23.6
% 56. idEye_Hand:idRight_Right         :           44.2
% ===== 10
% log(spikeRate_go_mintp) ~ 1 + t_p*idEye + t_p*idRight + idEye*idRight
% 57. (Intercept)         :           82.5
% 58. t_p         :           48.2
% 59. idEye_Hand         :           50.4
% 60. idRight_Left         :           29
% 61. t_p:idEye_Hand         :           44.3
% 62. t_p:idRight_Left         :           27.6
% 63. idEye_Hand:idRight_Left         :           47.9
% ===== 11
% log(spikeRate_reward_100) ~ 1 + rewarded + rewDur:rewarded
% 64. (Intercept)         :           89.2
% 65. rewarded_1         :           16
% 66. rewDur:rewarded_1         :           12.5
% ===== 12
% log(spikeRate_fix_m100) ~ 1 + idShort*idEye
% 67. (Intercept)         :           90.6
% 68. idShort_Long         :           21.3
% 69. idEye_Hand         :           20.6
% 70. idShort_Long:idEye_Hand         :           20

%% noSliding resutls using fitglm & 1000*mean (spike rate)
% ===== 
% log(spikeRate_fix_500) ~ 1 + idShort*idEye
% 1. (Intercept)         :           98.1
% 2. idShort_Long         :           55.8
% 3. idEye_Hand         :           59
% 4. idShort_Long:idEye_Hand         :           51.1
% ===== 2
% log(spikeRate_250_tOn) ~ 1 + idShort*idEye
% 5. (Intercept)         :           98.1
% 6. idShort_Long         :           59.5
% 7. idEye_Hand         :           61.7
% 8. idShort_Long:idEye_Hand         :           55.1
% ===== 3
% log(spikeRate_tOn_250) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 9. (Intercept)         :           97.3
% 10. idShort_Long         :           53.4
% 11. idEye_Hand         :           58.7
% 12. idRight_Left         :           48.3
% 13. idShort_Long:idEye_Hand         :           52.3
% 14. idShort_Long:idRight_Left         :           41.9
% 15. idEye_Hand:idRight_Left         :           46
% ===== 4
% log(spikeRate_ready_m250) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 16. (Intercept)         :           97.8
% 17. idShort_Long         :           54.8
% 18. idEye_Hand         :           62
% 19. idRight_Right         :           51.4
% 20. idShort_Long:idEye_Hand         :           53.6
% 21. idShort_Long:idRight_Right         :           43.4
% 22. idEye_Hand:idRight_Right         :           52.6
% ===== 5
% log(spikeRate_ready_480) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 23. (Intercept)         :           97.8
% 24. idShort_Long         :           62.7
% 25. idEye_Hand         :           67.7
% 26. idRight_Right         :           48.9
% 27. idShort_Long:idEye_Hand         :           55.9
% 28. idShort_Long:idRight_Right         :           40.5
% 29. idEye_Hand:idRight_Right         :           56.9
% ===== 6
% log(spikeRate_200_set) ~ 1 + t_s*idShort + t_s2*idShort
% 30. (Intercept)         :           97.9
% 31. t_s         :           72.3
% 32. t_s2         :           56.6
% 33. idShort_Long         :           68.9
% 34. t_s:idShort_Long         :           61.7
% 35. t_s2:idShort_Long         :           53.3
% ===== 7
% log(spikeRate_shortPrior) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 36. (Intercept)         :           94.8
% 37. idShort_Long         :           59.4
% 38. idEye_Hand         :           61
% 39. idRight_Right         :           45.8
% 40. idShort_Long:idEye_Hand         :           52.2
% 41. idShort_Long:idRight_Right         :           44
% 42. idEye_Hand:idRight_Right         :           60.1
% ===== 8
% log(spikeRate_200_set) ~ 1 + idShort*idEye + idShort*idRight + idEye*idRight
% 43. (Intercept)         :           97.9
% 44. idShort_Long         :           70.1
% 45. idEye_Hand         :           73.1
% 46. idRight_Right         :           56.6
% 47. idShort_Long:idEye_Hand         :           63.3
% 48. idShort_Long:idRight_Right         :           52.8
% 49. idEye_Hand:idRight_Right         :           65.8
% ===== 9
% log(spikeRate_set_mintp) ~ 1 + t_p*idEye + t_p*idRight + idEye*idRight
% 50. (Intercept)         :           95.6
% 51. t_p         :           63.9
% 52. idEye_Hand         :           66
% 53. idRight_Right         :           44.6
% 54. t_p:idEye_Hand         :           61
% 55. t_p:idRight_Right         :           42.5
% 56. idEye_Hand:idRight_Right         :           60.2
% ===== 10
% log(spikeRate_go_mintp) ~ 1 + t_p*idEye + t_p*idRight + idEye*idRight
% 57. (Intercept)         :           94.6
% 58. t_p         :           64.3
% 59. idEye_Hand         :           64.6
% 60. idRight_Left         :           48.7
% 61. t_p:idEye_Hand         :           61.2
% 62. t_p:idRight_Left         :           47.8
% 63. idEye_Hand:idRight_Left         :           65.1
% ===== 11
% log(spikeRate_reward_100) ~ 1 + rewarded + rewDur:rewarded
% 64. (Intercept)         :           97.8
% 65. rewarded_1         :           65.9
% 66. rewDur:rewarded_1         :           58.9
% ===== 12
% log(spikeRate_fix_m100) ~ 1 + idShort*idEye
% 67. (Intercept)         :           98
% 68. idShort_Long         :           68
% 69. idEye_Hand         :           68
%  70. idShort_Long:idEye_Hand         :           68.2

 %  log(spikeRate_200_set) ~ ts+t_s2+idShort+idEye+idRight+ idShort:ts+idEye:ts+idRight:ts+idShort:t_s2+idEye:t_s2+idRight:t_s2
% 1. (Intercept)         :           97.6
% 2. ts         :           69.2
% 3. t_s2         :           56
% 4. idShort         :           68.9
% 5. idEye         :           80
% 6. idRight         :           69.5
% 7. idShort:ts         :           63.1
% 8. idEye:ts         :           52.6
% 9. idRight:ts         :           61
% 10. idShort:t_s2         :           51.5
% 11. idEye:t_s2         :           53
% 12. idRight:t_s2         :           48.9

%% checking single-neuron stats
% % index for each stat
% % iAnalPerCoeff1=5; iCoeff=24;CoefficientNames1='idShort_Long';legName1='prior(early t_s)';
% iAnalPerCoeff2=9; jCoeff=52;CoefficientNames2='t_p';legName2='production(early t_p)';
% iAnalPerCoeff1=6; iCoeff=31;CoefficientNames1='t_s';legName1='measure(late t_s)';
% % iAnalPerCoeff2=6; jCoeff=31;CoefficientNames2='t_s';legName2='measure(late t_s)';
% 
% thP=0.05;
% h=figure; load pplot.mat;kColor=1;
% legName={legName1;legName2};
% for i=1:nCoeff
%     if i==1
%         disp(['===== ']);
%         disp(Formula{iAnalPerCoeff(i)});
%     else
%         if iAnalPerCoeff(i)-iAnalPerCoeff(i-1)==1
%             disp(['===== ']);
%             disp(Formula{iAnalPerCoeff(i)});
%         end
%     end
%     pSig=mean(pval(:,i)<thP)*100;    
%     disp([CoefficientNames{i} '         :           ' num2str(pSig,3)]);
%     % hist for regression coefficients; only main effect
%     if (iAnalPerCoeff(i)==iAnalPerCoeff1 & strcmp(CoefficientNames{i},CoefficientNames1)) ||...
%             (iAnalPerCoeff(i)==iAnalPerCoeff2 & strcmp(CoefficientNames{i},CoefficientNames2))
%         subplot(2,1,kColor);
%         htmp=histStairs(coeff(:,i),40,0,h,pplot.cmap{kColor});
%         axis tight; ylabel('# neurons');legend(legName{kColor},'location','best'); legend boxoff;plotVertical(gca);
%         kColor=kColor+1;
%     end
% end
%  xlabel('regression coefficient');
% % scatter plot
% idPriorEarlyTs=(iAnalPerCoeff==iAnalPerCoeff1 & strcmp(CoefficientNames,CoefficientNames1));
% idTsLateTs=(iAnalPerCoeff==iAnalPerCoeff2 & strcmp(CoefficientNames,CoefficientNames2));
% figure; set(gcf,'position',pplot.rect2_3);
% tmpX=coeff(:,idPriorEarlyTs);
% tmpY=coeff(:,idTsLateTs);
% plot(tmpX,tmpY,'k.'); axis tight; hold on; plotVertical(gca); plotHorizon(gca);plotIdentity(gca); 
% gmfit{1}=fitgmdist(coeff(:,[iCoeff jCoeff]),1);
% % multiple cluster?
% gmfit{2}=fitgmdist(coeff(:,[iCoeff jCoeff]),2);
% % ops=statset('MaxIter',1000);
% gmfit{3}=fitgmdist(coeff(:,[iCoeff jCoeff]),3); % ,'options',ops);
% disp(num2str([gmfit{1}.AIC gmfit{2}.AIC gmfit{3}.AIC;...
%     gmfit{1}.BIC gmfit{2}.BIC gmfit{3}.BIC]));
% disp('');
% idMinBic=min([gmfit{1}.BIC gmfit{2}.BIC gmfit{3}.BIC])==[gmfit{1}.BIC gmfit{2}.BIC gmfit{3}.BIC];
% ezcontour(@(x,y)pdf(gmfit{idMinBic},[x y]),[mean(tmpX)-[1 -1]*3*std(tmpX)],[mean(tmpY)-[1 -1]*3*std(tmpY)]);
% axis tight;title([]);
% xlabel(['\beta ' legName1]); ylabel(['\beta ' legName2]);
% 
% acosd((tmpX./norm(tmpX))'*(tmpY./norm(tmpY)))
% 
% % chi2 test for independence
% contMat=[nnz(pval(:,iCoeff)>0.05&pval(:,jCoeff)>0.05) nnz(pval(:,iCoeff)<0.05&pval(:,jCoeff)>0.05);...
% nnz(pval(:,iCoeff)>0.05&pval(:,jCoeff)<0.05) nnz(pval(:,iCoeff)<0.05&pval(:,jCoeff)<0.05)];
% 
% [chi,pChi]=chi2testIndep(contMat)










%% representative cell list; check STAT

% tmp=[161218 57 1;...
% 161221 4 1;...
% 161211 7 1;...
% 161220 16 1;...
% 161222 25 1;...
% 161214 39 1;...
% 161223 23 3;...
% 161222 33 4;...
% 161208 57 1;...
% 161206 51 4;...
% 161222 15 1;...
% 161223 41 1;...
% 161222 45 4;...
% 161220 47 4;...
% 161220 15 1;...
% 161210 67 1];
% 
% iMat=[];
% for i=1:length(tmp)
%     id=find(fidMat(:,1)==tmp(i,1) &...
%         fidMat(:,2)==tmp(i,2) &...
%         fidMat(:,3)==tmp(i,3));
%     iMat=[iMat; id];
%     table(statName,stat(:,id))
%     figure; waitforbuttonpress; close;
% end
% 
% % num2str(stat(:,iMat))
% 
% table(statName,round(mean(stat<0.05,2)*100),round(mean(stat<0.001,2)*100),'variablenames',{'stat','pCellp05','pCellp001'})

%% best sessions for single-trial analysis
%% pick up common trials across cells
% listSess=[161221;...
%       161220;...
%       161222;...
%       161223;...
%       161211];
%   load pplot.mat;
%   for i=1:length(listSess)
%       idSess=fidMat(:,1)==listSess(i);
%       nC=nnz(idSess);
%       nT=max(cnameMat(fidMat(:,1)==listSess(i),end));
%       tmp=zeros(nC,nT);
%       for j=1:nC
%           startT=cnameValid{strcmp(fname,['H_RSGprior_20' num2str(listSess(i))])}(j,end-1);
%           endT=cnameValid{strcmp(fname,['H_RSGprior_20' num2str(listSess(i))])}(j,end);
%           tmp(j,startT:endT)=1;
%       end
%       figure; set(gcf,'position',pplot.(['rect1_' num2str(i)]));
%       imagesc(tmp); colormap('gray');
%       xlabel('trials'); ylabel('cell');
%       
%       edges=[1:nT]-0.5;
%       nStartT=histc(fidMat(idSess,end-1),edges);
%       nTlength=histc(fidMat(idSess,end)-fidMat(idSess,end-1),edges);
%       nEndT=histc(fidMat(idSess,end),edges);
%       figure; set(gcf,'position',pplot.(['rect2_' num2str(i)]));
%       plot(1:nT,cumsum(nStartT),'r',1:nT,cumsum(nTlength),'g',1:nT,cumsum(nEndT),'b'); 
%       axis tight; xlabel('trials'); ylabel('cumulative # cells');
% %       figure; plot([1:nT]'.*(max(cumsum(nTlength))-cumsum(nTlength)));
%   end


%% cell list


%% exploratory analysis for example curvedx cells
% note from old>if not exist new sorting done
% last column for manual id for ramp(1)/curved(2)
cn=[161204        1042 0 
      161204        1083 0
      161204        1106 0
      161204        1107 0
      161204        1124 2 % curved
      161204        1131 2
      161204        1141 2
      161204        1156 0
      161205        1106 0
      161205        1136 0
      161205        1152 0
      161205        1154 2
      161205        1170 2
      161206        1011 0
      161206        1028 0
      161206        1097 0
      161206        1111 0
      161206        1202 0
      161206        2016 0
      161206        2087 0
      161206        2094 1
      161206        2120 0
      161210        3000 2
      161210        3009 2
      161210        3011 2
      161210        3040 1
      161210        3082 2
      161210        3104 2
      161211        1007 2
      161211        1012 2 % curved hand only
      161211        1019 2
      161211        1029 2
      161211        1044 2
      161211        1053 1
      161211        1102 2
      161211        1111 2
      161218        2002 2
      161218        2007 2
      161218        2022 2
      161219        2090 2
      161219        2109 2
      161220        2144 2
      161221        1116 2
      161221        1163 2
      161221        2114 2
      161221        2135 0
      161221        2150 2
      161222        1126 2
      170506        1004 2
      170506        1034 2
      170506        1111 2
      170506        2015 1
      170506        2058 2
      170507        1024 1
      170507        1097 0
      170817        1007 0
      170817        1147 1 %ramp down
      170818        1097 2
      170818        1136 2
      170818        2153 2
      170818        3024 1
      170818        3031 1
      170818        3088 1
      170821        1067 1
      170821        2058 1 %ramp down
      170822        1159  2
      170822        1179 2
      170822        2013 2
      170823        1003 2
      170823        1099 2
      170823        1101 1 % ramp
      170823        1104 2
      170823        1136 1
      170823        1148 2
      170823        2023 2
      170823        2107 2
      170823        2111 2
      170823        2147 2
      170823        3167 2
      ];

type={'N.A.';'ramp';'curved'};

% % from runPlotPSTH_CCN17_measOnly
%     [161210 3011;...%% prior
%     161210 3104;...%% prior
%     161211 1107;...%% prior
%     161211 1019;...%% prior
%     161211 1044;...%% prior
%     161211 1053;...%% ramp
%     161211 1102;...%% prior
%     161218 2002;...%% prior
%     161218 2022;...%% prior
%     161222 1126;...%% prior
%     161203 1034;...
%     161204 1003;...
%     161204 1009;...
%     161204 1042;...
%     161204 1083;...
%     161204 1106;...
%     161204 1141;... %% prior
%     161205 1136;...
%     161206 1011;...
%     161206 1028;...
%     161206 1097;...
%     161206 1111;...
%     161206 2016;...
%     161206 2087;...
%     161206 2120;...
%     % G
%     170506 2015;... % ramping
%     170506 1111;... % prior
%     170507 1127;... % motor
%     170507 1117;... % hand
%     170507 1097;...
%     170507 1024;... % ramp/prior
%     170817 1147;... % ramp down
%     170817 1007;... % ramp 
%     170818 3088;... % ramp 
%     170818 3031;... % ramp 
%     170818 2153;... % prior
%     170818 2018;... 
%     170818 1136;... 
%     170818 1097;... 
%     170821 2058;... 
%     170821 1067;... 
%     170822 2013;...
%     170822 1179;...
%     170822 1159;...
%     170822 2013;...
%     170823 3167;...
%     170823 2147;...
%     170823 2111;...
%     170823 2107;...
%     170823 2023;...
%     170823 1148;...
%     170823 1136;...
%     170823 1104;...
%     170823 1101;...
%     170823 1099;...
% ];
% 
% % red tagged in PSTH_kilosrotOld
% 161204 1106;...
% 161204 1107;...
% 161204 1124;...
% 161204 1131;...
% 161204 1141;...
% 161204 1156;...
% 161205 1106;...
% 161205 1152;...
% 161205 1154;...
% 161205 1170;...
% 161206 1111;...
% 161206 1202;...
% 161206 2094;...
% 161210 3000;...
% 161210 3009;...
% 161210 3040;...
% 161210 3082;...
% 161210 3104;...
% 161211 1007;...
% 161211 1012;...
% 161211 1019;...
% 161211 1029;...
% 161211 1044;...
% 161211 1102;...
% 161211 1111;...
% 161218 2002;...
% 161218 2007;...
% 161218 2022;...
% 161219 2090;...
% 161219 2109;...
% 161220 2111;...
% 161220 2144;...
% 161221 1116;...
% 161221 1163;...
% 161221 2114;...
% 161221 2135;...
% 161221 2150;...
% 161222 1126;...
% 170506 1004;...
% 170506 1034;...
% 170506 1111;...
% 170506 2058;...
% 170507 1024;...
% 170818 1097;...
% 170818 2037;...
% 170818 2153;...
% 170818 3024;...
% 170822 1159;...
% 170822 1179;...
% 170822 2013;...
% 170823 1003;...
% 170823 1099;...
% 170823 1104;...
% 170823 1148;...
% 170823 2107;...
% 170823 2111;...
% 170823 2147;...
% 161210 3011;...


%% colormap archives
function plotNormMaxColormap(ws)
% sorting by short prior
            idNorm=1; % 0; % 1;
            
% dS dL dSL

v2struct(ws);
% colormap(jet);

tSetS=tBuffer:80:size(dS,2);
tSetL=tBuffer:100:size(dL,2);

% normalize
nDS=bsxfun(@times,dS,1./max(dS,[],ndims(dS))); % making max frinig rate 1
nDL=bsxfun(@times,dL,1./max(dL,[],ndims(dL))); % making max frinig rate 1
nDSL=bsxfun(@times,dSL,1./max(dSL,[],ndims(dSL))); % making max frinig rate 1

% sorting for short prior
[Y,tMax]=max(nDS,[],ndims(nDS));
[Y,iSortTmax]=sort(tMax,1); % default: ascending

if idNorm
snDS=nDS(iSortTmax,:);
snDL=nDL(iSortTmax,:);
snDSL=nDSL(iSortTmax,:);
else
snDS=dS(iSortTmax,:);
snDL=dL(iSortTmax,:);
snDSL=dSL(iSortTmax,:);
    
end
% plot
hFig=figure;

if idNorm
title('normalized spike rate');
else
    title('spike rate');
end
subplot(2,2,1) % short
imagesc(snDS);
axis tight; drawnow; plotVertical(gca,tSetS,[]);

subplot(2,2,2);
tmpD=dS(iSortTmax(1:25:end),:); % dS(iSortTmax([-100:5:100]+round(length(iSortTmax)/2)),:);
plotmultiple([],tmpD,[],[]);
shadedErrorBar([],mean(dS,1),sem(dS,1),'k-');
% plotmultiple([],dS(iSortTmax(1:25:end),:),[],[]);
axis tight; drawnow; plotVertical(gca,tSetS,[]);
% subplot(2,2,2); % time course color coded by peak timing
% plotmultiple([],dS(iSortTmax,:),[],[]);
% axis tight; drawnow; plotVertical(gca,tSetS,[]);
% set(gca,'yscale','log');

subplot(2,2,3) % SL
imagesc(snDSL);
axis tight; drawnow; plotVertical(gca,tSetS,[]);
xlabel('time');
ylabel('neurons (sorted by peak timing of short)');

subplot(2,2,4) % long
imagesc(snDL);
axis tight; drawnow; plotVertical(gca,tSetL,[]);
xlabel('time');



% subplot(2,2,2);
%             hBar=colorbar; set(hBar.Label,'String','norm. firing rate (sp/s)');
% disp(ws.psthNm{ws.i});
% disp(ws.psthNm{ws.i}(1:(length(ws.psthNm{ws.i})-4)));

% if idSaveFig
%             exportfig(hFig,['spikeRateMap_sortShort_' ws.psthNm{ws.i}(1:(length(ws.psthNm{ws.i})-4)) '.' optsExpFig.Format],optsExpFig);
% 
% saveas(hFig,['spikeRateMap_sortShort_' ws.psthNm{ws.i}(1:(length(ws.psthNm{ws.i})-4)) '.fig']);
%             close(hFig);drawnow;
% end
            
            
            
            % sorting for long prior
[Y,tMax]=max(nDL,[],ndims(nDL));
[Y,iSortTmax]=sort(tMax,1); % default: ascending

if idNorm
snDS=nDS(iSortTmax,:);
snDL=nDL(iSortTmax,:);
snDSL=nDSL(iSortTmax,:);
else
    snDS=dS(iSortTmax,:);
snDL=dL(iSortTmax,:);
snDSL=dSL(iSortTmax,:);

end
% plot
hFig=figure;
if idNorm
title('normalized spike rate');
else
    title('spike rate');
end
subplot(2,2,1) % short
imagesc(snDS);
axis tight; drawnow; plotVertical(gca,tSetS,[]);

subplot(2,2,2); % time course color coded by peak timing
tmpD=dL(iSortTmax(1:25:end),:); % dL(iSortTmax([-100:5:100]+round(length(iSortTmax)/2)),:);
plotmultiple([],tmpD,[],[]);
shadedErrorBar([],mean(dL,1),sem(dL,1),'k-');
% plotmultiple([],dL(iSortTmax(1:25:end),:),[],[]);
axis tight; drawnow; plotVertical(gca,tSetS,[]);
% subplot(2,2,2); % time course color coded by peak timing
% plotmultiple([],dL(iSortTmax,:),[],[]);
% axis tight; drawnow; plotVertical(gca,tSetS,[]);
% set(gca,'yscale','log');

subplot(2,2,3) % SL
imagesc(snDSL);
axis tight; drawnow; plotVertical(gca,tSetS,[]);
xlabel('time');
ylabel('neurons (sorted by peak timing of long)');

subplot(2,2,4) % long
imagesc(snDL);
axis tight; drawnow; plotVertical(gca,tSetL,[]);
xlabel('time');



% subplot(2,2,2);
%             hBar=colorbar; set(hBar.Label,'String','norm. firing rate (sp/s)');
% disp(ws.psthNm{ws.i});
% disp(ws.psthNm{ws.i}(1:(length(ws.psthNm{ws.i})-4)));
if idSaveFig
            exportfig(hFig,['spikeRateMap_sortLong_' ws.psthNm{ws.i}(1:(length(ws.psthNm{ws.i})-4)) '.' optsExpFig.Format],optsExpFig);

saveas(hFig,['spikeRateMap_sortLong_' ws.psthNm{ws.i}(1:(length(ws.psthNm{ws.i})-4)) '.fig']);
            close(hFig);drawnow;
end
            

%
function plotRawColormap(ws)

v2struct(ws);
colormap(jet);
% sorting for each condition (EH/ts/direction)

rTmp=firingRatesAverage{idEpoch}; %[neuron x 10ts x 2EH x 2direction x time]
[Y,tMax]=max(rTmp,[],ndims(rTmp));
for i=1:nPr
    for j=1:nTspp % length(T{i})
        for k=1:nEH
            for m=1:nTarg
                clf;
                [Y,iSortTmax]=sort(tMax(:,(i-1)*nTspp+j,k,m),1); % default: ascending
%                 r(:,(i-1)*nTspp+j,k,m,:)=r(iSortTmax,(i-1)*nTspp+j,k,m,:);
                r=rTmp(iSortTmax,(i-1)*nTspp+j,k,m,:);

                r=permute(r,[5 1 4 2 3]); %[time x neuron x 2direction x 10ts x 2EH]
                r=r(:,:); % [time x sortNeurons]
                r=r(~isnan(nanmean(r,2)),:);
                nTpoint=size(r,1); % -nnz(isnan(r(:,1)));
%                 if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
%                 elseif idEpoch==2, x=-max(Tall)+[0:(nTpoint-1)]'; %[set-ts set+tp]
%                 else x=-(max(nCol3)--durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]
                if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
                elseif idEpoch==2, x=-T{i}(j)+[0:(nTpoint-1)]'; %[set-ts set+tp]
                else x=-(nTpoint-durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]
                
                h=imagesc(x(:),1:size(r,2),r'); % legend(leg{2},'location','best','Box','off');
                xlabel(['time from ' epochNm{idEpoch} ' (ms)']); ylabel('neurons (sorted by peak timing)'); axis tight; drawnow; plotVertical(gca);
                hBar=colorbar; set(hBar.Label,'String','firing rate (sp/s)');
                saveas(gcf,[dirName 'PSTHs/colormap_rawMean/' prNm{i}(1) num2str(T{i}(j)) ehNm{k}(1) targNm{m} '_' epochNm{idEpoch}],'fig');
                exportfig(gcf,[dirName 'PSTHs/colormap_rawMean/' prNm{i}(1) num2str(T{i}(j)) ehNm{k}(1)  targNm{m} '_' epochNm{idEpoch}],optsExpFig);
            end % m
        end % k
    end % j
end % i



function plotNormColormap(ws)

v2struct(ws);
colormap(jet);
normR=bsxfun(@times,firingRatesAverage{idEpoch},1./nanmean(firingRatesAverage{idEpoch},ndims(firingRatesAverage{idEpoch}))); % making mean frinig rate 1

% sorting for each condition (EH/ts/direction)
rTmp=normR; %[neuron x 10ts x 2EH x 2direction x time]
[Y,tMax]=max(rTmp,[],ndims(rTmp));
for i=1:nPr
    for j=1:nTspp % length(T{i})
        for k=1:nEH
            for m=1:nTarg
                clf;
                [Y,iSortTmax]=sort(tMax(:,(i-1)*nTspp+j,k,m),1); % default: ascending
%                 r(:,(i-1)*nTspp+j,k,m,:)=r(iSortTmax,(i-1)*nTspp+j,k,m,:);
                r=rTmp(iSortTmax,(i-1)*nTspp+j,k,m,:);
                
                r=permute(r,[5 1 4 2 3]); %[time x neuron x 2direction x 10ts x 2EH]
                r=r(:,:); % [time x sortNeurons]
                r=r(~isnan(nanmean(r,2)),:);
                nTpoint=size(r,1); % -nnz(isnan(r(:,1)));
%                 if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
%                 elseif idEpoch==2, x=-max(Tall)+[0:(nTpoint-1)]'; %[set-ts set+tp]
%                 else x=-(max(nCol3)-durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]
                if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
                elseif idEpoch==2, x=-T{i}(j)+[0:(nTpoint-1)]'; %[set-ts set+tp]
                else x=-(nTpoint-durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]
                
                h=imagesc(x(:),1:size(r,2),r'); % legend(leg{2},'location','best','Box','off');
                xlabel(['time from ' epochNm{idEpoch} ' (ms)']); ylabel('neurons (sorted by peak timing)'); axis tight; drawnow; plotVertical(gca);
                hBar=colorbar; set(hBar.Label,'String','norm. firing rate (sp/s)');
                
                saveas(gcf,[dirName 'PSTHs/colormap_normMean/' prNm{i}(1) num2str(T{i}(j)) ehNm{k}(1) targNm{m} '_' epochNm{idEpoch}],'fig');
                exportfig(gcf,[dirName 'PSTHs/colormap_normMean/' prNm{i}(1) num2str(T{i}(j)) ehNm{k}(1) targNm{m} '_' epochNm{idEpoch}],optsExpFig);

            end % m
        end % k
    end % j
end % i


%% populationAvgPSTH archives
function plotRaw(ws)

% firingRatesAverage{1}=nan(N,nTall,nEH,nTarg,nT(1)); % 10 ts x 2 eh x 2 direction [ready-250 ready+ts]
% firingRatesAverage{2}=nan(N,nTall,nEH,nTarg,nT(2)); % 10 ts x 2 eh x 2 direction [set-ts set+tp]
% firingRatesAverage{3}=nan(N,nTall,nEH,nTarg,nT(3)); % 10 ts x 2 eh x 2 direction [production-tp production+100]

v2struct(ws);

mR=squeeze(nanmean(nanmean(firingRatesAverage{idEpoch}(:,:,idEH,:,:),1),4)); % output [10ts  x time] % x 2direct
nTpoint=size(mR,ndims(mR));
mR=reshape(mR,[nTall nTpoint]); % [20 x time] *nTarg

% time axis
if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
elseif idEpoch==2, x=-max(Tall)+[0:(nTpoint-1)]'; %[set-ts set+tp] mR(mR==0)=NaN;
else x=-(max(nCol3)-durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]

% plot
h=plot(x(:),mR','linewidth',lw); % legend(leg{2},'location','best','Box','off'); 
% xlabel(['time from ' epochNm{idEpoch} ' (ms)']); ylabel('firing rate (sp/s)'); 
% if idEpoch==2, set(gca,'xlim',x([1 end-1])); end; %?
axis tight; drawnow; 

% xlim similar to individual neurons' SDF
if idEpoch==1, set(gca,'tickdir','out','xtick',unique([0 min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
            'xticklabel',[],'xlim',[-100 max(T{2})],'ytick',[5 6 7]); % ,'ylim',[4.5 7.5] {0,min(T{1}),[],max(T{1}),[],max(T{2})}
elseif idEpoch==2, 
    tmpXtick=unique([0 -min(T{1}) -mean(T{1}) -max(T{1}) -mean(T{2}) -max(T{2}) min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]);
    set(gca,'tickdir','out','xtick',tmpXtick,...
        'xticklabel',[],... % {[],[],-max(T{1}),[],-min(T{1}),0,min(T{1}),[],max(T{1}),[],[]}
        'ytick',[5 6 7]); % ,'ylim',[4.5 7.5]); % ,...
%                     'yticklabel',[],'xticklabel',[]); % {-max(T{2}),[],-max(T{1}),[],-min(T{1}),0,min(T{1}),[],max(T{1}),[],max(T{2})});
                xlim([-max(T{2}) min(T{1})]);
else         set(gca,'xtick',([-400 -200 0]),'tickdir','out','xlim',[-800 50],...
        'ytick',[5 6 7]); % ,'ylim',[4.5 7.5],'xticklabel',[]); % ,...
%             'yticklabel',[],'xticklabel',[]); 
end;
plotVertical(gca);
% setting color
for i=1:length(h),h(i).Color=tmpCmap{(i>nTspp)+1,idEH}(rem(i-1,nTspp)+1,:);h(i).LineStyle=linestyle{ceil(i./nTall)};
    if idEpoch==1,hold all; idLast=find(~isnan(h(i).YData),1,'last');
    plot(h(i).XData(idLast),h(i).YData(idLast),'o','markerfacecolor','w','color',tmpCmap{(i>nTspp)+1,idEH}(rem(i-1,nTspp)+1,:),'markersize',8); end;% plot marker for the last points     
if idEpoch==2,hold all;
    plot(h(i).XData(h(i).XData==0),h(i).YData(h(i).XData==0),'o','markerfacecolor','w','color',tmpCmap{(i>nTspp)+1,idEH}(rem(i-1,nTspp)+1,:),'markersize',8); end;% plot marker for the last points    
end
set(gca,'box','off');
saveas(gcf,[dirName '/rawMeanPopulation_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection_' fnmEnd],'fig');
exportfig(gcf,[dirName '/rawMeanPopulation_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection_' fnmEnd],optsExpFig);

%%
function plotNorm(ws)

v2struct(ws);

normR=bsxfun(@times,firingRatesAverage{idEpoch},1./nanmean(firingRatesAverage{idEpoch},ndims(firingRatesAverage{idEpoch}))); % making mean frinig rate 1

mR=squeeze(nanmean(nanmean(normR(:,:,idEH,:,:),1),4)); % output [10ts  x time] % x 2direct
nTpoint=size(mR,ndims(mR));
mR=reshape(mR,[nTall nTpoint]); % [20 x time] *nTarg

% time axis
if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
elseif idEpoch==2, x=-max(Tall)+[0:(nTpoint-1)]'; %[set-ts set+tp] mR(mR==0)=NaN;
else x=-(max(nCol3)-durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]

% plot
h=plot(x(:),mR','linewidth',lw); % legend(leg{2},'location','best','Box','off'); 
% xlabel(['time from ' epochNm{idEpoch} ' (ms)']); ylabel('norm. firing rate (sp/s)'); 

% if idEpoch==2, set(gca,'xlim',x([1 end-1])); end;
axis tight; drawnow; 

% xlim similar to individual neurons' SDF
if idEpoch==1, set(gca,'tickdir','out','xtick',unique([0 min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
            'xticklabel',{0,min(T{1}),[],max(T{1}),[],max(T{2})},'xlim',[-100 max(T{2})],'ytick',[1 1.2]); % ,'ylim',[0.85 1.25]);
elseif idEpoch==2, 
    tmpXtick=unique([0 -min(T{1}) -mean(T{1}) -max(T{1}) -mean(T{2}) -max(T{2}) min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]);
    set(gca,'tickdir','out','xtick',tmpXtick,...
        'xticklabel',{-max(T{1}),[],-min(T{1}),0,min(T{1}),[],max(T{1})},...
        'ytick',[1 1.2]); % ,'ylim',[0.85 1.25]); % ,...
%                     'yticklabel',[],'xticklabel',[]); % {-max(T{2}),[],-max(T{1}),[],-min(T{1}),0,min(T{1}),[],max(T{1}),[],max(T{2})});
                xlim([-max(T{1}) max(T{1})]);
else         set(gca,'xtick',([-400 -200 0]),'tickdir','out','xlim',[-400 50],...
       'ytick',[1 1.2]); % ,'ylim',[0.85 1.25]); % ,...
%             'yticklabel',[],'xticklabel',[]); 
end;
plotVertical(gca);
% setting color
for i=1:length(h),h(i).Color=tmpCmap{(i>nTspp)+1,idEH}(rem(i-1,nTspp)+1,:);h(i).LineStyle=linestyle{ceil(i./nTall)};end
set(gca,'box','off');
exportfig(gcf,[dirName '/normMeanPopulation_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection_' fnmEnd],optsExpFig);
saveas(gcf,[dirName '/normMeanPopulation_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection_' fnmEnd],'fig');

%%
function plotVar(ws)
% variance across neurons, not across trials
v2struct(ws);

mR=squeeze(nanvar(nanmean(firingRatesAverage{idEpoch}(:,:,idEH,:,:),4),1)); % output [10ts  x time] % x 2direct
nTpoint=size(mR,ndims(mR));
mR=reshape(mR,[nTall nTpoint]); % [20 x time] *nTarg

% time axis
if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
elseif idEpoch==2, x=-max(Tall)+[0:(nTpoint-1)]'; %[set-ts set+tp] mR(mR==0)=NaN;
else x=-(max(nCol3)-durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]

% plot
h=plot(x(:),mR','linewidth',lw); % legend(leg{2},'location','best','Box','off'); 
% xlabel(['time from ' epochNm{idEpoch} ' (ms)']); ylabel('firing rate (sp/s)'); 
% if idEpoch==2, set(gca,'xlim',x([1 end-1])); end; %?
axis tight; drawnow; 

% xlim similar to individual neurons' SDF
if idEpoch==1, set(gca,'tickdir','out','xtick',unique([0 min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
            'xticklabel',[],'xlim',[-100 max(T{2})]); % ,'ytick',[5 6 7]); % ,'ylim',[4.5 7.5] {0,min(T{1}),[],max(T{1}),[],max(T{2})}
elseif idEpoch==2, 
    tmpXtick=unique([0 -min(T{1}) -mean(T{1}) -max(T{1}) -mean(T{2}) -max(T{2}) min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]);
    set(gca,'tickdir','out','xtick',tmpXtick,...
        'xticklabel',[]); % ,... % {[],[],-max(T{1}),[],-min(T{1}),0,min(T{1}),[],max(T{1}),[],[]}
%         'ytick',[5 6 7],'ylim',[4.5 7.5]); % ,...
%                     'yticklabel',[],'xticklabel',[]); % {-max(T{2}),[],-max(T{1}),[],-min(T{1}),0,min(T{1}),[],max(T{1}),[],max(T{2})});
                xlim([-max(T{2}) min(T{1})]);
else         set(gca,'xtick',([-400 -200 0]),'tickdir','out','xlim',[-800 50]); % ,...
%         'ytick',[5 6 7],'ylim',[4.5 7.5],'xticklabel',[]); % ,...
%             'yticklabel',[],'xticklabel',[]); 
end;
plotVertical(gca);
% setting color
for i=1:length(h),h(i).Color=tmpCmap{(i>nTspp)+1,idEH}(rem(i-1,nTspp)+1,:);h(i).LineStyle=linestyle{ceil(i./nTall)};
    if idEpoch==1,hold all; idLast=find(~isnan(h(i).YData),1,'last');
    plot(h(i).XData(idLast),h(i).YData(idLast),'o','markerfacecolor','w','color',tmpCmap{(i>nTspp)+1,idEH}(rem(i-1,nTspp)+1,:),'markersize',8); end;% plot marker for the last points     
if idEpoch==2,hold all;
    plot(h(i).XData(h(i).XData==0),h(i).YData(h(i).XData==0),'o','markerfacecolor','w','color',tmpCmap{(i>nTspp)+1,idEH}(rem(i-1,nTspp)+1,:),'markersize',8); end;% plot marker for the last points    
end
set(gca,'box','off');
saveas(gcf,[dirName '/var_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection_' fnmEnd],'fig');
exportfig(gcf,[dirName '/var_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection_' fnmEnd],optsExpFig);

% nTpoint=size(mR,ndims(mR));
% mR=reshape(mR,[nTall nTpoint]); % [20 x time] *nTarg
% if idEpoch==1,x=-minTon(1)+[0:(nTpoint-1)]'; % [ready-250 ready+ts]
% elseif idEpoch==2, x=-max(Tall)+[0:(nTpoint-1)]'; %[set-ts set+tp]mR(mR==0)=NaN;
% else x=-(max(nCol3)-durAfterGo)+[0:(nTpoint-1)]'; end % [production-tp production+100]
% h=plot(x(:),mR'); % legend(leg{2},'location','best','Box','off'); 
% xlabel(['time from ' epochNm{idEpoch} ' (ms)']); ylabel('\sigma^2(firing rate) (sp/s)'); 
% if idEpoch==2, set(gca,'xlim',x([1 end-1])); end;
% axis tight; drawnow; plotVertical(gca);
% for i=1:length(h),h(i).Color=tmpCmap{idEH}(rem(i-1,nTall)+1,:);h(i).LineStyle=linestyle{ceil(i./nTall)};end
% saveas(gcf,[dirName '/var_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection'],'fig');
% exportfig(gcf,[dirName '/var_' ehNm{idEH} '_' epochNm{idEpoch} '_avgDirection'],optsExpFig);