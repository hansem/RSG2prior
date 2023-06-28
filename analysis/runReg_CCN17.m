function runReg_CCN17

% run regression for prior-span activity
% log(@spikes in 80ms bin before set) = a*idTs + b*(idTs^2) + c*idPrior + 
%                                                           d*(idTs*idPrior) + e*(idTs^2*idPrior)
%
% spec of regressors
%   - idTs: normalize([-2 -1 0 1 2])
%   - idTs2: normalize([-2 -1 0 1 2].^2)
%   - idPrior: 1 for short, 0 for long
%
% idPlot: plot raster, PSTH, fitted PSTH with GLM
%
% note
% - ignore modality, direciton, tp, reward for now (note 45%,30% of neurons were modulated by them)
% - now using kilosort results

%% init
initRSG2prior;
dirName='/Users/hansem/Dropbox (MIT)/CCN17/regression/';
cd([dirName]);

% model setting
% nReg=5;
b=[]; % beta regression coefficient [# neurons x nReg]
p=[]; % p value of regression [# neurons x nReg]
c=[]; % cell ID [session cell]
r=[];

% misc setting
idSUonly=0;
eventId=4; % set

idPlot=1;

load pplot.mat; % tmpCmap

for iAnimal=1:1 % length(animalNm):(-1):1
    % for now, analyze only H's data (G's kilosort TBD)
    % session
    for i=1:length(fn{iAnimal})
        fnm=fn{iAnimal}{i}; % H_RSGprior_20161203
        disp(['===== ' fnm ' =====']);
        fid=fnm(end-5:end); % '161203'
        load(fullfile(neuDir2,fid)); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo
        
        % extract behavior
        beh=load([behDir animalNm{iAnimal} '_RSGprior_DMFC.mat']);
        
        % remove outlier
        idOut2=false(size(tInfo,1),1);
        movingBallTrials=unique(tInfo(:,1)); % size same as idNoise
        idSess=str2num(fid)==beh.sessId;
        if strcmp(fid,'161207'),idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
        if nnz(idSess)~=length(movingBallTrials),disp('smth wrong for trial align'); end;
        iOutSess=find(beh.idOut(idSess)); disp(['outlier trials (incl. abort/wait): ' num2str(length(iOutSess)) '/' num2str(length(movingBallTrials)) ',' ...
            num2str(length(iOutSess)/length(movingBallTrials)*100) '%']);
        idValid=beh.t(idSess)>0&beh.t(idSess)<3*beh.T(idSess);
        disp(['outlier trials (w/o abort/wait): ' num2str(nnz(beh.idOut(idSess)&idValid)) '/' num2str(length(movingBallTrials)) ',' ...
            num2str(nnz(beh.idOut(idSess)&idValid)/length(movingBallTrials)*100) '%']);
        for iOut=1:length(iOutSess)
            idOut2(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
        end
        
        % loop through each unit
        nUnit=size(idSU,1); % both SU & MU
        for j=1:nUnit
            
            if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
                % SU only
            else % SU & MU
                tmpCellId=idSU(j,1);
                disp(['cluster#' num2str(tmpCellId) ': #sp=' num2str(nnz(tmpCellId==idClust))]);
                % consistent with plotSDFexport
                NEV=kilosort2NEV(sp,idClust,idSU); % sp [#spikeTiming x 1], idClust [cellId(XXXX) of # spikes x 1], idSU[cellID idSU]
                tmpId=[idSU(j,1) idSU(j,1) 1 length(movingBallTrials)]; % id(electrode,Unit,start trial,end trial)
                NEV.MetaTags.DateTime=fid;
                
                % run actual regression
%                 try % # trials too small
                    statTmp=reg_CCN17(NEV,tmpId,tInfo(~idOut2,:),eventId,idPlot); % {[2(b,p) x nReg]; meanFR[10 x1]}
                    if ~isempty(statTmp)
                        b=[b; statTmp{1}(1,:)]; % const, ts, ts2, prior tsXpr ts2Xpr
                        p=[p; statTmp{1}(2,:)];
                        c=[c; str2num(fid) tmpCellId];
                        r=cat(3,r,statTmp{2}); % [2priors x 5 ts x #cell]
                    end
                    close all;
%                     drawnow;
%                     exportfig(gcf,fullfile(dirName,'PSTHs',[fid '_' num2str(tmpCellId) '.png']),'color','rgb','Height',13,'Width',11,'format','png'); % ,'Reference',hTmp);
%                     close all;

%                 end
                
            end %                     if idSUonly & idSU(j,2)==0 % MU when idSUonly==1            
        end % j nUnit
    end % for i=1:length(fn{iAnimal})
end % animal

save('reg.mat','b','p','c','r','statTmp');

%% checking proportion of prior, ramp, curved neurons
load('reg.mat');

%
%% check p(prior), p(ramp), p(curve)
% # valid cells: 592
% # modulation cells: 450
% % modulation cells: 0.76014
% i(ts)
% # cells: 258
% % cells: 0.43581
% i(ts)^2
% # cells: 124
% % cells: 0.20946
% prior
% # cells: 358
% % cells: 0.60473

% setting
regNm={'i(ts)','i(ts)^2','prior'}; % ,'i(ts) x prior','i(ts)^2  x prior'};
idCorrectMC=1; % controling for FDR

nReg=length(regNm);

% controling for FDR
if idCorrectMC
    H=[];
    P=p;
    thP=nan(1,nReg);
%     P=nan(size(p));
    for i=2:size(p,2) % across regressors        
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p(:,i),0.05,'pdep','no'); % desired false discovery rate, method:assuming indep/pos dep, no report
        thP(i-1)=crit_p;
        H=[H h(:)];
    end    
else
    P=p;
    thP=repmat(0.05,1,nReg);
end
save('reg.mat','H','thP','-append');
%     figure;plot(sort(p(:,2)));
%     ha; plot([1:size(p,1)]*0.05/size(p,1),'r-');
%     legend('sorted p value','FDR criteria'); xlabel('id for multiple comparison'); ylabel('p value');
    
% window
winSize=80;
stepSize=[80 100]; % short/long
% attrition

% thresholding with # trials for the longest ts (attrition)
thNT=15;% minimum 15 trials
thFR=1; % [sp/s] if trial-averaged FR <thFR, discard

disp('check p(prior), p(ramp), p(curve)');

% total number of celss
% 364 single 623 multi
% -18 46;... % 12/7
% -17 20;...
%
% 364+623-18-46-17-20=886

% 162 cells too low FR or too small # trials
nValidCell=size(b,1); % 724
disp(['# valid cells: ' num2str(nValidCell)]);

iPr=4;iRamp=2;iCurv=3;iRampPr=5;iCurvPr=6;

% task modulation
idTaskMod=sum(P(:,2:end)<repmat(thP,size(P,1),1),2)>0;
disp(['# modulation cells: ' num2str(nnz(idTaskMod))]); % 534
disp(['% modulation cells: ' num2str(mean(idTaskMod))]); % 74

for iReg=2:(nReg+1) % const 1st columbn
    disp(regNm{iReg-1});
    
    idSig=P(:,iReg)<thP(iReg-1); % TBD: ramp& curv separately
    disp(['# cells: ' num2str(nnz(idSig))]);
    disp(['% cells: ' num2str(mean(idSig))]);    
    
end % iReg

% check mean of curved neurons
r3=permute(r,[3 2 1]);
r3=r3(:,:); % [neurons x 10]
idCurved=p(:,iCurv)<thP(iCurv-1);
idConcave=b(:,iCurv)>0;
figure; 
errorbar(1:size(r3(idCurved&idConcave,1:5),2),mean(r3(idCurved&idConcave,1:5),1),sem(r3(idCurved&idConcave,1:5),1),'color','r');ha;
errorbar(1:size(r3(idCurved&idConcave,6:10),2),mean(r3(idCurved&idConcave,6:10),1),sem(r3(idCurved&idConcave,6:10),1),'color','b');ha;
set(gca,'xtick',1:5,'xticklabel',[]);xlabel('t(set)'); ylabel('mean FR (sp/s)');
figure; 
errorbar(1:size(r3(idCurved&~idConcave,1:5),2),mean(r3(idCurved&~idConcave,1:5),1),sem(r3(idCurved&~idConcave,1:5),1),'color','r');ha;
errorbar(1:size(r3(idCurved&~idConcave,6:10),2),mean(r3(idCurved&~idConcave,6:10),1),sem(r3(idCurved&~idConcave,6:10),1),'color','b');
set(gca,'xtick',1:5,'xticklabel',[]);xlabel('t(set)'); ylabel('mean FR (sp/s)');


% % ramp: scatter for slope
% % figure; hist(b(:,iRamp),100);
% tmpX=normalize(([1:5]-3),2); % regressor: change into
% tmpX(4); % 0.3162
% 
% figure; ha;
% nBin=15;cmap=[.5 .5 .5];markersize=10;
% idSigRamp=P(:,iRamp)<thP;
% b2=b(idSigRamp,iRampPr)*tmpX(4)/(winSize/1000);
% [n,x]=hist(b2,nBin);
% bar(x,n,'FaceColor',cmap,'EdgeColor','none');
% axis tight; maxYlim=max(ylim);
% plot(mean(b2),1.1*maxYlim,'color',cmap,'markerfacecolor',cmap,'marker','v','markersize',markersize);
% plotVertical(gca,0,[]);
% xlabel('differential ramping slope (spike/s)'); ylabel('# neurons');
% 
%         set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%         'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
%         'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
%         'LineWidth', 1)   
% 
% % curv: scatter for curv
% % figure; hist(b(:,iCurv),100);
% tmpX=normalize(([1:5]-3).^2,2); % regressor: change into
% 
% figure; ha;
% nBin=15;cmap=[.5 .5 .5];markersize=10;
% idSigCurv=P(:,iCurv)<thP;
% b2=b(idSigCurv,iCurvPr)*tmpX(4)/(winSize/1000);
% [n,x]=hist(b2,nBin);
% bar(x,n,'FaceColor',cmap,'EdgeColor','none');
% axis tight; maxYlim=max(ylim);
% plot(mean(b2),1.1*maxYlim,'color',cmap,'markerfacecolor',cmap,'marker','v','markersize',markersize);
% plotVertical(gca,0,[]);
% xlabel('differential curvature (spike/s)'); ylabel('# neurons');
% 
%         set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%         'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
%         'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
%         'LineWidth', 1)   

% colormap of mean FR: plot with log
% figure; hist(log(r(:)),100);
r=permute(r,[2 1 3]);
r2=reshape(r,size(r,1)*size(r,2),size(r,3)); % [10(tsShort,tsLong) x #neuron]
% r2=log(r2); % transform Poisson to more of normal
%PCA
[coeff,score]=pca(r2');
figure; plot(coeff(:,1:3));
figure; hist(score(:,1),50);
% TBD: sorting & smoothing, 
figure; imagesc(r2');
colorbar;
% mean: concluding pretty noisy
figure; errorbar(1:10,mean(r2,2),std(r2,0,2));
% check representative prior cell > code okay but issue is log transform????
i=find(c(:,1)==161210 & c(:,2)==3011); % 104);
 figure; plot(r2(:,i))

%% 
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