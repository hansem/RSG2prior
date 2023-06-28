function runPlotPSTH

% plot PSTH using kilosort results (no manual start/end trials)
% main dependent function: plotPSTH
% almost identical runPlotSDF
ccc;

% note
% difference b/t runPlotSDFnew vs runPlotSDF: % new stat: poisson regression

% save data from around 800 trials in 12/7
% data transfer speed warning issue in 12/13? two cells in NS3 data

%% init
initRSG2prior;
% fn={fnameH fnameG}; H_RSGprior_20161203
nAnimal=1;
cd(neuDir2); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo

%% 
idPlot=1;
load pplot.mat;
optsExpFig.Height='13'; % 7;
optsExpFig.Width='11';
optsExpFig.Format='png'; % 'png';
optsExpFig.LineMode='scaled';
        
% new stat: poisson regression
nAnal=11; % 10;
nCoeff=67;
Formula=cell(nAnal,1);
CoefficientNames=[]; % =cell(nCoeff,1);
coeff=[];coeffTmp=[];
pval=[];pvalTmp=[];

%     '1)    ' prior x eyeHand 22GLM [fix fix+500]'                                  0.24361      0.06015
%     '2)    ' prior x eyeHand 22GLM [tOn-250 tOn]'                                  0.16692     0.021053
%     '3)    ' prior x eyeHand x target 222GLM [tOn tOn+250]'                        0.06015    0.0015038
%     '4)    ' prior x eyeHand x target 222GLM [ready-250 ready]'                   0.046617            0
%     '5) ' prior x eyeHand x target 222GLM   ' prior x eyeHand 22ANOVA [ready ready+480]'                              0.26617     0.082707
%     '6)   ' ts + eyeHand x target regression ' prior x eyeHand 22ANOVA [set-480 set]'                                  0.26917     0.096241
%     '7) prior x eyeHand x target 222GLM [ready+480 ready+800] for overlap ts inc. all long prior?                           0.29624      0.14436       
%     ?8) ' prior x eyeHand x target 222GLM   ' prior x eyeHand 22ANOVA [set-480 set]'
%     '9) intercept tp + eyeHand x directions regression [set set+min(tp)]'      0.8812      0.73835
%     '10) intercept tp + eyeHand x directions regression [production-min(tp) production]';
%     '11) intercept reward regression [rewardOnset onset+100]'                 0.99098      0.95489
%     ' rewarded regression [rewardOnset onset+100]'                             0.1203     0.013534
%     ' rewDur*rewarded regression [rewardOnset onset+100]'                    0.079699    0.0045113
% stat=[];

for iAnimal=1:nAnimal % length(animalNm)% :(-1):1
    for i=length(fn{iAnimal}):(-1):1 % session
        fnm=fn{iAnimal}{i}; % H_RSGprior_20161203
        disp(['===== ' fnm ' =====']);
        fid=fnm(end-5:end);
        load(fid); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo
        
        % extract behavior
        beh=load([behDir animalNm{iAnimal} '_20' fid '.mat']);
        
        % loop through each unit
        nUnit=size(idSU,1); % both SU & MU
        for j=1:nUnit
            disp(['cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust))]);
            
            % consistent with plotSDFexport
            NEV=kilosort2NEV(sp,idClust,idSU);
            tmpId=[idSU(j,1) idSU(j,1) 1 length(beh.t)]; % id(electrode,Unit,start trial,end trial)
            
            CoefficientNames=[];
            
            try
            statTmp=plotSDFexport(NEV,tmpId,tInfo,0,idPlot); % p values cell
%             for k=1:nAnal
%                 Formula{k}=statTmp{k}.Formula;
%                 CoefficientNames=[CoefficientNames(:); statTmp{k}.CoefficientNames(:)];
%                 coeffTmp=[coeffTmp table2array(statTmp{k}.Coefficients(:,1))']; % 1 x nCoeff
%                 pvalTmp=[pvalTmp table2array(statTmp{k}.Coefficients(:,end))'];
%             end
%             coeff=[coeff;coeffTmp];coeffTmp=[];
%             pval=[pval;pvalTmp ];pvalTmp=[];
            
        
        if idPlot
            set(gcf,'PaperPositionMode','auto');
            saveas(gcf,[neuDir 'PSTH_kilosort/' fid '_' num2str(idSU(j,1)) '_' sortQ{idSU(j,2)+1} '.png']);% sortQ={'MU','SU'};
        end
        
        close all;
            catch
                disp(['ERROR cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust))]);
            end
        end % for j=1:nUnit
        
    end % for i=1:length(fname)
    
end % for iAnimal=1:length(animalNm)

iAnalPerCoeff=[];n=[];for i=1:length(statTmp),n=n+statTmp{i}.NumCoefficients;iAnalPerCoeff=[iAnalPerCoeff;repmat(i,statTmp{i}.NumCoefficients,1)];end
%  save([neuDir 'singleNeuronStat_smallWindow.mat'],'CoefficientNames','Formula','coeff','pval','nCoeff','nAnal','iAnalPerCoeff','statTmp');
 


%     statName={'1) prior 22ANOVA [fix fix+500]';
%     ' eyeHand 22ANOVA [fix fix+500]';
%     ' prior x eyeHand 22ANOVA [fix fix+500]';
%     '2) prior 22ANOVA [tOn-250 tOn]';
%     ' eyeHand 22ANOVA [tOn-250 tOn]';
%     ' prior x eyeHand 22ANOVA [tOn-250 tOn]';
%     '3) prior 222ANOVA [tOn tOn+250]';
%     ' eyeHand 222ANOVA [tOn tOn+250]';
%     ' target 222ANOVA [tOn tOn+250]';
%     ' prior x eyeHand 222ANOVA [tOn tOn+250]';
%     ' prior x target 222ANOVA [tOn tOn+250]';
%     ' eyeHand x target 222ANOVA [tOn tOn+250]';
%     ' prior x eyeHand x target 222ANOVA [tOn tOn+250]';
%     '4) prior 222ANOVA [ready-250 ready]';
%     ' eyeHand 222ANOVA [ready-250 ready]';
%     ' target 222ANOVA [ready-250 ready]';
%     ' prior x eyeHand 222ANOVA [ready-250 ready]';
%     ' prior x target 222ANOVA [ready-250 ready]';
%     ' eyeHand x target 222ANOVA [ready-250 ready]';
%     ' prior x eyeHand x target 222ANOVA [ready-250 ready]';
%     '5) prior 22ANOVA [ready ready+400]';
%     ' eyeHand 22ANOVA [ready ready+400]';
%     ' prior x eyeHand 22ANOVA [ready ready+400]';
%     '6) prior 22ANOVA [set-400 set]';
%     ' eyeHand 22ANOVA [set-400 set]';
%     ' prior x eyeHand 22ANOVA [set-400 set]';
%     '7) prior 22ANOVA [set-320 set] for overlap ts';
%     ' eyeHand 22ANOVA [set-320 set] for overlap ts';
%     ' prior x eyeHand 22ANOVA [set-320 set] for overlap ts';
%     '8) intercept tp + eyeHand x directions regression [set set+min(tp)]';
%     ' tp regression [set set+min(tp)]';
%     ' eyeHand regression [set set+min(tp)]';
%     ' directions regression [set set+min(tp)]';
%     ' tp x eyeHand regression [set set+min(tp)]';
%     ' tp x directions regression [set set+min(tp)]';
%     ' eyeHand x directions regression [set set+min(tp)]';
%     ' tp x eyeHand x directions regression [set set+min(tp)]';
%     '9) intercept tp + eyeHand x directions regression [production-min(tp) production]';
%     ' tp regression [production-min(tp) production]';
%     ' eyeHand regression [production-min(tp) production]';
%     ' directions regression [production-min(tp) production]';
%     ' tp x eyeHand regression [production-min(tp) production]';
%     ' tp x directions regression [production-min(tp) production]';
%     ' eyeHand x directions regression [production-min(tp) production]';
%     ' tp x eyeHand x directions regression [production-min(tp) production]';
%     '10) intercept reward regression [rewardOnset onset+100]';
%     ' rewarded regression [rewardOnset onset+100]';
%     ' rewDur*rewarded regression [rewardOnset onset+100]';};
%     save([neuDir 'stat.mat'],'stat','statName'); 
    
    
%     Generalized Linear regression model:
%     log(numSpikes) ~ 1 + idShort*idEye
%     Distribution = Poisson
% 
% Estimated Coefficients:
%                                Estimate      SE        tStat       pValue  
%                                ________    _______    _______    __________
% 
%     (Intercept)                 -1.2879      0.125    -10.303    6.8415e-25
%     idShort_Long                 0.5025    0.16428     3.0588     0.0022225
%     idEye_Hand                  0.37914     0.1583     2.3951      0.016617
%     idShort_Long:idEye_Hand    -0.35537    0.21156    -1.6798      0.092998
% 
% 
% 945 observations, 941 error degrees of freedom
% Dispersion: 1
% Chi^2-statistic vs. constant model: 14.2, p-value = 0.00262

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
% 
% %% representative cell list; check STAT
% 
% % tmp=[161218 57 1;...
% % 161221 4 1;...
% % 161211 7 1;...
% % 161220 16 1;...
% % 161222 25 1;...
% % 161214 39 1;...
% % 161223 23 3;...
% % 161222 33 4;...
% % 161208 57 1;...
% % 161206 51 4;...
% % 161222 15 1;...
% % 161223 41 1;...
% % 161222 45 4;...
% % 161220 47 4;...
% % 161220 15 1;...
% % 161210 67 1];
% % 
% % iMat=[];
% % for i=1:length(tmp)
% %     id=find(fidMat(:,1)==tmp(i,1) &...
% %         fidMat(:,2)==tmp(i,2) &...
% %         fidMat(:,3)==tmp(i,3));
% %     iMat=[iMat; id];
% %     table(statName,stat(:,id))
% %     figure; waitforbuttonpress; close;
% % end
% % 
% % % num2str(stat(:,iMat))
% % 
% % table(statName,round(mean(stat<0.05,2)*100),round(mean(stat<0.001,2)*100),'variablenames',{'stat','pCellp05','pCellp001'})
% 
% %% best sessions for single-trial analysis
% %% pick up common trials across cells
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

% function NEV=kilosort2NEV(sp,idClust,idSU)
% 
% % convert kilosort result to NEV
% 
% NEV.MetaTags.SampleRes=30000;
% NEV.Data.Spikes.TimeStamp=sp;
% % in kilosort idClust is unique across electrodes & units
% NEV.Data.Spikes.Electrode=idClust;
% NEV.Data.Spikes.Unit=idClust;
