function behAnal_GRSGprior

% 2019/3/25
% show deltaTp for sigmoid now between two priors
% plot deltaTp histogram combining across animals (idPool)

% 2019/2/18: 
% model-based comparison b/t linear & BLS
% check var(tp|extreme ts)< E(var(tp|priorMean))
% full probabilistic linear model fit: idFitLlinear, estLinear.m

% 2019/1/10: model-free test of sigmodial

% 9/16/2018
% sFig1: distribution of weights
% sFig3: bias-variance (BLS-data)

% 8/31/2018
% rank sum test for overalp: ignoring EH/RL
% supplementary figure 1: copy of figure1 for each condition (ER,EL,HR,HL)
% supplementary figure 2: bias-var (BLS), data vs model

% 8/18/2018
% everythign condition-speciifc (animal, priorX, effector, direction)
% condition-specific BLS fitting

% 8/17/2018
% slope with weighted regression

% 8/3/2018
% stat for # trial/block
% plot tp histogram right after switch @ overlap ts

% 5/10/2018
% plot tp-ts, modified DnPlotResponses2P with PlotDiff

% 4/23/2018
% do stat on linear regression slope for tp-ts

% 3/5/2018
% switching across blocks (bias=f(trials after switch), only overlap?
%   clean up the code using initRSG2prior

% 2/8/2018
% compare bias b/t longest vs shortest ts for each prior

% save NEV, NS3

% 2017/8/24:
%  setting noSet trials as all outlier (tp in those trials are not <0  always)

%% init
initRSG2prior;

iAnimal=2;
fname=fn{iAnimal};
Tmat=T;

idExpFig=0;

idFitBLS=0; % 1; % 0;

idFitLlinear=0; % 1; % 0; % 1;

idPool=1;

%% merge data and save into one file

% % T 
% % t 
% % idShortTrial
% % idHandEye
% % theta
% % fixTimeDur
% % targetTimeDur
% % iti
% % reward
% % FeedbackTimeDur
% sessId=[];
% for i=1:nS
%     
%     disp(['===== ' fname{i} ' =====']);
% %     cd([dirName fname{i}]);
%     
%     % MAT file: RSG prior
%     fid=fname{i}(end-5:end);
%     
%     if i==1
%         load([behDir 'G_20' fid '.mat']);
%         Ttmp=T;
%         % dealing with noSet
%         acquisitionTime=[];
%         tStart=[];
%         nCorrectNoSet=[];
%         pNoSet=[];
%         idNoSet=[];
%         nNoSet=[];
%         idSuccessNoSet=[];
%     else
%         beh=load([behDir 'G_20' fid '.mat']);
%         fn=fieldnames(beh);
%         for j=1:length(fn)
%             if size(beh.([fn{j}]),1)>=length(beh.T)
%                 eval([fn{j} '=[' fn{j} '; beh.([fn{j}]) ];']);
%             end
%         end
%         Ttmp=beh.T;
%     end
%     sessId=[sessId; repmat(str2num(fid),length(Ttmp),1)];
% end
% % % save([behDir 'G_RSGprior_DMFC.mat'],'sessId','-append'); 
% save([behDir 'G_RSGprior_DMFC.mat'],'T','t','radius','ballAlpha','idShortTrial','winF','theta','bonusRewDur','rewardDur','fixTimeDur','targetTimeDur','idHandEye','dtProdRew','FeedbackTimeDur','iti','timeout','sessId');
% 
% % noSet: pNoSet idNoSet idSuccessNoSet
% idSessNoSet=false(size(T)); idSessNoSet((end-(length(idNoSet)-1)):end)=true;
% T=T(idSessNoSet); 
% t=t(idSessNoSet);   
% radius=radius(idSessNoSet);     
% ballAlpha= ballAlpha(idSessNoSet);     
% idShortTrial=idShortTrial(idSessNoSet);      
% winF= winF(idSessNoSet);     
% theta=  theta(idSessNoSet);    
% bonusRewDur= bonusRewDur(idSessNoSet);     
% rewardDur = rewardDur(idSessNoSet);    
% fixTimeDur =  fixTimeDur(idSessNoSet);   
% targetTimeDur  = targetTimeDur(idSessNoSet);   
% idHandEye= idHandEye(idSessNoSet);     
% dtProdRew = dtProdRew(idSessNoSet);    
% FeedbackTimeDur=  FeedbackTimeDur(idSessNoSet);    
% iti= iti(idSessNoSet);     
% timeout=timeout(idSessNoSet); 
% sessId=sessId(idSessNoSet); 
% % save([behDir 'G_RSGprior_DMFC_noSet.mat'],'sessId','-append'); 
% save([behDir 'G_RSGprior_DMFC_noSet.mat'],'T','t','radius','ballAlpha','idShortTrial','winF','theta','bonusRewDur','rewardDur','fixTimeDur','targetTimeDur','idHandEye','dtProdRew','FeedbackTimeDur','iti','timeout',...
%     'pNoSet','idNoSet','idSuccessNoSet');
% clear idNoSet;
% 
% % 
% % % save([behDir 'H_RSGprior_DMFC.mat'],fn{1});
% % % for j=2:length(fn)
% % %     save([behDir 'H_RSGprior_DMFC.mat'],fn{j},'-append');
% % % end
% % % %  idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)

%% check across-session variability of wm,wp,offset (wFit w/o outlier removal?)
% figure; nMeas=4; 
% for iS=1:length(fname)
%     load([behDir 'G_' fname{iS}(end-7:end) '.mat']); % T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)
%     subplot(nMeas,1,1);plot(iS,wFit.w_m,'o'); ylabel('w_m');hold all;
%     subplot(nMeas,1,2);plot(iS,wFit.w_p,'o'); ylabel('w_p');hold all;
%     subplot(nMeas,1,3);plot(iS,wFit.offset1,'o'); ylabel('offset1');hold all;
%     subplot(nMeas,1,4);plot(iS,wFit.offset2,'o'); ylabel('offset2');hold all;
% end

%% outlier removal for each session
% also prior-, modality-, direction-, ts-specific

% idOut=[]; rangeOut=[0 3]; 
% pOut=[];
% idPlot=0; % 1;
% for iS=1:nS
%     fid=fname{iS}(end-5:end);
%     load([behDir 'G_20' fid '.mat']);% T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)
%     idOutS=true(size(T)); % assumig all outliers (inc. t<0)
%     for i=1:nPr
%         for j=1:nEH
%             for k=1:nTarg
%                 for m=1:length(Tmat{i})        
% 
% %  setting noSet trials as all outlier (N.B. tp in those trials are not <0  always)
% if exist('idNoSet')
%     id=idShortTrial==(2-i) &... % idShortTrial
%         idHandEye==(2-j) &...  % idHandEye
%         theta==(k-1)*180 &... % target location
%         T==Tmat{i}(m) &...
%         t>0 & t<3*Tmat{i}(m) &...
%         ~idNoSet; % only valid trials for now
% else
%     id=idShortTrial==(2-i) &... % idShortTrial
%         idHandEye==(2-j) &...  % idHandEye
%         theta==(k-1)*180 &... % target location
%         T==Tmat{i}(m) &...
%         t>0 & t<3*Tmat{i}(m); % only valid trials for now
% end
%                     u0=rangeOut*Tmat{i}(m); %[0.1 1.9]*overlapT; % [0.5 1.5]*overlapT; %[0 2]*overlapT; % with least outlier removal; or max(t)
%                     [tClean,idClean,pOutTmp,tOut]=removeOutlier(t(id),0,u0); %nSDremoveOutlier);
%                     idOutS(id)=~idClean;
%                     pOut=[pOut;mean(~idClean)];
%                     if idPlot
%                         figure; set(gcf,'position',pplot.(['rect' num2str(k) '_' num2str(m)]));hold all;
%                         histfit(tClean,35); if ~isempty(tOut), plot(tOut,0,'rx'); end; axis tight;
%                         title([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', ts' num2str(Tmat{i}(m)) ', p(outlier): ' num2str(mean(~idClean)*100) '%']);
%                         disp([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', ts' num2str(Tmat{i}(m)) ', p(outlier): ' num2str(mean(~idClean)*100) '%']);
%                         if mean(~idClean)>10
%                             disp('');
%                         end
%                     end                                        
%                 end % m ts
%             end % k targert
%             if idPlot, waitforbuttonpress; close all; end;
%         end % j EH
%     end % prior
%     idOut=[idOut; idOutS];
% end % for iS=1:nS
% save([behDir 'G_RSGprior_DMFC.mat'],'idOut','rangeOut','pOut','-append');
% 
% % noSet: pNoSet idNoSet idSuccessNoSet
% idOut=idOut(idSessNoSet); 
% save([behDir 'G_RSGprior_DMFC_noSet.mat'],'idOut','rangeOut','pOut','-append');
% 
% % load([behDir 'G_RSGprior_DMFC.mat'],'T','t','idShortTrial','idHandEye','theta');

load([behDir 'G_RSGprior_DMFC.mat']); % T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)

%% 2019/3/25
% show deltaTp for sigmoid now between two priors
% 2019/1/10: show behavior is sigmoidal by mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
% 2019/1/10: model comparison with linear model fit
% getting BIC LogPtp_Meas

% 2019/2/18: 
% model-based comparison b/t linear & BLS
% check var(tp|extreme ts)< E(var(tp|priorMean))
% opacity~# trials for model-free (idUseOp)
% condition-specific marker; red for short, blue for long

idChk=0; % 1;
idChkModel= 0; % 1;
idSep=0; % 1; % separately for shorter and longer ts
idUseOp=1;
sidUni=unique(sessId); % 170823
nSess=length(sidUni);

tmpMarker={'o';'o'}; % fliplr({'^','o','d','s'}); % {'^','o','d','s'}; % ER>EL>HR>HL
cmapTmp=[1 0 0;0 0 1]; % [0 .6 0;1 .5 0;0 0 1;1 0 0];
markersize=2; % 4;
% optsExpFig.Width=3; % 10/2.54;
% optsExpFig.Height=3; % 7.3/2.54;

if idSep
%     dTp=nan(nSess,nPr,nEH,nTarg,2,2); % last for middle vs extreme
%     nTp=nan(nSess,nPr,nEH,nTarg,2); % last for mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
else
    dTp=nan(nSess,nPr,nEH,nTarg,2); % last for mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
    nTp=nan(nSess,nPr,nEH,nTarg); % last for mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
    cid=nan(nSess,nPr,nEH,nTarg); % 12345678 pr>EH>RL
end

RMSE=nan(nSess,nPr,nEH,nTarg,2); % BLS linear
vTp=nan(nSess,nPr,nEH,nTarg,nTspp);
mTp=nan(nSess,nPr,nEH,nTarg,nTspp);
resLM=nan(nSess,nPr,nEH,nTarg,nTspp);

if idUseOp
    minOp=0; % 0.8; % 9; % linearly scaled by white and original color (max#trials>original; min#trial>c+(1-c)*0.9
else
    minMS=1; % min marker size
    maxMS=10; % max marker size
end

for iS=length(sidUni):(-1):1 % 1:length(sidUni)
    disp(['===== ' num2str(sidUni(iS)) ' =====']);
    idS=sessId==sidUni(iS);
    
    for i=1:nPr % 1 short 2 long
        for j=1:nEH 
            for l=1:nTarg
                
                id=idShortTrial==(2-i) &... % idShortTrial
                    idHandEye==(2-j) &...            % eye hand
                    theta==180*(l-1) &... % right left
                    idS &... % session
                    ~idOut; % outlier
                
                if idChk
                    idPlot=1; plotTpTs(T(id),t(id),idPlot); setFigPos(1,1);
                    %                     waitforbuttonpress; close;
                end
                
                % get mean tp for each ts
                muTp=nan(length(Tmat{i}),1);
                varTp=nan(length(Tmat{i}),1);
                for k=1:length(Tmat{i})
                    tmpId=id & T==Tmat{i}(k);
                    muTp(k)=mean(t(tmpId));
                    varTp(k)=std(t(tmpId));
%                     disp(nnz(tmpId)); % # data points
                end
                
                % model-based comparison b/t linear & BLS
                wFitTmp=wFitSessCond(iS,j,l);
                nParam=0; % 2;
                % no fitting
                sqErr=@(p,q) sum((muTp(:)-(BLS(q,p(1),[min(q) max(q)],'uniform')+p(2))).^2); % p: wm, offset, ts
                SE0=sqErr([wFitTmp.w_m,wFitTmp.(['offset' num2str(3-i)])],Tmat{i}(:));

                RMSE(iS,i,j,l,1)=sqrt(SE0/(length(Tmat{i})-nParam)); % SE>RMSE (N.B. now taking into account degree of freedom)
                LM=fitlm(Tmat{i},muTp,'linear');
                resLM(iS,i,j,l,:)=table2array(LM.Residuals(:,1));
                RMSE(iS,i,j,l,2)=sqrt(mean((muTp-LM.predict).^2)); % LM.RMSE; N.B. LinearModel's RMSE taking into account degree of freedom
%                 disp(['RMSE(BLS vs linear): ' num2str(RMSE(iS,i,j,l,:))]);
                
                % check var(tp|extreme ts)< E(var(tp|priorMean))
                mTp(iS,i,j,l,:)=muTp;
                vTp(iS,i,j,l,:)=varTp;
                
                if idSep
%                     dTp(iS,i,j,l,1,1)=diff(muTp(2:3));
%                     dTp(iS,i,j,l,2,1)=diff(muTp(3:4));
%                     dTp(iS,i,j,l,1,2)=diff(muTp(1:2));
%                     dTp(iS,i,j,l,2,2)=diff(muTp(4:5));
%                     nTp(iS,i,j,l,1)=nnz(id & T<=median(Tmat{i}));
%                     nTp(iS,i,j,l,2)=nnz(id & T>=median(Tmat{i}));
                else
                    dTp(iS,i,j,l,1)=mean(diff(muTp(2:4)));
                    dTp(iS,i,j,l,2)=mean([diff(muTp(1:2)) diff(muTp(4:5))]);
                    nTp(iS,i,j,l)=nnz(id);
                    cid(iS,i,j,l)=l+(j-1)*nTarg+(i-1)*nEH*nTarg;
                end
                
                if idChk
                    figure; setFigPos(2,1); ha;
                    plot(diff(muTp(3:4)),diff(muTp(4:5)),'o'); 
                    plot(diff(muTp(2:3)),diff(muTp(1:2)),'o'); 
%                     axis tight; 
                    plotIdentity(gca);
                    waitforbuttonpress; close all;
                end
                
            end % l targ
        end % j EH
    end % i pr
end % sess

nTp0=nTp; nTp=nTp(:); tmpMS=(nTp-min(nTp))/(max(nTp)-min(nTp)); % [0 1]
cid=cid(:);
disp(['# trials: ' num2str(max(nTp)) ', ' num2str(min(nTp))]);

figure; ha; setFigPos(1,2); % scatter plot for model free
if idSep
%     tmpX=dTp(:,:,:,:,:,1); tmpY=dTp(:,:,:,:,:,2);
%     tmpX=tmpX(:);tmpY=tmpY(:);
%     nTp=nTp(:); tmpMS=(nTp-min(nTp))/(max(nTp)-min(nTp)); % [0 1]
else
    tmpX=dTp(:,:,:,:,1); tmpX=tmpX(:);
    tmpY=dTp(:,:,:,:,2); tmpY=tmpY(:);
%     nTp=nTp(:); tmpMS=(nTp-min(nTp))/(max(nTp)-min(nTp)); % [0 1]
%     cid=cid(:);
end
% disp(['# trials: ' num2str(max(nTp)) ', ' num2str(min(nTp))]);
for i=1:length(nTp)
    if idUseOp
        % condition<prior-specific marker; red for short, blue for long
        idShortTmp=cid(i)<=(nEH*nTarg); % 1 for short, 0 for long
        cidTmp=cid(i)-(1-idShortTmp)*nEH*nTarg; % 1234 for ER EL HR HL
        tmpC=cmapTmp(2-idShortTmp,:);
%         tmpC=pplot.cmap{3-2*idShortTmp}; % red for short, blue for long
        tmpC=tmpC+(1-tmpC)*minOp*(1-tmpMS(i));
        plot(tmpX(i),tmpY(i),tmpMarker{2-idShortTmp},... % fliplr({'^','o','d','s'});
            'color',tmpC,'markersize',markersize,'linewidth',0.5); 
    else % idUseOp
        plot(tmpX(i),tmpY(i),'ko','markersize',tmpMS(i)*(maxMS-minMS)+minMS);
    end
end
axis tight; 
xlim([20 125]);ylim([20 125]); set(gca,'xtick',30:30:120,'ytick',30:30:120,'tickdir','out'); 
plotIdentity(gca);
ylabel('difference in mean t_p for boundary t_s');
xlabel('difference in mean t_p for middle t_s');
remTickLabel;remLabel;
applytofig(gcf,optsExpFig);

[p,h]=signrank(tmpX,tmpY); disp(['sigmoid: pooled across sessions,prior,effector,direction: ' num2str(p,3)]);

% stat for each condition; pooled across prior
for i=1:(length(unique(cid))/nPr)
    idCondTmp=i==(mod(cid-1,4)+1); % 12341234
    [p,h]=signrank(tmpX(idCondTmp),tmpY(idCondTmp)); disp(['cond' num2str(i) ': ' num2str(p,3)]);
end
% for i=1:length(unique(cid))
%     idCondTmp=i==cid;
%     [p,h]=signrank(tmpX(idCondTmp),tmpY(idCondTmp)); disp(['cond' num2str(i) ': ' num2str(p,3)]);
% end

% stat separately for short/long
idShortTmp=cid<=(nEH*nTarg);
[p,h]=signrank(tmpX(idShortTmp),tmpY(idShortTmp)); disp(['short: ' num2str(p,3)]);
[p,h]=signrank(tmpX(~idShortTmp),tmpY(~idShortTmp)); disp(['long: ' num2str(p,3)]);

% histogram: now condition-specific<prior-specific
figure; ha; setFigPos(1,3);
% hHist=histogram(tmpY-tmpX,25);
% set(hHist,'Facecolor','w');
% plotVertical(gca,0,[]);
% ylabel('# data sets');
% xlabel('degree of being sigmoidal');
nBin=15;barWidth=0.4;markersize=10;

if idPool
    dataH=load('/Users/hansem/Dropbox (MIT)/figuresRSG2prior/revision/deltaTp_H.mat'); % close;
    tmpY=[dataH.tmpY; tmpY];
    tmpX=[dataH.tmpX; tmpX];
    cid=[dataH.cid; cid];    
end
xHist=linspace(min(tmpY-tmpX),max(tmpY-tmpX),nBin); % getting bin locations
idShortTmp=cid<=(nEH*nTarg); % 1 for short, 0 for long

% plot
diffYX=cell(nPr,1);
for i=1:nPr   
%      idCondTmp=i==(mod(cid-1,4)+1); % 12341234
%     diffXY{i}=-(tmpY(idCondTmp)-tmpX(idCondTmp));
    diffYX{i}=tmpY(i==(2-idShortTmp))-tmpX(i==(2-idShortTmp)); % 2-idShortTmp: 1 for short, 2 for long
    hHist=histogram(-diffYX{i},xHist);
     Xval    = hHist.BinEdges + hHist.BinWidth*0.5 +hHist.BinWidth*(i-(nPr+1)/2)/2; % (i-1)/(nEH*nTarg)+1/(nEH*nTarg*2); % (i-(nPr+1)/2)/2; % binWidth*0.25 for short, 0.75 for long
     Xval    = Xval(1:end-1);
     Yval    = hHist.Values;
     delete(hHist);

        tmpC=cmapTmp(i,:); % 1234 for ER EL HR HL
     hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',tmpC,'EdgeColor','none'); %  pplot.cmap{2*i-1}
end
axis tight; 
maxYlim=max(ylim);
% plot mean
for i=1:nPr
         tmpC=cmapTmp(i,:); % 1234 for ER EL HR HL
    plot(mean(-diffYX{i}),0,'color',tmpC,'markerfacecolor',tmpC,'marker','o','markersize',markersize); %  pplot.cmap{2*i-1}
end % for i=1:nPr

% plot model prediction
 tmpShort=diff(BLS(480:80:800,0.05,[480 800],'uniform'));
tmpLong=diff(BLS(800:100:1200,0.05,[800 1200],'uniform'));
tmpLong=-mean(tmpLong([1 end]))+mean(tmpLong([2 3]));
tmpShort=-mean(tmpShort([1 end]))+mean(tmpShort([2 3]));
plot(tmpShort,1.1*maxYlim,'color',cmapTmp(1,:),'markerfacecolor','w','marker','v','markersize',markersize); %  pplot.cmap{2*i-1}
plot(tmpLong,1.1*maxYlim,'color',cmapTmp(2,:),'markerfacecolor','w','marker','v','markersize',markersize); %  pplot.cmap{2*i-1}

ylabel('# data sets');
xlabel('degree of being sigmoidal');
plotVertical(gca,0,[]);
set(gca,'tickdir','out','xtick',-60:30:60,'ytick',0:10:50); xlim([-73 73]);
remTickLabel;remLabel;
applytofig(gcf,optsExpFig);

% stat 19/4/12
signrank(diffYX{1},diffYX{2})

% return;

%% full probabilistic linear model fit: idFitLlinear, estLinear.m: everythign condition-speciifc (animal, priorX, effector, direction)
% redo for 161206 161207 161210 161221
% change estWmWpOld: nInitW: 6>15, nInit5>10, also increase offset10, offset20 (1.5>3, 0.5>1.5)

if idFitLlinear
    
    sidUni=unique(sessId);
    bicLinearSessCond=nan(length(sidUni),nEH,nTarg); % ,nPr
    if exist('wFitSessCond') % redo
        wFitSessCond0=wFitSessCond;
    end
    wFitLinearSessCond=[];
    for iS=length(sidUni):(-1):1 % 1:length(sidUni)
        disp(['===== ' num2str(sidUni(iS)) ' =====']);
        idS=sessId==sidUni(iS);
        
        % redo
%         if sidUni(iS)==161206 | sidUni(iS)==161207 | sidUni(iS)==161210 | sidUni(iS)==161221 % redo for 161206 161207 161210 161221
            
            %     for i=1:nPr
            for j=1:nEH
                for l=1:nTarg
                    
                    % redo
                    disp([wFitSessCond0(iS,j,l).w_m(:) wFitSessCond0(iS,j,l).w_p(:) wFitSessCond0(iS,j,l).offset1(:) wFitSessCond0(iS,j,l).offset2(:)]);
                    
                    id=... % idShortTrial==(2-i) &... % idShortTrial
                        idHandEye==(2-j) &...            % eye hand
                        theta==180*(l-1) &... % right left
                        idS &... % session
                        ~idOut; % outlier
                    
                    try
                        %                     [wFitTmp,bicTmp]=estWmWpOld(t(id),T(id),true(nnz(id),1),'mmse1Wm',[],-1,true); % display: off
                        [wFitTmp,bicTmp]=estLinear(t(id),T(id),idShortTrial(id),'mmse2offset',[],-1,true); % display: off
                    catch
                        disp('');
                    end
                    fieldnm=fieldnames(wFitTmp);
                    for iField=1:length(fieldnm)
                        wFitLinearSessCond(iS,j,l).(fieldnm{iField})=wFitTmp.(fieldnm{iField}); % mmse2offset: w_m w_p offset1 offset2 ,i
                    end
                    try
                        bicLinearSessCond(iS,j,l)=bicTmp;
                    catch
                        disp(bicTmp);
                    end
                    disp([wFitTmp.w_m(:) wFitTmp.w_p(:) wFitTmp.sS(:) wFitTmp.sL(:) wFitTmp.b0S(:) wFitTmp.b0L(:)]);
                end % l targ
            end % j EH
            %     end % i pr
            
            % redo
%         end % if sidUni(iS)==161206 | sidUni(iS)==161207 | sidUni(iS)==161210 | sidUni(iS)==161221
        
    end % serss
    save([behDir 'G_RSGprior_DMFC.mat'],'wFitLinearSessCond','bicLinearSessCond','-append');

    % 2019/2/25: stop in the middle
%     ===== 170823 =====
%     0.0392    0.0857  -70.5369  -22.2766
%     0.0750    0.0028    0.8905    0.8806   50.8370   67.3359
%     0.0488    0.0567  -63.4318  -18.8874
%     0.0777    0.0040    0.8694    0.7864   64.3520  151.2879
%     0.0493    0.0493    0.0582    0.0582
%     0.0851    0.0036    0.8539    0.8387   77.9189  128.0329
%     0.0617    0.0520  -45.6227   -9.2026
%     0.0549    0.0512    0.8357    0.7042   96.8445  246.6639
% 
% ===== 170822 =====
%     0.0488    0.0488    0.0572    0.0572
%     0.0240    0.0656    0.8701    0.7930   62.8462  148.3389
%     0.0400    0.0400    0.0590    0.0590
%     0.0100    0.0662    0.9200    0.7957   25.1885  136.3587
%     0.0499    0.0636  -13.5581  -17.4881
%     0.0521    0.0588    0.8912    0.7408   51.3676  243.5534
%     0.0489    0.0599  -18.0817   -3.4540
%     0.0620    0.0452    0.9352    0.7122   37.8967  266.6718
% 
% ===== 170821 =====
%     0.0295    0.0636  -39.3385  -15.6528
%     0.0664    0.0054    0.8924    0.9129   53.6085   48.9749
%     0.0365    0.0630  -38.2195  -16.7220
%     0.0777    0.0042    0.9100    0.8691   37.9093   91.5339
%     0.0568    0.0657  -48.1563  -15.1653
%     0.0481    0.0663    0.8565    0.7187   74.4572  232.0722
%     0.0372    0.0372    0.0710    0.0710
%     0.0843    0.0032    0.8709    0.8722   58.8835   98.4421
% 
% ===== 170818 =====
%     0.0652    0.0895  -40.0621  -26.3505
%     0.0957    0.0037    0.8066    0.8244  100.1485  154.6531
%     0.0459    0.0601  -24.6619  -30.6169
%     0.0789    0.0035    0.8329    0.8492   75.5104  122.8491
%     0.0613    0.0605   -7.8150   13.3173
%     0.0105    0.0719    0.8800    0.6624   91.6247  327.3944
%     0.0654    0.0551  -24.4207    2.7754
    
% check result
figure; setFigPos(1,1); h=histogram(bicLinearSessCond(:)-bicSessCond(:),20); set(h,'facecolor','w');
xlabel('BIC(linear model)-BIC(Bayesian model)'); ylabel('# data sets'); plotVertical(gca,0,[]);
figure; setFigPos(2,1); hist([wFitLinearSessCond.w_m],25); xlabel('w_m');
figure; setFigPos(2,2); hist([wFitLinearSessCond.w_p],25); xlabel('w_p');
figure; setFigPos(2,3); hist([wFitLinearSessCond.sS],25); xlabel('slope(short)');
figure; setFigPos(2,4); hist([wFitLinearSessCond.sL],25); xlabel('slope(long)');
figure; setFigPos(2,5); hist([wFitLinearSessCond.b0S],25); xlabel('intercept(short)');
figure; setFigPos(2,6); hist([wFitLinearSessCond.b0L],25); xlabel('intercept(long)');


return;


end % if idFitLlinear
%% checking condition-speciifc BLS fit: mean +- SD across sessions
% wFitSessCond [17 x 2EH x 2RL]. w_m w_p offset1 offset2
for iEH=1:nEH
    for iRL=1:nTarg
        [m,s]=meanSD([wFitSessCond(:,iEH,iRL).w_m]); % wm
        disp([num2str(m,3) '  ' num2str(s,3)]);
        [m,s]=meanSD([wFitSessCond(:,iEH,iRL).w_p]); % wp
        disp([num2str(m,3) '  ' num2str(s,3)]);
        [m,s]=meanSD([wFitSessCond(:,iEH,iRL).offset2]); % offset2
        disp([num2str(m,3) '  ' num2str(s,3)]);
        [m,s]=meanSD([wFitSessCond(:,iEH,iRL).offset1]); % offset1 for long
        disp([num2str(m,3) '  ' num2str(s,3)]);
        disp('-----');
    end
end


%% rank sum test for overalp: ignoring EH/RL
idS=T==800 & idShortTrial==1 & t>0 & t<3*T & ~idOut;
idL=T==800 & idShortTrial==0 & t>0 & t<3*T & ~idOut;
[p,h,stat]=ranksum(t(idS),t(idL),'tail','left')

% p =
% 
%    5.9384e-76
% 
% 
% h =
% 
%   logical
% 
%    1
% 
% 
% stat = 
% 
%   struct with fields:
% 
%        zval: -18.4055
%     ranksum: 6952806
% 
% [mean(t(idS)) mean(t(idL))]
% 
% ans =
% 
%   775.0541  808.8846
  
%% 2019/1/10: show behavior is sigmoidal by mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
% 2019/1/10: model comparison with linear model fit
% getting BIC LogPtp_Meas

% 2019/2/18: 
% model-based comparison b/t linear & BLS
% check var(tp|extreme ts)< E(var(tp|priorMean))
% opacity~# trials for model-free (idUseOp)
% condition-specific marker; red for short, blue for long

idChk=0; % 1;
idChkModel= 0; % 1;
idSep=0; % 1; % separately for shorter and longer ts
idUseOp=1;
sidUni=unique(sessId); % 170823
nSess=length(sidUni);

tmpMarker=fliplr({'^','o','d','s'}); % {'^','o','d','s'}; % ER>EL>HR>HL
cmapTmp=[0 .6 0;1 .5 0;0 0 1;1 0 0];
markersize=2; % 4;
% optsExpFig.Width=3; % 10/2.54;
% optsExpFig.Height=3; % 7.3/2.54;

if idSep
%     dTp=nan(nSess,nPr,nEH,nTarg,2,2); % last for middle vs extreme
%     nTp=nan(nSess,nPr,nEH,nTarg,2); % last for mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
else
    dTp=nan(nSess,nPr,nEH,nTarg,2); % last for mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
    nTp=nan(nSess,nPr,nEH,nTarg); % last for mean(ts2-ts1,ts5-ts4) vs mean(ts4-ts3,ts3-ts2)
    cid=nan(nSess,nPr,nEH,nTarg); % 12345678 pr>EH>RL
end

RMSE=nan(nSess,nPr,nEH,nTarg,2); % BLS linear
vTp=nan(nSess,nPr,nEH,nTarg,nTspp);
mTp=nan(nSess,nPr,nEH,nTarg,nTspp);
resLM=nan(nSess,nPr,nEH,nTarg,nTspp);

if idUseOp
    minOp=0; % 0.8; % 9; % linearly scaled by white and original color (max#trials>original; min#trial>c+(1-c)*0.9
else
    minMS=1; % min marker size
    maxMS=10; % max marker size
end

for iS=length(sidUni):(-1):1 % 1:length(sidUni)
    disp(['===== ' num2str(sidUni(iS)) ' =====']);
    idS=sessId==sidUni(iS);
    
    for i=1:nPr % 1 short 2 long
        for j=1:nEH 
            for l=1:nTarg
                
                id=idShortTrial==(2-i) &... % idShortTrial
                    idHandEye==(2-j) &...            % eye hand
                    theta==180*(l-1) &... % right left
                    idS &... % session
                    ~idOut; % outlier
                
                if idChk
                    idPlot=1; plotTpTs(T(id),t(id),idPlot); setFigPos(1,1);
                    %                     waitforbuttonpress; close;
                end
                
                % get mean tp for each ts
                muTp=nan(length(Tmat{i}),1);
                varTp=nan(length(Tmat{i}),1);
                for k=1:length(Tmat{i})
                    tmpId=id & T==Tmat{i}(k);
                    muTp(k)=mean(t(tmpId));
                    varTp(k)=std(t(tmpId));
%                     disp(nnz(tmpId)); % # data points
                end
                
                % model-based comparison b/t linear & BLS
                wFitTmp=wFitSessCond(iS,j,l);
                nParam=0; % 2;
                % no fitting
                sqErr=@(p,q) sum((muTp(:)-(BLS(q,p(1),[min(q) max(q)],'uniform')+p(2))).^2); % p: wm, offset, ts
                SE0=sqErr([wFitTmp.w_m,wFitTmp.(['offset' num2str(3-i)])],Tmat{i}(:));
%                 [wm,offset,SE]=LSEfitBLS(muTp,Tmat{i},wFitTmp.w_m,wFitTmp.(['offset' num2str(3-i)])); % offset1 for long
%                 disp(['wm: ' num2str(wFitTmp.w_m) ' vs ' num2str(wm)]);
%                 disp(['offset: ' num2str(wFitTmp.(['offset' num2str(3-i)])) ' vs ' num2str(offset)]);
%                 if idChkModel
%                     figure; setFigPos(2,3); ha;  plot(Tmat{i},BLS(Tmat{i},wFitTmp.w_m,[min(Tmat{i}) max(Tmat{i})],'uniform')+wFitTmp.(['offset' num2str(3-i)]),'bo-'); plot(Tmat{i},muTp,'ro-'); 
%                     plot(Tmat{i},BLS(Tmat{i},wm,[min(Tmat{i}) max(Tmat{i})],'uniform')+offset,'go-'); legend('old BLS','raw data','new LSE BLS fit','location','northwest'); plotIdentity(gca);
%                     applytofig4keynote;
%                     waitforbuttonpress; close;
%                 end
                RMSE(iS,i,j,l,1)=sqrt(SE0/(length(Tmat{i})-nParam)); % SE>RMSE (N.B. now taking into account degree of freedom)
                LM=fitlm(Tmat{i},muTp,'linear');
                resLM(iS,i,j,l,:)=table2array(LM.Residuals(:,1));
                RMSE(iS,i,j,l,2)=sqrt(mean((muTp-LM.predict).^2)); % LM.RMSE; N.B. LinearModel's RMSE taking into account degree of freedom
%                 disp(['RMSE(BLS vs linear): ' num2str(RMSE(iS,i,j,l,:))]);
                
                % check var(tp|extreme ts)< E(var(tp|priorMean))
                mTp(iS,i,j,l,:)=muTp;
                vTp(iS,i,j,l,:)=varTp;
                
                if idSep
%                     dTp(iS,i,j,l,1,1)=diff(muTp(2:3));
%                     dTp(iS,i,j,l,2,1)=diff(muTp(3:4));
%                     dTp(iS,i,j,l,1,2)=diff(muTp(1:2));
%                     dTp(iS,i,j,l,2,2)=diff(muTp(4:5));
%                     nTp(iS,i,j,l,1)=nnz(id & T<=median(Tmat{i}));
%                     nTp(iS,i,j,l,2)=nnz(id & T>=median(Tmat{i}));
                else
                    dTp(iS,i,j,l,1)=mean(diff(muTp(2:4)));
                    dTp(iS,i,j,l,2)=mean([diff(muTp(1:2)) diff(muTp(4:5))]);
                    nTp(iS,i,j,l)=nnz(id);
                    cid(iS,i,j,l)=l+(j-1)*nTarg+(i-1)*nEH*nTarg;
                end
                
                if idChk
                    figure; setFigPos(2,1); ha;
                    plot(diff(muTp(3:4)),diff(muTp(4:5)),'o'); 
                    plot(diff(muTp(2:3)),diff(muTp(1:2)),'o'); 
%                     axis tight; 
                    plotIdentity(gca);
                    waitforbuttonpress; close all;
                end
                
            end % l targ
        end % j EH
    end % i pr
end % sess

nTp0=nTp; nTp=nTp(:); tmpMS=(nTp-min(nTp))/(max(nTp)-min(nTp)); % [0 1]
cid=cid(:);
disp(['# trials: ' num2str(max(nTp)) ', ' num2str(min(nTp))]);

% % model-based comparison b/t linear & BLS
% figure; ha; setFigPos(2,1); % scatter plot for model free
% tmpX=RMSE(:,:,:,:,1); tmpX=tmpX(:);
% tmpY=RMSE(:,:,:,:,2); tmpY=tmpY(:);
% for i=1:length(nTp)
%     if idUseOp
%         % condition-specific marker; red for short, blue for long
%         idShortTmp=cid(i)<=(nEH*nTarg); % 1 for short, 0 for long
%         tmpC=pplot.cmap{3-2*idShortTmp}; % red for short, blue for long
%         tmpC=tmpC+(1-tmpC)*minOp*(1-tmpMS(i));
%         plot(tmpX(i),tmpY(i),tmpMarker{cid(i)-(1-idShortTmp)*nEH*nTarg},...
%             'color',tmpC,'markersize',markersize,'linewidth',0.5);
%     else % idUseOp
%         plot(tmpX(i),tmpY(i),'ko','markersize',tmpMS(i)*(maxMS-minMS)+minMS);
%     end
% end
% axisEqual; % axis tight; 
% % xlim([20 125]);ylim([20 125]); set(gca,'xtick',30:30:120,'ytick',30:30:120,'tickdir','out'); 
% plotIdentity(gca);
% ylabel('RMSE(linear model)');
% xlabel('RMSE(Bayesian model)');
% remTickLabel;remLabel;
% applytofig(gcf,optsExpFig);
% % figure; h=histogram(tmpX-tmpY,25); set(h,'facecolor','w'); applytofig4keynote

% % check linear model's residual
% figure; ha; setFigPos(2,1); % scatter plot for model free
% for i=1:size(resLM,1)
%     for j=1:size(resLM,3) % EH
%         for k=1:size(resLM,4) % targ
%             for iPr=1:size(resLM,2)
%             plot(Tmat{iPr}(:),squeeze(resLM(i,iPr,j,k,:)),'-','color',pplot.cmap{2*(iPr-1)+1});
%             end
%         end
%     end
% end
% for iPr=1:size(resLM,2) % plot mean
%     tmpY=squeeze(resLM(:,iPr,:,:,:));
%     tmpY=reshape(tmpY,size(resLM,1)*size(resLM,3)*size(resLM,4),size(resLM,5));
%     shadedErrorBar(Tmat{iPr}(:),mean(tmpY,1),sem(tmpY,1),{'-','color',pplot.cmap{2*(iPr-1)+1},'linewidth',2});
% %     plot(Tmat{iPr}(:),squeeze(mean(mean(mean(tmpY,1),3),4)),'-','color',pplot.cmap{2*(iPr-1)+1},'linewidth',2);
% end
% axis tight; 
% plotHorizon(gca,0,[]);
% xlabel('t_s (ms)');
% ylabel('Residual of linear model');
% % figure; histogram(tmpY(:,1),15);

% % check var(tp|extreme ts)< E(var(tp|priorMean))
% tmpX=mTp(:); tmpY=vTp(:); nTp0=(nTp0-min(nTp))./(max(nTp)-min(nTp));
% % figure 1: sigma vs mu (only extreme ts & prior mean)
% figure; ha; setFigPos(2,2);
% for i=1:nPr
%     for j=1:(nEH*nTarg)
%         iEH=floor((j-1)/nTarg)+1; iTarg=mod(j-1,nTarg)+1;
%         for k=1:2:5 % ts
%             for iS=1:size(mTp,1)
%                 if idUseOp
%                     tmpC=tmpCmap{i,1}(k,:); %
%                     tmpC=tmpC+(1-tmpC)*minOp*(1-nTp0(iS,i,iEH,iTarg));
%                     plot(mTp(:,i,iEH,iTarg,k),vTp(:,i,iEH,iTarg,k),tmpMarker{j},'color',tmpC,'markersize',markersize,'linewidth',0.5);
%                 else
%                     plot(mTp(:,i,iEH,iTarg,k),vTp(:,i,iEH,iTarg,k),tmpMarker{j},'color',tmpCmap{i,1}(k,:),'markersize',markersize,'linewidth',0.5);
%                 end
%             end % iS
%         end % k
%     end % j
% end % i
% for i=1:nPr % plot mean
%     for k=5:(-2):1 % 1:2:5 % ts
%             tmpX=squeeze(mTp(:,i,:,:,k)); tmpX=tmpX(:);
%             tmpY=squeeze(vTp(:,i,:,:,k)); tmpY=tmpY(:);
%             if k==3 % plot line for reference in prior mean
%                 plot(Tmat{i},Tmat{i}*mean(tmpY)/mean(tmpX),':','color',pplot.cmap{2*(i-1)+1});
%             end
%             plotErrorbar(mean(tmpX),std(tmpX),mean(tmpY),std(tmpY),tmpCmap{i,1}(k,:),2); % last for linewdith
%             plot(mean(tmpX),mean(tmpY),'o','markersize',12,'linewidth',2,'color',tmpCmap{i,1}(k,:),'markerfacecolor','w');
%     end
% end
% % axisEqual; % axis tight; 
% % xlim([20 125]);ylim([20 125]); 
% set(gca,'xtick',unique([Tmat{1}(1:2:end)' Tmat{2}(1:2:5)']),'tickdir','out');  % ,'ytick',30:30:120,'tickdir','out'); 
% xlabel('\mu(t_p|t_s)');
% ylabel('\sigma(t_p|t_s)');
% remTickLabel;remLabel;
% applytofig(gcf,optsExpFig);
% % figure 2: observed sigma vs expected sigma
% figure; ha; setFigPos(2,3);
% tmpV=nan(nPr,nTspp-3,2,size(vTp,1)*size(vTp,3)*size(vTp,4)); % last for obs/exp
% iTmpV=1;
% for iS=1:size(vTp,1)
%     for i=1:size(vTp,2)
%         for j=1:size(vTp,3)
%             for k=1:size(vTp,4)
%                 cv=vTp(iS,i,j,k,3)/Tmat{i}(3); % mTp(iS,i,j,k,3); % based on ts, not mu(tp)
%                 tmpY=cv*Tmat{i}(1); tmpX=vTp(iS,i,j,k,1); tmpC=tmpCmap{i,1}(1,:); tmpC=tmpC+(1-tmpC)*minOp*(1-nTp0(iS,i,j,k)); % ts1
%                 tmpV(i,1,:,iTmpV)=[tmpX tmpY];
%                 plot(tmpX,tmpY,tmpMarker{(j-1)*nTarg+k},'color',tmpC,'markersize',markersize,'linewidth',0.5);
%                 tmpY=cv*Tmat{i}(5); tmpX=vTp(iS,i,j,k,5); tmpC=tmpCmap{i,1}(5,:); tmpC=tmpC+(1-tmpC)*minOp*(1-nTp0(iS,i,j,k)); % ts5
%                 tmpV(i,2,:,iTmpV)=[tmpX tmpY];
%                 plot(tmpX,tmpY,tmpMarker{(j-1)*nTarg+k},'color',tmpC,'markersize',markersize,'linewidth',0.5);
%                 iTmpV=iTmpV+1;
%             end
%         end
%     end
% end
% axisEqual; % axis tight; 
% % xlim([20 125]);ylim([20 125]); set(gca,'xtick',30:30:120,'ytick',30:30:120,'tickdir','out'); 
% plotIdentity(gca);
% xlabel('observed \sigma(t_p|t_s)');
% ylabel('expected \sigma(t_p|t_s)');
% remTickLabel;remLabel;
% applytofig(gcf,optsExpFig);
% % stat
% [p,h]=signrank(squeeze(tmpV(1,1,1,:)),squeeze(tmpV(1,1,2,:))); disp(p); % 480
% [p,h]=signrank(squeeze(tmpV(1,2,1,:)),squeeze(tmpV(1,2,2,:))); disp(p); % short800
% [p,h]=signrank(squeeze(tmpV(2,1,1,:)),squeeze(tmpV(2,1,2,:))); disp(p);
% [p,h]=signrank(squeeze(tmpV(2,2,1,:)),squeeze(tmpV(2,2,2,:))); disp(p);

figure; ha; setFigPos(1,2); % scatter plot for model free
if idSep
%     tmpX=dTp(:,:,:,:,:,1); tmpY=dTp(:,:,:,:,:,2);
%     tmpX=tmpX(:);tmpY=tmpY(:);
%     nTp=nTp(:); tmpMS=(nTp-min(nTp))/(max(nTp)-min(nTp)); % [0 1]
else
    tmpX=dTp(:,:,:,:,1); tmpX=tmpX(:);
    tmpY=dTp(:,:,:,:,2); tmpY=tmpY(:);
%     nTp=nTp(:); tmpMS=(nTp-min(nTp))/(max(nTp)-min(nTp)); % [0 1]
%     cid=cid(:);
end
% disp(['# trials: ' num2str(max(nTp)) ', ' num2str(min(nTp))]);
for i=1:length(nTp)
    if idUseOp
        % condition<prior-specific marker; red for short, blue for long
        idShortTmp=cid(i)<=(nEH*nTarg); % 1 for short, 0 for long
        cidTmp=cid(i)-(1-idShortTmp)*nEH*nTarg; % 1234 for ER EL HR HL
        tmpC=cmapTmp(cidTmp,:);
%         tmpC=pplot.cmap{3-2*idShortTmp}; % red for short, blue for long
        tmpC=tmpC+(1-tmpC)*minOp*(1-tmpMS(i));
        plot(tmpX(i),tmpY(i),tmpMarker{cidTmp},... % fliplr({'^','o','d','s'});
            'color',tmpC,'markersize',markersize,'linewidth',0.5); 
    else % idUseOp
        plot(tmpX(i),tmpY(i),'ko','markersize',tmpMS(i)*(maxMS-minMS)+minMS);
    end
end
axis tight; 
xlim([20 125]);ylim([20 125]); set(gca,'xtick',30:30:120,'ytick',30:30:120,'tickdir','out'); 
plotIdentity(gca);
ylabel('difference in mean t_p for boundary t_s');
xlabel('difference in mean t_p for middle t_s');
remTickLabel;remLabel;
applytofig(gcf,optsExpFig);

% stat for each condition; pooled across prior
for i=1:(length(unique(cid))/nPr)
    idCondTmp=i==(mod(cid-1,4)+1); % 12341234
    [p,h]=signrank(tmpX(idCondTmp),tmpY(idCondTmp)); disp(['cond' num2str(i) ': ' num2str(p,3)]);
end
% for i=1:length(unique(cid))
%     idCondTmp=i==cid;
%     [p,h]=signrank(tmpX(idCondTmp),tmpY(idCondTmp)); disp(['cond' num2str(i) ': ' num2str(p,3)]);
% end

% % stat separately for short/long
% idShortTmp=cid<=(nEH*nTarg);
% [p,h]=signrank(tmpX(idShortTmp),tmpY(idShortTmp)); disp(['short: ' num2str(p,3)]);
% [p,h]=signrank(tmpX(~idShortTmp),tmpY(~idShortTmp)); disp(['long: ' num2str(p,3)]);

% histogram: now condition-specific<prior-specific
figure; ha; setFigPos(1,3);
% hHist=histogram(tmpY-tmpX,25);
% set(hHist,'Facecolor','w');
% plotVertical(gca,0,[]);
% ylabel('# data sets');
% xlabel('degree of being sigmoidal');
nBin=15;barWidth=0.2;markersize=10;
xHist=linspace(min(tmpY-tmpX),max(tmpY-tmpX),nBin); % getting bin locations
idShortTmp=cid<=(nEH*nTarg); % 1 for short, 0 for long
% plot
diffXY=cell(nPr,1);
for i=1:(length(unique(cid))/nPr) % nPr   
     idCondTmp=i==(mod(cid-1,4)+1); % 12341234
    diffXY{i}=-(tmpY(idCondTmp)-tmpX(idCondTmp));
%     diffYX{i}=tmpY(i==(2-idShortTmp))-tmpX(i==(2-idShortTmp)); % 2-idShortTmp: 1 for short, 2 for long
    hHist=histogram(diffXY{i},xHist);
     Xval    = hHist.BinEdges + hHist.BinWidth*0.5 +hHist.BinWidth*(i-1)/(nEH*nTarg)+1/(nEH*nTarg*2); % (i-(nPr+1)/2)/2; % binWidth*0.25 for short, 0.75 for long
     Xval    = Xval(1:end-1);
     Yval    = hHist.Values;
     delete(hHist);

        tmpC=cmapTmp(i,:); % 1234 for ER EL HR HL
     hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',tmpC,'EdgeColor','none'); %  pplot.cmap{2*i-1}
end
axis tight; 
maxYlim=max(ylim);
% plot mean
for i=1:(length(unique(cid))/nPr) % nPr
         tmpC=cmapTmp(i,:); % 1234 for ER EL HR HL
    plot(mean(diffXY{i}),1.1*maxYlim,'color',tmpC,'markerfacecolor',tmpC,'marker','v','markersize',markersize); %  pplot.cmap{2*i-1}
end % for i=1:nPr
ylabel('# data sets');
xlabel('degree of being sigmoidal');
plotVertical(gca,0,[]);
set(gca,'tickdir','out','xtick',-60:30:60,'ytick',0:5:10); xlim([-50 50]);
remTickLabel;remLabel;
applytofig(gcf,optsExpFig);

% % check BLS predict higher sigmoid for long prior
% teS=BLS(480:80:800,0.0508,[480 800],'uniform');
% teL=BLS(800:100:1200,0.0508,[800 1200],'uniform');
% x=diff(teS);
% y=diff(teL);
% mean([x(1) x(4)])-mean([x(2) x(3)]) % -23.7240
% mean([y(1) y(4)])-mean([y(2) y(3)]) % -33.1630

% % scatter hist diff
% figure; ha; setFigPos(1,4);
% scatterHistDiff(tmpX,tmpY,[],[],[]);

% return;
  
%% BLS fit for each conditoin/sessions: everythign condition-speciifc (animal, prior, effector, direction)

if idFitBLS

    % conditon, pooled across sessions
    bicCond=nan(nEH,nTarg); % ,nPr
    wFitCond=[];
        
        %     for i=1:nPr
        for j=1:nEH
            for l=1:nTarg
                        disp(['===== ' ehNm{j} targNm{l}  ' =====']);

                id=... % idShortTrial==(2-i) &... % idShortTrial
                    idHandEye==(2-j) &...            % eye hand
                    theta==180*(l-1) &... % right left
                    ~idOut; % outlier % t<0 & t>3*T already removed
                
                try
                    %                     [wFitTmp,bicTmp]=estWmWpOld(t(id),T(id),true(nnz(id),1),'mmse1Wm',[],-1,true); % display: off
                    [wFitTmp,bicTmp]=estWmWpOld(t(id),T(id),idShortTrial(id),'mmse2offset',[],-1,true); % display: off
                catch
                    disp('');
                end
                fieldnm=fieldnames(wFitTmp);
                for iField=1:length(fieldnm)
                    wFitCond(j,l).(fieldnm{iField})=wFitTmp.(fieldnm{iField}); % mmse2offset: w_m w_p offset1 offset2 ,i
                end
                try
                    bicCond(j,l)=bicTmp;
                catch
                    disp(bicTmp);
                end
                disp([wFitTmp.w_m(:) wFitTmp.w_p(:) wFitTmp.offset1(:) wFitTmp.offset2(:)]);
            end % l targ
        end % j EH
        %     end % i pr
    save([behDir 'G_RSGprior_DMFC.mat'],'wFitCond','bicCond','-append');

    
% % Session
% sidUni=unique(sessId);
% bicSessCond=nan(length(sidUni),nEH,nTarg); % ,nPr
% wFitSessCond=[];
% for iS=1:length(sidUni)
%     disp(['===== ' num2str(sidUni(iS)) ' =====']);
%     idS=sessId==sidUni(iS);
%     
% %     for i=1:nPr
%         for j=1:nEH
%             for l=1:nTarg
%                 id=... % idShortTrial==(2-i) &... % idShortTrial
%                     idHandEye==(2-j) &...            % eye hand
%                     theta==180*(l-1) &... % right left
%                     idS &... % session
%                     ~idOut; % outlier
%                 
%                 try
% %                     [wFitTmp,bicTmp]=estWmWpOld(t(id),T(id),true(nnz(id),1),'mmse1Wm',[],-1,true); % display: off
%                     [wFitTmp,bicTmp]=estWmWpOld(t(id),T(id),idShortTrial(id),'mmse2offset',[],-1,true); % display: off
%                 catch
%                     disp('');
%                 end
%                 fieldnm=fieldnames(wFitTmp);
%                 for iField=1:length(fieldnm)
%                     wFitSessCond(iS,j,l).(fieldnm{iField})=wFitTmp.(fieldnm{iField}); % mmse2offset: w_m w_p offset1 offset2 ,i
%                 end
%                 try
%                     bicSessCond(iS,j,l)=bicTmp;
%                 catch
%                     disp(bicTmp);
%                 end
%                 disp([wFitTmp.w_m(:) wFitTmp.w_p(:) wFitTmp.offset1(:) wFitTmp.offset2(:)]);
%             end % l targ
%         end % j EH
% %     end % i pr
% end % serss
% save([behDir 'G_RSGprior_DMFC.mat'],'wFitSessCond','bicSessCond','-append');

% return;


end % if idFitBLS
%% check # trials
% remove no set
d=load('G_RSGprior_DMFC_noSet.mat'); % to retrieve noSet

tRemove1=find(d.sessId==170822,1); % see chkTInfo
tRemove2=find(d.sessId==170823,1);
d.idNoSet([tRemove1 tRemove2])=nan;
d.idNoSet=d.idNoSet(~isnan(d.idNoSet));

iNoSet=find(sessId>170512,1)-1+find(d.idNoSet); % 1902 among 40317
idNoSet=false(size(t));
idNoSet(iNoSet)=true;
% save('G_RSGprior_DMFC.mat','idNoSet','-append');

idComplete=(t(~idNoSet)<3*T(~idNoSet) & t(~idNoSet)>0);
disp('# total trials:');
disp(nnz(t(~idNoSet))); % 38415

disp('# completed trials:');
disp(nnz(idComplete)); % 30777

disp('# abort trials:');
disp(nnz(~idComplete)); % 7638

disp('# completed & outlier trials:');
disp(nnz(idComplete & idOut(~idNoSet))); % 1755

disp('% outlier among completed trials:');
disp(100*nnz(idComplete & idOut(~idNoSet))/nnz(idComplete)); % 5.7023 % 100*mean(idOut(idComplete)));

%% check block size
idC=(2-idShortTrial-1)*nEH+(1-idHandEye)+1; % [32839 x 1] 1234:shorteye/shorthand/longeye/longhand
idComplete=(t<3*T & t>0);

idTSw=[true;[(abs(diff(idC(idComplete)))>0)|(diff(sessId(idComplete))~=0)]]; % 1 for the trial after switchn
nTrialBlock=diff(find(idTSw));
% figure; hist(nTrialBlock,100);
[m,sd]=meanSD(nTrialBlock)


%% larger bias for long prior @ 20171117
% @2018/8/17: slope estimated with weighted regression
% do linear regression between mean(tp|ts) and ts for each prior > compare slope
% two prior x two EH slopes from each session

sessUni=unique(sessId);
slope=nan(nPr,nEH,length(sessUni));
sFull=nan(nPr,nEH,nTarg,length(sessUni)); % for all comb (pr, EH, targ)
for iSess=1:length(sessUni)
    % session data
    idSess=sessId==sessUni(iSess);        
    ts=T(~idOut & idSess);
    tp=t(~idOut & idSess);
    ctidx=2-idShortTrial(~idOut & idSess); % 1 for short, 2 for long
    citdx2=2-idHandEye(~idOut & idSess); % 1 for eye, 2 for hand
    citdx3=2-(theta(~idOut & idSess)==0); % 1 for right, 2 for left
    
    % weighted regression with precision of tp as weights
    % getting std(tp|ts)
    stp=nan(size(ts));
    stp2=nan(size(ts));
    for i=1:nPr
        for k=1:nTspp
            for j=1:nEH
                id=ctidx==i &... % idShortTrial
                    ts==Tmat{i}(k) &... % ts
                    citdx2==j;            % eye hand
                stp(id)=std(tp(id));
                for l=1:nTarg
                    id2=ctidx==i &... % idShortTrial
                        ts==Tmat{i}(k) &... % ts
                        citdx2==j &...            % eye hand
                        citdx3==l; % right left
                    stp2(id2)=std(tp(id2));
                end
            end
        end
    end
    for i=1:nPr
        for j=1:nEH
            id=ctidx==i &... % idShortTrial
                citdx2==j;            % eye hand
            LM=fitlm(ts(id),tp(id),'linear',...
                'Weights',1./stp(id));
            slope(i,j,iSess)=table2array(LM.Coefficients(2,1));
%             disp(slope(i,j,iSess)); % debug
%             figure; plot(ts(id),tp(id),'o'); ha; plot(ts(id),predict(LM),'-'); waitforbuttonpress; close;
            for l=1:nTarg
                id2=ctidx==i &... % idShortTrial
                    citdx2==j &...            % eye hand
                    citdx3==l; % right left
                LM2=fitlm(ts(id2),tp(id2),'linear',...
                    'Weights',1./stp2(id2));
                sFull(i,j,l,iSess)=table2array(LM2.Coefficients(2,1));
%                 disp(sFull(i,j,k,l,iSess)); % debug
%                 figure; plot(ts(id),tp(id),'o'); ha; plot(ts(id),predict(LM),'-'); waitforbuttonpress; close;
            end % l targ
        end % j EH
    end % i Pr
                 
     % old slope estimation with mean tp
%     % getting mean(tp|ts)
%     mtp=nan(nPr,nEH,nTspp);
%     mtp2=nan(nPr,nEH,nTarg,nTspp);
%     for i=1:nPr
%         for k=1:nTspp
%             for j=1:nEH
%                 id=ctidx==i &... % idShortTrial
%                     ts==Tmat{i}(k) &... % ts
%                     citdx2==j;            % eye hand
%                 mtp(i,j,k)=mean(tp(id));
%                 for l=1:nTarg
%                     id2=ctidx==i &... % idShortTrial
%                         ts==Tmat{i}(k) &... % ts
%                         citdx2==j &...            % eye hand
%                         citdx3==l; % right left
%                     mtp2(i,j,l,k)=mean(tp(id2));
%                 end
%             end
%         end
%     end
%     
%     % slope
%     for i=1:nPr
%         for j=1:nEH
%             tmp=regstats(squeeze(mtp(i,j,:)),Tmat{i}(:),'linear',{'beta'}); % beta
%             slope(i,j,iSess)=tmp.beta(2);
%             for k=1:nTarg
%                 tmp=regstats(squeeze(mtp2(i,j,k,:)),Tmat{i}(:),'linear',{'beta'}); % beta
%                 sFull(i,j,k,iSess)=tmp.beta(2);
%             end
%         end
%     end
end

% % plot
% hFig=figure;H=[];
% for i=1:size(slope,1)
%     for j=1:size(slope,2)
%         x=squeeze(slope(i,j,:));
%         tmptmpcmap=pplot.cmap{(i-1)*size(slope,2)+j};
%         [~,~,hTmp]=histStairs(x,15,0,hFig,tmptmpcmap);H=[H;hTmp];ha;
%         plot(mean(x),max(get(gca,'ylim')),'v','color',tmptmpcmap,'markerfacecolor',tmptmpcmap);
%     end
% end
% legend(H,'ShortEye','ShortHand','LongEye','LongHand','location','best'); legend boxoff;

% figure for slopes
sS=squeeze(sFull(1,:,:,:));
sL=squeeze(sFull(2,:,:,:));
msize=5;

figure; setFigPos(2,1); ha; % no condition info
plot([1;2],[sS(:)'; sL(:)'],'-','color',[.5 .5 .5],'linewidth',.5); % connecting lines
plot(1,sS(:),'.','color',tmpCmap{1,1}(1,:)); % individual points
plot(2,sL(:),'.','color',tmpCmap{2,1}(1,:));
plot([1;2],mean([sS(:)'; sL(:)'],2),'k-','linewidth',2); % connecting lines
plot(1,mean(sS(:)),'o','markersize',msize,'markerfacecolor','w','color',tmpCmap{1,1}(1,:)); % individual points
plot(2,mean(sL(:)),'o','markersize',msize,'markerfacecolor','w','color',tmpCmap{2,1}(1,:));
set(gca,'xtick',1:2,'xticklabel',{[];[]},'ytick',0:.5:1,'xlim',[0.9 2.1],'tickdir','out','ticklength',[0.03 0.03]);

disp('signrank test for all conditions');
[p,h,stats]=signrank(sS(:),1,'tail','left');
disp(['H0(slope(short)<1): ' num2str([p stats.signedrank])]);
[p,h,stats]=signrank(sL(:),1,'tail','left');
disp(['H0(slope(long)<1): ' num2str([p stats.signedrank])]);
[p,h,stats]=signrank(sS(:),0,'tail','right');
disp(['H0(slope(short)>0): ' num2str([p stats.signedrank])]);
[p,h,stats]=signrank(sL(:),0,'tail','right');
disp(['H0(slope(long)>0): ' num2str([p stats.signedrank])]);
[p,h,stats]=signrank(sS(:),sL(:),'tail','right');
disp(['H0(slope(long)<slope(short)): ' num2str([p stats.signedrank])]);

% 
% figure; setFigPos(1,1); ha; % plot condition-specific but avgAcrossSessions
% for iEH=1:nEH
%     for iTarg=1:nTarg
%         sS=squeeze(sFull(1,iEH,iTarg,:));
%         sL=squeeze(sFull(2,iEH,iTarg,:));
%         plot([1;2],mean([sS(:)'; sL(:)'],2),'-','color',[.5 .5 .5],'linewidth',.5); % connecting lines
%         plot(1,mean(sS(:)),'.','color',tmpCmap{1,1}(1,:)); % individual points
%         plot(2,mean(sL(:)),'.','color',tmpCmap{2,1}(1,:));
%     end % for iTarg=1:nTarg
% end % for iEH=1:nEH
% sS=squeeze(sFull(1,:,:,:));
% sL=squeeze(sFull(2,:,:,:));
% plot([1;2],mean([sS(:)'; sL(:)'],2),'k-','linewidth',2); % connecting lines
% plot(1,mean(sS(:)),'o','markersize',msize,'markerfacecolor','w','color',tmpCmap{1,1}(1,:)); % individual points
% plot(2,mean(sL(:)),'o','markersize',msize,'markerfacecolor','w','color',tmpCmap{2,1}(1,:));
% set(gca,'xtick',1:2,'xticklabel',{[];[]},'ytick',0:.5:1,'xlim',[0.9 2.1],'tickdir','out','ticklength',[0.03 0.03]);

figure; setFigPos(1,1); ha; % plot session-specific but avgAcrossConditions
for iSess=1:length(sessUni)
    sS=squeeze(sFull(1,:,:,iSess));
    sL=squeeze(sFull(2,:,:,iSess));
    plot([1;2],mean([sS(:)'; sL(:)'],2),'-','color',[.5 .5 .5],'linewidth',.5); % connecting lines
    plot(1,mean(sS(:)),'.','color',tmpCmap{1,1}(1,:)); % individual points
    plot(2,mean(sL(:)),'.','color',tmpCmap{2,1}(1,:));
end %  for iSess=1:length(sessUni)
sS=squeeze(sFull(1,:,:,:));
sL=squeeze(sFull(2,:,:,:));
plot([1;2],mean([sS(:)'; sL(:)'],2),'k-','linewidth',2); % connecting lines
plot(1,mean(sS(:)),'o','markersize',msize,'markerfacecolor','w','color',tmpCmap{1,1}(1,:)); % individual points
plot(2,mean(sL(:)),'o','markersize',msize,'markerfacecolor','w','color',tmpCmap{2,1}(1,:));
set(gca,'xtick',1:2,'xticklabel',{[];[]},'ytick',0:.5:1,'xlim',[0.9 2.1],'tickdir','out','ticklength',[0.03 0.03]);

disp('signrank test ignoring conditions');
disp(['H0(slope(short)<1): ' num2str(signrank(sS(:),1,'tail','left'))]);
disp(['H0(slope(long)<1): ' num2str(signrank(sL(:),1,'tail','left'))]);
disp(['H0(slope(short)>0): ' num2str(signrank(sS(:),0,'tail','right'))]);
disp(['H0(slope(long)>0): ' num2str(signrank(sL(:),0,'tail','right'))]);
disp(['H0(slope(long)<slope(short)): ' num2str(signrank(sS(:),sL(:)))]);

% anova
% matPr=repmat([1;2],[1 nEH length(sessUni)]);
% matEH=repmat([1 2;1 2],[1 1 length(sessUni)]);
% stat=anovan(slope(:),{matPr(:) matEH(:)},'model','interaction','display','on','varname',{'idPr','idEH'});

matPr=repmat([1;2],[1 nEH nTarg length(sessUni)]);
matEH=repmat([1 2],[nPr 1 nTarg length(sessUni)]);
matTarg=repmat(cat(3,1,2),[nPr nEH 1 length(sessUni)]);

stat=anovan(sFull(:),{matPr(:) matEH(:) matTarg(:)},'model','interaction','display','on','varname',{'idPr','idEH','idRL'});


% new supplementary figure 1 (distribution of slope)
% sS=squeeze(sFull(1,:,:,:)); % EH RL session
tmpMarker={'s','d';'o','^'};

idHist=1; % 0; % 1;

figure; setFigPos(2,1); ha;
nBin=7;barWidth=0.4;markersize=4;
optsExpFig.Width=3; % 10/2.54;
optsExpFig.Height=3; % 7.3/2.54;
pause(0.1);
% for iEH=1:nEH
%     for iRL=1:nTarg
% %         figure; setFigPos(iEH,iRL); ha;
%         x=sS(iEH,iRL,:);        y=sL(iEH,iRL,:);
                x=squeeze(mean(mean(sFull(1,:,:,:),2),3)); % sessions by 1
        y=squeeze(mean(mean(sFull(2,:,:,:),2),3));
        
        if idHist
            b=linspace(min([x(:);y(:)]),max([x(:);y(:)]),nBin);
            
            h=histogram(x,b); i=1; % short
            Xval    = h.BinEdges + h.BinWidth*0.5 +h.BinWidth*(i-(nPr+1)/2)/2; % binWidth*0.25 for short, 0.75 for long
            Xval    = Xval(1:end-1);Yval    = h.Values; delete(h);
            hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',pplot.cmap{(i-1)*nPr+1},'EdgeColor','none');
            
            h=histogram(y,b); i=2; % long
            Xval    = h.BinEdges + h.BinWidth*0.5 +h.BinWidth*(i-(nPr+1)/2)/2; % binWidth*0.25 for short, 0.75 for long
            Xval    = Xval(1:end-1);Yval    = h.Values; delete(h);
            hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',pplot.cmap{(i-1)*nPr+1},'EdgeColor','none');
            
            axis tight; % plot mean
            maxYlim=max(ylim);
            plot(mean(x),1.1*maxYlim,'color',pplot.cmap{(1-1)*nPr+1},'markerfacecolor',pplot.cmap{(i-1)*nPr+1},'marker','v','markersize',markersize);
            plot(mean(y),1.1*maxYlim,'color',pplot.cmap{(i-1)*nPr+1},'markerfacecolor',pplot.cmap{(i-1)*nPr+1},'marker','v','markersize',markersize);
            
            plotVertical(gca,0,[]); plotVertical(gca,1,[]);
            
            xlim([-0.2 1.3]);
            set(gca,'xtick',0:0.5:1,'ytick',0:3:20,'tickdir','out','ticklength',[0.02 0.02]);
            pause(0.1);
            
        else % scatter
            plot(x(:),y(:),['k' tmpMarker{iEH,iRL}],'markerfacecolor','w','markersize',markersize,'linewidth',0.5,'color',pplot.cmap{(iEH-1)*nTarg+iRL});
            if iEH==1 & iRL==2 % EL
                axis([-0.1 1.2 -0.1 1.3]);
            else
                axis([-0.1 1.2 -0.1 1.2]);
            end
            plotIdentity(gca);  plotVertical(gca,0,[]); plotVertical(gca,1,[]);plotHorizon(gca,0,[]); plotHorizon(gca,1,[]);
            set(gca,'xtick',0:0.25:2,'ytick',0:0.25:2,'tickdir','out','ticklength',[0.02 0.02]); remTickLabel;
            pause(0.1);
            applytofig(gcf,optsExpFig);
        end
        
        % state
        disp('===== G =====');
        disp([ehNm{iEH} targNm{iRL}]);
        disp('signrank test ignoring conditions');
        disp(['H0(slope(short)<1): ' num2str(signrank(x(:),1,'tail','left'))]);
        disp(['H0(slope(long)<1): ' num2str(signrank(y(:),1,'tail','left'))]);
        disp(['H0(slope(short)>0): ' num2str(signrank(x(:),0,'tail','right'))]);
        disp(['H0(slope(long)>0): ' num2str(signrank(x(:),0,'tail','right'))]);
        disp(['H0(slope(long)<slope(short)): ' num2str(signrank(y(:),x(:)))]);
        
%     end % for iRL=1:nTarg
% end % for iEH=1:nEH
disp([]);
%% 3/5/2018: switching across blocks (bias=f(trials after switch), only overlap?

% 12/6, 12/7, part of 12/10 [1 3 15] used one-trial block
sessTBT=[161206; 161207; 161210];
for iSess=1:length(sessTBT)
    idOut(sessId==sessTBT(iSess))=true;
end % for iSess=1:length(sessTBT)

% identifying condition: prior(S/L) by effector (E/H) 
% 1234: shortE>shortH>longE>longH
idC=(2-idShortTrial-1)*nEH+(1-idHandEye)+1; % [32839 x 1] 1234:shorteye/shorthand/longeye/longhand
idTSw=[true;[(abs(diff(idC))>0)|(diff(sessId)~=0)]]; % 1 for the trial after switchn
iTSw=idC(idTSw); % [6507 blocks x 1]
idBlock=cumsum(idTSw); % [11112222...]
idComplete=(t<3*T & t>0);

normTp=(t-T)./T;
idOverlap=(T==Tmat{1}(end));

nBlock=max(idBlock);
normTp2=[]; % cell(1,nBlock);
idOverlap2=[]; % cell(1,nBlock);
idOut2=[]; % cell(1,nBlock);
idSwPr=[]; % -1 long2short, 1 for short2long, ignoring EH switch
iSwPr=1;
for iB=1:nBlock
    % only considering # trials(complete)
%     normTp2{iB}=normTp(idBlock==iB & idComplete); %[trial/block x 1]
%     idOverlap2{iB}=idOverlap(idBlock==iB & idComplete);
%     idOut2{iB}=idOut(idBlock==iB & idComplete);
    if iB==1 % assume transition
        if iTSw(iB)<=2 % short
            idSwPr=-1; % -1 long2short
        else
            idSwPr=1;
        end
        normTp2{iB}=normTp(idBlock==iB & idComplete); %[trial/block x 1]
        idOverlap2{iB}=idOverlap(idBlock==iB & idComplete);
        idOut2{iB}=idOut(idBlock==iB & idComplete);
    else
        if iTSw(iB-1)<=2 & iTSw(iB)>2 % short2long
            idSwPr=[idSwPr;1];iSwPr=iSwPr+1;
            normTp2{iSwPr}=normTp(idBlock==iB & idComplete); %[trial/block x 1]
            idOverlap2{iSwPr}=idOverlap(idBlock==iB & idComplete);
            idOut2{iSwPr}=idOut(idBlock==iB & idComplete);
        elseif iTSw(iB-1)>2 & iTSw(iB)<=2 % long2short
            idSwPr=[idSwPr;-1];iSwPr=iSwPr+1;
            normTp2{iSwPr}=normTp(idBlock==iB & idComplete); %[trial/block x 1]
            idOverlap2{iSwPr}=idOverlap(idBlock==iB & idComplete);
            idOut2{iSwPr}=idOut(idBlock==iB & idComplete);
        else % appending to previous blocks (changing EH/session w/o changing prior)
            normTp2{iSwPr}=[normTp2{iSwPr};normTp(idBlock==iB & idComplete)]; %[trial/block x 1]
            idOverlap2{iSwPr}=[idOverlap2{iSwPr};idOverlap(idBlock==iB & idComplete)];
            idOut2{iSwPr}=[idOut2{iSwPr};idOut(idBlock==iB & idComplete)];
        end
    end
end % iB=1:max(idBlock)

nTpB=cellfun(@length,normTp2);

normTp2=fillNanCell(normTp2); % fillNanCell accepts cells of vectors [6507 blocks x 26 trials]
idOverlap2=fillNanCell(idOverlap2);
idOverlap2(idOverlap2~=1)=nan;idOverlap2(idOverlap2==1)=0;
idOut2=fillNanCell(idOut2); % 
idOut2(idOut2==1)=nan; % nan if outlier

% TBD: plot w/ avg
trial=0:(size(normTp2,2)-1);
% 1) normTp2 w/o outlier: all trials
figure; setFigPos(1,2); ha;
iSwPr=unique(idSwPr); % -1 1
for iSw=1:length(iSwPr) % short2long vs long2short
    normTp2woOut=normTp2(idSwPr==iSwPr(iSw),:)+idOut2(idSwPr==iSwPr(iSw),:); % nan if outlier
    mNormTp2=nanmean(normTp2woOut,1);
    semNormTp2=sem(normTp2woOut,0,1,'omitnan');
    shadedErrorBar(trial,mNormTp2,semNormTp2,{'-','color',tmpCmap{iSw,1}(1,:),'linewidth',2},1); drawnow;
end
nTrialEnd=30;
set(gca,'xlim',[0 nTrialEnd],'ylim',[-0.1 0.1],'ytick',[-0.1:0.05:0.1],'tickdir','out','ticklength',[0.03 0.03]);
xlabel('Trials after prior switch'); ylabel('(t_p - t_s)/t_s');

% 2) normTp2 w/o outlier: only overlap trials
figure; setFigPos(1,3); ha;
for iSw=1:length(iSwPr) % short2long vs long2short
    normTp2woOut=normTp2(idSwPr==iSwPr(iSw),:)+idOut2(idSwPr==iSwPr(iSw),:)+idOverlap2(idSwPr==iSwPr(iSw),:); % nan if outlier
    mNormTp2=nanmean(normTp2woOut,1);
    semNormTp2=sem(normTp2woOut,0,1,'omitnan');
    shadedErrorBar(trial,mNormTp2,semNormTp2,{'-','color',tmpCmap{iSw,1}(nTspp-(nTspp-1)*(iSw-1),:),'linewidth',2},1); drawnow;
end
nTrialEnd=20;
set(gca,'xlim',[0 nTrialEnd],'ylim',[-0.1 0.1],'ytick',[-0.1:0.05:0.1],'xtick',[0:5:nTrialEnd],'tickdir','out','ticklength',[0.03 0.03]);
xlabel('Trials after prior switch'); ylabel('(t_p-800)/800 (overlap t_s)');

%% compare bias b/t longest vs shortest ts for each prior
% % initially pooled across all sessions/effectors/directions
% % : issue>for short, bias(shortest ts)>bias(longest ts)
% % offset? correct by bias at prior mean? didn't work (small like 6)
% % separate by sessions,effectors
% 
% idCompBias=1; % 0; % 0; % 1;
% 
% figure;
% tmpMarker={'o','d';'^','v'}; % eye/hand x right/left
% 
% sessUni=unique(sessId);
% 
% mtp=nan(length(sessUni),nEH,nTarg,nPr,2);stp=nan(length(sessUni),nEH,nTarg,nPr,2);
% offsetMuPr=nan(length(sessUni),nEH,nTarg,nPr);
% 
% for iSess=1:length(sessUni)
%     % session data
%     disp(sessUni(iSess));
%     idSess=sessId==sessUni(iSess);
%     
%     ts=T(~idOut&idSess);
%     tp=t(~idOut&idSess);
%     ctidx=2-idShortTrial(~idOut & idSess); % 1 for short, 2 for long
%     citdx2=2-idHandEye(~idOut & idSess); % 1 for eye, 2 for hand
%     citdx3=2-(theta(~idOut & idSess)==0); % 1 for right, 2 for left
%     
%     for i=1:nPr
%         for j=1:nEH
%             for k=1:nTarg
%             % shortest ts
%             id1=ctidx==i &... % idShortTrial
%                 ts==Tmat{i}(1) &... % ts
%                 citdx2==j&... %;            % eye hand
%                 citdx3==k; % left right
%             mtp(iSess,j,k,i,1)=mean(tp(id1));
%             stp(iSess,j,k,i,1)=sem(tp(id1));
%             % longest ts
%             id2=ctidx==i &... % idShortTrial
%                 ts==Tmat{i}(end) &... % ts
%                 citdx2==j&... %;            % eye hand
%                 citdx3==k; % left right
%             mtp(iSess,j,k,i,2)=mean(tp(id2));
%             stp(iSess,j,k,i,2)=sem(tp(id2));
%             % bias @ prior mean
%             id=ctidx==i &... % idShortTrial
%                 ts==median(Tmat{i}) &... % ts
%                 citdx2==j&... %;            % eye hand
%                 citdx3==k; % left right
%             offsetMuPr(iSess,j,k,i)=mean(tp(id))-ts(find(id,1));
%             
%             %         plot(tp(id1)-ts(id1),tp(id2)-ts(id2),'.','color',pplot.cmap{(i-1)*2+1},'markersize',10);ha;
%             mx=abs(mtp(iSess,j,k,i,1)-offsetMuPr(iSess,j,k,i)*idCompBias-ts(find(id1,1)));
%             my=abs(mtp(iSess,j,k,i,2)-offsetMuPr(iSess,j,k,i)*idCompBias-ts(find(id2,1)));
%             errorbarXY_meanSem(mx,my,...
%                 stp(iSess,j,i,1),stp(iSess,j,i,2),pplot.cmap{(i-1)*2+1},tmpMarker{j,k},10); ha;
%             end
%         end
%     end
% end
% axis tight;
% tmpX=get(gca,'xlim');
% tmpY=get(gca,'ylim');
% 
% axis([min([tmpX tmpY]) max([tmpX tmpY]) min([tmpX tmpY]) max([tmpX tmpY])]);
% plotIdentity(gca);
% xlabel('mean tp-shortest ts');ylabel('mean tp-longest ts');
% 
% % % check correlation b/t offsetMuPr and BLS offset (not separate for EH)
% % % conclusion: overall correlated
% % figure; setFigPos(1,1);ha;
% % plot(offsetMuPr(:,1,1),[wFitSess.offset1],'ro'); % eye short
% % plotIdentity(gca);
% % figure; setFigPos(2,1);ha;
% % plot(offsetMuPr(:,2,1),[wFitSess.offset1],'r^'); % hand
% % plotIdentity(gca);
% % figure; setFigPos(1,2);ha;
% % plot(offsetMuPr(:,1,2),[wFitSess.offset2],'bo'); % eye long
% % plotIdentity(gca);
% % figure; setFigPos(2,2);ha;
% % plot(offsetMuPr(:,2,2),[wFitSess.offset2],'b^'); % hand
% % plotIdentity(gca);

%% SELECT BEST SESSION FOR VISUALIZATION
idExSess=1;
iSess=170512; % 818;
if idExSess
    idSess=sessId==iSess;
else
    idSess=true(size(sessId));
end

d.ts=T(~idOut & idSess);
d.tp=t(~idOut & idSess);
ctidx=2-idShortTrial(~idOut & idSess); % 1 for short, 2 for long
ipidx=nan(size(ctidx));

%% figure 1B: tp vs ts
close all;
PlotDiff=0; % 1;

% data formatting: [# interval x maxRepetition]
% find max # trials 1st
for i=1:nPr
    for k=1:nTspp
        id=idShortTrial(~idOut  & idSess)==(2-i) &... % idShortTrial
            d.ts==Tmat{i}(k); % ts
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
        id=idShortTrial(~idOut  & idSess)==(2-i) &... % idShortTrial
            d.ts==Tmat{i}(k); % ts
        Ts((i-1)*nTspp+k,1:nnz(id))=d.ts(id);
        Tp((i-1)*nTspp+k,1:nnz(id))=d.tp(id);
    end
end
d.Ts=Ts;
d.Tp=Tp;
% theory
if idExSess
    idSessTmp=find(unique(sessId)==iSess);
    wm=wFitSess(idSessTmp).w_m;
    wp=wFitSess(idSessTmp).w_p;    
    offset=[wFitSess(idSessTmp).offset2; wFitSess(idSessTmp).offset1];
else
    wm=wFit.w_m;
    wp=wFit.w_p;
    offset=[wFit.offset2; wFit.offset1];
end
for i=1:nPr
    theory{i}.tm=[min(Tmat{i}):max(Tmat{i})]'; % theory{x}.tm - Tx1 - sample corresponding to above
    theory{i}.te=BLS(theory{i}.tm,wm,[min(Tmat{i}) max(Tmat{i})],'uniform')+offset(i);theory{i}.te=theory{i}.te(:);
    theory{i}.range=1:length(theory{i}.tm);
    theory{i}.wm=wm;
end
DnPlotResponses2P(d,theory,ctidx,ipidx,PlotDiff); % PlotDiff for the last input

% paper
if PlotDiff
    axis tight;
    xlim([450 1250]);
    set(gca,'xtick',[480 640 800 1000 1200],'ytick',-400:50:300,'xticklabel',[],'yticklabel',[]);
%     ylim([350 1400]);
    xlabel([]);ylabel([]);
    title([])
    load pplot.mat;
    
    % optsExpFig.Width=5; % 2.6/2.54*2;
    % optsExpFig.Height=5; % 1.9/2.54*2;
    optsExpFig.Format='eps';
    % exportfig(gcf,'tp_ts_H_small.eps',optsExpFig);
    % poster
    optsExpFig.Width=4/2.54; %6.5; % 10/2.54;
    optsExpFig.Height=4/2.54; %6.5; % 7.3/2.54;
    if idExpFig
        savefig(gcf,'tp_ts_H.fig');
        exportfig(gcf,'tp_ts_H.eps',optsExpFig);
    end
else
xlim([450 1250]);
set(gca,'xtick',[480 640 800 1000 1200],'ytick',[480 640 800 1000 1200],'xticklabel',[],'yticklabel',[]);
ylim([350 1400]);
xlabel([]);ylabel([]);
title([])
load pplot.mat;
% optsExpFig.Width=5; %2.6/2.54*2;
% optsExpFig.Height=5; % 1.9/2.54*2;
optsExpFig.Format='eps';
% optsExpFig.Format='eps';exportfig(gcf,'tp_ts_G_small.eps',optsExpFig);
% poster
optsExpFig.Width=4/2.54; % 6.5; % 10/2.54;
optsExpFig.Height=4/2.54; % 6.5; % 7.3/2.54;
if idExpFig
exportfig(gcf,'tp_ts_G.eps',optsExpFig);
end
end

%% supplementary figure 1: copy of figure1 for each condition (ER,EL,HR,HL)
close all;
v=v2struct;
plotSF1(v);


%% figure 1C: hist

HsHistTpOverlap(d,ctidx);
xlim([800-350 800+350]);
set(gca,'xtick',600:200:1000,'ytick',0:20:40,'tickdir','out'); xlabel([]);ylabel([]);
set(gca,'xticklabel',[],'yticklabel',[]);
load pplot.mat;
% optsExpFig.Width=2; % 2.6/2.54*2;
% optsExpFig.Height=2; % 1.9/2.54*2;
optsExpFig.Format='eps';
% exportfig(gcf,'histTpOverlap_G_small.eps',optsExpFig);
% poster
optsExpFig.Width=3; % 10/2.54;
optsExpFig.Height=2; % 7.3/2.54;
if idExpFig
    savefig(gcf,'histTpOverlap_G.fig');

exportfig(gcf,'histTpOverlap_G.eps',optsExpFig);
end

%% supplementary figure 2: bias-var (BLS), data vs model
% now session-specific (line data-model, filled for model, circle for data)
% wFitSessCond: 17sess x 2EH x 2RL: w_m, w_p, offset1(long), offset2
HsPlotBiasVar2(sessId,idHandEye,theta,T,t,idShortTrial,wFitSessCond,idOut);


%% figure 1D: bias variance


dBLS.wm=wm; % wFit.w_m;
dBLS.wp=wp; % wFit.w_p;
dBLS.offset1=offset(2); % wFit.offset1; % for long
dBLS.offset2=offset(1); % wFit.offset2;
% with offset, discrepancy is lower; TBD: model comparison with prior-dependent wm,wp
HsPlotBiasVar(d,ctidx,dBLS);

set(gca,'xtick',0:50:150,'ytick',0:50:150,'tickdir','out','ticklength',[0.03 0.03]); xlabel([]);ylabel([]);

load pplot.mat;
optsExpFig.Width=2.6/2.54;optsExpFig.Height=1.9/2.54;
optsExpFig.Format='eps';
if idExpFig
exportfig(gcf,'biasVar_G.eps',optsExpFig);
end
% % poster
% optsExpFig.Width=10/2.54;optsExpFig.Height=7.3/2.54;
% exportfig(gcf,'biasVar_G_large.eps',optsExpFig);



%% stat for overlap


idAbort=t<=0; disp(['abort:' num2str([nnz(idAbort) nnz(idAbort)/nnz(t)])]);
idNoResp=t>=3*T; disp(['noResp:' num2str([nnz(idNoResp) nnz(idNoResp)/nnz(t)])]);
idValid=t>0&t<3*T; disp(['outlier response:' num2str([nnz(idOut&idValid) nnz(idOut&idValid)/nnz(idValid)])]);

% anova with new outlier selected
overlapT=800; 
tmpId=~idOut&T==overlapT;
stat=anovan(t(tmpId),{idShortTrial(tmpId) idHandEye(tmpId) theta(tmpId)/180},'model',3,'display','on','varnames',{'idShort' 'idEye' 'idRight'});

idPlot=0; % 1;
for i=1:nPr
    for j=1:nEH
        for k=1:nTarg
            id=idShortTrial==(2-i) &... % idShortTrial
                idHandEye==(2-j) &...  % idHandEye
                theta==(k-1)*180 &... % target location
                T==overlapT &...
                ~idOut;
            
                        disp([prNm{i} ',' ehNm{j} ',' targNm{k} num2str(mean(t(id)))]);
            
            if idPlot
                figure; set(gcf,'position',pplot.(['rect' num2str(i) '_' num2str((j-1)*nTarg+k)]));hold all;
                histfit(t(id),35); axis tight; plotHorizon(gca,mean(t(id)),[]);
                title([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, \mu&SD: ' num2str(meanSD(t(id))) ' ms']);
%                 disp([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, \mu&SD: ' num2str(meanSD(t(id))) ' ms']);
                xlim([400 1200]);
            end
            
        end % for k=1:nTarg
    end % for j=1:nEH
    id=idShortTrial==(2-i) &... % idShortTrial
        T==overlapT &...
        ~idOut;
    disp([prNm{i} ', overlap, \mu&SD: ' num2str(meanSD(t(id))) ' ms']);
end % for i=1:nPr

%% BLS fit for each sessions
% for now, wm, wp common for all condition (prior, ts, modality, directions) & offset separately for prior
% model comparion to be done: promising (offset,wp separate for modality)

% sidUni=unique(sessId);
% bicSess=[];wFitSess=[];
% for iS=1:length(sidUni)
%     disp(['===== ' num2str(sidUni(iS)) ' =====']);
%     idS=sessId==sidUni(iS);
%     id=idS&~idOut;
%     
%     [wFitTmp,bicTmp]=estWmWpOld(t(id),T(id),idShortTrial(id),'mmse2offset',[],-1,true); % display: off
%     fieldnm=fieldnames(wFitTmp);
%     for j=1:length(fieldnm)
%         wFitSess(iS).(fieldnm{j})=wFitTmp.(fieldnm{j}); % mmse2offset: w_m w_p offset1 offset2
%     end
%     bicSess=[bicSess; bicSess];
%     disp([wFitTmp.w_m(:) wFitTmp.w_p(:) wFitTmp.offset1(:) wFitTmp.offset2(:)]);
%     save([behDir 'G_RSGprior_DMFC.mat'],'wFitSess','bicSess','-append');
% 
% end
% 
% % also BLS fit for pooled data
% analHandEye2PriorRSG('G_RSGprior_DMFC.mwk'); 

%% check across-session variability of wm,wp,offset (wFit w/ outlier removal)
figure; nMeas=4; 
load([behDir 'G_RSGprior_DMFC.mat']); % T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)

for iS=1:length(wFitSess)
    subplot(nMeas,1,1);plot(iS,wFitSess(iS).w_m,'o'); ylabel('w_m');hold all;
    subplot(nMeas,1,2);plot(iS,wFitSess(iS).w_p,'o'); ylabel('w_p');hold all;
    subplot(nMeas,1,3);plot(iS,wFitSess(iS).offset1,'o'); ylabel('offset1');hold all;
    subplot(nMeas,1,4);plot(iS,wFitSess(iS).offset2,'o'); ylabel('offset2');hold all;
end
disp('check offset is needed (long,short)')
p=signrank([wFitSess.offset1])
p=signrank([wFitSess.offset2])

% output

% abort:4283     0.25882282
% noResp:30   0.00181291
% outlier response:840     0.0686555
% Short, overlap, \mu&SD: 787.0071 ms
% Long, overlap, \mu&SD: 799.6886 ms
% ===== 170506 =====
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 4382.398474 seconds.
%     0.0208    0.0208    0.1409    0.1409
% 
% ===== 170507 =====
% 
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% plotTraj('traj_periSet_H_attrition.mat');cc;plotTraj('traj_periSet_G_attrition.mat');
% Elapsed time is 7312.201795 seconds.
%     0.0114    0.1522  -43.3054   14.1497
% 
% ===== 170508 =====
% 
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 4302.464814 seconds.
%     0.0471    0.1438  -23.9119   27.3071
% 
% ===== 170510 =====
% 
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 4349.493720 seconds.
%     0.0225    0.1142  -45.9624    1.9498
% 
% ===== 170511 =====
% 
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 5880.654135 seconds.
%     0.0605    0.0865  -49.0636  -10.3220
% 
% ===== 170512 =====
% 
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 3416.583438 seconds.
%     0.0637    0.0947  -14.5295    0.3105

%% define D data structure
% tSort=sort(t(t>0));tB=tSort([1 round([(1/nTp):(1/nTp):1]*length(tSort))]); % 0 0.2 0.4 ... 1
% tLong=sort(t(idShortTrial==0&t>0));tLB=tLong([1 round([(1/nTp):(1/nTp):1]*length(tLong))]); % 0 0.2 0.4 ... 1
% tShort=sort(t(idShortTrial==1&t>0));tSB=tShort([1 round([(1/nTp):(1/nTp):1]*length(tShort))]); % 0 0.2 0.4 ... 1

%% session-pooled outlier removal
% % 2 prior x 2 modality x 2 direction ANOVA for overlap ts
% % after remove outliers (2 SD) for each combination
% idPlot=1;
% nSDremoveOutlier=2;
% overlapT=800; 
% tMat=[];iPr=[];iEH=[];iRL=[];
% for i=1:nPr
%     for j=1:nEH
%         for k=1:nTarg
%             id=idShortTrial==(2-i) &... % idShortTrial
%                 idHandEye==(2-j) &...  % idHandEye
%                 theta==(k-1)*180 &... % target location
%                 T==overlapT &...
%                 t>0 & t<3*overlapT; % only valid trials for now
%             u0=[0.5 1.5]*overlapT; %[0.1 1.9]*overlapT; % [0.5 1.5]*overlapT; %[0 2]*overlapT; % with least outlier removal; or max(t)
%             [tClean,idClean,pOut,tOut]=removeOutlier(t(id),0,u0); %nSDremoveOutlier);
%             tMat=[tMat;tClean(:)];
%             iPr=[iPr; repmat(prNm(i),[nnz(tClean) 1])];iEH=[iEH; repmat(ehNm(j),[nnz(tClean) 1])];iRL=[iRL; repmat(targNm(k),[nnz(tClean) 1])];
%             
%             if idPlot
%                 figure; set(gcf,'position',pplot.(['rect' num2str(i) '_' num2str((j-1)*nTarg+k)]));hold all;
%                 histfit(tClean,35); if ~isempty(tOut), plot(tOut,0,'rx'); end; axis tight;
%                 title([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, p(outlier): ' num2str(pOut) '%']);
%                 disp([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, p(outlier): ' num2str(pOut) '%']);
%                 xlim([400 1200]);
%             end
%             
%         end % for k=1:nTarg
%     end % for j=1:nEH
% end % for i=1:nPr
% stat=anovan(tMat,{iPr iEH iRL},'model',3,'display','on','varnames',{'idShort' 'idEye' 'idRight'})

return;

%% old init
try
    cd('/Users/hansem/Dropbox (MIT)/figuresRSG2prior/1beh'); % SFN17/1beh'); %
    behDir='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/dataMat/';
    neuDir='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/neuralData/';
    psthDir='/Users/hansem/Dropbox (MIT)/psthDataHigh/';
    dirName='/Users/hansem/Documents/Recording/'; % cd(dirName);
catch
    cd('/Users/seonminahn/Dropbox (MIT)/figuresRSG2prior/1beh'); %
    behDir='/Users/seonminahn/Dropbox (MIT)/fileFromNHPrig/dataMat/';
    neuDir='/Users/seonminahn/Dropbox (MIT)/fileFromNHPrig/neuralData/';
    psthDir='/Users/seonminahn/Dropbox (MIT)/psthDataHigh/';
    dirName='/Users/seonminahn/Documents/Recording/'; % cd(dirName);
end

%%

fname={...
    'G_RSGprior_20170506';... % t2080 (c1491)
    'G_RSGprior_20170507';... % t2506 (c1907)
    'G_RSGprior_20170508';... % t2388 (c1700)
%     'G_RSGprior_20170509';... % blocksize [100+exp(-n/20) 200]
    'G_RSGprior_20170510';... % t2003 (c1687) % block seq: LH>SH>SE>LE (randTmp0=0)
    'G_RSGprior_20170511';... % t4338 (c3104)
    'G_RSGprior_20170512'; % t3028 (c2349)
    
    % no set
    'G_RSGprior_20170803';...% 8/3: t2448 (c1649)
    'G_RSGprior_20170817';...% 8/17: t4214 (c3455) 2 probes
    'G_RSGprior_20170818';...%8/18: t4527 (c3883)
    'G_RSGprior_20170821';...%8/21: t3703 (c2982)
    'G_RSGprior_20170822';...%8/22: t4752 (c3893)
    'G_RSGprior_20170823';...%8/23: t4127 (c3251)
    };



load pplot.mat;

nPr=2;    prNm={'Short','Long'};
nEH=2;ehNm={'Eye','Hand'};
nTarg=2; targNm={'Right','Left'};
nTspp=5;Tmat=[]; Tmat{1}=linspace(480,800,nTspp); Tmat{2}=linspace(800,1200,nTspp); Tall=[Tmat{1} Tmat{2}];
nTp=5; 

%% %% supplementary figure 1: copy of figure1 for each condition (ER,EL,HR,HL)
function plotSF1(v)

v2struct(v);

%% SELECT BEST SESSION FOR VISUALIZATION
% iSess=161218;
% d.ts=T(~idOut & idSess);
% d.tp=t(~idOut & idSess);
% ctidx=2-idShortTrial(~idOut & idSess); % 1 for short, 2 for long
% ipidx=nan(size(ctidx));

%% figure 1B: tp vs ts

for iEH=1:nEH
    for iTarg=1:nTarg

        disp(['---- ' ehNm{iEH} ' ' targNm{iTarg} ' -----']);
        
        d.ts=T(~idOut & idSess&...
            idHandEye==(2-iEH) &...
            theta==(iTarg-1)*180);
        d.tp=t(~idOut & idSess&...
            idHandEye==(2-iEH) &...
            theta==(iTarg-1)*180);
        ctidx=2-idShortTrial(~idOut & idSess&...
            idHandEye==(2-iEH) &...
            theta==(iTarg-1)*180); % 1 for short, 2 for long
        ipidx=nan(size(ctidx));
        nMaxTrial=[];
        disp([num2str(nnz(d.tp)) ' trials']);
        
        PlotDiff=0; % 1;
        % data formatting: [# interval x maxRepetition]
        % find max # trials 1st
        for i=1:nPr
            for k=1:nTspp
                tmpId=idHandEye==(2-iEH) &...
                    theta==(iTarg-1)*180;
                id=idShortTrial(~idOut & idSess & tmpId)==(2-i) &... % idShortTrial
                    T(~idOut & idSess& tmpId)==Tmat{i}(k); % ts
                nMaxTrial=max([nMaxTrial; nnz(id)]);
                ipidx(id)=k; % for each trial in d.ts specifies which time bin or ts interval it belongs to. Values: 1-5
            end
        end
        % filling matrix
        Ts=nan(length(Tall),nMaxTrial); % [K x R]
        Tp=nan(length(Tall),nMaxTrial);
        for i=1:nPr
            for k=1:nTspp
                tmpId=idHandEye==(2-iEH) &...
                    theta==(iTarg-1)*180;
                id=idShortTrial(~idOut & idSess & tmpId)==(2-i) &... % idShortTrial
                    d.ts==Tmat{i}(k); % ts
                Ts((i-1)*nTspp+k,1:nnz(id))=d.ts(id);
                Tp((i-1)*nTspp+k,1:nnz(id))=d.tp(id);
            end
        end
        d.Ts=Ts;
        d.Tp=Tp;
        % theory
        if idExSess
            idSessTmp=find(unique(sessId)==iSess);
%             wm=wFitSess(idSessTmp).w_m;
%             wp=wFitSess(idSessTmp).w_p;
%             offset=[wFitSess(idSessTmp).offset2; wFitSess(idSessTmp).offset1];
            wm=wFitSessCond(idSessTmp,iEH,iTarg).w_m;
            wp=wFitSessCond(idSessTmp,iEH,iTarg).w_p;
            offset=[wFitSessCond(idSessTmp,iEH,iTarg).offset2; wFitSessCond(idSessTmp,iEH,iTarg).offset1]; % offset1 for long prior
        else
            wm=wFit.w_m;
            wp=wFit.w_p;
            offset=[wFit.offset2; wFit.offset1];
        end
        for i=1:nPr
            theory{i}.tm=[min(Tmat{i}):max(Tmat{i})]'; % theory{x}.tm - Tx1 - sample corresponding to above
            theory{i}.te=BLS(theory{i}.tm,wm,[min(Tmat{i}) max(Tmat{i})],'uniform')+offset(i);theory{i}.te=theory{i}.te(:);
            theory{i}.range=1:length(theory{i}.tm);
            theory{i}.wm=wm;
        end
        DnPlotResponses2P(d,theory,ctidx,ipidx,PlotDiff); % PlotDiff for the last input
        setFigPos(iEH,iTarg);
        
        xlim([450 1250]);
        set(gca,'xtick',[480 640 800 1000 1200],'ytick',[480 640 800 1000 1200],'xticklabel',[],'yticklabel',[]);
        ylim([350 1400]);
        xlabel([]);ylabel([]);
        title([]); remTickLabel;
        load pplot.mat;
        optsExpFig.Width=3.5/2.54;optsExpFig.Height=3.5/2.54;
        applytofig(gcf,optsExpFig);
        
        tmpD='/Users/hansem/Dropbox (MIT)/figuresRSG2prior/suppleFinal/sF1';
        savefig(gcf,['tp_ts_' num2str(iSess) '_' ehNm{iEH} '_' targNm{iTarg} '.fig']);
%         saveas(gcf,['tp_ts_' num2str(iSess) '_' ehNm{iEH} '_' targNm{iTarg} '.eps']);

        PlotDiff=1; % for inset
        DnPlotResponses2P(d,theory,ctidx,ipidx,PlotDiff); % PlotDiff for the last input
        setFigPos(iEH,iTarg+2);
        
        axis tight;
        xlim([450 1250]);
        set(gca,'xtick',[480 640 800 1000 1200],'ytick',-300:50:300,'xticklabel',[],'yticklabel',[]);
            ylim([-120 120]);
        xlabel([]);ylabel([]);
        title([])
        load pplot.mat;
        applytofig(gcf,optsExpFig);
        savefig(gcf,['tpMinusTs_' num2str(iSess) '_' ehNm{iEH} '_' targNm{iTarg} '.fig']);
%         saveas(gcf,['tpMinusTs_' num2str(iSess) '_' ehNm{iEH} '_' targNm{iTarg} '.eps']);
        
        %% figure 1C: hist
        
        HsHistTpOverlap(d,ctidx);
        setFigPos(iEH,iTarg+4);
        
        set(gca,'ytick',0:5:40,'xtick',[480 640 800 1000 1200],'xticklabel',[],'yticklabel',[],'tickdir','out', 'TickLength', [.02 .02]);
        xlim([350 1400]); ylim([0 40]); remTickLabel; remAxLabel;
        load pplot.mat;
        optsExpFig.Width=3.5/2.54;optsExpFig.Height=3.5/2.54;
        applytofig(gcf,optsExpFig);
        savefig(gcf,['histTpOverlap_' num2str(iSess) '_' ehNm{iEH} '_' targNm{iTarg} '.fig']);
%         saveas(gcf,['histTpOverlap_' num2str(iSess) '_' ehNm{iEH} '_' targNm{iTarg} '.eps']);

    end % for iEH=1:nEH
end % for iTarg=1:nTarg

% ---- Eye Right -----
% 547 trials
% ---- Eye Left -----
% 539 trials
% ---- Hand Right -----
% 571 trials
% ---- Hand Left -----
% 526 trials