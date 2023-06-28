function behAnal_GRSGprior_figure1

% save NEV, NS3

% note
% save data from around 800 trials in 12/7
% data transfer speed warning issue in 12/13? two cells in NS3 data

% 2017/8/24:
%  setting noSet trials as all outlier (tp in those trials are not <0  always)

% 2019/1/10: model comparison with linear model fit

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


idExportFig=0;
load pplot.mat;

nPr=2;    prNm={'Short','Long'};
nEH=2;ehNm={'Eye','Hand'};
nTarg=2; targNm={'Right','Left'};
nTspp=5;Tmat=[]; Tmat{1}=linspace(480,800,nTspp); Tmat{2}=linspace(800,1200,nTspp); Tall=[Tmat{1} Tmat{2}];
nTp=5; 
%% merge data and save into one file
nS=length(fname);

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

%% larger bias for long prior @ 20171117
% do linear regression between mean(tp|ts) and ts for each prior > compare slope
% two prior x two EH slopes from each session

sessUni=unique(sessId);
slope=nan(nPr,nEH,length(sessUni));
for iSess=1:length(sessUni)
    % session data
    idSess=sessId==sessUni(iSess);        
    ts=T(~idOut & idSess);
    tp=t(~idOut & idSess);
    ctidx=2-idShortTrial(~idOut & idSess); % 1 for short, 2 for long
    citdx2=2-idHandEye(~idOut & idSess); % 1 for eye, 2 for hand
    
    % getting mean(tp|ts)
    mtp=nan(nPr,nEH,nTspp);
    for i=1:nPr
        for k=1:nTspp
            for j=1:nEH
            id=ctidx==i &... % idShortTrial
                ts==Tmat{i}(k) &... % ts
                citdx2==j;            % eye hand
            mtp(i,j,k)=mean(tp(id));
            end
        end
    end
    
    % slope
    for i=1:nPr
        for j=1:nEH
            tmp=regstats(squeeze(mtp(i,j,:)),Tmat{i}(:),'linear',{'beta'}); % beta
            slope(i,j,iSess)=tmp.beta(2);
        end
    end
end

% plot
hFig=figure;H=[];
for i=1:size(slope,1)
    for j=1:size(slope,2)
        x=squeeze(slope(i,j,:));
        tmptmpcmap=pplot.cmap{(i-1)*size(slope,2)+j};
        [~,~,hTmp]=histStairs(x,15,0,hFig,tmptmpcmap);H=[H;hTmp];ha;
        plot(mean(x),max(get(gca,'ylim')),'v','color',tmptmpcmap,'markerfacecolor',tmptmpcmap);
    end
end
legend(H,'ShortEye','ShortHand','LongEye','LongHand','location','best'); legend boxoff;

% anova
matPr=repmat([1;2],[1 nEH length(sessUni)]);
matEH=repmat([1 2;1 2],[1 1 length(sessUni)]);

stat=anovan(slope(:),{matPr(:) matEH(:)},'model','interaction','display','on','varname',{'idPr','idEH'});

close all;

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
plotBLS=1; plotTpMTs=1;
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
DnPlotResponses2P(d,theory,ctidx,ipidx,plotBLS,plotTpMTs);

% paper
if plotTpMTs
    set(gca,'xtick',[480 640 800 1000 1200],'ytick',-500:100:500,'xticklabel',[],'yticklabel',[]);
    axis tight; %     ylim([350 1400]);
else
    set(gca,'xtick',[480 640 800 1000 1200],'ytick',[480 640 800 1000 1200],'xticklabel',[],'yticklabel',[]);
    ylim([350 1400]);
end
xlim([450 1250]);
xlabel([]);ylabel([]);
title([])
savefig(gcf,'tp_ts_G.fig');
load pplot.mat;
% optsExpFig.Width=5; %2.6/2.54*2;
% optsExpFig.Height=5; % 1.9/2.54*2;
optsExpFig.Format='eps';
% optsExpFig.Format='eps';exportfig(gcf,'tp_ts_G_small.eps',optsExpFig);
% poster
optsExpFig.Width=4/2.54; % 6.5; % 10/2.54;
optsExpFig.Height=4/2.54; % 6.5; % 7.3/2.54;
if idExportFig
exportfig(gcf,'tp_ts_G.eps',optsExpFig);
end

%% figure 1C: hist

HsHistTpOverlap(d,ctidx);
xlim([800-350 800+350]);
set(gca,'xtick',600:200:1000,'ytick',0:20:40,'tickdir','out'); xlabel([]);ylabel([]);
set(gca,'xticklabel',[],'yticklabel',[]);
savefig(gcf,'histTpOverlap_G.fig');
load pplot.mat;
% optsExpFig.Width=2; % 2.6/2.54*2;
% optsExpFig.Height=2; % 1.9/2.54*2;
optsExpFig.Format='eps';
% exportfig(gcf,'histTpOverlap_G_small.eps',optsExpFig);
% poster
optsExpFig.Width=3; % 10/2.54;
optsExpFig.Height=2; % 7.3/2.54;
if idExportFig
exportfig(gcf,'histTpOverlap_G.eps',optsExpFig);
end

%% figure 1D: bias variance


dBLS.wm=wm; % wFit.w_m;
dBLS.wp=wp; % wFit.w_p;
dBLS.offset1=offset(2); % wFit.offset1; % for long
dBLS.offset2=offset(1); % wFit.offset2;
% with offset, discrepancy is lower; TBD: model comparison with prior-dependent wm,wp
HsPlotBiasVar(d,ctidx,dBLS);

set(gca,'xtick',0:80:160,'ytick',0:80:160); xlabel([]);ylabel([]);

load pplot.mat;
optsExpFig.Width=2.6/2.54;optsExpFig.Height=1.9/2.54;
optsExpFig.Format='eps';
if idExportFig
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
% Parallel pool using the 'local' profile is shutting down.
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% plotTraj('traj_periSet_H_attrition.mat');cc;plotTraj('traj_periSet_G_attrition.mat');
% Elapsed time is 7312.201795 seconds.
%     0.0114    0.1522  -43.3054   14.1497
% 
% ===== 170508 =====
% Parallel pool using the 'local' profile is shutting down.
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 4302.464814 seconds.
%     0.0471    0.1438  -23.9119   27.3071
% 
% ===== 170510 =====
% Parallel pool using the 'local' profile is shutting down.
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 4349.493720 seconds.
%     0.0225    0.1142  -45.9624    1.9498
% 
% ===== 170511 =====
% Parallel pool using the 'local' profile is shutting down.
% Starting parallel pool (parpool) using the 'local' profile ... connected to 4 workers.
% Elapsed time is 5880.654135 seconds.
%     0.0605    0.0865  -49.0636  -10.3220
% 
% ===== 170512 =====
% Parallel pool using the 'local' profile is shutting down.
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

%%

