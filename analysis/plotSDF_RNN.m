function stat=plotSDF_RNN(ws)
v2struct(ws);
if ~exist(dirSingleRNN,'dir')
    mkdir(dirSingleRNN);
end
% x r r2[time x trial x neurons] y idPr ts tp readyOnset2

% input is on 50ms always
% ready is given randomly [100 200] after input on

% trial stage codes (2,3,4,5,6,8/9 for now)
% 1. Fix On
% 2. Fix % photodiode sync
% 3. target On % photodiode sync
% 4. Ready % photodiode sync
% 5. Set % photodiode sync
% 6. Production 
% 7. target acquired 
% 8. reward delay % photodiode sync
% 9. reward: bonusRewDur*1000+rewardDur*1000*max(0.00001,1-abs(productionInterval - interval)/interval/win_fraction ) % photodiode sync
% 10. trial end
% 0. ITI
% 8. incorrect
% 0. Bad

% tInfo: trial#, stageCode, t(blackRock), 
%           idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)]

% idPlot=1;
% output: p value for ANOVA
nAnal=11; % 10;
stat=cell(nAnal,1);

% new stat: poisson regression
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

% old stat
% 1) 2 prior x 2 eyeHand ANOVA [fix fix+500]: 3
%     2) ~ [tOn-250 tOn]
% 3) 2 prior x 2 eyeHand x 2 target ANOVA [tOn tOn+250]
%     4) ~ [ready-250 ready]
% 5) 2 prior x 2 eyeHand [ready ready+400]
%     6) ~ [set-400 set]
% 7) 2 prior x 2 eyeHand ANOVA [ready+480 ready+800] for overlap ts
% 8) tp regression + 2 eyeHand x 2 directions [set set+min(tp)]
%     9) ~ [production-min(tp) production]
% 10) reward regression [rewardOnset onset+100]



%% init
load pplot.mat;
nPr=2;    prNm={'Short','Long'};
nTspp=5;
T=[];
T{1}=linspace(500,820,nTspp); T{2}=linspace(820,1220,nTspp); Tall=[T{1} T{2}];        
        
if idPlot
figure; set(gcf,'position',[0 0 1920 1200],'color','w');
    idNoTrun=1; % truncated by min of trial duration across ts
end
mPlot=1;
nPlot=3;
yLimMat=NaN(mPlot*nPlot,2);iYlim=0;
try
%% 1. Ready 
% 10 ts(short/long) autumn/winter
dur1=149;
if idNoTrun
dur2=max(ts);
else
    dur2=min(ts);
end


tmpR=nan(dur1+dur2+1,ntrials);
for iT=1:ntrials
    if idNoTrun
        tmpR(1:length(-dur1:round(ts(iT))),iT)=r(readyOnset2(iT)+[-dur1:round(ts(iT))],iT,iN);    % time x trials
    else
        tmpR(:,iT)=r(readyOnset2(iT)+[-dur1:dur2],iT,iN);    % time x trials
    end
end


tmpCmap={flipud(autumn(nTspp));... % , flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
            winter(nTspp)}; % , tmp2(1:nTspp,:)}; % longEye, longHand
        legNmTmp=cell(nPr,length(T{1}));hMat=[];
for i=1:nPr
    for j=1:length(T{i})
        if idPlot, hTmp=subplot(mPlot,nPlot,1);hold all;end
        tmpId=T{i}(j)==ts & idPr==(i-1);
        if idPlot, 
            h=shadedErrorBar([-dur1:(dur2)]',mean(tmpR(:,tmpId),2,'omitnan'),sem(tmpR(:,tmpId),2,'omitnan'),{'-','color',tmpCmap{i}(j,:),'linewidth',1},1); 
% h.mainLine=plot([-dur1:dur2]',mean(tmpR(:,tmpId),2),'-','color',tmpCmap{i}(j,:),'linewidth',1); 
            hMat=[hMat;h.mainLine];
            legNmTmp{i,j}= [prNm{i} num2str(T{i}(j))];
        end
    end % j
end % i        
axis tight; xlabel('time from ready (ms)'); %ylabel('spike rate (/s)');
% plotVertical(gca);   
iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;

legNmTmp=legNmTmp';
legend(hMat,legNmTmp(:),'location','southeast'); legend boxoff;drawnow;



%% 2. Set % photodiode sync [-480/800/1200 min(tp)]
if idNoTrun
dur1=max(ts);
dur2=max(tp);
else
dur1=min(ts);
dur2=min(tp);
end
tmpR=nan(dur1+dur2+1,ntrials);
for iT=1:ntrials
    if idNoTrun
        tmpR(dur1-ts(iT)+[1:length([-round(ts(iT)):round(tp(iT))])],iT)=r(readyOnset2(iT)+ts(iT)+[-round(ts(iT)):round(tp(iT))],iT,iN);    % time x trials
%         disp(length([-round(ts(iT)):round(tp(iT))]))
    else
        tmpR(:,iT)=r(readyOnset2(iT)+ts(iT)+[-dur1:dur2],iT,iN);    % time x trials
    end
end

tmpCmap={flipud(autumn(nTspp));... % , flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
            winter(nTspp)}; % , tmp2(1:nTspp,:)}; % longEye, longHand
        legNmTmp=cell(nPr,length(T{1}));hMat=[];
for i=1:nPr
    for j=1:length(T{i})
        if idPlot, subplot(mPlot,nPlot,2);hold all;end
        tmpId=T{i}(j)==ts & idPr==(i-1);
        if idPlot, 
            h=shadedErrorBar([-dur1:(dur2)]',mean(tmpR(:,tmpId),2,'omitnan'),sem(tmpR(:,tmpId),2,'omitnan'),{'-','color',tmpCmap{i}(j,:),'linewidth',1},1); 
% h.mainLine=plot([-dur1:dur2]',mean(tmpR(:,tmpId),2),'-','color',tmpCmap{i}(j,:),'linewidth',1); 
            hMat=[hMat;h.mainLine];
            legNmTmp{i,j}= [prNm{i} num2str(T{i}(j))];
        end
    end % j
end % i        
axis tight; xlabel('time from set (ms)'); %ylabel('spike rate (/s)');
% plotVertical(gca);   
iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;

% hMat=hMat';legNmTmp=legNmTmp';
% legend(hMat(:),legNmTmp(:),'location','best'); legend boxoff;drawnow;


%% 3. Production [-tp 50]
if idNoTrun
    dur1=max(tp);
else
    dur1=min(tp);
end
dur2=round(.5*(time-(max(readyOnset2)+max(ts)+max(tp))));

tmpR=nan(dur1+dur2+1,ntrials);
for iT=1:ntrials
    if idNoTrun
        tmpR(dur1-tp(iT)+[1:length([-round(tp(iT)):(dur2)])],iT)=r(readyOnset2(iT)+ts(iT)+tp(iT)+[-round(tp(iT)):dur2],iT,iN);    % time x trials
    else
        tmpR(:,iT)=r(readyOnset2(iT)+ts(iT)+tp(iT)+[-dur1:dur2],iT,iN);    % time x trials
    end
end

tmpCmap={flipud(autumn(nTspp));... % , flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
            winter(nTspp)}; % , tmp2(1:nTspp,:)}; % longEye, longHand
        legNmTmp=cell(nPr,length(T{1}));hMat=[];
for i=1:nPr
    for j=1:length(T{i})
        if idPlot, subplot(mPlot,nPlot,3);hold all;end
        tmpId=T{i}(j)==ts & idPr==(i-1);
        if idPlot, 
            h=shadedErrorBar([-dur1:dur2]',mean(tmpR(:,tmpId),2,'omitnan'),sem(tmpR(:,tmpId),2,'omitnan'),{'-','color',tmpCmap{i}(j,:),'linewidth',1},1); 
% h.mainLine=plot([-dur1:dur2]',mean(tmpR(:,tmpId),2),'-','color',tmpCmap{i}(j,:),'linewidth',1); 
            hMat=[hMat;h.mainLine];
            legNmTmp{i,j}= [prNm{i} num2str(T{i}(j))];
        end
    end % j
end % i        
axis tight; xlabel('time from go (ms)'); %ylabel('spike rate (/s)');
% plotVertical(gca);   
iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;


%%
%matching ylim
if idPlot, 
for i=1:(mPlot*nPlot)
    h=subplot(mPlot,nPlot,i);
%     set(h,'position',[0.2*rem(i-1,5)+0.02*(i==1||i==nPlot+1) 0.5*(i<6)+0.03 0.18-0.02*(i==1||i==nPlot+1) 0.42]);
    ylim([min(yLimMat(:,1)) max(yLimMat(:,2))]);
    plotVertical(gca);
    if i==1, plotVertical(gca,T{1}(3),[1 0 0]);plotVertical(gca,T{2}(3),[0 0 1]); end;
    if i~=1 %&& i~=(nPlot+1), 
        set(gca,'yticklabel',[]); end;
end
end
% disp('');
optsExpFig=rmfield(optsExpFig,'Width');optsExpFig=rmfield(optsExpFig,'Height');
optsExpFig.Format='png';optsExpFig.Width=9;optsExpFig.Height=5;
% optsExpFig.Reference=hTmp;
optsExpFig.Bounds='loose';
exportfig(gcf,[dirSingleRNN '/' num2str(iN) '.png'],optsExpFig); % 'color','rgb','format','png','Reference',hTmp); % ,'Height',3,'Width',6); 

catch
    disp([]);

end
%% back up
