function stat=plotPSTH(NEV,id,tInfo,eventId,idPlot,varargin)

% input: NEV data, id(electrode,Unit,start trial,end trial), eventID, idPlot, idReg(varargin)
% use Gaussian kernel for now (note also photodiode sync is not used forg now)

% TBD 20170511: raster for each time points (also for more informative selection of start/end trials)
%       spike waveform
% ts tuning
% tp correlation

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
%           idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12), idNoset(13, if exist), idSuccessNoSet(13, if exist)]

%% FINAL regression
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

% % idPlot=1;
% % output: p value for ANOVA
% nAnal=16; % 10;
% stat=cell(nAnal,1);

% new stat: poisson regression
%     '1)    ' prior x eyeHand 22GLM [fix fix+500]'                                  0.24361      0.06015
%     '2)    ' prior x eyeHand 22GLM [tOn-250 tOn]'                                  0.16692     0.021053
%     '3)    ' prior x eyeHand x target 222GLM [tOn tOn+250]'                        0.06015    0.0015038
%     '4)    ' prior x eyeHand x target 222GLM [ready-250 ready]'                   0.046617            0
%     '5) ' prior x eyeHand x target 222GLM   '  [ready ready+480]'                              0.26617     0.082707
%     '6)   ' ts + eyeHand x target regression ' [set-480 set]'                                  0.26917     0.096241
%     '7) prior x eyeHand x target 222GLM [ready+480 ready+800] for overlap ts inc. all long prior?                           0.29624      0.14436       
%     ?8) ' prior x eyeHand x target 222GLM   [set-480 set]'
%     '9) intercept tp + eyeHand x directions regression [set set+min(tp)]'      0.8812      0.73835
%     '10) intercept tp + eyeHand x directions regression [production-min(tp) production]';
%     '11) intercept reward regression [rewardOnset onset+100]'                 0.99098      0.95489
%     ' rewarded regression [rewardOnset onset+100]'                             0.1203     0.013534
%     ' rewDur*rewarded regression [rewardOnset onset+100]'                    0.079699    0.0045113
% 12)    ' prior x eyeHand 22GLM [fix-100 fix]'  
% 13) 100 sliding window  ' prior x eyeHand 22GLM [fix-100 fix+500]'
% 14) 100 sliding window   prior x eyeHand x target 222GLM [tOn-250 tOn+250]' 
% 15) 100 sliding window   prior x eyeHand x target 222GLM [ready-250 ready+480]'    
% 16 ) 100 sliding window  [set-480 set+100]
%         % log(spikeRate_100_set) ~ b1 + b2*t_s +  b3*idShort + b4*idEye + b5*idRight ...
%         %                                         + b6*t_s*idShort + b7*t_s*idEye + b8*t_s*idRight ...
%         %                                           + b9*t_s^2 + b10*t_s^2* *idShort + b11*t_s^2*idEye + b12*t_s^2*idRight

% % % old new stat: poisson regression
% % %     '1)    ' prior x eyeHand 22GLM [fix fix+500]'                                  0.24361      0.06015
% % %     '2)    ' prior x eyeHand 22GLM [tOn-250 tOn]'                                  0.16692     0.021053
% % %     '3)    ' prior x eyeHand x target 222GLM [tOn tOn+250]'                        0.06015    0.0015038
% % %     '4)    ' prior x eyeHand x target 222GLM [ready-250 ready]'                   0.046617            0
% % %     '5) ' prior x eyeHand x target 222GLM   ' prior x eyeHand 22ANOVA [ready ready+480]'                              0.26617     0.082707
% % %     '6)   ' ts + eyeHand x target regression ' prior x eyeHand 22ANOVA [set-480 set]'                                  0.26917     0.096241
% % %     '7) prior x eyeHand x target 222GLM [ready+480 ready+800] for overlap ts inc. all long prior?                           0.29624      0.14436       
% % %     ?8) ' prior x eyeHand x target 222GLM   ' prior x eyeHand 22ANOVA [set-480 set]'
% % %     '9) intercept tp + eyeHand x directions regression [set set+min(tp)]'      0.8812      0.73835
% % %     '10) intercept tp + eyeHand x directions regression [production-min(tp) production]';
% % %     '11) intercept reward regression [rewardOnset onset+100]'                 0.99098      0.95489
% % %     ' rewarded regression [rewardOnset onset+100]'                             0.1203     0.013534
% % %     ' rewDur*rewarded regression [rewardOnset onset+100]'                    0.079699    0.0045113
% % 
% % % old stat
% % % 1) 2 prior x 2 eyeHand ANOVA [fix fix+500]: 3
% % %     2) ~ [tOn-250 tOn]
% % % 3) 2 prior x 2 eyeHand x 2 target ANOVA [tOn tOn+250]
% % %     4) ~ [ready-250 ready]
% % % 5) 2 prior x 2 eyeHand [ready ready+400]
% % %     6) ~ [set-400 set]
% % % 7) 2 prior x 2 eyeHand ANOVA [ready+480 ready+800] for overlap ts
% % % 8) tp regression + 2 eyeHand x 2 directions [set set+min(tp)]
% % %     9) ~ [production-min(tp) production]
% % % 10) reward regression [rewardOnset onset+100]



%% init
stat=[];
if ~isempty(varargin)
    idReg=varargin{1};
else
idReg=0; % wheter to do regression analysis
end

% sliding window: step size 1ms
winSize=100; 

idE=id(1); idU=id(2); % unitMap=[0 1 2 3 255];  idU=unitMap(id(2));
load pplot.mat;
sampF=double(NEV.MetaTags.SampleRes); % 30000 /s
sp=double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==idE&NEV.Data.Spikes.Unit==idU))/sampF*1000; % blackrock sample unit into ms
if isempty(sp), disp('no channel/unit'); return; end;
% wf=double(NEV.Data.Spikes.Waveform(:,NEV.Data.Spikes.Electrode==idE&NEV.Data.Spikes.Unit==idU)); % [48 x spikes]
% strdate=datestr(NEV.MetaTags.DateTime,'yymmdd');
clear NEV;
% unit 0 1 2 3 255>1 2 3 4 5?

% removing trials without no neural data (given uncertain end trials ID)
tLastSp=max(sp); % [ms]
tInfo=tInfo(tInfo(:,1)>=id(3)&tInfo(:,1)<=id(4)&tInfo(:,3)*1000<tLastSp,:); % start end trials
% tInfo=tInfo(tInfo(:,1)>=id(3)&tInfo(:,1)<=id(4),:); % start end trials

% treating noSet experiment
if size(tInfo,2)>12
    noSet=1;idSuccessNoSet=tInfo(:,14);idNoSet=tInfo(:,13);waitDur=400;
else
    noSet=0;
end

nPr=2;    prNm={'Short','Long'};
nEH=2;ehNm={'Eye','Hand'};
nTarg=2; targNm={'Right','Left'};
nTspp=5; T{1}=unique(tInfo(tInfo(:,4)==1,7)); T{2}=unique(tInfo(tInfo(:,4)==0,7));
nSplitTp=5; t=tInfo(tInfo(:,8)>0,8); tSort=sort(t); tBound=tSort([1 round([(1/nSplitTp):(1/nSplitTp):1]*length(t))]); % 0 0.2 0.4 ... 1
nSplitRew=5; rewValid=tInfo(tInfo(:,8)>0&tInfo(:,12)>0,end); rSort=sort(rewValid); rBound=rSort([1 round([(1/nSplitRew):(1/nSplitRew):1]*length(rewValid))]); % 0 0.33 0.66 1

if idPlot
figure; set(gcf,'position',[0 0 1920 1200],'color','w');
end
mPlot=2;
nPlot=5;
yLimMat=NaN(mPlot*nPlot,2);iYlim=0;

%% Gaussian kernel for now (note also photodiode sync is not used forg now)
td=30; nTd=3; xFiltG=-td*nTd:1:td*nTd; gaussampFilter=exp(-(xFiltG.^2)/(2*(td^2))); gaussampFilter=gaussampFilter/sum(gaussampFilter); % 40
% tg=1; td=20; xFilt=0:1:100; % smoothing as in Hanes et al., JNP98
% epspFilter=(1-exp(-xFilt/tg)).*exp(-xFilt/td); epspFilter=[zeros(length(epspFilter)-1,1);epspFilter(:)]; epspFilter=epspFilter/sum(epspFilter);

%% different categorization for different trials events
% % waveform
% figure(21); set(gcf,'position',pplot.rect2_1,'color','w'); hold all;
% nSWF=size(wf,1);
% plot(1000*[0:(nSWF-1)]'/sampF,wf(:,round(linspace(1,size(wf,2),min([size(wf,2) 100])))));
% h=shadedErrorBar(1000*[0:(nSWF-1)]'/sampF,mean(wf,2),sem(wf,2),{'-','color',[0 0 0],'linewidth',2},1); % std(wf,0,2)
% legend(h.mainLine,[strdate '-ch' num2str(id(1)) '-' num2str(id(2))],'location','SouthEast'); legend boxoff;
% axis tight; xlabel('time (ms)'); ylabel('amplitude');drawnow;

%% 2. Fix % photodiode sync [-100 500]
% ShortEye, ShortHand, LongEye, LongHand (Red, Orange, Blue, Green)
if eventId==2 || eventId==0 % 0 for plot all
    
    [idState,idSw,iSw,tNB,tmp]=sortSwitch(tInfo);%idState% 0123 for shortEye, shortHand, LongEye, LongHand; tmp=unique(tInfo(:,[1 4 5 8]),'rows','stable'); % t#, idLS, idHE
    tSw=tmp(find(diff(idState)~=0)+1,1);
    tmpIdSw=true(size(tInfo,1),1); for iT=1:length(tSw), tmpIdSw(tSw(iT)==tInfo(:,1))=false; end; % no false for switching trials
    
    dur1=200; % 100; % ms
    dur2=500; tDur=(-dur1+1):dur2;
    nr=[]; iPr=[]; iEH=[]; nr2=[]; nr3=[];% for stat
    rectWin=rect_filter(winSize,dur1+dur2,false); % [time x time] row NOT normalized with last argin to make just sum(spikes)
    % note winSize/2 at start and end is not useful
                    iPr4=[]; iEH4=[]; iRL4=[]; iSw=[];
    if idPlot, hTmp=subplot(mPlot,nPlot,1); end
%     figure(11); set(gcf,'position',pplot.rect1_1,'color','w'); 
    hold all; tmpCmap={[1 0 0] [1 .5 0];[0 0 1] [0 1 0]};hMat=[];
    legName={'ShortEye','ShortHand';'LongEye','LongHand'};% hMat=cell(nPr,nEH);
    
    r=cell(nPr,nEH);
    for i=1:nPr,
        for j=1:nEH,
            tmpId=tInfo(:,2)==1 & ... % 2 &... % event id 2.fix; 1 for fixOn (not photodioide) but goal is to check context effect before fixOn
                tInfo(:,4)==(2-i) &... % idShortTrial
                tInfo(:,5)==(2-j) &...  % idHandEye
                tInfo(:,8)>0; % only valid trials for now
            tE=tInfo(tmpId,3)*1000; % into ms
            t1=tE-dur1;            t2=tE+dur2;
            r{i,j}=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
            tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
            for k=1:nnz(tmpId)
                tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
                if ~isempty(tSp),r{i,j}(k,ceil(tSp))=1;end;
                tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
                if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
            end % k for trials
            nr=[nr;sum(r{i,j}(:,(dur1+1):end),2)]; % mean instead of sum for spike rate
            nr2=[nr2;sum(r{i,j}(:,1:dur1),2)];  % 12)    ' prior x eyeHand 22GLM [fix-100 fix]'  
            nr3=[nr3;[rectWin*r{i,j}']'];% [nr3;[1000*rectWin*r{i,j}']'];  % [trial x winId] 13) 100 sliding window  ' prior x eyeHand 22GLM [fix-100 fix+500]'
            iPr=[iPr; repmat(prNm(i),[nnz(tmpId) 1])];iEH=[iEH; repmat(ehNm(j),[nnz(tmpId) 1])];
                                            iPr4=[iPr4; repmat(strcmp(prNm{i},'Short'),[nnz(tmpId) 1])];iEH4=[iEH4; repmat(strcmp(ehNm{j},'Eye'),[nnz(tmpId) 1])];
                                            iSw=[iSw;tmpIdSw(tmpId)];
            rConv=conv2(tmpR,gaussampFilter,'valid')*1000;
            if idPlot, 
            h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap{i,j},'linewidth',2},1); hMat(i,j)=h.mainLine;
            end
        end % j for EyeHand
    end % i for ShortLong
    if idPlot, 
    axis tight; xlabel('time from fixation (ms)'); ylabel('spike rate (/s)');%plotVertical(gca);
    legend(hMat(:),legName(:),'location','best'); legend boxoff;    drawnow;
    end
    iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;
    
    if idReg
    % 1) 1)    ' prior x eyeHand 22GLM [fix fix+500]
    nrT=table(nr(:,1),iPr,iEH,'VariableNames',{'spikeRate_fix_500' 'idShort' 'idEye'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{1}=fitglm(nrT,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye'},'ResponseVar','spikeRate_fix_500'); % 'Exclude' for outlier
%     stat{1}=anovan(nr,{iPr iEH},'model',2,'display','off','varnames',{'idShort' 'idEye'});
% %     nrT=table(nr(:,1),logical(nr(:,2)),logical(nr(:,3)),'VariableNames',{'spikeRate' 'idShort' 'idEye'});

    % GLM: Coefficients (table: regressors x estimate/SE/tStat/pValue), ModelCriterion.AIC/BIC, Deviance, Rsquared.Ordinary/Adjusted/Deviance, ResponseName, PredictorNames, Variables, SSE/SST/SSR, plotResiduals

    % 12)    ' prior x eyeHand 22GLM [fix-100 fix]'  
    nrT2=table(nr2(:,1),iPr,iEH,'VariableNames',{'spikeRate_fix_m100' 'idShort' 'idEye'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{12}=fitglm(nrT2,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye'},'ResponseVar','spikeRate_fix_m100'); % 'Exclude' for outlier
    
    % 13) 100 sliding window  ' prior x eyeHand 22GLM [fix-100 fix+500]'
%     stat{13}=[];
%     X=x2fx([iPr4 iEH4],[1  0;0 1 ;... % 1st level
%         1 1 ]); % interaction with ts
%     for iWin=1:size(nr3,2) % (dur1+dur2)
%         %     nrT3=table(nr3(:,iWin),iPr,iEH,'VariableNames',{'spikeRate_perifix' 'idShort' 'idEye'}); % logical(nr(:,2)),logical(nr(:,3))
%         %     glmTmp=fitglm(nrT3,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye'},'ResponseVar','spikeRate_perifix'); % 'Exclude' for outlier
%         %
%         %     stat{13}= cat(3,stat{13},[table2array(glmTmp.Coefficients(:,1))';... % 1 x nCoeff
%         %    table2array(glmTmp.Coefficients(:,end))']);
%         
%         [B,DEV,glmStat]=glmfit(X,nr3(:,iWin),'poisson'); % glmStat.p
%         stat{13}= cat(3,stat{13},[B(:)';... % 1 x nCoeff
%             glmStat.p(:)']);
% 
%     % for later ROC analysis
% % [FA,Hit,AUC,pBoot]=rocCurve([lgG(tTmp,:) lgC(tTmp,:)],[ones(size(lgG(tTmp,:))) zeros(size(lgC(tTmp,:)))],0,flagBoot); % AreaUnderROC/ChoiceProbability: switch sign for identical threshold
%     end
    
    % FINAL:log(nSpikesWin100_fixOn[m200_500]) ~ 1 + idLong + idHand + idSwitch +idLong*idHand (5)
    stat{13}=[];
    X=x2fx([iPr4 iEH4 iSw],[1  0 0;0 1 0;0 0 1;... % 1st level
        1 1 0],... % idShort*idEye
        [1 2 3]); % categorial
    for iWin=1:size(nr3,2) % (dur1+dur2)
        [B,DEV,glmStat]=glmfit(X,nr3(:,iWin),'poisson'); % glmStat.p
        stat{13}= cat(3,stat{13},[B(:)';... % 1 x nCoeff
            glmStat.p(:)']);
        
        % for later ROC analysis
        % [FA,Hit,AUC,pBoot]=rocCurve([lgG(tTmp,:) lgC(tTmp,:)],[ones(size(lgG(tTmp,:))) zeros(size(lgC(tTmp,:)))],0,flagBoot); % AreaUnderROC/ChoiceProbability: switch sign for identical threshold
    end
    end % idReg
end

%% SWITCH 2. Fix % photodiode sync [-100 500]
% ShortEye, ShortHand, LongEye, LongHand (Red, Orange, Blue, Green)
if eventId==2 || eventId==0 % 0 for plot all
    % indexing trial number in blocks
    [idState,idSw,iSw,tNB,tmp]=sortSwitch(tInfo);
    iSwEnd=[7;4;10;11;1;8;2;12;5;6;9;3]; % to shortEye, shortHand, LongEye, LongHand x  acrPr, acrMod, acrBot
    dur1=100; % ms
    dur2=500; tDur=(-dur1+1):dur2;
    if idPlot, subplot(mPlot,nPlot,10);end
%         figure(22); set(gcf,'position',pplot.rect2_2,'color','w'); 
        hold all; tmpCmap={[1 0 0] [1 0 0]*2/3 [1 0 0]/3 [1 .5 0] [1 .5 0]*2/3 [1 .5 0]/3 [0 0 1] [0 0 1]*2/3 [0 0 1]/3 [0 1 0] [0 1 0]*2/3 [0 1 0]/3};hMat=[];
        legName={'acrPr2ShE','acrMod2ShE','acrBoth2ShE',...
            'acrPr2ShH','acrMod2ShH','acrBoth2ShH',...
            'acrPr2LoE','acrMod2LoE','acrBoth2LoE',...
            'acrPr2LoH','acrMod2LoH','acrBoth2LoH'};% hMat=cell(nPr,nEH); 

    r=cell(length(iSwEnd)); 
    for i=1:length(iSwEnd), 
        tmpIdTN=tmp(iSwEnd(i)==iSw & tmp(:,end)>0,1); % trial # with valid tp
        r{i}=zeros(length(tmpIdTN),dur1+dur2); % [trial x time]
        tmpR=zeros(length(tmpIdTN),dur1+dur2+td*nTd*2); % for convoluted plot
        
        for j=1:length(tmpIdTN)
            tmpId=tInfo(:,2)==2 &... % event id 2.fix
                tInfo(:,1)==tmpIdTN(j) &... % trial #
                tInfo(:,8)>0; % only valid trials for now
%             if ~isempty(tmpId)
            tE=tInfo(tmpId,3)*1000; % into ms
            t1=tE-dur1;            t2=tE+dur2; % disp(size(t1));
                tSp=sp(t1<sp & sp<=t2)-t1;
                if ~isempty(tSp),r{i}(j,ceil(tSp))=1;end;
                tSp2=sp(t1-td*nTd<sp & sp<=t2+td*nTd)-(t1-td*nTd);
                if ~isempty(tSp2),tmpR(j,ceil(tSp2))=1;end;
%             end
        end % j=1:length(tmpIdTN)   
        rConv=conv2(tmpR,gaussampFilter,'valid')*1000;
                if idPlot, hMat(i)=plot(tDur(:)',mean(rConv,1),'-','color',tmpCmap{i},'linewidth',2); end
%                 h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-'
%                 ,'color',tmpCmap{i},'linewidth',2},1); hMat(i,j)=h.mainLine;
        
    end % i for all switch types
    if idPlot, 
    axis tight; xlabel('time from fixation (ms)'); %ylabel('spike rate (/s)');plotVertical(gca);
    legend(hMat(:),legName(:),'location','best'); legend boxoff;    drawnow;
    end
    yLimMat(end,:)=ylim;
end

%% 3. target On % photodiode sync [-250 250]
% ShortEye, ShortHand, LongEye, LongHand x target (0/180:solid:dashed)
if eventId==3 || eventId==0 % 0 for plot all
    dur1=250; % ms
    dur2=250; tDur=(-dur1+1):dur2;
        nr=[]; nr2=[];iPr=[]; iEH=[]; iRL=[]; nr3=[];% for stat
rectWin=rect_filter(winSize,dur1+dur2,false); % [time x time] row NOT normalized
                iPr4=[]; iEH4=[]; iRL4=[];
    if idPlot, subplot(mPlot,nPlot,6);end
%     figure(12); set(gcf,'position',pplot.rect1_2,'color','w'); 
    hold all; tmpCmap={[1 0 0] [1 .5 0];[0 0 1] [0 1 0]}; linestyle={'-','--'};hMat=[];
    legName=cell(nPr,nEH,nTarg);
    legName{1,1,1}='ShortEye0';
    legName{1,1,2}='ShortEye180';
    legName{1,2,1}='ShortHand0';
    legName{1,2,2}='ShortHand180';
    legName{2,1,1}='LongEye0';
    legName{2,1,2}='LongEye180';
    legName{2,2,1}='LongHand0';
    legName{2,2,2}='LongHand180'; %hMat=cell(nPr,nEH,nTarg);
    
    r=cell(nPr,nEH,nTarg);
    for i=1:nPr,
        for j=1:nEH,
            for l=1:nTarg
                tmpId=tInfo(:,2)==3 &...  % event id 3.tOn
                    tInfo(:,4)==(2-i) &... % idShortTrial
                    tInfo(:,5)==(2-j) &... % idHandEye
                    tInfo(:,6)==(l-1)*180 &... % target location
                    tInfo(:,8)>0; % only valid trials for now
                tE=tInfo(tmpId,3)*1000; % into ms
                t1=tE-dur1;            t2=tE+dur2;
                
                r{i,j,l}=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
                tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
                for k=1:nnz(tmpId)
                    tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
                    if ~isempty(tSp),r{i,j,l}(k,ceil(tSp))=1;end;
                    tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
                if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
                end % k for trials
                rConv=conv2(tmpR,gaussampFilter,'valid')*1000;
                if idPlot, 
                if j==1 & l==1
                    h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{linestyle{l},'color',tmpCmap{i,j},'linewidth',2},1); hMat(i,j,l)=h.mainLine;
                else
                    hMat(i,j,l)=plot(tDur(:)',mean(rConv,1),linestyle{l},'color',tmpCmap{i,j},'linewidth',2);
                end
                end
                nr=[nr;sum(r{i,j,l}(:,1:dur1),2)]; nr2=[nr2;sum(r{i,j,l}(:,(dur1+1):end),2)]; 
                nr3=[nr3;[rectWin*r{i,j,l}']'];%[nr3;[1000*rectWin*r{i,j,l}']'];  % [trial x winId] % 14) 100 sliding window   prior x eyeHand x target 222GLM [tOn-250 tOn+250]' 
                iPr=[iPr; repmat(prNm(i),[nnz(tmpId) 1])];iEH=[iEH; repmat(ehNm(j),[nnz(tmpId) 1])];iRL=[iRL; repmat(targNm(l),[nnz(tmpId) 1])];
                                iPr4=[iPr4; repmat(strcmp(prNm{i},'Short'),[nnz(tmpId) 1])];iEH4=[iEH4; repmat(strcmp(ehNm{j},'Eye'),[nnz(tmpId) 1])];iRL4=[iRL4; strcmp(targNm(1+(tInfo(tmpId,6)==180))','Right')];
            end % l for target
        end % j for EyeHand
    end % i for ShortLong
    if idPlot, 
    axis tight; xlabel('time from target on (ms)'); ylabel('spike rate (/s)');%plotVertical(gca);
    legend(hMat(:),legName(:),'location','best'); legend boxoff;drawnow;
    end
    iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;

    if idReg
    %     2)    ' prior x eyeHand 22GLM [tOn-250 tOn]
    nrT=table(nr(:,1),iPr,iEH,'VariableNames',{'spikeRate_250_tOn' 'idShort' 'idEye'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{2}=fitglm(nrT,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye'},'ResponseVar','spikeRate_250_tOn'); % 'Exclude' for outlier
    %             stat{2}=anovan(nr,{iPr iEH},'model',2,'display','off','varnames',{'idShort' 'idEye'});
    % 3)    ' prior x eyeHand x target 222GLM [tOn tOn+250]
    nrT2=table(nr2(:,1),iPr,iEH,iRL,'VariableNames',{'spikeRate_tOn_250' 'idShort' 'idEye' 'idRight'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{3}=fitglm(nrT2,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar','spikeRate_tOn_250'); % 'Exclude' for outlier
    %         stat{3}=anovan(nr2,{iPr iEH iRL},'model',3,'display','off','varnames',{'idShort' 'idEye' 'idRight'});
    
    % 14) 100 sliding window   prior x eyeHand x target 222GLM [tOn-250 tOn+250]'
%     stat{14}=[];
%     X=x2fx([iPr4 iEH4 iRL4],[1 0 0;0 1 0 ;0 0 1;... % 1st level
%         1 1 0;1 0 1; 0 1 1]); % interaction with ts
%     for iWin=1:size(nr3,2) % (dur1+dur2)
%         %     nrT3=table(nr3(:,iWin),iPr,iEH,iRL,'VariableNames',{'spikeRate_periTon' 'idShort' 'idEye' 'idRight'});  % logical(nr(:,2)),logical(nr(:,3))
%         %         glmTmp=fitglm(nrT3,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar','spikeRate_periTon'); % 'Exclude' for outlier
%         %
%         %     stat{14}= cat(3,stat{14},[table2array(glmTmp.Coefficients(:,1))';... % 1 x nCoeff
%         %    table2array(glmTmp.Coefficients(:,end))']);
%         [B,DEV,glmStat]=glmfit(X,nr3(:,iWin),'poisson'); % glmStat.p
%         stat{14}= cat(3,stat{14},[B(:)';... % 1 x nCoeff
%             glmStat.p(:)']);
%     end
    
    % FINAL: log(nSpikesWin100_tOn[m250_250]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
    stat{14}=[];
    X=x2fx([iPr4 iEH4 iRL4],[1 0 0;0 1 0 ;0 0 1;... % 1st level
        1 1 0;1 0 1; 0 1 1],... % interaction with ts
        [1 2 3]); % categorical
    for iWin=1:size(nr3,2) % (dur1+dur2)
        [B,DEV,glmStat]=glmfit(X,nr3(:,iWin),'poisson'); % glmStat.p
        stat{14}= cat(3,stat{14},[B(:)';... % 1 x nCoeff
            glmStat.p(:)']);
    end
    end % idReg
end

%% 4. Ready % photodiode sync [-200 480/800/1200]
% 10 ts(short/long) for each eye, hand (autumn/winter, spring/summer; use only 200 to avoid yellow)     nTspp=5;
if eventId==4 || eventId==0 % 0 for plot all
    dur1=250; % ms
    periReadyWindow=480; % 240; %T{1}(1); % 400; %%%%%%
    periReadyWindow2=T{1}(end);
    periReadyWindow3=T{1}(1);
                nr=[]; nr2=[];nr3=[];nr4=[];iPr=[]; iEH=[]; iRL=[]; ts=[];% for stat
                iPr3=[]; iEH3=[]; iRL3=[];
                iPr4=[]; iEH4=[]; iRL4=[];
rectWin=rect_filter(winSize,dur1+periReadyWindow,false); % [time x time] row normalized

%     figure(101); set(gcf,'position',pplot.rect1_3,'color','w'); hold all; % eye
%     figure(102); set(gcf,'position',pplot.rect2_3,'color','w'); hold all; % hand
    tmp1=spring(round(nTspp*1.5));
    tmp2=summer(round(nTspp*1.5));
    tmpCmap={flipud(autumn(nTspp)), flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
        winter(nTspp), tmp2(1:nTspp,:)}; % longEye, longHand
    legName=cell(nPr,nTspp,nEH);
hMat=[];
    %hMat=cell(nPr,nTspp,nEH);
    
    r=cell(nPr,nTspp,nEH);
    for l=1:nEH        
%         figure(100+l);
    if idPlot, subplot(mPlot,nPlot,2+(l-1)*nPlot);hold all;end
        for i=1:nPr,
            for j=1:nTspp,
                dur2=T{i}(j); tDur=(-dur1+1):dur2;
                legName{i,j,l}=[prNm{i} num2str(T{i}(j))];

                if noSet
                    tmpId=tInfo(:,2)==4 &...  % event id 4.ready
                        tInfo(:,4)==(2-i) &... % idShortTrial
                        tInfo(:,5)==(2-l) &... % idHandEye
                        tInfo(:,7)==T{i}(j) &... % ts
                        tInfo(:,8)>0 &... % only valid trials for now
                        ~idNoSet; % NoSet trials has max ts of each prior
                else % no NoSet trials
                    tmpId=tInfo(:,2)==4 &...  % event id 4.ready
                        tInfo(:,4)==(2-i) &... % idShortTrial
                        tInfo(:,5)==(2-l) &... % idHandEye
                        tInfo(:,7)==T{i}(j) &... % ts
                        tInfo(:,8)>0; % only valid trials for now
                end
                tE=tInfo(tmpId,3)*1000; % into ms
                t1=tE-dur1;            t2=tE+dur2;
                
                r{i,j,l}=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
                tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
                for k=1:nnz(tmpId)
                    tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
                    if ~isempty(tSp),r{i,j,l}(k,ceil(tSp))=1;end;
                    tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
                    if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
                end % k for trials
                rConv=conv2(tmpR,gaussampFilter,'valid')*1000;            
                if (i==1 & j==nTspp) || (i==nPr & j==1), lw=4; % overlapped interval
                    if idPlot, h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap{i,l}(j,:),'linewidth',lw},1); hMat(i,j,l)=h.mainLine;end
                    nr3=[nr3;sum(r{i,j,l}(:,(dur1+periReadyWindow3):(dur1+periReadyWindow2)),2)]; % 250+480:250+800
                    iPr3=[iPr3; repmat(prNm(i),[nnz(tmpId) 1])];iEH3=[iEH3; repmat(ehNm(l),[nnz(tmpId) 1])];iRL3=[iRL3; targNm(1+(tInfo(tmpId,6)==180))'];
                else lw=1.5;
                    if idPlot, hMat(i,j,l)=plot(tDur(:)',mean(rConv,1),'color',tmpCmap{i,l}(j,:),'linewidth',lw);end
                    if i==nPr % inc. all long prior
                        nr3=[nr3;sum(r{i,j,l}(:,(dur1+periReadyWindow3):(dur1+periReadyWindow2)),2)];
                        iPr3=[iPr3; repmat(prNm(i),[nnz(tmpId) 1])];iEH3=[iEH3; repmat(ehNm(l),[nnz(tmpId) 1])];iRL3=[iRL3; targNm(1+(tInfo(tmpId,6)==180))'];
                    end
                end;            
                nr=[nr;sum(r{i,j,l}(:,1:dur1),2)]; nr2=[nr2;sum(r{i,j,l}(:,(dur1+1):(dur1+periReadyWindow)),2)];
                iPr=[iPr; repmat(prNm(i),[nnz(tmpId) 1])];iEH=[iEH; repmat(ehNm(l),[nnz(tmpId) 1])];iRL=[iRL; targNm(1+(tInfo(tmpId,6)==180))'];
                ts=[ts; repmat(T{i}(j),[nnz(tmpId) 1])];
                                iPr4=[iPr4; repmat(strcmp(prNm{i},'Short'),[nnz(tmpId) 1])];iEH4=[iEH4; repmat(strcmp(ehNm{l},'Eye'),[nnz(tmpId) 1])];iRL4=[iRL4; strcmp(targNm(1+(tInfo(tmpId,6)==180))','Right')];
                                nr4=[nr4;[rectWin*r{i,j,l}(:,1:(dur1+periReadyWindow))']'];% [nr4;[1000*rectWin*r{i,j,l}(:,1:(dur1+periReadyWindow))']'];  % [trial x winId] % 15) 100 sliding window   prior x eyeHand x target 222GLM [ready-250 ready+480]'    
            end % j for different ts  
            % noSet trials
            if noSet
                dur2=max(T{i})+waitDur; tDur=(-dur1+1):dur2;
                tmpId=tInfo(:,2)==4 &...  % event id 4.ready
                    tInfo(:,4)==(2-i) &... % idShortTrial
                    tInfo(:,5)==(2-l) &... % idHandEye
                    idSuccessNoSet==1 &... % if animals wait >400+max(ts|prior)
                    idNoSet==1; % NoSet trials has max ts of each prior
                tE=tInfo(tmpId,3)*1000; % into ms
                t1=tE-dur1;            t2=tE+dur2;
                
                rNoSet=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
                tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
                for k=1:nnz(tmpId)
                    tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
                    if ~isempty(tSp),rNoSet(k,ceil(tSp))=1;end;
                    tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
                    if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
                end % k for trials
                rConv=conv2(tmpR,gaussampFilter,'valid')*1000;
                lw=4;
                 if idPlot, h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap{i,l}((i-1)*(1-nTspp)+nTspp,:)./2,'linewidth',lw},1); end % darker than overlap
            end % noSet
            
        end % i for ShortLong    
        if idPlot, 
        title(ehNm{l});
        axis tight; xlabel('time from ready (ms)'); %ylabel('spike rate (/s)');plotVertical(gca);
        hMatTmp=hMat(:,:,l); legNmTmp=legName(:,:,l);
        legend(hMatTmp(:),legNmTmp(:),'location','best'); legend boxoff;drawnow;
        end
        iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;
    end % l for EyeHand
    
    if idReg
    %     4)    ' prior x eyeHand x target 222GLM [ready-250 ready
    nrT=table(nr(:,1),iPr,iEH,iRL,'VariableNames',{'spikeRate_ready_m250' 'idShort' 'idEye' 'idRight'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{4}=fitglm(nrT,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar','spikeRate_ready_m250'); % 'Exclude' for outlier
    % stat{4}=anovan(nr,{iPr iEH iRL},'model',3,'display','off','varnames',{'idShort' 'idEye' 'idRight'});
    
    % 5) ' prior x eyeHand x target 222GLM   ' prior x eyeHand 22ANOVA [ready ready+480]
    nrT2=table(nr2(:,1),iPr,iEH,iRL,'VariableNames',{['spikeRate_ready_' num2str(periReadyWindow)] 'idShort' 'idEye' 'idRight'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{5}=fitglm(nrT2,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar',['spikeRate_ready_' num2str(periReadyWindow)]); % 'Exclude' for outlier
    %     stat{5}=anovan(nr2,{iPr iEH},'model',2,'display','off','varnames',{'idShort' 'idEye'});
    
    % 7) prior x eyeHand x target 222GLM [ready+480 ready+800] for overlap ts inc. all long prior
    nrT3=table(nr3(:,1),iPr3,iEH3,iRL3,'VariableNames',{'spikeRate_shortPrior' 'idShort' 'idEye' 'idRight'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{7}=fitglm(nrT3,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar','spikeRate_shortPrior'); % 'Exclude' for outlier
    %     stat{7}=anovan(nr3,{iPr2 iEH2},'model',2,'display','off','varnames',{'idShort' 'idEye'});
    % %     stat{7}=ranksum(nr2(strcmp(prNm(1),iShortOverlap)),nr2(strcmp(prNm(2),iShortOverlap)));
    
    % 15) 100 sliding window   prior x eyeHand x target 222GLM [ready-250 ready+480]'
%     stat{15}=[];
%     X=x2fx([iPr4 iEH4 iRL4],[1 0 0;0 1 0 ;0 0 1;... % 1st level
%         1 1 0;1 0 1; 0 1 1]); % interaction with ts
%     for iWin=1:size(nr4,2) % (dur1+dur2)
%         %     nrT4=table(nr4(:,iWin),iPr,iEH,iRL,'VariableNames',{'spikeRate_250_ready_480' 'idShort' 'idEye' 'idRight'});  % logical(nr(:,2)),logical(nr(:,3))
%         %       glmTmp=fitglm(nrT4,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar','spikeRate_250_ready_480'); % 'Exclude' for outlier
%         %
%         %     stat{15}=cat(3,stat{15}, [table2array(glmTmp.Coefficients(:,1))';... % 1 x nCoeff
%         %    table2array(glmTmp.Coefficients(:,end))']);
%         
%         [B,DEV,glmStat]=glmfit(X,nr4(:,iWin),'poisson'); % glmStat.p
%         stat{15}= cat(3,stat{15},[B(:)';... % 1 x nCoeff
%             glmStat.p(:)']);
%         
%     end
    
    % FINAL:  log(nSpikesWin100_ready[m250_480]) ~ 1 + idLong + idHand + idLeft + idLong*idHand + idLong*idLeft + idHand*idLeft (7)
    stat{15}=[];
    X=x2fx([iPr4 iEH4 iRL4],[1 0 0;0 1 0 ;0 0 1;... % 1st level
        1 1 0;1 0 1; 0 1 1],... % interaction with ts
        [1 2 3]);% categorical
    for iWin=1:size(nr4,2) % (dur1+periReadyWindow)
        [B,DEV,glmStat]=glmfit(X,nr4(:,iWin),'poisson'); % glmStat.p
        stat{15}= cat(3,stat{15},[B(:)';... % 1 x nCoeff
            glmStat.p(:)']);
        
    end
    end % idReg
end

%% 5. Set % photodiode sync [-480/800/1200 min(tp)]
% 10 ts(short/long) for each eye, hand (autumn/winter, spring/summer; use only 200 to avoid yellow)
if eventId==5 || eventId==0 % 0 for plot all
%     dur2=floor(min(t));
periSetWindow=200; % 240;%T{1}(1); % 400;
periSetWindow2=300; % 100; % 240; % 320
%     figure(201); set(gcf,'position',pplot.rect1_4,'color','w'); hold all; % eye
%     figure(202); set(gcf,'position',pplot.rect2_4,'color','w'); hold all; % hand
            nr=[]; nr2=[];nr3=[];iPr=[]; iEH=[];iPr2=[]; iEH2=[]; iRL=[];ts=[];tp=[]; iShortOverlap=[];% for stat
            iPr3=[]; iEH3=[]; iRL3=[];
            nr4=[];rectWin=rect_filter(winSize,T{1}(1)+periSetWindow2,false); % [time x time] row NOT normalized

            % find valid min tp
            tmp=unique(tInfo(:,[1 7 8]),'rows'); tpMinTs=tmp(tmp(:,2)==T{1}(1)&tmp(:,3)>0,3); % assuming mintp from shortest Ts
            minTpMinTs=mean(tpMinTs)-2*std(tpMinTs);
            dur2min=floor(minTpMinTs);
            
    tmp1=spring(round(nTspp*1.5));
    tmp2=summer(round(nTspp*1.5));
    tmpCmap={flipud(autumn(nTspp)), flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
        winter(nTspp), tmp2(1:nTspp,:)}; % longEye, longHand
    legName=cell(nPr,nTspp,nEH);
    prNm={'Short','Long'};
    ehNm={'Eye','Hand'};hMat=[];
%     hMat=cell(nPr,nTspp,nEH);
    
    r=cell(nPr,nTspp,nEH);
    for l=1:nEH        
        if idPlot, subplot(mPlot,nPlot,3+(l-1)*nPlot);hold all;end
%         figure(200+l);
        for i=1:nPr,
            for j=1:nTspp,
                dur1=T{i}(j); % ms
                dur2=max(periSetWindow2,floor(min(tInfo(tInfo(:,5)==(2-l) &tInfo(:,4)==(2-i) &tInfo(:,7)==T{i}(j)&tInfo(:,8)>0,8)))); %%%%%%
%                 dur2=max(dur2min,floor(min(tInfo(tInfo(:,5)==(2-l) &tInfo(:,4)==(2-i) &tInfo(:,7)==T{i}(j)&tInfo(:,8)>0,8)))); %%%%%%
                tDur=(-dur1+1):dur2;
                legName{i,j,l}=[prNm{i} num2str(T{i}(j))];
                    
                tmpId=tInfo(:,2)==5 &...  % event id 5.set
                tInfo(:,4)==(2-i) &... % idShortTrial
                tInfo(:,5)==(2-l) &... % idHandEye
                tInfo(:,7)==T{i}(j)&... % ; %  &... % Ts
                tInfo(:,8)>periSetWindow2; % dur2min; % only valid trials for now; 0 %%%%%% % assuming removeOuliter works fine
                tE=tInfo(tmpId,3)*1000; % into ms
                t1=tE-dur1;            t2=tE+dur2;
                
                r{i,j,l}=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
                tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
                for k=1:nnz(tmpId)
                    tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
                    if ~isempty(tSp),r{i,j,l}(k,ceil(tSp))=1;end;
                    tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
                    if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
                end % k for trials
                rConv=conv2(tmpR,gaussampFilter,'valid')*1000;                
                if (i==1 & j==nTspp) || (i==nPr & j==1), lw=4; % overlapped interval
                    if idPlot, 
                    h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap{i,l}(j,:),'linewidth',lw},1); hMat(i,j,l)=h.mainLine;end
                    nr2=[nr2;sum(r{i,j,l}(:,(dur1-periSetWindow2):(dur1)),2)];iShortOverlap=[iShortOverlap;repmat(prNm(i),[nnz(tmpId) 1])];
                    iPr2=[iPr2; repmat(prNm(i),[nnz(tmpId) 1])];iEH2=[iEH2; repmat(ehNm(l),[nnz(tmpId) 1])];
                else lw=1.5; 
                    if idPlot, 
                    hMat(i,j,l)=plot(tDur(:)',mean(rConv,1),'color',tmpCmap{i,l}(j,:),'linewidth',lw); end
                end;
                try
                nr=[nr;sum(r{i,j,l}(:,(dur1-periSetWindow+1):dur1),2)]; 
                catch
                    disp('');
                end
                nr3=[nr3;sum(r{i,j,l}(:,(dur1+1):(dur1+periSetWindow2)),2)]; %  dur2min
                 nr4=[nr4;[rectWin*r{i,j,l}(:,(dur1-T{1}(1)+1):(dur1+periSetWindow2))']'];  %[nr4;[1000*rectWin*r{i,j,l}(:,(dur1-T{1}(1)+1):(dur1+periSetWindow2))']'];  % [trial x winId] % 15) 100 sliding window   prior x eyeHand x target 222GLM [ready-250 ready+480]'    
                iPr=[iPr; repmat(prNm(i),[nnz(tmpId) 1])];iEH=[iEH; repmat(ehNm(l),[nnz(tmpId) 1])];iRL=[iRL; targNm(1+(tInfo(tmpId,6)==180))'];  tp=[tp;tInfo(tmpId,8)];
                iPr3=[iPr3; repmat(strcmp(prNm{i},'Short'),[nnz(tmpId) 1])];iEH3=[iEH3; repmat(strcmp(ehNm{l},'Eye'),[nnz(tmpId) 1])];iRL3=[iRL3; strcmp(targNm(1+(tInfo(tmpId,6)==180))','Right')];
                ts=[ts; repmat(T{i}(j),[nnz(tmpId) 1])];
            end % j for different ts
        end % i for ShortLong        
        if idPlot, 
        title(ehNm{l});
        axis tight; xlabel('time from set (ms)'); %ylabel('spike rate (/s)');plotVertical(gca);
        hMatTmp=hMat(:,:,l); legNmTmp=legName(:,:,l);
        legend(hMatTmp(:),legNmTmp(:),'location','best'); legend boxoff;drawnow;
        end
    end % l for EyeHand
    
    if idReg
    % formatting regressors appropriate for regression (centering for continuous/0&1 for discrete dummy)
    muPr=strcmp(iPr,'Short')*mean(T{1})+strcmp(iPr,'Long')*mean(T{2});
    ts2=(ts-muPr).^2; 
%     ts2=ts2-mean(ts2);
%     ts=ts-mean(ts);
ctp=tp-mean(tp);
    
    %     '6)   ' ts + ts^2 + idPrior regression ' [set-200 set]'
    nrT=table(nr(:,1),ts-mean(ts),ts2-mean(ts2),iPr,'VariableNames',{['spikeRate_' num2str(periSetWindow) '_set'] 't_s' 't_s2' 'idShort'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{6}=fitglm(nrT,['spikeRate_' num2str(periSetWindow) '_set ~ idShort*(t_s+t_s2)'],'Distribution','poisson','CategoricalVars',{'idShort'},'ResponseVar',['spikeRate_' num2str(periSetWindow) '_set']); % 'Exclude' for outlier
    %     stat{6}=anovan(nr,{iPr iEH},'model',2,'display','off','varnames',{'idShort' 'idEye'});
    
    % 8) ' prior x eyeHand x target 222GLM    [set-200 set]
    nrT2=table(nr(:,1),iPr,iEH,iRL,'VariableNames',{['spikeRate_' num2str(periSetWindow) '_set'] 'idShort' 'idEye' 'idRight'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{8}=fitglm(nrT2,'interactions','Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar',['spikeRate_' num2str(periSetWindow) '_set']); % 'Exclude' for outlier
    
    % 9) intercept tp + eyeHand x directions regression [set set+min(tp)]
    nrT3=table(nr3(:,1),tp,iEH,iRL,'VariableNames',{'spikeRate_set_mintp' 't_p' 'idEye' 'idRight'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{9}=fitglm(nrT3,'interactions','Distribution','poisson','CategoricalVars',{'idEye' 'idRight'},'ResponseVar','spikeRate_set_mintp'); % 'Exclude' for outlier
    %     tbl=table(nr3,tp,iEH,iRL,'VariableNames',{'spikeRate','t_p','idEye','idRight'});
    %     lm=fitlm(tbl,'spikeRate~t_p*idEye*idRight','CategoricalVar',{'idEye','idRight'});
    %     stat{8}=table2array(lm.Coefficients(:,end)); % [8 x 1]
    
    % % 16 ) 100 sliding window  [set-480 set+100]
    % %         % log(nSpikes_100_set) ~ b1 + b2*t_s +  b3*idShort + b4*idEye + b5*idRight ...
    % %         %                                         + b6*t_s*idShort + b7*t_s*idEye + b8*t_s*idRight ...
    % %         %                                           + b9*t_s^2 + b10*t_s^2* *idShort + b11*t_s^2*idEye + b12*t_s^2*idRightfor iWin=1:size(nr4,2) % (dur1+dur2)
    % stat{16}=[];
    % for iWin=1:size(nr4,2) % (dur1+dur2)
    %     nrT4=table(nr4(:,iWin),ts-mean(ts),ts2-mean(ts2),iPr,iEH,iRL,'VariableNames',{'spikeRate_480_set_100' 'ts' 't_s2' 'idShort' 'idEye' 'idRight'});  % logical(nr(:,2)),logical(nr(:,3))
    %     glmTmp=fitglm(nrT4,['spikeRate_480_set_100 ~ ts+t_s2+idShort+idEye+idRight+ idShort:ts+idEye:ts+idRight:ts+idShort:t_s2+idEye:t_s2+idRight:t_s2'],...
    %         'Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar','spikeRate_480_set_100'); % 'Exclude' for outlier
    %
    %
    %      stat{16}= cat(3,stat{16},[table2array(glmTmp.Coefficients(:,1))';... % 1 x nCoeff
    %    table2array(glmTmp.Coefficients(:,end))']);
    %
    %
    % end
    %
    % % nrT4=table(nr(:,1),ts-mean(ts),ts2-mean(ts2),iPr,iEH,iRL,'VariableNames',{['spikeRate_' num2str(periSetWindow) '_set'] 'ts' 't_s2' 'idShort' 'idEye' 'idRight'});  % logical(nr(:,2)),logical(nr(:,3))
    % % glmTmp=fitglm(nrT4,['spikeRate_' num2str(periSetWindow) '_set ~ ts+t_s2+idShort+idEye+idRight+ idShort:ts+idEye:ts+idRight:ts+idShort:t_s2+idEye:t_s2+idRight:t_s2'],...
    % %         'Distribution','poisson','CategoricalVars',{'idShort' 'idEye' 'idRight'},'ResponseVar',['spikeRate_' num2str(periSetWindow) '_set']); % 'Exclude' for outlier
    % %
    % %      stat{16}= [table2array(glmTmp.Coefficients(:,1))';... % 1 x nCoeff
    % %    table2array(glmTmp.Coefficients(:,end))'];
    
%     % use glmfit (hopefully save analysis time)
%     stat{17}=[];stat{16}=[];
%     X=x2fx([ts-mean(ts) ts2-mean(ts2) iPr3],[1 0 0;0 1 0;0 0 1;1 0 1;0 1 1]); % idShort*(t_s+t_s2)'
%     X2=x2fx([ts-mean(ts) ts2-mean(ts2) iPr3 iEH3 iRL3],[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;... % 1st level
%         1 0 1 0 0; 1 0 0 1 0;1 0 0 0 1;... % interaction with ts
%         0 1 1 0 0; 0 1 0 1 0;0 1 0 0 1]); % interaction with ts2
%     for iWin=1:size(nr4,2) % (dur1+dur2)
%         %          nrT5=table(nr4(:,iWin),ts-mean(ts),ts2-mean(ts2),iPr,'VariableNames',{['spikeRate_' num2str(periSetWindow) '_set'] 't_s' 't_s2' 'idShort'}); % logical(nr(:,2)),logical(nr(:,3))
%         %     stat{17}=fitglm(nrT5,['spikeRate_' num2str(periSetWindow) '_set ~ idShort*(t_s+t_s2)'],'Distribution','poisson','CategoricalVars',{'idShort'},'ResponseVar',['spikeRate_' num2str(periSetWindow) '_set']); % 'Exclude' for outlier
%         
%         
%         %     '6)   ' ts + ts^2 + idPrior regression ' [set-200 set]'
%         [B,DEV,glmStat]=glmfit(X,nr4(:,iWin),'poisson'); % glmStat.p
%         stat{17}= cat(3,stat{17},[B(:)';... % 1 x nCoeff
%             glmStat.p(:)']);
%         
%         [B,DEV,glmStat]=glmfit(X2,nr4(:,iWin),'poisson'); % glmStat.p
%         stat{16}= cat(3,stat{16},[B(:)';... % 1 x nCoeff
%             glmStat.p(:)']);
%     end % iWin
    
    % FINAL
    % log(nSpikesWin100_set[m480_300]) ~ 1 + catTs (8)+ idLong + tp + idHand + idLeft + catTs*tp (8) (24) 

    % use glmfit (hopefully save analysis time)
    stat{16}=[];
    X=x2fx([ts iPr3 ctp iEH3 iRL3],[1 0 0 0 0 ;0 1 0 0 0 ;0 0 1 0 0 ; 0 0 0 1 0 ; 0 0 0 0 1 ;... % 1st level+catTs(800)*idLong
        1 0 1 0 0],...% ;0 0 1 1 0 ;... % 2nd level interaction: ccatTs*idLong (8) + catTs*tp (8)+ tp*idHand; %         0 0 1 1 1 0;],... % tp*idHand*idLeft
        [1 2 4 5]); % categorical
% log(nSpikesWin100_set[m480_300]) ~ 1 + catTs (8)+ idLong + tp + idHand + idLeft + catTs(800)*idLong (1) + catTs*tp (8)+ tp*idHand + tp*idHand*idLeft (24) 
    % ill-conditioned X: (ts==800).*(~iPr3) can be constructed from constatnt+ts+iPr3?
    % also high-level interaction removed
%     X=x2fx([ts iPr3 ctp iEH3 iRL3 (ts==800).*(~iPr3)],[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0;0 0 0 0 0 1;... % 1st level+catTs(800)*idLong
%         1 0 1 0 0 0;0 0 1 1 0 0;... % 2nd level interaction: ccatTs*idLong (8) + catTs*tp (8)+ tp*idHand
%         0 0 1 1 1 0;],... % tp*idHand*idLeft
%         [1 2 4 5 6]); % categorical
    for iWin=1:size(nr4,2) % (dur1+dur2)
        [B,DEV,glmStat]=glmfit(X,nr4(:,iWin),'poisson'); % glmStat.p
        stat{16}= cat(3,stat{16},[B(:)';... % 1 x nCoeff
            glmStat.p(:)']);
        % TBD: f stat for all catTs related variables
        
    end % iWin
    end % idReg
end

%% 6. Production [-tp 50]
% 5 split of tp (jet)
% nSplitTp=5; t=tInfo(tInfo(:,8)>0); tSort=sort(t); tBound=tSort(round([(1/nSplitTp):(1/nSplitTp):1]*length(t))); % 0 0.2 0.4 ... 1
if eventId==6 || eventId==0 % 0 for plot all
    dur2=50;
    periGoWindow=300;
    nr=[]; iEH=[]; iRL=[];tp=[];% for stat
    % find valid min tp
    tmp=unique(tInfo(:,[1 7 8]),'rows'); tpMinTs=tmp(tmp(:,2)==T{1}(1)&tmp(:,3)>0,3); % assuming mintp from shortest Ts
    minTpMinTs=mean(tpMinTs)-2*std(tpMinTs);
    dur2min=floor(minTpMinTs);
            
    iEH3=[]; iRL3=[];
    nr4=[];rectWin=rect_filter(winSize,periGoWindow+dur2,false); % [time x time] row NOT normalized

%     figure(301); set(gcf,'position',pplot.rect1_5,'color','w'); hold all; % eye
%     figure(302); set(gcf,'position',pplot.rect2_5,'color','w'); hold all; % hand
        tmpCmap=(bone(round(1.5*(nTspp+1)))); tmpCmap=tmpCmap(1:(nTspp+1),:); % last for incorrect
    legName=cell(nSplitTp,nEH,nTarg);
    ehNm={'Eye','Hand'};linestyle={'-','--'};hMat=[];
%     hMat=cell(nSplitTp,nEH,nTarg);
    
    r=cell(nSplitTp,nEH,nTarg);
    for l=1:nEH        
        if idPlot, subplot(mPlot,nPlot,4+(l-1)*nPlot);hold all;end
%         figure(300+l);
        for i=1:nSplitTp,
            for j=1:nTarg,
                dur1=max(periGoWindow,round(tBound(i)));
%                 dur1=max(dur2min,round(tBound(i))); % ms %%%%%%
                tDur=(-dur1+1):dur2;
                legName{i,l,j}=['t_p:[' num2str(round(tBound(i))/1000) ' ' num2str(round(tBound(i+1))/1000) ']'];
                    
                tmpId=tInfo(:,2)==6 &...  % event id 6.production
                tInfo(:,5)==(2-l) &... % idHandEye
                tInfo(:,6)==(j-1)*180 &... % target location
                tInfo(:,8)<=tBound(i+1) &... % tp range                
                tInfo(:,8)>dur1; % tBound(i); % only valid trials for now %%%%%
                tE=tInfo(tmpId,3)*1000; % into ms
                t1=tE-dur1;            t2=tE+dur2;
                
                r{i,l,j}=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
                tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
                for k=1:nnz(tmpId)
                    tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
                    if ~isempty(tSp),r{i,l,j}(k,ceil(tSp))=1;end;
                    tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
                    if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
                end % k for trials
                rConv=conv2(tmpR,gaussampFilter,'valid')*1000;   
                if idPlot, 
                if i==1 || i==nSplitTp
                    h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{linestyle{j},'color',tmpCmap(i,:),'linewidth',2},1); hMat(i,l,j)=h.mainLine;
                else
                    hMat(i,l,j)=plot(tDur(:)',mean(rConv,1),linestyle{j},'color',tmpCmap(i,:),'linewidth',2);
                end
                end
                nr=[nr;sum(r{i,l,j}(:,(dur1-periGoWindow+1):dur1),2)]; %  dur2min
                iEH=[iEH; repmat(ehNm(l),[nnz(tmpId) 1])];iRL=[iRL; targNm(1+(tInfo(tmpId,6)==180))'];  tp=[tp;tInfo(tmpId,8)];
                
                nr4=[nr4;[rectWin*r{i,l,j}(:,(dur1-periGoWindow+1):(dur1+dur2))']']; 
                iEH3=[iEH3; repmat(strcmp(ehNm{l},'Eye'),[nnz(tmpId) 1])];iRL3=[iRL3; strcmp(targNm(1+(tInfo(tmpId,6)==180))','Right')];
            end %j=1:nTarg,
        end % i=1:nSplitTp, 
        if idPlot, 
        title(ehNm{l});
        axis tight; xlabel('time from production (ms)'); %ylabel('spike rate (/s)');plotVertical(gca);
        hMatTmp=hMat(:,l,:); legNmTmp=legName(:,l,:);
        legend(hMatTmp(:),legNmTmp(:),'location','best'); legend boxoff;drawnow;
        end
    iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;
    end % l for EyeHand
    
    if idReg
%     10) intercept tp + eyeHand x directions regression [production-min(tp) production]
    nrT3=table(nr(:,1),tp,iEH,iRL,'VariableNames',{'spikeRate_go_mintp' 't_p' 'idEye' 'idRight'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{10}=fitglm(nrT3,'interactions','Distribution','poisson','CategoricalVars',{'idEye' 'idRight'},'ResponseVar','spikeRate_go_mintp'); % 'Exclude' for outlier
%     tbl=table(nr,tp,iEH,iRL,'VariableNames',{'spikeRate','t_p','idEye','idRight'});
%     lm=fitlm(tbl,'spikeRate~t_p*idEye*idRight','CategoricalVar',{'idEye','idRight'});    
%     stat{9}=table2array(lm.Coefficients(:,end)); % [8 x 1]

% FINAL
% log(nSpikesWin100_go[m300_50]) ~ 1 + t_p + idHand + idLeft + t_p*idHand + t_p*idHand*idLeft (6)
ctp=tp-mean(tp); % centering t_p

    stat{17}=[];
    X=x2fx([ctp iEH3 iRL3],[1 0 0;0 1 0;0 0 1;...
        1 1 0;1 1 1],... % t_p*idHand + t_p*idHand*idLeft (6)
        [2 3]); 

    for iWin=1:size(nr4,2) % (dur1+dur2)
        [B,DEV,glmStat]=glmfit(X,nr4(:,iWin),'poisson'); % glmStat.p
        stat{17}= cat(3,stat{17},[B(:)';... % 1 x nCoeff
            glmStat.p(:)']);
    end % iWin
    end % idReg
end

%% 8. reward delay % photodiode sync [-tp 100] rewdur:90
% 9. reward: bonusRewDur*1000+rewardDur*1000*max(0.00001,1-abs(productionInterval - interval)/interval/win_fraction ) % photodiode sync
% 3 split for reward (bone): reward/incorrect (solid/dashed)
% nSplitRew=5; rewValid=tInfo(tInfo(:,8)>0&tInfo(:,12)>0,end); rSort=sort(rewValid); rBound=rSort(ceil([0:(1/nSplitRew):1]*length(rewValid))); % 0 0.33 0.66 1
if eventId==8 || eventId==9 || eventId==0 % 0 for plot all
    td=15; nTd=4; xFiltG=-td*nTd:1:td*nTd; gaussampFilter2=exp(-(xFiltG.^2)/(2*(td^2))); gaussampFilter2=gaussampFilter2/sum(gaussampFilter2); % 40

    dur2=100;dur1=100;
    nr4=[];rectWin=rect_filter(winSize,dur1+dur2,false); % [time x time] row NOT normalized
     nr=[]; iRew=[];rewDur=[]; % for stat

    if idPlot, subplot(mPlot,nPlot,nPlot);hold all;end
%     figure(16); set(gcf,'position',pplot.rect1_6,'color','w'); hold all; 
    tmpCmap=(copper(round(1.2*(nSplitRew+1)))); tmpCmap=tmpCmap(1:(nSplitRew+1),:); % last for incorrect
    legName=cell(nSplitRew+1,1);hMat=[];
%     hMat=cell(nSplitRew+1,1);    
    r=cell(nSplitRew+1,1);
    for i=1:nSplitRew,                
%         dur1=floor(min(tSort)); % ms
        tDur=(-dur1+1):dur2;        
        legName{i}=['reward:[' num2str(round(rBound(i))) ' ' num2str(round(rBound(i+1))) ']'];        
        tmpId=tInfo(:,2)==9 &...  % event id 8 for incorret, 9 for reward
            tInfo(:,12)<=rBound(i+1) &... % reward range
            tInfo(:,12)>rBound(i); % only valid trials for now
        tE=tInfo(tmpId,3)*1000; % into ms
        t1=tE-dur1;            t2=tE+dur2;        
        r{i}=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
        tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
        for k=1:nnz(tmpId)
            tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
            if ~isempty(tSp),r{i}(k,ceil(tSp))=1;end;
            tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
            if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
        end % k for trials
        nr=[nr;sum(r{i}(:,(dur1+1):end),2)]; rewDur=[rewDur; tInfo(tmpId,end)];iRew=[iRew; ones(nnz(tmpId),1)];

        nr4=[nr4;[rectWin*r{i}']']; 
        
        rConv=conv2(tmpR,gaussampFilter2,'valid')*1000;
        if idPlot, 
        h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap(i,:),'linewidth',2},1); hMat(i)=h.mainLine;
        end
    end % i for ShortLong
    % for no reward (incorrect)
    i=1+nSplitRew;
    legName{i}='no reward';
    tmpId=tInfo(:,2)==8 &...  % event id 8 for incorret, 9 for reward
        tInfo(:,12)<=0.001 &... % reward range
        tInfo(:,8)>0; % only valid trials for now
    tE=tInfo(tmpId,3)*1000; % into ms
    t1=tE-dur1;            t2=tE+dur2;
    r{i}=zeros(nnz(tmpId),dur1+dur2); % [trial x time]
    tmpR=zeros(nnz(tmpId),dur1+dur2+td*nTd*2); % for convoluted plot
    for k=1:nnz(tmpId)
        tSp=sp(t1(k)<sp & sp<=t2(k))-t1(k);
        if ~isempty(tSp),r{i}(k,ceil(tSp))=1;end;
        tSp2=sp(t1(k)-td*nTd<sp & sp<=t2(k)+td*nTd)-(t1(k)-td*nTd);
        if ~isempty(tSp2),tmpR(k,ceil(tSp2))=1;end;
    end % k for trials
            nr=[nr;sum(r{i}(:,(dur1+1):end),2)]; rewDur=[rewDur; zeros(nnz(tmpId),1)];iRew=[iRew; zeros(nnz(tmpId),1)];
            nr4=[nr4;[rectWin*r{i}']'];

    rConv=conv2(tmpR,gaussampFilter2,'valid')*1000;
    if idPlot, 
    h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap(i,:),'linewidth',2},1); hMat(i)=h.mainLine;
    axis tight; xlabel('time from reward onset (ms)'); %ylabel('spike rate (/s)');plotVertical(gca);
    legend(hMat(:),legName(:),'location','best'); legend boxoff;drawnow;
    iYlim=iYlim+1; yLimMat(iYlim,:)=ylim;
    end
    
    if idReg
    % 10) reward regression [rewardOnset onset+100]
    nrT3=table(nr(:,1),rewDur,iRew,'VariableNames',{'spikeRate_reward_100' 'rewDur' 'rewarded'}); % logical(nr(:,2)),logical(nr(:,3))
    stat{11}=fitglm(nrT3,'spikeRate_reward_100~rewDur*rewarded-rewDur','Distribution','poisson','CategoricalVars',{'rewarded'}); % 'Exclude' for outlier
%     tbl=table(nr,rewDur,iRew,'VariableNames',{'spikeRate','rewDur','rewarded'});
%     lm=fitlm(tbl,'spikeRate~rewDur*rewarded-rewDur','CategoricalVar',{'rewarded'});    
%     stat{10}=table2array(lm.Coefficients(:,end)); % [3 x 1]

% FINAL
% log(nSpikesWin100_reward[m100_100]) ~ 1 + rewarded + rewDur:rewarded (3)
cRewDur=rewDur-mean(rewDur);% centering rewDur

stat{18}=[];
X=x2fx([iRew cRewDur],[1 0 ;... % 1st level %17/7/27 main effect of cRewDur removed
    1 1],... % 2nd level interaction: rewDur:rewarded
    [1]); % categorical
for iWin=1:size(nr4,2) % (dur1+dur2)
    [B,DEV,glmStat]=glmfit(X,nr4(:,iWin),'poisson'); % glmStat.p
    stat{18}= cat(3,stat{18},[B(:)';... % 1 x nCoeff
        glmStat.p(:)']);
end % iWin

    end % idReg
end

%%
%matching ylim
if idPlot, 
for i=1:(mPlot*nPlot)
    h=subplot(mPlot,nPlot,i);
    set(h,'position',[0.2*rem(i-1,5)+0.02*(i==1||i==nPlot+1) 0.5*(i<6)+0.03 0.18-0.02*(i==1||i==nPlot+1) 0.42]);
    ylim([min(yLimMat(:,1)) max(yLimMat(:,2))]);
    plotVertical(gca);
    if i~=1 && i~=(nPlot+1), set(gca,'yticklabel',[]); end;
end
end
% disp('');
% exportfig(gcf,[neuDir fid '_' num2str(cname{i}(j,1)) '_' num2str(cname{i}(j,2)) 'new.png'],'color','rgb','Height',13,'Width',11,'format','png','Reference',hTmp);

function [idState,idSw,iSw,tNB,tmp]=sortSwitch(tInfo)
tmp=unique(tInfo(:,[1 4 5 8]),'rows','stable'); % t#, idLS, idHE
idState=(1-tmp(:,2))*2+(1-tmp(:,3)); % 0123 for shortEye, shortHand, LongEye, LongHand
idSw=0; % 123 acrPr, acrMod, acrBoth for shortEye, shortHand, LongEye, LongHand; all 12 types
    iSw=0; % for indexling all 12 types
    tNB=1; % trial number in block
    for i=1:size(tmp,1)-1
        if tmp(i,2)==tmp(i+1,2) && tmp(i,3)==tmp(i+1,3)
            tNB=[tNB; tNB(end)+1];
            idSw=[idSw; 0];
            iSw=[iSw; 0];
        else % switched
            tNB=[tNB; 1];
            switch idState(i)
                case 0 % shortEye
                    switch idState(i+1)
                        case 1 % acrMod
                            idSw=[idSw; 2];
                            iSw=[iSw; 1];
                        case 2 % acrPr
                            idSw=[idSw; 1];
                            iSw=[iSw; 2];
                        case 3 % acrBoth
                            idSw=[idSw; 3];
                            iSw=[iSw; 3];
                    end
                case 1 % shortHand
                    switch idState(i+1)
                        case 0 % acrMod
                            idSw=[idSw; 2];
                            iSw=[iSw;4];
                        case 2 % acrBoth
                            idSw=[idSw; 3];
                            iSw=[iSw; 5];
                        case 3 % acrPr
                            idSw=[idSw; 1];
                            iSw=[iSw; 6];
                    end
                case 2 % LongEye
                    switch idState(i+1)
                        case 0 % acrPr
                            idSw=[idSw; 1];
                            iSw=[iSw; 7];
                        case 1 % acrBoth
                            idSw=[idSw; 3];
                            iSw=[iSw; 8];
                        case 3 % acrMod
                            idSw=[idSw; 2];
                            iSw=[iSw; 9];
                    end
                case 3 % LongHand
                    switch idState(i+1)
                        case 0 % acrBoth
                            idSw=[idSw; 3];
                            iSw=[iSw; 10];
                        case 1 % acrPr
                            idSw=[idSw; 1];
                            iSw=[iSw; 11];
                        case 2 % acrMod
                            idSw=[idSw; 2];
                            iSw=[iSw; 12];
                    end
            end
        end
    end
    
%% back up
% for i=1:length(fname)
%     disp(['===== ' fname{i} ' =====']);
%     cd([dirName fname{i}]);
% %     chkSDF{i}=[];
%     % MAT file: RSG prior
%     fid=fname{i}(end-5:end);
% beh=load([behDir 'H_20' fid '.mat']);
% figure; hist(beh.t(beh.T==480&beh.t>0),1000);
% set(gca,'xtick',[-2*std(beh.t(beh.T==480&beh.t>0)) 0]+mean(beh.t(beh.T==480&beh.t>0)));
% xlim([0 max(xlim)]);
% waitforbuttonpress; close all;
% end