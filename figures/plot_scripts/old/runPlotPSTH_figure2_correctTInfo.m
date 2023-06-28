function runPlotPSTH_figure2_correctTInfo

% 2018/3/26
% plot PSTH for problematic sessions of 8/23 (8/22, 12/7)

% plot PSTH using kilosort results (no manual start/end trials)
% main dependent function: plotSDFexport_inFunc
% almost identical runPlotSDF
% ccc;

% note
% difference b/t runPlotSDFnew vs runPlotSDF: % new stat: poisson regression

% save data from around 800 trials in 12/7
% data transfer speed warning issue in 12/13? two cells in NS3 data

%% init
S=initRSG2prior;
% fn={fnameH fnameG}; H_RSGprior_20161203

v2struct(S);
cd(neuDir2); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo

%% specifying example cells
cname=... % sessId idClust    
[... % new selection
% H
% 161207 2055;... % new addition
% 161207 2058;... % new addition
% 161207 2079;... % new addition
% 161207 2096;... % new addition
161210 3011;...	%%%%%
% 161210 3097;...	%%%%%
% 161221 1117;...	
% 161221 1128;...	%%%%%
% 161221 1163;...	
161214 2159;...	%%%%%
% 161206 1011;...	%%%%%
% 161206 2016;...	
% 161206 2120;...	
% 161206 2201;...	
% 161206 2227;...	
% 161222 1126;...	
161218 2022;...	
% 161211 1029;...	
% 161211 1007;...	
% 161211 1003;...	
% 161211 1028;...	
% 161211 1030;...	
% 161211 1044;...	
% 161211 1049;...	
% 161211 1053;...	
% 161211 1115;...
% G
% 170822 2008;... % new addition
% 170822 2013;... $ new addition
% 170822 2123;... % new additioln
% 170822 2140;... % new addition
% 170822 3053;... % new addition
% 170822 3106;... % new addition
% 170822 3137;... % new addition
% 170823 1099;... % new addition
% 170823 1101;... % new addition
% 170823 1104;... % new addition
% 170823 1148;... % new addition
% 170818 3088;... 
% 170818 2153;... %%%%%
% 170823 1104;...
170506 1111;... %%%%%
% 170507 1127;...	%%%%%
% 170506 1002;...	
% 170506 1042;...	
% 170506 3035;...	
% 170506 2100;...	
% 170506 2098;...	
170506 2096;...	
% 170506 2019;...	
% 170506 2015;...	
% 170506 1095;...	
% 170507 3144;...	
% 170507 3136;...	%%%%%
% 170507 1117;...	
% 170507 1097;...	
% 170507 1024;...	
% 170508 3103;...	
% 170508 3009;...	%%%%%
% 170508 1111;...	
% 170511 3149;...	
% 170511 2117;...	
% 170817 1192;...	
% 170817 1007;...	
% 170818 3088;...	
% 170818 2097;...	
% 170818 2018;...	
% 170818 1167;...	
% 170818 1039;...	
% 170823 3167;...	
% 170823 2147;...	
170823 2111;...	
% 170823 1136;...
]; % 

% % old selection
%     [...% 161210 3000;... % 11;...%% prior
% %     161210 3104;...%% prior
%     161211 1007;...%% prior
% %     161211 1017;...%%%%%
% %     161211 1029;...%% prior
% % %     161211 1107;...%% prior
% % %     161211 1019;...%% prior
% % %     161211 1044;...%% prior
%     161211 1053;...%% ramp %%%%%%
% % %     161211 1102;...%% prior
% %     161218 2007;...%% prior
%     161218 2022;...%% prior
%     161218 2036;...%%%%%
% %     161222 1126;...%% prior
% %     161203 1034;...
% %     161221 2114;...
% %     161222 1126;...
%     161221 1163;... % vis
%     161221 1117;... % set align
%     161220 1214;...% set align
% %     161204 1003;...
%     161204 1009;... %%%%%%
% %     161204 1042;...
% %     161204 1083;...
% %     161204 1106;...
% %     161204 1141;... %% prior
% %     161205 1136;...
% %     161206 1011;...
% %     161206 1028;...
% %     161206 1097;...
% %     161206 1111;...
% %     161206 2016;...
% %     161206 2087;...
%     161206 2120;... %%%%%%
%     161206 2221;... %%%%%%
%     % G
% %     170506 2015;... % ramping
%     170506 1111;... % prior
%     170507 1127;... % motor %%%%%
% %     170507 1117;... % hand
% %     170507 1097;...
% %     170507 1024;... % ramp/prior
% %     170817 1147;... % ramp down
% %     170817 1007;... % ramp 
% %     170818 3088;... % ramp 
% %     170818 3031;... % ramp 
%     170818 2153;... % prior
% %     170818 2018;... 
% %     170818 1136;... 
% %     170818 1097;... 
% %     170821 2058;... 
% %     170821 1067;... 
% %     170822 2013;...
% %     170822 1179;...
% %     170822 1159;...
% %     170822 2013;...
% %     170823 3167;...
% %     170823 2147;...
% %     170823 2111;...
% %     170823 2107;...
% %     170823 2023;...
% %     170823 1148;...
% %     170823 1136;...
% %     170823 1104;...
% %     170823 1101;...
% %     170823 1099;...
%     ]; % override initRSG2prior

if strcmp(version('-release'),'2017a') % MBP
    dirFig2singleNeuron='/Users/seonminahn/Dropbox (MIT)/figuresRSG2prior/2singleNeuron/';
else
    dirFig2singleNeuron='/Users/hansem/Dropbox (MIT)/figuresRSG2prior/2singleNeuron/';
end

eventId=[4;5;6];
% prior, periReady (priors/modality, meanTsAttrition): 161210_71_1;     161218 57 1;       161210_67_1          161221 8 2(not peak at priorMean)
% visual, periReady, periSet (priors/modality, meanTsAttrition): 161211_81_1            161223 23 3            161222_12_1            161208_57_1
% eye, periGo (hand/eye, direction-dashed): 161206 51 4
% hand, periGo (hand/eye, direction-dashed): 161222 45 4

    %% archive from helpMJgrant (only measurement) 
%     'H_RSGprior_20161211';... [13 1 105 1774];
%     'H_RSGprior_20161221';...[4 2 1654 2139;42 4 1 2139];
%     'H_RSGprior_20161222';... [49 4 1 2158]%
%     'G_RSGprior_20170506';...  [12 2 1 2085;...    38 2  1 2085;...
    % & CSNE17
%     'H_RSGprior_20161206';... [51 4 1 1698;19 1 1 2797];..
%     'H_RSGprior_20161210';...[67 1 1 1420;71 1 1 1546];..
%     'H_RSGprior_20161211';... [7 1 1 1774;81 1 1 1774];... 
%     'H_RSGprior_20161218';...  [55 1 1 2097;57 1 1 2097];...%     [1
%     'H_RSGprior_20161221';...[21 2 138 1411];...%
%     'H_RSGprior_20161222';... [13 2 139 2158]%   
%     'G_RSGprior_20170506';...    [12 2 1 2085;...    38 2  1 2085;...
%     'G_RSGprior_20170508';...    [47 1 1 2388;    65 1 1 2388;    67 1 1 2388];...
    
    
%% 
idPlot=1;
load pplot.mat;
% figure's PaperUnits: 'inches'
optsExpFig.Height=2.1/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
optsExpFig.Width=3.05/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
optsExpFig.Width2=2.34/2.54; % 1.7*2*1.2; %
optsExpFig.Width3=2.11/2.54; % 1.3*2*1.2;
optsExpFig.FontSize='6';
optsExpFig.FontMode='fixed';
% optsExpFig.FontName='Helvetica';
% optsExpFig.FontWeight='normal';
%         rmfield(optsExpFig,'Width');
%         rmfield(optsExpFig,'Height');
                optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
%         optsExpFig.Color='cmyk';
optsExpFig.LockAxes=1;
optsExpFig.LineMode='scaled';
optsExpFig.LineWidthMin=0.5;
optsExpFig.LineWidthMin=1.5;

optsExpFig.Renderer='painters';

% optsExpFig.Height='5.5'; % 7;
% optsExpFig.Width='5.5';
% optsExpFig.Format='pdf'; % 'png';
% optsExpFig.LineMode='scaled';
        
% new stat: poisson regression
nAnal=11; % 10;
nCoeff=67;
Formula=cell(nAnal,1);
CoefficientNames=[]; % =cell(nCoeff,1);
coeff=[];coeffTmp=[];
pval=[];pvalTmp=[];


for iAnimal=1:length(animalNm)
    for i=1:length(fn{iAnimal}) % session
        fnm=fn{iAnimal}{i}; % H_RSGprior_20161203
        fid=fnm(end-5:end); % 161203
        if sum(cname(:,1)==str2num(fid)) % session is in example cell list
            disp(['===== ' fnm ' =====']);
            load(fid); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo
            
%             % tmp: check if tInfo
%             iTrial=unique(tInfo(:,1));
%             for iTrialId=length(iTrial):(-1):1
%                 if nnz(tInfo(:,1)==iTrial(iTrialId))>5
%                 break;
%             end
%             
%             tmpTInfo=diff(tInfo(tInfo(:,1)==100,3))';
%             disp([tmpTInfo(2:5)*1000; tInfo(find(tInfo(:,1)==100,1,'first'),[9 10 7 8])]);
% %         end
%     end
% end
            
            % extract behavior
            beh=load([behDir '/' animalNm{iAnimal} '_RSGprior_DMFC.mat']);
            
            % remove outlier
             idOut2=makeIdOut(tInfo,beh,fid);
            
%             idOut2=false(size(tInfo,1),1);
%             movingBallTrials=unique(tInfo(:,1)); % size same as idNoise
%             idSess=str2num(fid)==beh.sessId;
%             if strcmp(fid,'161207'),
%                 idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;
%             end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
%             if nnz(idSess)~=length(movingBallTrials),disp('smth wrong for trial align'); end;
%             iOutSess=find(beh.idOut(idSess)); disp(['outlier trials (incl. abort/wait): ' num2str(length(iOutSess)) '/' num2str(length(movingBallTrials)) ',' ...
%                 num2str(length(iOutSess)/length(movingBallTrials)*100) '%']);
%             idValid=beh.t(idSess)>0&beh.t(idSess)<3*beh.T(idSess);
%             disp(['outlier trials (w/o abort/wait): ' num2str(nnz(beh.idOut(idSess)&idValid)) '/' num2str(length(movingBallTrials)) ',' ...
%                 num2str(nnz(beh.idOut(idSess)&idValid)/length(movingBallTrials)*100) '%']);
%             for iOut=1:length(iOutSess)
%                 idOut2(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
%             end
            
            % loop through each unit
            nUnit=size(idSU,1); % both SU & MU
            for j=1:nUnit
                if sum(idSU(j,1)==cname(cname(:,1)==str2num(fid),2)) % cell matched
                    disp(['cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust))]);
                    
                    % consistent with plotSDFexport
                    NEV=kilosort2NEV(sp,idClust,idSU);
                    tmpId=[idSU(j,1) idSU(j,1) 1 length(unique(tInfo(:,1)))]; % id(electrode,Unit,start trial,end trial)
                    NEV.MetaTags.DateTime=fid;
                    CoefficientNames=[];
                    
%                     try
                        statTmp=plotSDFexport_inFunc(NEV,tmpId,tInfo(~idOut2,:),eventId,fid,idPlot,optsExpFig,dirFig2singleNeuron); % p values cell
                        if idPlot
%                             set(gcf,'PaperPositionMode','auto');
%                             saveas(gcf,[neuDir 'PSTH_kilosort/' fid '_' num2str(idSU(j,1)) '_' sortQ{idSU(j,2)+1} '.png']);% sortQ={'MU','SU'};
                        end               
%                         waitforbuttonpress;
%                         close all;
%                     catch
%                         close all;
%                         disp(['ERROR cluster#' num2str(idSU(j,1)) ': #sp=' num2str(nnz(idSU(j,1)==idClust))]);
%                     end
                    
                end % cell matched                % if sum(num2str(idSU(j,1))==cname(cname(:,1)==num2str(fid),2)) % cell matched
            end % for j=1:nUnit
            
        end % if sum(cname(:,1)==num2str(fid))
    end % for i=1:length(fname)
    
end % for iAnimal=1:length(animalNm)






%%

function stat=plotSDFexport_inFunc(NEV,id,tInfo,eventId,strdate,idPlot,optsExpFig,dirFig2singleNeuron)

% input: NEV data, id(electrode,Unit,start trial,end trial), eventID, idPlot
% use Gaussian kernel for now (note also photodiode sync is not used forg now)

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

% dirFig2singleNeuron='/Users/hansem/Dropbox (MIT)/SFN17/2singleNeurons/';
idEyeOnly=0; % 1;

%% init
idE=id(1); idU=id(2); % unitMap=[0 1 2 3 255];  idU=unitMap(id(2));
load pplot.mat pplot tmpCmap;
sampF=double(NEV.MetaTags.SampleRes); % 30000 /s
sp=double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==idE&NEV.Data.Spikes.Unit==idU))/sampF*1000; % blackrock sample unit into ms
if isempty(sp), disp('no channel/unit'); return; end;
% wf=double(NEV.Data.Spikes.Waveform(:,NEV.Data.Spikes.Electrode==idE&NEV.Data.Spikes.Unit==idU)); % [48 x spikes]
% strdate=datestr(NEV.MetaTags.DateTime,'yymmdd');
clear NEV;
% unit 0 1 2 3 255>1 2 3 4 5?

% removing trials without no neural data (given uncertain end trials ID)
% update 8/30/17 remove initial trials too (kilosort)
tInitSp=min(sp);
tLastSp=max(sp); % [ms]
tInfo=tInfo(tInfo(:,1)>=id(3)&tInfo(:,1)<=id(4)&tInfo(:,3)*1000<tLastSp&tInfo(:,3)*1000>tInitSp,:); % start end trials
% tInfo=tInfo(tInfo(:,1)>=id(3)&tInfo(:,1)<=id(4),:); % start end trials

% treating noSet experiment
if size(tInfo,2)>12
    noSet=1;idSuccessNoSet=tInfo(:,14);idNoSet=tInfo(:,13);waitDur=400;
else
    noSet=0;
end

nPr=2;    prNm={'Short','Long'};
if idEyeOnly
    nEH=1;
else
    nEH=2;
    tmpCmap{1,2}=tmpCmap{1,1};tmpCmap{2,2}=tmpCmap{2,1};
end
ehNm={'Eye','Hand'};
nTarg=2; targNm={'Right','Left'};
nTspp=5; T{1}=unique(tInfo(tInfo(:,4)==1,7)); T{2}=unique(tInfo(tInfo(:,4)==0,7));
nSplitTp=1; % 5; 
t=tInfo(tInfo(:,8)>0,8); tSort=sort(t); tBound=tSort([1 round([(1/nSplitTp):(1/nSplitTp):1]*length(t))]); % 0 0.2 0.4 ... 1
nSplitRew=5; rewValid=tInfo(tInfo(:,8)>0&tInfo(:,12)>0,end); rSort=sort(rewValid); rBound=rSort([1 round([(1/nSplitRew):(1/nSplitRew):1]*length(rewValid))]); % 0 0.33 0.66 1

if idPlot
% figure; % set(gcf,'position',[0 0 1920 1200],'color','w');
end

nPlot=3;
yLimMat=NaN(nPlot,nEH,2);iYlim=0;

lw=1.5;
msize=4;

%% Gaussian kernel for now (note also photodiode sync is not used forg now)
td=25; nTd=3; xFiltG=-td*nTd:1:td*nTd; gaussampFilter=exp(-(xFiltG.^2)/(2*(td^2))); gaussampFilter=gaussampFilter/sum(gaussampFilter); % 40
% tg=1; td=20; xFilt=0:1:100; % smoothing as in Hanes et al., JNP98
% epspFilter=(1-exp(-xFilt/tg)).*exp(-xFilt/td); epspFilter=[zeros(length(epspFilter)-1,1);epspFilter(:)]; epspFilter=epspFilter/sum(epspFilter);

%% 4. Ready % photodiode sync [-200 480/800/1200]
% 10 ts(short/long) for each eye, hand (autumn/winter, spring/summer; use only 200 to avoid yellow)     nTspp=5;
if sum(eventId==4 | eventId==0) % 0 for plot all
    
    dur1=100; % ms
    periReadyWindow=240; %T{1}(1); % 400; %%%%%%
    periReadyWindow2=T{1}(end);
    periReadyWindow3=T{1}(1);
                nr=[]; nr2=[];nr3=[];iPr=[]; iEH=[]; iRL=[]; ts=[];% for stat
                iPr3=[]; iEH3=[]; iRL3=[];

%     figure(101); set(gcf,'position',pplot.rect1_3,'color','w'); hold all; % eye
%     figure(102); set(gcf,'position',pplot.rect2_3,'color','w'); hold all; % hand
%     tmp1=spring(round(nTspp*1.5));
%     tmp2=summer(round(nTspp*1.5));
%     tmpCmap={flipud(autumn(nTspp)), flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
%         winter(nTspp), tmp2(1:nTspp,:)}; % longEye, longHand
    legName=cell(nPr,nTspp,nEH);
hMat=[];%hold all;
    %hMat=cell(nPr,nTspp,nEH);
    
    r=cell(nPr,nTspp,nEH);
    for l=1:nEH   
        h1{l}=figure; setFigPos(l,1);ha;
%         figure(100+l);
%     if idPlot, subplot(mPlot,nPlot,2+(l-1)*nPlot);hold all;end
        for i=1:nPr,
            for j=nTspp:(-1):1 % 1:nTspp, % attrition
                dur2=T{i}(j); tDur=(-dur1+1):dur2;
                legName{i,j,l}=[prNm{i} num2str(T{i}(j))];
                    
                if noSet
                    tmpId=tInfo(:,2)==4 &...  % event id 4.ready
                        tInfo(:,4)==(2-i) &... % idShortTrial
                        tInfo(:,5)==(2-l) &... % idHandEye
                        tInfo(:,7)==T{i}(j) &... % ts % attrition %                 tInfo(:,7)>=T{i}(j) &... % ts
                        tInfo(:,8)>0 &... % only valid trials for now
                        ~idNoSet; % NoSet trials has max ts of each prior
                else
                    tmpId=tInfo(:,2)==4 &...  % event id 4.ready
                        tInfo(:,4)==(2-i) &... % idShortTrial
                        tInfo(:,5)==(2-l) &... % idHandEye
                        tInfo(:,7)==T{i}(j) &... % ts % attrition %                 tInfo(:,7)>=T{i}(j) &... % ts
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
%                 if (i==1 & j==nTspp) || (i==nPr & j==1), lw=4; % overlapped interval
%                     if idPlot, h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap{i,l}(j,:),'linewidth',lw},1); hMat(i,j,l)=h.mainLine;end
%                     nr3=[nr3;sum(r{i,j,l}(:,(dur1+periReadyWindow3):(dur1+periReadyWindow2)),2)]; % 250+480:250+800
%                     iPr3=[iPr3; repmat(prNm(i),[nnz(tmpId) 1])];iEH3=[iEH3; repmat(ehNm(l),[nnz(tmpId) 1])];iRL3=[iRL3; targNm(1+(tInfo(tmpId,6)==180))'];
%                 else
%                     lw=3; % 1.5;
                    if idPlot, 
%                         if j~=1 % attrition
%                             tmpX=[T{i}(j-1):dur2]';
%                             hMat(i,j,l)=plot(tmpX,mean(rConv(:,(end-length(tmpX)+1:end)),1),'color',tmpCmap{i,l}(j,:),'linewidth',lw);
%                         else
                            hMat(i,j,l)=plot(tDur(:)',mean(rConv,1),'color',tmpCmap{i,l}(j,:),'linewidth',lw); hold on;                            
%                         end
                    end
                    if i==nPr % inc. all long prior
                        nr3=[nr3;sum(r{i,j,l}(:,(dur1+periReadyWindow3):(dur1+periReadyWindow2)),2)];
                        iPr3=[iPr3; repmat(prNm(i),[nnz(tmpId) 1])];iEH3=[iEH3; repmat(ehNm(l),[nnz(tmpId) 1])];iRL3=[iRL3; targNm(1+(tInfo(tmpId,6)==180))'];
                    end
%                 end;            
                nr=[nr;sum(r{i,j,l}(:,1:dur1),2)]; nr2=[nr2;sum(r{i,j,l}(:,(dur1+1):(dur1+periReadyWindow)),2)];
                iPr=[iPr; repmat(prNm(i),[nnz(tmpId) 1])];iEH=[iEH; repmat(ehNm(l),[nnz(tmpId) 1])];iRL=[iRL; targNm(1+(tInfo(tmpId,6)==180))'];
                ts=[ts; repmat(T{i}(j),[nnz(tmpId) 1])];
                
            end % j for different ts            
        end % i for ShortLong    
        if idPlot, 
%         title(ehNm{l});
        axis tight; 
%         xlabel('Time from Ready (ms)'); ylabel('Spike rate (/s)');%plotVertical(gca);
        hMatTmp=hMat(:,:,l); legNmTmp=legName(:,:,l);
%         legend(hMatTmp(:),legNmTmp(:),'location','best'); legend boxoff;
        drawnow;
        set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',unique([0 min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
            'xticklabel',{0,min(T{1})/1000,[],max(T{1})/1000,[],max(T{2})/1000});
        end
        iYlim=iYlim+1; yLimMat(1,l,:)=ylim;
    end % l for EyeHand
    
    % plot dot @ set
%     if idEyeOnly
    for l=1:nEH
        figure(h1{l});
        for i=1:nPr,
            for j=nTspp:(-1):1 % 1:nTspp, % attrition
                
                idLast=find(~isnan(get(hMat(i,j,l),'YData')),1,'last'); % hMat(i,j,l)
                Xdata=get(hMat(i,j,l),'XData');
                Ydata=get(hMat(i,j,l),'YData');                
                plot(Xdata(idLast),Ydata(idLast),'o','markerfacecolor','w','color',tmpCmap{i,l}(j,:),'markersize',msize); % plot marker for the last points                    
            end % j for different ts
        end % i for ShortLong
    end % l for EyeHand
%     end
    

end

%% 5. Set % photodiode sync [-480/800/1200 min(tp)]
% 10 ts(short/long) for each eye, hand (autumn/winter, spring/summer; use only 200 to avoid yellow)
if sum(eventId==5 | eventId==0) % 0 for plot all
%     h2=figure; set(gcf,'position',pplot.rect1_2);
%     dur2=floor(min(t));
periSetWindow=240;%T{1}(1); % 400;
periSetWindow2=240; % 320
%     figure(201); set(gcf,'position',pplot.rect1_4,'color','w'); hold all; % eye
%     figure(202); set(gcf,'position',pplot.rect2_4,'color','w'); hold all; % hand
            nr=[]; nr2=[];nr3=[];iPr=[]; iEH=[];iPr2=[]; iEH2=[]; iRL=[];ts=[];tp=[]; iShortOverlap=[];% for stat
            % find valid min tp
            tmp=unique(tInfo(:,[1 7 8]),'rows'); tpMinTs=tmp(tmp(:,2)==T{1}(1)&tmp(:,3)>0,3); minTpMinTs=mean(tpMinTs)-2*std(tpMinTs);
            dur2min=floor(minTpMinTs);
            
%     tmp1=spring(round(nTspp*1.5));
%     tmp2=summer(round(nTspp*1.5));
%     tmpCmap={flipud(autumn(nTspp)), flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
%         winter(nTspp), tmp2(1:nTspp,:)}; % longEye, longHand
    legName=cell(nPr,nTspp,nEH);
    prNm={'Short','Long'};
    ehNm={'Eye','Hand'};hMat=[];
%     hMat=cell(nPr,nTspp,nEH);
    
    r=cell(nPr,nTspp,nEH);
    for l=1:nEH        
%         if idPlot, subplot(mPlot,nPlot,3+(l-1)*nPlot);hold all;end
% %         figure(200+l);
        h2{l}=figure; setFigPos(l,2);  ha;
        for i=1:nPr,
            for j=nTspp:(-1):1 % 1:nTspp,
                dur1=T{i}(j); % ms
                dur2=max(dur2min,floor(min(tInfo(tInfo(:,5)==(2-l) &tInfo(:,4)==(2-i) &tInfo(:,7)==T{i}(j)&tInfo(:,8)>0,8)))); %%%%%%
                tDur=(-dur1+1):dur2;
                legName{i,j,l}=[prNm{i} num2str(T{i}(j))];
                    
                if noSet
                    tmpId=tInfo(:,2)==5 &...  % event id 5.set
                        tInfo(:,4)==(2-i) &... % idShortTrial
                        tInfo(:,5)==(2-l) &... % idHandEye
                        tInfo(:,7)==T{i}(j) &... % Ts % tInfo(:,7)==T{i}(j) &... % Ts
                        tInfo(:,8)>dur2min &... % only valid trials for now; 0 %%%%%%
                        ~idNoSet; % NoSet trials has max ts of each prior
                else
                    tmpId=tInfo(:,2)==5 &...  % event id 5.set
                        tInfo(:,4)==(2-i) &... % idShortTrial
                        tInfo(:,5)==(2-l) &... % idHandEye
                        tInfo(:,7)==T{i}(j) &... % Ts % tInfo(:,7)==T{i}(j) &... % Ts
                        tInfo(:,8)>dur2min; % only valid trials for now; 0 %%%%%%
                    
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
%                 if (i==1 & j==nTspp) || (i==nPr & j==1), lw=4; % overlapped interval
%                     if idPlot, 
%                     h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap{i,l}(j,:),'linewidth',lw},1); hMat(i,j,l)=h.mainLine;end
%                     nr2=[nr2;sum(r{i,j,l}(:,(dur1-periSetWindow2):(dur1)),2)];iShortOverlap=[iShortOverlap;repmat(prNm(i),[nnz(tmpId) 1])];
%                     iPr2=[iPr2; repmat(prNm(i),[nnz(tmpId) 1])];iEH2=[iEH2; repmat(ehNm(l),[nnz(tmpId) 1])];
%                 else
%                     lw=3; % 1.5;
                    if idPlot, 
                    hMat(i,j,l)=plot(tDur(:)',mean(rConv,1),'color',tmpCmap{i,l}(j,:),'linewidth',lw); end
%                 end;
                nr=[nr;sum(r{i,j,l}(:,(dur1-periSetWindow+1):dur1),2)]; 
                nr3=[nr3;sum(r{i,j,l}(:,(dur1+1):(dur1+dur2min)),2)];
                iPr=[iPr; repmat(prNm(i),[nnz(tmpId) 1])];iEH=[iEH; repmat(ehNm(l),[nnz(tmpId) 1])];iRL=[iRL; targNm(1+(tInfo(tmpId,6)==180))'];  tp=[tp;tInfo(tmpId,8)];
                ts=[ts; repmat(T{i}(j),[nnz(tmpId) 1])];
            end % j for different ts
        end % i for ShortLong        
        
        if idPlot,
            %         title(ehNm{l});
            axis tight;
            %         xlabel('Time from Set (ms)'); %ylabel('spike rate (/s)');plotVertical(gca);
            xlim([-200 max(T{1})]); %%%%%
%             xlim([-max(T{2})-1 max(T{1})]);
            hMatTmp=hMat(:,:,l); legNmTmp=legName(:,:,l);
            %         legend(hMatTmp(:),legNmTmp(:),'location','best'); legend boxoff;
            drawnow;
            set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',[0 max(T{1})/2 max(T{1})], ...% unique([0 -min(T{1}) -mean(T{1}) -max(T{1}) -mean(T{2}) -max(T{2}) min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
                'yticklabel',[],'xticklabel',[0 max(T{1})/2/1000 max(T{1})/1000]); % {-max(T{2}),[],-max(T{1}),[],-min(T{1}),0,min(T{1}),[],max(T{1}),[],max(T{2})});
%             xlim([-max(T{2}) max(T{1})]);
        end
        iYlim=iYlim+1; yLimMat(2,l,:)=ylim;
        
    end % l for EyeHand
    
    % plot dot @ set
%     if idEyeOnly
    for l=1:nEH
        figure(h2{l});
        for i=1:nPr,
            for j=nTspp:(-1):1 % 1:nTspp, % attrition
                Xdata=get(hMat(i,j,l),'XData');
                Ydata=get(hMat(i,j,l),'YData');                
                plot(Xdata(Xdata==0),Ydata(Xdata==0),'o','markerfacecolor','w','color',tmpCmap{i,l}(j,:),'markersize',msize); % plot marker for the last points                    
            end % j for different ts
        end % i for ShortLong
    end % l for EyeHand
%     end

end

%% 6. Production [-tp 50]
% 5 split of tp (jet)
% nSplitTp=5; t=tInfo(tInfo(:,8)>0); tSort=sort(t); tBound=tSort(round([(1/nSplitTp):(1/nSplitTp):1]*length(t))); % 0 0.2 0.4 ... 1
% td=10; nTd=3; xFiltG=-td*nTd:1:td*nTd; gaussampFilter=exp(-(xFiltG.^2)/(2*(td^2))); gaussampFilter=gaussampFilter/sum(gaussampFilter); % 40

if sum(eventId==6 | eventId==0) % 0 for plot all
    
    dur2=100; % 50;
    hold all;
    nr=[]; iEH=[]; iRL=[];tp=[];% for stat
    % find valid min tp
    tmp=unique(tInfo(:,[1 7 8]),'rows'); tpMinTs=tmp(tmp(:,2)==T{1}(1)&tmp(:,3)>0,3); minTpMinTs=mean(tpMinTs)-2*std(tpMinTs);
    dur2min=400; % floor(minTpMinTs);
            
% prior/ts
r=cell(nPr,nTspp,nEH);
for l=1:nEH
    h3{l}=figure; setFigPos(l,3); ha;
    
    for i=1:nPr,
        for j=nTspp:(-1):1 % 1:nTspp,

            if noSet
                tmpId=tInfo(:,2)==6 &...  % event id 6.production
                    tInfo(:,4)==(2-i) &... % idShortTrial
                    tInfo(:,5)==(2-l) &... % idHandEye
                    tInfo(:,7)==T{i}(j)&... %  &... % Ts % tInfo(:,7)==T{i}(j) &... % Ts
                    ~idNoSet; % NoSet trials has max ts of each prior
                %                 tInfo(:,8)>dur2min; % tBound(i); % only valid trials for now %%%%%
            else
                tmpId=tInfo(:,2)==6 &...  % event id 6.production
                    tInfo(:,4)==(2-i) &... % idShortTrial
                    tInfo(:,5)==(2-l) &... % idHandEye
                    tInfo(:,7)==T{i}(j); %  &... % Ts % tInfo(:,7)==T{i}(j) &... % Ts
            end

tmpTp=tInfo(tmpId,8);
            dur1=round(min(tmpTp)); % dur2min; % max(dur2min,round(tBound(i))); % ms %%%%%%
            tDur=(-dur1+1):dur2;
            
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
            %                 if (i==1 & j==nTspp) || (i==nPr & j==1), lw=4; % overlapped interval
            %                     if idPlot,
            %                     h=shadedErrorBar(tDur(:)',mean(rConv,1),sem(rConv,1),{'-','color',tmpCmap{i,l}(j,:),'linewidth',lw},1); hMat(i,j,l)=h.mainLine;end
            %                     nr2=[nr2;sum(r{i,j,l}(:,(dur1-periSetWindow2):(dur1)),2)];iShortOverlap=[iShortOverlap;repmat(prNm(i),[nnz(tmpId) 1])];
            %                     iPr2=[iPr2; repmat(prNm(i),[nnz(tmpId) 1])];iEH2=[iEH2; repmat(ehNm(l),[nnz(tmpId) 1])];
            %                 else
%             lw=3; % 1.5;
            if idPlot,
                try
                    hMat(i,j,l)=plot(tDur(:)',mean(rConv,1),'color',tmpCmap{i,l}(j,:),'linewidth',lw);
                catch
                    disp([]);
                end
            end
            %                 end;
%             nr=[nr;sum(r{i,j,l}(:,(dur1-periSetWindow+1):dur1),2)];
%             nr3=[nr3;sum(r{i,j,l}(:,(dur1+1):(dur1+dur2min)),2)];
%             iPr=[iPr; repmat(prNm(i),[nnz(tmpId) 1])];iEH=[iEH; repmat(ehNm(l),[nnz(tmpId) 1])];iRL=[iRL; targNm(1+(tInfo(tmpId,6)==180))'];  tp=[tp;tInfo(tmpId,8)];
%             ts=[ts; repmat(T{i}(j),[nnz(tmpId) 1])];
        end % j for different ts
    end % i for ShortLong
    if idPlot,
        axis tight;
        drawnow;
        set(gca,'xtick',([-max(T{1}) -max(T{1})/2 0]),'tickdir','out','TickLength', [.02 .02],...
            'yticklabel',[],'xticklabel',[-max(T{1})/1000 -max(T{1})/2/1000 0]);
        xlim([-max(T{1}) dur2]); %%%%%
        
    end
    iYlim=iYlim+1; yLimMat(3,l,:)=ylim;
end % l for EyeHand
    
    
    
end


%%
% optsExpFig.Height=2*1.2; % 0.8; % '1.2'; % '2'; % 7;
% optsExpFig.Width=1.5*1.2*2; % '2.4'; % '4';
% optsExpFig.FontSize='16';
% %         rmfield(optsExpFig,'Width');
% %         rmfield(optsExpFig,'Height');
%                 optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
% %         optsExpFig.Color='cmyk';
% optsExpFig.LockAxes=1;
% optsExpFig.LineMode='scaled';

Width2=optsExpFig.Width2;
 Width3=optsExpFig.Width3;
 optsExpFig=rmfield(optsExpFig,'Width2');
 optsExpFig=rmfield(optsExpFig,'Width3');
 
%matching ylim
if idPlot, 
    mkdir([dirFig2singleNeuron strdate '_' num2str(id(1))]);
    
    for l=1:nEH
    
        figure(h1{l});set(gca,'ylim',[min(yLimMat(:,l,1)) max(yLimMat(:,l,2))]);plotVertical(gca);
%         ytickTmp=get(gca,'ytick');    %set(gca,'ytick',unique([floor(min(ytickTmp)/10)*10 floor(max(ytickTmp)/10)*10]));
        set(gca,'ytick',unique(round(min(yLimMat(:,l,1))+(max(yLimMat(:,l,2))-min(yLimMat(:,l,1)))*[1/3 2/3])));
        applytofig(gcf,optsExpFig);
%         saveas(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/' strdate '_' num2str(id(1)) '_' ehNm{l} '_periReady.fig']);
%         %     axis off;
%         exportfig(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/' strdate '_' num2str(id(1)) '_' ehNm{l} '_periReady.eps'],optsExpFig); % pdf tiff
%         %     saveas(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/periReady.tiff']);
        
        figure(h2{l});set(gca,'ylim',[min(yLimMat(:,l,1)) max(yLimMat(:,l,2))]);plotVertical(gca);
        set(gca,'ytick',unique(round(min(yLimMat(:,l,1))+(max(yLimMat(:,l,2))-min(yLimMat(:,l,1)))*[1/3 2/3])));
%         set(gca,'ytick',unique([ceil(min(ytickTmp)/10)*10 floor(max(ytickTmp)/10)*10]));
        optsExpFig.Width=Width2; % 1.7*2*1.2; % '2.4';% 4.5; % optsExpFig.Width*length(tDur)/1300;
        applytofig(gcf,optsExpFig);
%         saveas(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/' strdate '_' num2str(id(1)) '_' ehNm{l} '_periSet' num2str(1) '.fig']);
%         %     axis off;
%         exportfig(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/' strdate '_' num2str(id(1)) '_' ehNm{l} '_periSet' num2str(1) '.eps'],optsExpFig);% pdf tiff
%         %     saveas(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/periSet' num2str(1) '.tiff']);
        
        figure(h3{l});set(gca,'ylim',[min(yLimMat(:,l,1)) max(yLimMat(:,l,2))]);plotVertical(gca);
        set(gca,'ytick',unique(round(min(yLimMat(:,l,1))+(max(yLimMat(:,l,2))-min(yLimMat(:,l,1)))*[1/3 2/3])));
%         set(gca,'ytick',unique([ceil(min(ytickTmp)/10)*10 floor(max(ytickTmp)/10)*10]));
        optsExpFig.Width=Width3; % 1.3*2*1.2; % '2'; % 1.2; % optsExpFig.Width*length(tDur)/1300;
        applytofig(gcf,optsExpFig);
%         optsExpFig.LockAxes=0; % don't know why xticklabel flipped
%         saveas(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/' strdate '_' num2str(id(1)) '_' ehNm{l} '_periGo.fig']);
%         %     axis off;
%         exportfig(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/' strdate '_' num2str(id(1)) '_' ehNm{l} '_periGo.eps'],optsExpFig);% pdf tiff
%         %     saveas(gcf,[dirFig2singleNeuron strdate '_' num2str(id(1)) '/periGo.tiff']);
    
    end
    
% for i=1:(mPlot*nPlot)
%     h=subplot(mPlot,nPlot,i);
%     set(h,'position',[0.2*rem(i-1,5)+0.02*(i==1||i==nPlot+1) 0.5*(i<6)+0.03 0.18-0.02*(i==1||i==nPlot+1) 0.42]);
%     ylim([min(yLimMat(:,1)) max(yLimMat(:,2))]);
%     plotVertical(gca);
%     if i~=1 && i~=(nPlot+1), set(gca,'yticklabel',[]); end;
% end
end
% disp('');
% exportfig(gcf,[neuDir fid '_' num2str(cname{i}(j,1)) '_' num2str(cname{i}(j,2)) 'new.png'],'color','rgb','Height',13,'Width',11,'format','png','Reference',hTmp);


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
% end_

function idOut=makeIdOut(tInfo,beh,fid)

% function idOut=makeIdOut(tInfo,beh,fid)
% make idOut [size(tInfo,1) x 1] for a given session's tInfo
% tInfo
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
%
% tInfo: trial#, stageCode, t(blackRock), 
%           idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)]
%
% fid: string e.g. 161218

%% remove outlier

% initial
idOut=false(size(tInfo,1),1); % 42602
movingBallTrials=unique(tInfo(:,1)); % ~4000
idSess=str2num(fid)==beh.sessId; % ~40000

% special treat for individual sessions
if strcmp(fid,'161207'),
    idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;
end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
if nnz(idSess)~=length(movingBallTrials),disp('smth wrong for trial align'); end;

iOutSess=find(beh.idOut(idSess)); % 1157

% report
disp(['outlier trials (incl. abort/wait): ' num2str(length(iOutSess)) '/' num2str(length(movingBallTrials)) ',' ...
    num2str(length(iOutSess)/length(movingBallTrials)*100) '%']);
idValid=beh.t(idSess)>0&beh.t(idSess)<3*beh.T(idSess);
disp(['outlier trials (w/o abort/wait): ' num2str(nnz(beh.idOut(idSess)&idValid)) '/' num2str(length(movingBallTrials)) ',' ...
    num2str(nnz(beh.idOut(idSess)&idValid)/length(movingBallTrials)*100) '%']);

% main
for iOut=1:length(iOutSess)
    idOut(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
end

function NEV=kilosort2NEV(sp,idClust,idSU)

% convert kilosort result to NEV

NEV.MetaTags.SampleRes=30000;
NEV.Data.Spikes.TimeStamp=sp;
% in kilosort idClust is unique across electrodes & units
NEV.Data.Spikes.Electrode=idClust;
NEV.Data.Spikes.Unit=idClust;