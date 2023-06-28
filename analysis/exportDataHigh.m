function D=exportDataHigh(NEV,id,tInfo,eventId,varargin)
%
% input: NEV data, id(electrode,Unit,start trial,end trial), tInfo,
% eventID, trialId (if empty, trial-averaged)

% tInfo: trial#(movingBall_trials; not complete), stageCode, t(blackRock), 
%           idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)]

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

% output: trial-avaraged PSTH % 40 cond x time

%% exclusion crteria
nTrialExc=5; % if # trials per condition<nTrialExc, cell is removed from analysis

%% epoch
nEpoch=5;
% 1) [fix-100 fix+500]
% 2) [tOn-500 tOn+250]
% 3) [ready-250 ready+ts]
% 4) [set-ts set+tp]
% 5) [production-tp production+100]

%% init
T{1}=unique(tInfo(tInfo(:,4)==1,7)); T{2}=unique(tInfo(tInfo(:,4)==0,7)); % setting ts before selecting subset trials

idE=id(1); idU=id(2); % unitMap=[0 1 2 3 255];  idU=unitMap(id(2));
sampF=double(NEV.MetaTags.SampleRes); % 30000 /s
sp=double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Electrode==idE&NEV.Data.Spikes.Unit==idU))/sampF*1000; % blackrock sample unit into ms

if ~isempty(varargin) % restrict subset of trials
    tmpTrialId=varargin{1};
    if numel(tmpTrialId)==1 % only one trial selected
        tInfo=tInfo(tInfo(:,1)==tmpTrialId(1),:);
    elseif numel(tmpTrialId)==2  % trial range
        tInfo=tInfo(tInfo(:,1)>=tmpTrialId(1)&tInfo(:,1)<=tmpTrialId(2),:);
    else % assuming id for selected trials are provided
        tInfo=tInfo(tmpTrialId,:);
    end
else % all trials from the cell    
    % removing trials without no neural data (given uncertain end trials ID)
    % update 8/30/17 remove initial trials too (kilosort)
    tInitSp=min(sp);
    tLastSp=max(sp); % [ms]
    tInfo=tInfo(tInfo(:,1)>=id(3)&tInfo(:,1)<=id(4)&tInfo(:,3)*1000<tLastSp&tInfo(:,3)*1000>tInitSp,:); % start end trials
    % tInfo=tInfo(tInfo(:,1)>=id(3)&tInfo(:,1)<=id(4),:); % start end trials
end
load pplot.mat;
if isempty(sp), disp('no channel/unit'); D=[]; return; end;
% wf=double(NEV.Data.Spikes.Waveform(:,NEV.Data.Spikes.Electrode==idE&NEV.Data.Spikes.Unit==idU)); % [48 x spikes]
% strdate=datestr(NEV.MetaTags.DateTime,'yymmdd');
clear NEV;
% unit 0 1 2 3 4 5 255>unclassified 1 2 3 4 5 noise

nPr=2;    prNm={'Short','Long'};
nEH=2;ehNm={'Eye','Hand'};
nTarg=2; targNm={'Right','Left'};
nTspp=5; 
nSplitTp=5; t=tInfo(tInfo(:,8)>0,8); tSort=sort(t); tBound=tSort([1 round([(1/nSplitTp):(1/nSplitTp):1]*length(t))]); % 0 0.2 0.4 ... 1
nSplitRew=5; rewValid=tInfo(tInfo(:,8)>0&tInfo(:,12)>0,12); rSort=sort(rewValid); rBound=rSort([1 round([(1/nSplitRew):(1/nSplitRew):1]*length(rewValid))]); % 0 0.33 0.66 1

nCond=nPr*nTspp*nEH*nTarg;
D=cell(nCond,1);

%% Gaussian kernel for now (note also photodiode sync is not used forg now)
td=30; nTd=4; xFiltG=-td*nTd:1:td*nTd; gaussampFilter=exp(-(xFiltG.^2)/(2*(td^2))); gaussampFilter=gaussampFilter/sum(gaussampFilter); % 40
% tg=1; td=20; xFilt=0:1:100; % smoothing as in Hanes et al., JNP98
% epspFilter=(1-exp(-xFilt/tg)).*exp(-xFilt/td); epspFilter=[zeros(length(epspFilter)-1,1);epspFilter(:)]; epspFilter=epspFilter/sum(epspFilter);

%% 
%% 1) [fix-100 fix+500]
% 2. Fix % photodiode sync [-100 500]
% ShortEye, ShortHand, LongEye, LongHand (Red, Orange, Blue, Green)

%% 2) [tOn-500 tOn+250]
% 3. target On % photodiode sync [-500 250]

%% 3) [ready-250 ready+ts]
% 4. Ready % photodiode sync [-200 480/800/1200]
% 10 ts(short/long) for each eye, hand (autumn/winter, spring/summer; use
% only 200 to avoid yellow)     nTspp=5;

%% 4) [set-ts set+tp]
% 5. Set % photodiode sync [-480/800/1200 min(tp)]
% 10 ts(short/long) for each eye, hand (autumn/winter, spring/summer; use
% only 200 to avoid yellow)

%% 5) [production-tp production+100]
% 6. Production [-tp 100]

iD=1;
for i=1:nPr
    for j=1:length(T{i})
        for k=1:nEH
            for m=1:nTarg
                if eventId==1 % fixation
                    dur1=100; dur2=500; tDur=(-dur1+1):dur2;durmin=0;
                elseif eventId==2 % targetON
                    dur1=500; dur2=250; tDur=(-dur1+1):dur2;durmin=0;
                elseif eventId==3 % periReady
                    dur1=250; dur2=T{i}(j); tDur=(-dur1+1):dur2;durmin=0;
                elseif eventId==4 % periSet
                    tmp=unique(tInfo(:,[1 4 5 6 7 8]),'rows'); % movingBall_trials, idShortTrial(4), idHandEye(5), theta(6), T(7), t(8)
                    tpMinTs=tmp(tmp(:,2)==(2-i) &tmp(:,3)==(2-k) &tmp(:,4)==(m-1)*180 &tmp(:,5)==T{i}(j)&tmp(:,6)>0,6);
                    if isempty(tpMinTs) % use max(tp|minTs) regardless of type
                        if size(tmp,1)==1 % single trial
                            durmin=floor(tmp(end));
                        else
                            durmin=floor(max(tmp(tmp(:,5)==T{1}(1)&tmp(:,6)>0,6)));
                        end
                    else
                        minTpMinTs=max([min(tpMinTs) mean(tpMinTs)-2*std(tpMinTs)]);durmin=floor(minTpMinTs);
                    end
                    dur2=durmin;
                    dur1=T{i}(j);tDur=(-dur1+1):dur2;
                elseif eventId==5 % periGo
                    tmp=unique(tInfo(:,[1 4 5 6 7 8]),'rows'); % movingBall_trials, idShortTrial(4), idHandEye(5), theta(6), T(7), t(8)
                    tpMinTs=tmp(tmp(:,2)==(2-i) &tmp(:,3)==(2-k) &tmp(:,4)==(m-1)*180 &tmp(:,5)==T{i}(j)&tmp(:,6)>0,6); 
                    if isempty(tpMinTs) % use max(tp|minTs) regardless of type
                        durmin=floor(max(tmp(tmp(:,5)==T{1}(1)&tmp(:,6)>0,6)));
                    else
                        minTpMinTs=max([min(tpMinTs) mean(tpMinTs)-2*std(tpMinTs)]);durmin=floor(minTpMinTs);
                    end                    
                    dur1=durmin;
                    dur2=100;tDur=(-dur1+1):dur2;
                end
                
                
                % selecting trials
                tmpId=tInfo(:,2)==(1+eventId) &... % event id 2.fix
                    tInfo(:,4)==(2-i) &... % idShortTrial
                    tInfo(:,5)==(2-k) &...  % idHandEye
                    tInfo(:,6)==(m-1)*180 &... % target location
                    tInfo(:,7)==T{i}(j) &... % ts
                    tInfo(:,8)>=durmin; % tp>mean(tp|ts,EH,Targ)-2*std to prevent data too short
                tE=tInfo(tmpId,3)*1000; % into ms
                t1=tE-dur1;            t2=tE+dur2;
                
                if nnz(tmpId)<= nTrialExc % ==0
                    D=[];
                    return;
%                     tmpR=zeros(1,dur1+dur2); % all zeros if empty>excluded in DataHigh later if FR too low                
                else
                    tmpR=zeros(nnz(tmpId),dur1+dur2); % [# trials x timePoints]
                    for kT=1:nnz(tmpId)
                        tSp=sp(t1(kT)<sp & sp<=t2(kT))-t1(kT);
                        if ~isempty(tSp),tmpR(kT,ceil(tSp))=1;end;
                    end % k for trials
                end
%                 if nnz(tmpR)==0 % for debug
%                     disp('');
%                 end

                % trial average
                if isempty(varargin) % all trials from the cell
                    D{iD}=mean(tmpR,1); % mean across trials
                else
                    if numel(tmpTrialId)==1 % only one trial
                        tmp=unique(tInfo(:,[1 4 5 6 7 8]),'rows');% movingBall_trials, idShortTrial(4), idHandEye(5), theta(6), T(7), t(8)
                        if tmp(:,2)==(2-i) &tmp(:,3)==(2-k) &tmp(:,4)==(m-1)*180 &tmp(:,5)==T{i}(j)&tmp(:,6)>0 % nnz(tmpR)>0 % 
                            D=tmpR; return; % regardless of spikes or not
                        end
                    elseif numel(tmpTrialId)==2 % trial range
                        D{iD}=tmpR;
                    else % assuming id for selected trials are provided
                        D{iD}=mean(tmpR,1); % mean across trials
                    end
                end
%                 D(iD).condition=[prNm{i} num2str(T{i}(j)) '/' ehNm{k} '/' targNm{m}];
%                 D(iD).data=[];
%                 D(iD).epochColors=tmpCmap{i,k}(j,:)./m;
                iD=iD+1;
            end % m nTarg
        end % k EH
    end % j ts
end % i prior
%     D=mean(tmpR,1);
    
