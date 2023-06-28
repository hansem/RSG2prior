function fig3de

% 2019/3/6
% fig3d: speed vs Xu
% use new tp data set: % 2019/1/20% re-run after manual update of behavioral idOut% remove 'bin20_smth40' at file name end
% idNewTraj=1
% initT=1; % now speed measured from set

% 2018/6/11
% fig3f: speed vs IC
% fig3g: tp vs speed

% speed: v(idAnimal,iEH,iTh,idPr,iTs,1:nV)=tmpV; % [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
% tp: mtp(idAnimal,iEH,iTarg,iPr,iTs)=mean(beh.t(idTp));
% IC: sICa (axIC=normalize(mean(diff(sSetTmp,1,1),1)',1); % [dim x 1] % for now across prior)

% isSepPr: if 1, use separate readout vector for sICa

% 2018/5/16
% change te vs ts into te-ts vs ts (idTpMTs)
% also trimmed

% 2018/4/9: direct link to bias
% apply kinet to tp>bias in speed, IC
% plot S(IC), 1/speed = f(ts) across all data sets (trial-avg PSTH): tSp, spl
% linear regression to S(IC), 1/speed to match it to BLS
% prior-specific IC axis, speed after IC (200; initT)

% figures
% 1: tp & speed
% 2: speed & s(IC)
% 3
% 4: t(nearestS) & tRef (short)
% 5: t(nearestS) & tRef (long)
% 6: distance(nearestS) & tRef (short)
% 7
% 8: distance(nearestS) & tRef (long)
% 1000: te(speed) & ts
% 1001: te(IC) & ts

% 2018/2/23: cosyne
% 1) correlation b/t tp (avgAcrSess) and speed (avgAcrTimeAfterIC) across all data sets
% 2) 1D projection on IC axis vs speed

% 2018/2/5 Q: why bigger bias for long?
% 1) IC vector size (ref to each prior mean) or their relative angle (more curved?)
% 2) angle b/t overlap vector & time/ts vector: common time-warped axis?
% 3) angle/magnitude of time vectors (ref: prior mean)

% 2018/1/29
% 1) testing nontrivial prediction: all/many curved neuron respond to set
%     (or in the same subspace as set-responsive neurons), revealing funcitonal role of curved trajectory
%    - angle between ts vector & time vector: orthogonal?
%    - assuming curved trajectory is in 2D (check: scree?), projecting time vector to the 2D plane>angle b/t the projected time vector and rotating vector
%    - do PCA separately for set and IC: angle between proj_matrix (or CCA?)
%    - regression? how to choose set-responsive neuron?
% 2) check angle/magnitude b/t readout vectors across ts, prior: 2D cmap, if <90, parellel-not state dependent
%    - kinet: nearest distance coalesce during set transient?
% 3) curvature b/t Set vs IC (200ms): scatter, separately for short and long

% available from analSpeed_ts: speed, |2nd derivative|, curvature, kinet
% + angle b/t velocity vectors of short vs long as a function of time: if <90, parallel
% + shortest path length between short and long: constant?
% + angle b/t tangent and shortest path vector: if <90, on the same manifold

% 2018/1/19
% estimating instantaneous speed during tp

% OLD note
% plot state space trajectory for each condition (of EH and targets) & two epoches
% data in DataHigh format, 40 conditions (2 pr x 5 ts x 2 EH x 2 targets)
% input: filename (assuming in /Users/hansem/Dropbox (MIT)/psthDataHigh/)
%   e.g. plotTraj('traj_periSet_G_attrition.mat'); plotTraj('traj_periSet_H_attrition.mat')
% if nm contains '_attrition', attrition
% assuming two epoches exist (epochStart(2)) as in periset
% save both in png and fig formats
% e.g. plotTraj('traj_periSet_ER_H.mat')

%% initial
idNewTraj=1;
idTpMTs=1;  % plot te-ts vs ts

idDebug=0; %1;

idCosyne=0; % 1;

dim=inf; % 3; % 5; % 4; % 2; % 3; %%%%%
useOptimD=1; % 0; % if 1, use optimD
maxDim=10; % tp-specific; H/G, ER/EL/HR/HL: 6 5 7 6 8 9 7 8

isSepPr=0; % : if 1, use separate readout vector for sICa

%% figure setting
try
    load('/Users/hansem/Dropbox (MIT)/commonCode/pplot.mat'); % tmpCmap
catch
    try
        load('/Users/seonminahn/Dropbox (MIT)/commonCode/pplot.mat'); % tmpCmap
    catch
        load pplot.mat; % tmpCmap
    end
end

cmap={pplot.cmap{1}; pplot.cmap{3}; pplot.cmap{5}}; % for short/long/SL
cmap3=[tmpCmap{1,1};tmpCmap{2,1}];

tmpMarker={'^','v','<','>','o','s','d','h'};

if idCosyne
    % figure's PaperUnits: 'inches'
    optsExpFig.Height=2.1/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
    optsExpFig.Width=3.15/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
    Width2=2.34/2.54; % 1.7*2*1.2; %
    Width3=2.1/2.54; % 1.3*2*1.2;
    optsExpFig.FontSize='6';
    optsExpFig.FontMode='fixed';
    optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
    optsExpFig.LockAxes=1;
    optsExpFig.LineMode='scaled';
    optsExpFig.LineWidthMin=0.5;
    optsExpFig.LineWidthMin=1.2;
    optsExpFig.Renderer='painters';
    
    msize=10; % 4; % 6; % 2;
    lw=3; % 1.5;
    lw2=.5; % for each data set
    % H/G, E/H, R/L
    cmapHG=[0 0 0; .5 .5 .5]; % animal & E/H
    cmapCond={[0 0 0; 1 1 1];
        [.5 .5 .5; 1 1 1]};
    markers={'o','^','o','^'};
    
else
    % figure's PaperUnits: 'inches'
    optsExpFig.Height=4.2/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
    optsExpFig.Width=4.2/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
    optsExpFig.FontSize='6';
    optsExpFig.FontMode='fixed';
    optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
    optsExpFig.LockAxes=1;
    optsExpFig.LineMode='fixed'; % 'scaled';
    %     optsExpFig.LineWidthMin=0.5;
    %     optsExpFig.LineWidthMax=1.5;
    optsExpFig.LineWidth=.5; %1;
    optsExpFig.Renderer='painters';
    
    lw=1.5;
    lw2=1;
    msize=4; % 5; % 10; % 4;
    msize2=2;
    
    % H/G, E/H, R/L
    cmapHG=[0 0 0; .5 .5 .5]; % animal & E/H
    %     cmapCond={[0 0 0; 1 1 1];
    %         [.5 .5 .5; 1 1 1]};
    markers={'o','^','o','^'};
    
    nCmap=8; % nAnimal*nEH*nTarg;
    cmapCond=parula(nCmap); % hsv(nCmap);
    icmap=1;
    
    cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
    cmapPr2=[rgb('IndianRed'); rgb('DodgerBlue'); rgb('ForestGreen')]; % for each data set
end
idExpFig=0; % 0; % 1;
idPlot=0; % for plotTraj (projected onto meas space)

%% initial
S=initRSG2prior;
v2struct(S);

nTs=nTspp;
try
    cd(psthDir);
catch
    psthDir='/Users/seonminahn/Dropbox (MIT)/psthDataHigh';
    cd(psthDir);
end

binSize=20;
smthWidth=40;
nBins=round([T{1}(:); T{2}(:)]/binSize); % 24 28 ... 60
maxNbin=44; % tp-specific

% speed & arc length (whole meas/overlap)
v=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % speed [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
vv=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin,maxDim); % velocity vector vv(idAnimal,iEH,iTh,idPr,iTs,1:nV,1:dim)[2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints x dim]
v2=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % ||2nd derivative|| [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
c=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % curvature [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
nEpoch=1; % 3; % wholeMeas/overlap, short/long, SIL
x=nan(nAnimal,nEH,nTarg,nPr,nTs,nEpoch); % arc length [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]

%% kinet (distance between nearest neighbors)
% less noisy (distance is noisy)
% not highly dependent on binSize
% assuming same manifold

% use mean of each prior
refPr=[1 2]; % 2; % 1; % 2;
refTs= round(cellfun(@length,T)/2); % length(T{refPr}); %1; % length(T{refPr}); %1; % length(T{refPr}); % longest of long (overlap)

% speed & distance
sp=nan(nAnimal,nEH,nTarg,nPr,nTs); % [2HG x 2EH x 2RL x 2SL x 5ts]
dist=nan(nAnimal,nEH,nTarg,nPr,nTs); % [2HG x 2EH x 2RL x 2SL x 5ts]

% 1) angle b/t velocity vectors of short vs long as a function of time: if <90, parallel
% 2) shortest path length between short and long: constant?
% 3) angle b/t tangent and shortest path vector: if <90, on the same manifold

% shortest path legnth between short and long: after reaching final state for ts, meaningless
sSp=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin,maxDim); % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints]
tSp=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints] in refD's bin unit
spv=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin,maxDim); % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints x dim]

% 1) angle b/t velocity vectors of short vs long as a function of time: if <90, parallel
av=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin,nPr,nTs); % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints x (2SL x 5ts)]
av2=cell(nAnimal,nEH,nTarg,nPr,maxNbin); % [2HG x 2EH x 2RL x within/betweenPriors x timePoints]
% 2) shortest path length between short and long: constant?
spl=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints]
% 3) angle b/t tangent and shortest path vector: if <90, on the same manifold
aSpvTv=nan(nAnimal,nEH,nTarg,nPr,nTs,2,maxNbin); % [2HG x 2EH x 2RL x 2SL x 5ts x 2tv(ts,ref) x timePoints]


%% tp-specific
tInitCond=200;%400; % 200;
initT=1; % tInitCond/binSize; % 10th bin
t0sv=40/binSize; % given speed profiles: time vector [40 160]
t1sv=160/binSize;

% 1) testing nontrivial prediction: all/many curved neuron respond to set
%     (or in the same subspace as set-responsive neurons), revealing funcitonal role of curved trajectory
%    - angle between ts vector & time vector: orthogonal?                % extract vv from analSpeed_ts?
aTsvSv2=nan(nAnimal,nEH,nTarg,nPr,nTs); % prod-space % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints]; for short800, tsv to long800, for 1200, tsv from 1100
%    - assuming curved trajectory is in 2D (check: scree?), projecting time vector to the 2D plane>angle b/t the projected time vector and rotating vector
aTsvSv1=nan(nAnimal,nEH,nTarg,nPr,nTs); % meas-space
%    - do PCA separately for set and IC: angle between proj_matrix (or CCA?)
%    - regression? how to choose set-responsive neuron?
% 2) check angle/magnitude b/t readout vectors across ts, prior: 2D cmap, if <90, parellel-not state dependent
%    - kinet: nearest distance coalesce during set transient?
asv=nan(nAnimal,nEH,nTarg,nPr*nTs,nPr*nTs);
asv2=cell(nAnimal,nEH,nTarg,nPr); % within/betweenPriors
msv=nan(nAnimal,nEH,nTarg,nPr,nTs);
% 3) curvature b/t Set vs IC (200ms): scatter, separately for short and long
%   - 2D @ set, 1D @ set
% only in meas-subspace as rotation is not clear in prod-subspace
cAcrTs=nan(nAnimal,nEH,nTarg,nPr,2,nTspp-2); % [2HG x 2EH x 2RL x 2SL x @set/@IC x 3data points]

% 1) IC vector size (ref to each prior mean) or their relative angle (more curved?): without finding out nearest neighbors
tsv=nan(nAnimal,nEH,nTarg,nPr*nTs-1,maxDim,maxNbin); % [2HG x 2EH x 2RL x 9 diffTs x time x dim]: for short800, tsv to long800
mtsv=nan(nAnimal,nEH,nTarg,nPr*nTs-2,maxNbin); % magnitude of tsv [2HG x 2EH x 2RL x 9 diffTs x time]: for short800, tsv to long800
atsv=nan(nAnimal,nEH,nTarg,nPr*nTs-1,maxNbin); % angle b/t nearby tsv [2HG x 2EH x 2RL x 9 diffTs x time]: for short800, tsv to long800
% 2) angle b/t overlap vector & time/ts vector: common time-warped axis?
aTsvOv=nan(nAnimal,nEH,nTarg,nPr*nTs-2,maxNbin);
aVvOv=nan(nAnimal,nEH,nTarg,nPr*nTs,maxNbin);% vv=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % velocity vector [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints x dim]
% 3) angle/magnitude of time vectors (ref: prior mean)
% use asv, msv
% 4) angle b/t common IC axis
aCommIC=nan(nAnimal,nEH,nTarg,maxNbin);

% available from analSpeed_ts: speed, |2nd derivative|, curvature, kinet
% + angle b/t velocity vectors of short vs long as a function of time: if <90, parallel
% + shortest path length between short and long: constant?
% + angle b/t tangent and shortest path vector: if <90, on the same manifold

%% cosyne
mtp=nan(nAnimal,nEH,nTarg,nPr,nTs); % speed [2HG x 2EH x 2RL x 2SL x 5ts]
sICa=nan(nAnimal,nEH,nTarg,nPr,nTs); % projection of states onto IC axis

xuAll=nan(nAnimal,nEH,nTarg,nPr,nTs);

if idNewTraj
    d=dir('trajKS_tp_*.mat'); % bin20_smth40>use new tp data set
else
    d=dir('trajKS_tp_*bin20_smth40.mat'); % bin20_smth40>use new tp data set
end

for iFile=1:length(d)
    %% for all conditions
    nm=d(iFile).name;
    load(nm); % D keep_neurons binSize optimD proj_matrix smthWidth
    %     size(D(10).data) 44
    disp(nm);
    
    if useOptimD % isinf(dim)
        load(nm,'optimD'); dim=optimD;
    end
    
    % finding animal, condition > use that projection
    if idNewTraj
        if ~isempty(strfind(d(iFile).name,'_H.mat'))
            idAnimal=1; iAnimalNm='H';
        else
            idAnimal=2; iAnimalNm='G';
        end
    else
        if ~isempty(strfind(d(iFile).name,'_H_'))
            idAnimal=1; iAnimalNm='H';
        else
            idAnimal=2; iAnimalNm='G';
        end
    end
    if ~isempty(strfind(d(iFile).name,'_ER_'))
        iEH=1; iTarg=1;
    elseif ~isempty(strfind(d(iFile).name,'_EL_'))
        iEH=1; iTarg=2;
    elseif ~isempty(strfind(d(iFile).name,'_HR_'))
        iEH=2; iTarg=1;
    elseif ~isempty(strfind(d(iFile).name,'_HL_'))
        iEH=2; iTarg=2;
    end
    tmpId=[iEH; iTarg];
    
    % load Xu
    load(['trajKS_prior_' nm(11:14) '_avgAttB4PCA.mat'],'Xu'); %{1x3pr} [1x17 1x21 1x17]
    for iPr=1:nPr
        xuAll(idAnimal,iEH,iTarg,iPr,:)=Xu{iPr}(linspace(1,length(Xu{iPr}),nTspp)); % short
    end
    
    epNm={'fixation';'targetOn';'periReady';'periSet';'production';'ts';'tp';'_priorS_';'_priorL_';'_priorSL_';'t_s';'t_p'};
    epNm2={'preFix','postFix';'preTOn','postTOn';'preReady','postReady';'','';'prod','postProd';'','';'','';'','';'','';'','';'','';'',''};
    for i=1:length(epNm)
        if ~isempty(strfind(nm,epNm{i})), break; end;
    end
    idEp=i;
    
    if strfind(nm,'periSet') % trajKS_periSet_ER_H
        idProd=1;
        idMeas=1;
        idAtt=0;
        
        curv=cell(nPr,1);
    elseif strfind(nm,'ts')
        idProd=0; % 1; % plotting together
        idMeas=1; % 0; % 1;
        idAtt=0;
        
        curv=cell(nPr+1,1);
    elseif strfind(nm,'tp')
        idProd=1; % 1; % plotting together
        idMeas=0; % 1; % 0; % 1;
        idAtt=0;
    elseif strfind(nm,'prior')
        idProd=0; % 1; % plotting together
        idMeas=1; % 0; % 1;
        idAtt=0; % 1; % 0; % 1;
        %     if strfind(nm,'priorSL')
        %         idAtt=0;
        %     else
        %     idAtt=1;
        %     end
    end
    
    
    idCondSpecific=true; % false;
    conditionSpecificPCANm={'condSpecPCA'};
    
    %% main
    % condId=1; % eyeRight
    % for iEH=1:nEH
    %     for iTh=1:nTarg
    %         if idCondSpecific
    tmpD=D;iEH=tmpId(1);iTh=tmpId(2);
    %             if iEH>2, iEH=2; end; % E2, H2 avgDir
    %         else
    %             tmpD=D(condId:4:length(D));
    %         end
    
    tmpMin=[]; % for setting measure/produce in the same range
    tmpMax=[];
    
    sv=nan(length(tmpD),dim); % time vector [10ts x optimD]
    
    if idProd % idMeas
        %% production only
        %         hFig=figure;setFigPos(1,1);% set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+1)]));hold all;
        %         hFig2=figure;setFigPos(1,2);
        %         hFig3=figure;setFigPos(1,3);
        
        %         figure;
        if idDebug
            figure; grid on; % debug
        end
        
        for i=1:length(tmpD) % 10 ts
            if length(tmpD(i).epochStarts)>1
                iSet=tmpD(i).epochStarts(2);
            else
                iSet=size(tmpD(i).data,2);
            end
            idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
            iTs=rem(i-1,5)+1;
            
            endPre=1;
            tmpX=tmpD(i).data(1,1:iSet);
            tmpY=tmpD(i).data(2,1:iSet);
            tmpZ=tmpD(i).data(3,1:iSet);
            % PC456
            %                 tmpX2=tmpD(i).data(4,1:iSet);
            %                 tmpY2=tmpD(i).data(5,1:iSet);
            %                 tmpZ2=tmpD(i).data(6,1:iSet);
            
            %                 plot3(tmpX,tmpY,tmpZ);ha;
            
            %             if i<=length(tmpD)/2 % bias+ % for split tp
            
            %% speed
            tmp=tmpD(i).data(1:dim,1:iSet); % [dim x time]
            if idDebug,            ha; plot3(tmp(1,:),tmp(2,:),tmp(3,:),'o-','color',tmpCmap{idPr}(iTs,:)); end; % debug \
            tmpV=sqrt(sum(diff(tmp,1,2).^2,1));       % [1 x (time-1)]
            nV=length(tmpV);
            
            v(idAnimal,iEH,iTh,idPr,iTs,1:nV)=tmpV; % [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
            
%             figure; histogram(tmpV(initT:end),15); title([iAnimalNm ' ' ehNm{iEH} ' ' targNm{iTh} ' ' prNm{idPr} ' ts#' num2str(iTs)]); waitforbuttonpress; close;
            
            vv(idAnimal,iEH,iTh,idPr,iTs,1:nV,1:dim)=diff(tmp,1,2)';
            
            % norm of 2nd derivative
            tmpV2=sqrt(sum(diff(tmp,2,2).^2,1));       % [1 x (time-2)]
            nV=length(tmpV2);
            v2(idAnimal,iEH,iTh,idPr,iTs,1:nV)=tmpV2;
            
            %% arc length
            % wholeMeas/overlap, short, long, SIL
            x(idAnimal,iEH,iTh,idPr,iTs,1)=sum(tmpV); % wholeMeas(longest ts/prior)/overlap(overlap ts/prior)
            
            % 2) check angle/magnitude b/t readout vectors across ts, prior: 2D cmap, if <90, parellel-not state dependent
            %    - kinet: nearest distance coalesce during set transient?
            vTmp=tmp(:,t1sv)-tmp(:,t0sv); % tmp(:,initT)-tmp(:,1);
            sv(i,:)=normalize(vTmp(:),1);
            msv(idAnimal,iEH,iTh,idPr,iTs)=sqrt(sum(vTmp.^2,1));
            
            %             tStartS=length(tmpD(1).data(1,:)); % shortest tp for ts=480
            %             tEndS=length(tmpD(5).data(1,:))-1; % shortest tp for ts=800 short
            %             tStartPr=length(tmpD((idPr-1)*nTs+1).data(1,:)); % 24 for short, 40 for long
            %             if i~=1 & i~=(nTs+1)
            %                 x(idAnimal,iEH,iTh,idPr,iTs,2)=sum(tmpV(tStartPr:end)); %
            %             end % if i~=1 & i~=(nTs+1)
            %             if idPr>1 % SIL
            %                 x(idAnimal,iEH,iTh,idPr,iTs,3)=sum(tmpV(tStartS:tEndS)); %
            %             end
            
            %% kinet (distance between nearest neighbors)
            % less noisy (distance is noisy)
            % not highly dependent on binSize
            % assuming same manifold
            
            % 1) angle b/t velocity vectors of short vs long as a function of time: if <90, parallel
            % 2) shortest path length between short and long: constant?
            % 3) angle b/t tangent and shortest path vector: if <90, on the same manifold
            
            %%%%% prior mean as refs
            if i>nTs % long prior
                refD=tmpD((refPr(2)-1)*nTs+refTs(2)).data(1:dim,:); % prior mean
            else % short
                refD=tmpD((refPr(1)-1)*nTs+refTs(1)).data(1:dim,:);
            end
            nT_ts=size(tmp,2);
            
            for iRefT=1:size(refD,2) % for each time points of ref
                tmpDist=tmp-repmat(refD(:,iRefT),1,nT_ts); % [dim x time(comparison ts)]
                [nearstS,nearstT]=min(sqrt(sum(tmpDist.^2,1)));
                tSp(idAnimal,iEH,iTh,idPr,iTs,iRefT)=nearstT;
                sSp(idAnimal,iEH,iTh,idPr,iTs,iRefT,1:dim)=tmp(1:dim,nearstT);
                spv(idAnimal,iEH,iTh,idPr,iTs,iRefT,1:dim)=tmp(1:dim,nearstT)-refD(:,iRefT);
                spl(idAnimal,iEH,iTh,idPr,iTs,iRefT)=nearstS;
                
                % 3) angle b/t tangent and shortest path vector: if <90, on the same manifold
                vvTmp=squeeze(vv(idAnimal,iEH,iTh,idPr,iTs,nearstT,1:dim));
                spvTmp=squeeze(spv(idAnimal,iEH,iTh,idPr,iTs,iRefT,1:dim));
                
                if prod(~isnan(vvTmp))
                    aSpvTv(idAnimal,iEH,iTh,idPr,iTs,1,iRefT)=acosd(dot(normalize(vvTmp,1),normalize(spvTmp,1)));
                end
                try
                    vvRef=refD(:,iRefT+1)-refD(:,iRefT);
                    aSpvTv(idAnimal,iEH,iTh,idPr,iTs,2,iRefT)=acosd(dot(normalize(vvRef,1),normalize(spvTmp,1)));
                catch
                end
                
                if nearstT==nT_ts % after reaching final state for ts, meaningless
                    break;
                end
            end % for iRefT=1:size(refD,2) % for each time points of ref
            % speed, distance %%%%%
            tmpTsp=tSp(idAnimal,iEH,iTh,idPr,iTs,:); tmpTsp=tmpTsp(~isnan(tmpTsp));
            tmpTsp=(tmpTsp(initT:end)-initT)*binSize; % speed after IC
            sp(idAnimal,iEH,iTh,idPr,iTs)=squeeze(tmpTsp)\(binSize*[0:(length(tmpTsp)-1)]');
            %             if idAnimal==2 & iEH==1 & iTh==1 & idPr==1 & iTs==5 % debug
            %                 disp([]);
            %             end
            
            % or S=regstats((binSize*[0:(length(tmpTsp)-1)]'),squeeze(tmpTsp),[1],{'beta','tstat','fstat','yhat'})
            if iTs>refTs(idPr) % flip sign for ts longer than prior mean
                dist(idAnimal,iEH,iTh,idPr,iTs)=-spl(idAnimal,iEH,iTh,idPr,iTs,initT); % state(IC)
            else
                dist(idAnimal,iEH,iTh,idPr,iTs)=spl(idAnimal,iEH,iTh,idPr,iTs,initT); % state(IC)
            end
            
            % 1) angle b/t velocity vectors of short vs long as a function of time: if <90, parallel
            for jTmpD=1:length(tmpD) % for all other ts
                
                jdPr=floor((jTmpD-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
                jTs=rem(jTmpD-1,5)+1;
                
                for iRefT=1:size(refD,2) % for each time points of ref
                    t1=tSp(idAnimal,iEH,iTh,idPr,iTs,iRefT);
                    t2=tSp(idAnimal,iEH,iTh,jdPr,jTs,iRefT);
                    if ~isnan(t1) & ~isnan(t2)
                        vvTmp1=squeeze(vv(idAnimal,iEH,iTh,idPr,iTs,t1,1:dim));
                        vvTmp2=squeeze(vv(idAnimal,iEH,iTh,jdPr,jTs,t2,1:dim));
                        if ~isnan(vvTmp1) & ~isnan(vvTmp2)
                            tmpAngle=acosd(dot(normalize(vvTmp1,1),normalize(vvTmp2,1)));
                            av(idAnimal,iEH,iTh,idPr,iTs,iRefT,jdPr,jTs)=tmpAngle;
                            
                            if jTmpD~=i % only for different data
                                if idPr==jdPr % same prior
                                    av2{idAnimal,iEH,iTh,1,iRefT}=[av2{idAnimal,iEH,iTh,1,iRefT}; tmpAngle]; % [2HG x 2EH x 2RL x within/betweenPriors x timePoints]
                                else
                                    av2{idAnimal,iEH,iTh,2,iRefT}=[av2{idAnimal,iEH,iTh,2,iRefT}; tmpAngle]; % [2HG x 2EH x 2RL x within/betweenPriors x timePoints]
                                end
                            end % vvTmp1
                        end % if ~isnan(vvTmp1) & ~isnan(vvTmp2)
                    end % if ~isnan(t1) & ~isnan(t2)
                end % for iRefT=1:size(refD,2) % for each time points of ref
            end
            
            %% curvature ~= angle b/t successive vectors?
            tmpForC=tmpD(i).data(1:3,1:iSet); % [dim x time]
            tmpCurv=LineCurvature2D(tmpForC');
            tmpCurv=tmpCurv(2:(length(tmpCurv)-1));
            
            c(idAnimal,iEH,iTh,idPr,iTs,1:length(tmpCurv))=tmpCurv;
            
        end % length(tmpD)
        
        %% tp-specific
        % prep data: truncate [set IC]
        tmp=nan(nPr*nTs,dim,initT); % [10ts x dim x time]
        for i=1:length(tmpD)
            %                 idPr=floor((i-1)/nTs)+1; % 1 1 1 1 1 2 ...
            %                 iTs=rem(i-1,nTs)+1;
            tmpN=min([size(tmpD(i).data,2) initT]);
            tmp(i,:,1:tmpN)=tmpD(i).data(1:dim,1:tmpN);
        end %             for i=1:length(tmpD)
        
        % finding  tsv & sv
        tsv2=diff(tmp(:,:,1),1,1); % [9 diffTs x dim x tim] %             for short800, tsv to long800
        tsv2=cat(1,tsv2,tsv2(end,:,:)); % for 1200, tsv2 from 1100
        tsv2=squeeze(tsv2); % [10ts x dim]
        
        sv2=squeeze(diff(tmp(:,:,[1 end]),1,3)); % [10ts x dim]
        for i=1:length(tmpD)
            idPr=floor((i-1)/nTs)+1; % 1 1 1 1 1 2 ...
            iTs=rem(i-1,nTs)+1;
            aTsvSv2(idAnimal,iEH,iTh,idPr,iTs)=acosd(dot(normalize(tsv2(i,1:dim)',1),normalize(sv2(i,1:dim)',1)));
        end
        
        % 2) sICa=nan(nAnimal,nEH,nTarg,nPr,nTs); % projection of states onto IC axis
        sSetTmp=squeeze(tmp(:,:,end)); % [10ts x dim] initT tInitCond=200;%400; % 200;
        if isSepPr==0
            axIC=normalize(mean(diff(sSetTmp,1,1),1)',1); % [dim x 1] % for now across prior
        end
        for iPr=1:nPr
            idTsTmp=(iPr-1)*nTspp+[1:nTspp];
            if isSepPr
                axIC=normalize(mean(diff(sSetTmp(idTsTmp,:),1,1),1)',1); % [dim x 1] % for now across prior
            end
            sICa(idAnimal,iEH,iTh,iPr,:)=sSetTmp(idTsTmp,:)*axIC;
            %                      disp(squeeze(sICa(idAnimal,iEH,iTh,iPr,:)));
        end
        
        
    end %  if idMeas
end % for i=1:length(d)

%              if idCosyne
% 1) correlation b/t tp (avgAcrSess) and speed (avgAcrTimeAfterIC) across all data sets
% 2) 1D projection on IC axis vs speed

% v=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % speed [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
% mtp=nan(nAnimal,nEH,nTarg,nPr,nTs); % speed [2HG x 2EH x 2RL x 2SL x 5ts]
for idAnimal=1:nAnimal
    % 1) tp (avgAcrSess)
    beh=load([animalNm{idAnimal} '_RSGprior_DMFC.mat']); % T idOut t idShortTrial theta idHandEye
    for iEH=1:nEH
        for iTarg=1:nTarg
            for iPr=1:nPr
                for iTs=1:nTspp
                    idTp=~beh.idOut & beh.idHandEye==(2-iEH) & beh.theta==(180*(iTarg-1))...
                        & beh.idShortTrial==(2-iPr) & beh.T==T{iPr}(iTs);
                    mtp(idAnimal,iEH,iTarg,iPr,iTs)=mean(beh.t(idTp));
                    %                                  % plot
                    %                                  figure; hist(beh.t(idTp),100);
                    %                                  title([animalNm{idAnimal} ' ' ehNm{iEH} ' ' targNm{iTarg} ' ' prNm{iPr} num2str(T{iPr}(iTs))]);
                    %                                  waitforbuttonpress;close;
                end
            end
        end
    end
end
%              end % if idCosyne

%% plot
% v=nan(nAnimal,nEH,nTarg,nPr,nTs,round(T{end}(end)/binSize)); % [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
% nEpoch=3; % wholeMeas/overlap, short/long, SIL
% x=nan(nAnimal,nEH,nTarg,nPr,nTs,nEpoch); % [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]

% if idCosyne

h1=figure; ha;setFigPos(1,1); % 1); % speed vs time (avg. across ts) two lines for short vs long
h2=figure; ha;setFigPos(1,2); % 2); % avv
h3=figure; ha;setFigPos(1,3);
% figure(3); setFigPos(1,3); % avv colormap for each data set: short
% figure(4); setFigPos(2,3); % smae for long

condNm={'ER','EL','HR','HL'};
nCmap=length(animalNm)*length(condNm);
cmap3=parula(nCmap); % hsv(nCmap);
icmap=1;

X=[];Y=[];
X2=[];Y2=[];
X3=[];Y3=[];


for i=1:nAnimal
    for j=1:nEH
        for k=1:nTarg
            iCond= (j-1)*nTarg+k;
            
            mv=squeeze(nanmean(v(i,j,k,:,:,initT:end),6))'; % avgAcrTimeAfterIC [5ts x 2pr]
            %% fig3g: tp vs speed

            %1) correlation b/t tp (avgAcrSess) and speed (avgAcrTimeAfterIC) across all data sets
            figure(h1); 
            tmpX=mv(:); % [10 ts x 1]
            tmpY=[squeeze(mtp(i,j,k,1,:)); squeeze(mtp(i,j,k,2,:))];
            %                 plot(tmpX,tmpY,'-','color',cmapHG(i,:),'linewidth',lw2); 
            % for each data set
            plot(tmpX(1:nTs),tmpY(1:nTs),'-','color',cmapPr2(1,:),'linewidth',lw2); %short
            plot(tmpX(nTs+[1:nTs]),tmpY(nTs+[1:nTs]),'-','color',cmapPr2(2,:),'linewidth',lw2); %long
            X=[X; tmpX(:)']; % [8dataSet x 10ts]
            Y=[Y; tmpY(:)']; % [8dataSet x 10ts]
            if icmap==nCmap
                % avg across data sets
                for iPr=1:nPr % lines
                    idTsTmp=(iPr-1)*nTspp+[1:nTspp];
                    plot(nanmean(X(:,idTsTmp),1),...
                        nanmean(Y(:,idTsTmp),1),'-','color',cmapPr(iPr,:),'linewidth',lw); %long
                end
                for l=1:size(X,2) % for each ts
                    idPr=floor((l-1)/nTspp)+1; %  1     1     1     1     1     2     2     2     2     2
                    iTs=rem(l-1,5)+1;
                    tmpX2=nanmean(X(:,l),1);tmpXS=sem(X(:,l),1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    tmpY2=nanmean(Y(:,l),1);tmpYS=sem(Y(:,l),1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    errorbarXY_meanSem(tmpX2,tmpY2,tmpXS,tmpYS,...
                        tmpCmap{idPr,1}(iTs,:),'o',msize,'w',tmpCmap{idPr,1}(iTs,:)); % tmpCmap{idPr,1}(iTs,:),'k');
                end
                axis tight;
                set(gca,'tickdir','out','xtick',0.1:0.1:0.4,'xticklabel',{'0.1',[],'0.3',[]},'xlim',[0.08 0.4],... % ,'xlim',[T{1}(1)-binSize T{end}(end)+binSize],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)}
                    'ylim',[450 1250],'ytick',[T{1}(1);T{1}(3);T{1}(end);T{2}(3);T{2}(end)],'yticklabel',[T{1}(1);[];T{1}(end);[];T{2}(end)]./1000);
            end %  if icmap==nCmap
            
            %% fig3f: speed vs IC
            % 2) sICa=nan(nAnimal,nEH,nTarg,nPr,nTs); % projection of states onto IC axis
            figure(h2); 
            tmpY=mv(:); % [10 ts x 1]
            tmpX=[squeeze(sICa(i,j,k,1,:)); squeeze(sICa(i,j,k,2,:))];
            %             plot(tmpX,tmpY,'-','color',cmapHG(i,:),'linewidth',lw2); % for each data set
            plot(tmpX(1:nTs),tmpY(1:nTs),'-','color',cmapPr2(1,:),'linewidth',lw2); %short
            plot(tmpX(nTs+[1:nTs]),tmpY(nTs+[1:nTs]),'-','color',cmapPr2(2,:),'linewidth',lw2); %long
            X2=[X2; tmpX(:)']; % [8dataSet x 10ts]
            Y2=[Y2; tmpY(:)']; % [8dataSet x 10ts]
            if icmap==nCmap
                % avg across data sets
                for iPr=1:nPr
                     idTsTmp=(iPr-1)*nTspp+[1:nTspp];
                    plot(nanmean(X2(:,idTsTmp),1),...
                        nanmean(Y2(:,idTsTmp),1),'-','color',cmapPr(iPr,:),'linewidth',lw); %long
                end
                for l=1:size(X2,2) % for each ts
                    idPr=floor((l-1)/nTspp)+1; %  1     1     1     1     1     2     2     2     2     2
                    iTs=rem(l-1,5)+1;
                    tmpX2=nanmean(X2(:,l),1);tmpXS=sem(X2(:,l),1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    tmpY2=nanmean(Y2(:,l),1);tmpYS=sem(Y2(:,l),1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    errorbarXY_meanSem(tmpX2,tmpY2,tmpXS,tmpYS,...
                        tmpCmap{idPr,1}(iTs,:),'o',msize,'w',tmpCmap{idPr,1}(iTs,:)); % ,tmpCmap{idPr,1}(iTs,:),'k');
                end
                axis tight;
                set(gca,'tickdir','out','ytick',0.1:0.1:0.4,'yticklabel',{'0.1',[],'0.3',[]},'ylim',[0.08 0.4],... %); % ,'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xlim',[T{1}(1)-binSize T{end}(end)+binSize],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                    'xtick',-2:1);
                %                     'ylim',[0 0.2],'ytick',0:0.1:0.2);
            end %  if icmap==nCmap
            
            %% fig3d: speed vs Xu
            figure(h3); 
            tmpY=mv(:); % [10 ts x 1]
            tmpX=[squeeze(xuAll(i,j,k,1,:)); squeeze(xuAll(i,j,k,2,:))];
            %             plot(tmpX,tmpY,'-','color',cmapHG(i,:),'linewidth',lw2); % for each data set
            plot(tmpX(1:nTs),tmpY(1:nTs),'-','color',cmapPr2(1,:),'linewidth',lw2); %short
            plot(tmpX(nTs+[1:nTs]),tmpY(nTs+[1:nTs]),'-','color',cmapPr2(2,:),'linewidth',lw2); %long
            X3=[X3; tmpX(:)']; % [8dataSet x 10ts]
            Y3=[Y3; tmpY(:)']; % [8dataSet x 10ts]
            if icmap==nCmap
                % avg across data sets
                for iPr=1:nPr
                     idTsTmp=(iPr-1)*nTspp+[1:nTspp];
                    plot(nanmean(X3(:,idTsTmp),1),...
                        nanmean(Y3(:,idTsTmp),1),'-','color',cmapPr(iPr,:),'linewidth',lw); %long
                end
                for l=1:size(X3,2) % for each ts
                    idPr=floor((l-1)/nTspp)+1; %  1     1     1     1     1     2     2     2     2     2
                    iTs=rem(l-1,5)+1;
                    tmpX2=nanmean(X3(:,l),1);tmpXS=sem(X3(:,l),1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    tmpY3=nanmean(Y3(:,l),1);tmpYS=sem(Y3(:,l),1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    errorbarXY_meanSem(tmpX2,tmpY3,tmpXS,tmpYS,...
                        tmpCmap{idPr,1}(iTs,:),'o',msize,'w',tmpCmap{idPr,1}(iTs,:)); % ,tmpCmap{idPr,1}(iTs,:),'k');
                end
                axis tight;
                set(gca,'tickdir','out','ytick',0.1:0.1:0.4,'yticklabel',{'0.1',[],'0.3',[]},'ylim',[0.08 0.4],... %); % ,'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xlim',[T{1}(1)-binSize T{end}(end)+binSize],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                    'xtick',-2:1);
                %                     'ylim',[0 0.2],'ytick',0:0.1:0.2);
            end %  if icmap==nCmap
            
            icmap=icmap+1;
            
        end % k targ
    end % j EH
end % i animal

Xs=X(:,1:nTs);
Xl=X(:,(nTs+1):end);
Ys=Y(:,1:nTs);
Yl=Y(:,(nTs+1):end);
[r,p]=corr(Xs(:),Ys(:))
[r,p]=corr(Xl(:),Yl(:))

% X2s=X2(:,1:nTs);
% X2l=X2(:,(nTs+1):end);
% Y2s=Y2(:,1:nTs);
% Y2l=Y2(:,(nTs+1):end);
% [r,p]=corr(X2s(:),Y2s(:))
% [r,p]=corr(X2l(:),Y2l(:))

X3s=X3(:,1:nTs);
X3l=X3(:,(nTs+1):end);
Y3s=Y3(:,1:nTs);
Y3l=Y3(:,(nTs+1):end);
[r,p]=corr(X3s(:),Y3s(:))
[r,p]=corr(X3l(:),Y3l(:))

% else % cosyne
%
% close all;


% [r,p]=corr(X(:),Y(:))
% 
% r =
% 
%    -0.7390
% 
% 
% p =
% 
%    5.0253e-15
% 
% [r,p]=corr(X2(:),Y2(:))
% 
% r =
% 
%    -0.7812
% 
% 
% p =
% 
%    1.2615e-17


% 

return;

%%
idPlotEx=1;
if idPlotEx
    ii=2; % G
    jj=2; % hand
    kk=1; % right
    ll= 1; % short prior
end

% figure(1); setFigPos(1,1); % speed vs time
% figure(2); setFigPos(1,2); % arc length vs all ts (whole meas)
% figure(3); setFigPos(1,6); % 3); % speed from kinet
% % curvature
% figure(4); setFigPos(1,4); % norm of 2nd derivagtive
% figure(5); setFigPos(1,5); % curvature
% % kinet
% figure(6); setFigPos(2,1); % 1) angle b/t velocity vectors of short vs long as a function of time: if <90, parallel
% figure(7); setFigPos(2,6); % 1,4); % 2,2); % 2) shortest path length between short and long: constant?  - kinet: nearest distance coalesce during set transient?: spl=nan(nAnimal,nEH,nTarg,nPr,nTs,maxNbin); % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints]
% figure(8); setFigPos(2,3); % 3) angle b/t tangent and shortest path vector: if <90, on the same manifold
% % tp-specific
% figure(9); setFigPos(2,4); % %    - angle between ts vector & time vector: orthogonal?  % aTsvSv2=nan(nAnimal,nEH,nTarg,nPr,nTs); % prod-space % [2HG x 2EH x 2RL x 2SL x 5ts x timePoints]; for short800, tsv to long800, for 1200, tsv from 1100
% figure(10); setFigPos(2,5); % aTsvSv1=nan(nAnimal,nEH,nTarg,nPr,nTs); % meas-space
% figure(11); setFigPos(1,1); % asv=nan(nAnimal,nEH,nTarg,nPr,nTs,nPr,nTs); angle b/t readout vectors across ts, prior: 2D cmap, if <90, parellel-not state dependent
% % asv2=cell(nAnimal,nEH,nTarg,nPr); % within/betweenPriors
% figure(12); setFigPos(1,2); % msv=nan(nAnimal,nEH,nTarg,nPr,nTs);
% figure(13); setFigPos(1,3); % % 3) curvature b/t Set vs IC (200ms): scatter, separately for short and long % cAcrTs=nan(nAnimal,nEH,nTarg,nPr,nTs,2); % [2HG x 2EH x 2RL x 2SL x 5ts x @set/@IC]
%
% figure(14); setFigPos(1,1);% mtsv [8ts x time(initT)]
% figure(19); setFigPos(1,3);% mtsv [8ts x time(initT)]
% figure(15); setFigPos(1,2); % atsv [9ts x time(initT)]
% figure(16); setFigPos(1,3);% aTsvOv=nan(nAnimal,nEH,nTarg,nPr*nTs-2,maxNbin);
% figure(17); setFigPos(2,1);% aVvOv=nan(nAnimal,nEH,nTarg,nPr*nTs,maxNbin);
% figure(18); setFigPos(2,2); % asv=nan(nAnimal,nEH,nTarg,nPr*nTs,nPr*nTs);
% figure(20); setFigPos(2,3); % aCommIC=nan(nAnimal,nEH,nTarg,maxNbin);

figure(1000); setFigPos(2,4); % 1/speed = f(ts) with BLS
teSp=nan(nAnimal*nEH*nTarg,nPr,nTs); % 1/speed fit to BLS to avg
figure(1001); setFigPos(2,5); % distance = f(ts) with BLS
teDist=nan(nAnimal*nEH*nTarg,nPr,nTs); % dist fit to BLS

condNm={'ER','EL','HR','HL'};
nCmap=length(animalNm)*length(condNm);
cmap4=parula(nCmap); % hsv(nCmap);
icmap=1;

tmpMS=[];tmpML=[];

for i=1:nAnimal
    for j=1:nEH
        for k=1:nTarg
            iCond= (j-1)*nTarg+k;
            %%%%%%
            % BLS
            wm=0.05; %%%%% 0.048 for G, 0.05 for H
            
            figure(1000);  ha; % 1/speed = f(ts) with BLS
            for iPr=1:nPr
                if i==1 & j==1 & k==1
                    tmpTs=min(T{iPr}):max(T{iPr});
                    teBLS=BLS(tmpTs,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                    if idTpMTs
                        plot(tmpTs,teBLS-tmpTs,'-','linewidth',lw,'color',tmpCmap{iPr,1}(1,:));
                        if iPr==nPr, plotHorizon(gca,0,[]); end; %  axis tight;                        
                    else
                        plot(tmpTs,teBLS,'-','linewidth',lw,'color',tmpCmap{iPr,1}(1,:));
                        if iPr==nPr, plotIdentity(gca); end; %  axis tight;
                    end
                end
                teSpTmp=1./sp(i,j,k,iPr,:); % [5 x1]
                stats=regstats(BLS(T{iPr},wm,[min(T{iPr}) max(T{iPr})],'uniform')',...
                    teSpTmp(:),'linear',{'beta','rsquare','yhat'});
                teSp(icmap,iPr,:)=stats.yhat(:); % to avg
%                 disp(['R^2(speed): ' num2str(stats.rsquare,2)]);
                
%                 for iTs=1:nTs
%                     if idTpMTs
%                         plot(T{iPr}(iTs),stats.yhat(iTs)-T{iPr}(iTs),'.','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw2,'markersize',msize); %cmap4(icmap,:),'markerfacecolor','w');
%                     else
%                         plot(T{iPr}(iTs),stats.yhat(iTs),'.','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw2,'markersize',msize); %cmap4(icmap,:),'markerfacecolor','w');
%                     end
%                 end
            end % for iPr=1:nPr
            if icmap==nCmap % avg across data sets
                for iPr=1:nPr
                    mTeSp=squeeze(nanmean(teSp(:,iPr,:),1)); % ,'omitnan'));
                    sTeSp=squeeze(sem(teSp(:,iPr,:),1,'omitnan'));
                    if idTpMTs
                        h=shadedErrorBar(T{iPr}(:),mTeSp-T{iPr}(:),sTeSp,{'.','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize},1); drawnow; % 2
                         c=corr(mTeSp-T{iPr}(:),BLS(T{iPr}(:),0.05,[min(T{iPr}(:)) max(T{iPr}(:))],'uniform')-T{iPr}(:));
                        disp(['R^2(BLS-ts): ' num2str(c.^2,2)]);
                    else
                        h=shadedErrorBar(T{iPr}(:),mTeSp,sTeSp,{'.','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize},1); drawnow; % 2
                        c=corr(mTeSp(:),BLS(T{iPr}(:),0.05,[min(T{iPr}(:)) max(T{iPr}(:))],'uniform'));
                        disp(['R^2(BLS): ' num2str(c.^2,2)]);
                    end
                    %                     plot(T{iPr}(:),mTeSp,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw);                     % ,'markersize',msize); % line
                    for iTs=1:nTs % ts-speciifc colored dot
                        if idTpMTs
                            plot(T{iPr}(iTs),mTeSp(iTs)-T{iPr}(iTs),'o','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw,'markersize',msize2); % ,'markerfacecolor','w'); % ,'markersize',msize);
                        else
                            plot(T{iPr}(iTs),mTeSp(iTs),'o','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw,'markersize',msize2); % ,'markerfacecolor','w'); % ,'markersize',msize);
                        end
                    end
                end % for iPr=1:nPr
                if idTpMTs
                    set(gca,'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
                    xlabel('t_s (s)'); ylabel('transformed 1/speed - t_s (s)');
                else
                    set(gca,'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                        'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
                    xlabel('t_s (s)'); ylabel('transformed 1/speed (s)');
                end
                %                 'ytick',0.1:0.2:0.3,'ylim',[0.08 0.4],... %); % ,'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xlim',[T{1}(1)-binSize T{end}(end)+binSize],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                %                     'xtick',-2:1);
                %                     'ylim',[0 0.2],'ytick',0:0.1:0.2);
                
            end %  if icmap==nCmap
            
            
            figure(1001); ha; % distance = f(ts) with BLS
            for iPr=1:nPr
                if i==1 & j==1 & k==1
                    tmpTs=min(T{iPr}):max(T{iPr});
                    teBLS=BLS(tmpTs,wm,[min(T{iPr}) max(T{iPr})],'uniform');
                    if idTpMTs
                        plot(tmpTs,teBLS-tmpTs,'-','linewidth',lw,'color',tmpCmap{iPr,1}(1,:));
                        if iPr==nPr,  plotHorizon(gca,0,[]); end; %  axis tight;
                    else
                        plot(tmpTs,teBLS,'-','linewidth',lw,'color',tmpCmap{iPr,1}(1,:));
                        if iPr==nPr,  plotIdentity(gca); end; %  axis tight;
                    end
                    
                end
                teDistTmp=dist(i,j,k,iPr,:); % [5 x1]
                stats=regstats(BLS(T{iPr},wm,[min(T{iPr}) max(T{iPr})],'uniform')',...
                    teDistTmp(:),'linear',{'beta','rsquare','yhat'});
%                 disp(['R^2(distance): ' num2str(stats.rsquare,2)]);
                teDist(icmap,iPr,:)=stats.yhat(:); % to avg
%                 for iTs=1:nTs
%                     if idTpMTs
%                         plot(T{iPr}(iTs),stats.yhat(iTs)-T{iPr}(iTs),'.','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw2,'markersize',msize); %  cmap4(icmap,:) ,'markerfacecolor','w');
%                     else
%                         plot(T{iPr}(iTs),stats.yhat(iTs),'.','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw2,'markersize',msize); %  cmap4(icmap,:) ,'markerfacecolor','w');
%                     end
%                 end
            end % for iPr=1:nPr
            if icmap==nCmap % avg across data sets
                for iPr=1:nPr
                    mTeDist=squeeze(nanmean(teDist(:,iPr,:),1)); % ,'omitnan'));
                    sTeDist=squeeze(sem(teDist(:,iPr,:),1,'omitnan'));
                    if idTpMTs
                        h=shadedErrorBar(T{iPr}(:),mTeDist-T{iPr}(:),sTeDist,{'.','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize},1); drawnow; % 2
                        c=corr(mTeSp(:)-T{iPr}(:),BLS(T{iPr}(:),0.05,[min(T{iPr}(:)) max(T{iPr}(:))],'uniform')-T{iPr}(:));
                        disp(['R^2(BLS-ts): ' num2str(c.^2,2)]);
                    else
                        h=shadedErrorBar(T{iPr}(:),mTeDist,sTeDist,{'.','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize},1); drawnow; % 2
                        c=corr(mTeSp(:),BLS(T{iPr}(:),0.05,[min(T{iPr}(:)) max(T{iPr}(:))],'uniform'));
                        disp(['R^2(BLS): ' num2str(c.^2,2)]);
                    end
                    %                     plot(T{iPr}(:),mTeDist,'-','color',tmpCmap{iPr,1}(1,:),'linewidth',lw);                     % ,'markersize',msize); % line
                    for iTs=1:nTs
                        if idTpMTs
                            plot(T{iPr}(iTs),mTeDist(iTs)-T{iPr}(iTs),'o','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw,'markersize',msize2); % ,'markerfacecolor','w'); % ,'markersize',msize);
                        else
                            plot(T{iPr}(iTs),mTeDist(iTs),'o','color',tmpCmap{iPr,1}(iTs,:),'linewidth',lw,'markersize',msize2); % ,'markerfacecolor','w'); % ,'markersize',msize);
                        end
                    end
                end % for iPr=1:nPr
                if idTpMTs
                    set(gca,'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
                    xlabel('t_s (s)'); ylabel('transformed S(IC) - t_s (s)');
                else
                set(gca,'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                    'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
                xlabel('t_s (s)'); ylabel('transformed S(IC) (s)');
                end
                %                 'ytick',0.1:0.2:0.3,'ylim',[0.08 0.4],... %); % ,'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xlim',[T{1}(1)-binSize T{end}(end)+binSize],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                %                     'xtick',-2:1);
                %                     'ylim',[0 0.2],'ytick',0:0.1:0.2);
                
            end %  if icmap==nCmap
            
%             figure(3);
            if idPlotEx
                ha;
                if i==ii & j==jj & k==kk
                    for ll=1:nPr
                        figure; setFigPos(1,ll);ha;
                        for m=1:nTspp
                            tmpM=squeeze(tSp(i,j,k,ll,m,:))*binSize/1000; % for each ts
                            tmpT=([1:length(tmpM)]-1)*binSize/1000;
                            plot(tmpT,tmpM(:)','-','color',tmpCmap{ll,1}(m,:),'linewidth',lw);
                        end
                    end
                    if ll==1
                        set(gca,'tickdir','out','xtick',0:0.1:0.4,... % [T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                            'ytick',0:0.1:0.6);
                        %                         'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000})
                    else
                        set(gca,'tickdir','out','xtick',0:0.2:0.8,... % [T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                            'ytick',0:0.2:0.8);
                    end
                    xlabel('t_s in ref.(s)'); ylabel('time of nearest state (s)');
                    drawnow;
                end
            else
                subplot(nEH*nTarg,nAnimal,iCond*nAnimal-(nAnimal-i)); ha;
                title([animalNm{i} ':' ehNm{j} targNm{k}]);
                for l=1:nPr
                    for m=1:nTspp
                        tmpM=squeeze(tSp(i,j,k,l,m,:))*binSize; % for each ts
                        tmpT=([1:length(tmpM)]-1)*binSize;
                        plot(tmpT,tmpM(:)','-','color',tmpCmap{l,1}(m,:),'linewidth',lw);
                    end
                end
                axis tight;
                if iCond~=4
                    set(gca,'xticklabel',[],'ticklength',[0.03 0.03]);
                else
                    set(gca,'ticklength',[0.03 0.03]);
                end
            end
            
            
            %% kinet
            cmap2=[0 0 0;.5 .5 .5];
            
            
%             figure(7); % 2) shortest path length between short and long: constant?
            if idPlotEx
                ha;
                if i==ii & j==jj & k==kk
                    for ll=1:nPr
                        figure; setFigPos(2,ll);ha;
                        for m=1:nTspp
                            if m>refTs(ll) % flip sign for ts longer than prior mean
                                tmpM=squeeze(spl(i,j,k,ll,m,:));
                            else
                                tmpM=-squeeze(spl(i,j,k,ll,m,:));
                            end
                            tmpT=([1:length(tmpM)]-1)*binSize/1000;
                            plot(tmpT,tmpM(:)','-','color',tmpCmap{ll,1}(m,:),'linewidth',lw);
                        end % m ts
                        if ll==1
                            set(gca,'tickdir','out','xtick',0:0.1:0.4,...[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                                'ytick',-1:0.5:1);
                        else
                            set(gca,'tickdir','out','xtick',0:0.2:0.8,...[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
                                'ytick',-1:0.5:1);
                        end
                        %             'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000})
                        xlabel('time after Set (s)'); ylabel('distance to nearest state');
                    end % ll pr
                end
            else
                subplot(nEH*nTarg,nAnimal,iCond*nAnimal-(nAnimal-i)); ha;
                title([animalNm{i} ':' ehNm{j} targNm{k}]);
                for l=1:nPr
                    for m=1:nTspp
                        %                     if l==refPr % >1
                        tmpM=squeeze(spl(i,j,k,l,m,:)); % avg across ts
                        %                     else % flip sign for non-ref prior
                        %                         tmpM=-squeeze(spl(i,j,k,l,m,:)); % avg across ts
                        %                     end
                        tmpT=([1:length(tmpM)]-1)*binSize;
                        plot(tmpT,tmpM(:)','-','color',tmpCmap{l,1}(m,:),'linewidth',lw);
                    end
                end
                axis tight;
                if iCond~=4
                    set(gca,'xticklabel',[],'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'ticklength',[0.03 0.03]);
                else
                    set(gca,'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'ticklength',[0.03 0.03]);
                end
            end
            
            if iCond==3 & i==1
                % %                 figure(1);                ylabel('speed (a.u.)');
                % %                 figure(4);                ylabel('norm of 2nd derivative (a.u.)');
                % %                 figure(5);                ylabel('curvature (a.u.)');
                % %                 figure(2);                ylabel('travelling distance from ready (a.u.)');
                %                 figure(3);                ylabel('time to nearest states');
                %
                % %                 figure(6); ylabel('angle b/t time vectors');
                %                 figure(7); ylabel('nearest distance');
                % %                 figure(8); ylabel('angle b/t time and nearest vectors');
                % %
                % %                 figure(9); ylabel('angle b/t ts and time vectors in prodSpace');
                % %                 figure(10); ylabel('angle b/t ts and time vectors in measSpace');
                % %                 figure(11); ylabel('angle b/t time vectors');
                % %                 figure(12); ylabel(['|time vectors [' num2str(t0sv*binSize) ',' num2str(t1sv*binSize)  ']|']);
                % %
                % %                 figure(14); ylabel('distance from prior mean trajectory'); %%%%%
                % %                 figure(15); ylabel('angle b/t successive across-ts vectors');
                % %                 figure(16); ylabel('angle b/t across-ts vectors vs overlap vector');
                % %                 figure(17); ylabel(['angle b/t time vectors [' num2str(t0sv*binSize) ',' num2str(t1sv*binSize)  '] vs overlap vector']);
                % %                 figure(18); ylabel(['angle of time vectors [' num2str(t0sv*binSize) ',' num2str(t1sv*binSize)  ']']);
                % %                 figure(19); ylabel('distance from prior mean trajectory');
                % %                 figure(20); ylabel('angle b/t avg. ts vectors of short vs long');
            end
            
            
            icmap=icmap+1;
            
        end % k targ
    end % j EH
end % i animal

