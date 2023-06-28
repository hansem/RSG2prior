function analSpeed_ts_supportOnly_fig3b

% 2018/7/9
% curvature: using stencil methods, rather than diff (estCurvature.m in psthDataHigh)
% picking up peaks across [short in long]+[long], rather than [long],
% rather than using fitGauss
% use trajKS_prior_*_stitchSiL.mat

% 2018/6/29
% peak time as a measure of curvature (after separate PCA 'idSepPC')
% check PVAF S vs L vs SiL (across dataSets, PCs): idPVAF=1
% correct centering (using meanPSTH in trajKS_XXX.mat)

% 2018/6/25
% scatter plot for curvature
% use seconds as the unit (unit of PC: spk/s): isUseMs=0
% PVAF estimated same for S/L/SL: idPVAF=0

% 2018/6/15
% (TBD) idAvgAtt: 1 if after curvature estimation; 0 if before

% 2018/6/11
% isUseMs: unit for x axis

% 2018/6/7
% fig3b1: example PC time course (idDataFig)
% fig3b2: curvature vs % var explained
%   commenting other figures

% summary of parameters
% 1) dimension: use optimal (PVAF>75%) - useOptimD
% 2) use median for summary: curvature distribution across time for each data set & PCs is skewed - idMed
% 3) PCA done for each data set(an/HE/RL)>PVAF & curv measured for S/L/SL
% 4) averaging across ts with attrition before estimating curvature
% 5) remove outlier: G hand Left PC4 short (thOut)

% 2018/4/13
% curvature =f(varExplained), for each PC, colored by S/L/SL
% - %VarExplained calculated for each data set & S/L/SL
% - curvature: fit quadratic function (idUseQuadFit)

% 2018/4/12
% new analysis to show no rotation in shortInLong: angle b/t readout vector, avrv
% (shortest to longest ts) and difference vector
% dAngle to use low D for old angle estimation (idIncShortInLong)

% 2018/4/4
% paper figure 3i: angle with SL too
% update with correct tInfo
% condition-specific; idIncShortInLong=1>trajKS_prior_XXX_bin20_smth40.mat vs _supportOnly
% as 60ms truncated, no use of idExcB4pr(=0)
% use whole meas epoch if idIncShortInLong=2 (TBD: tBuffer after set not done)

% 2018/2/22 COSYNE18
% speed profile & angle
% from PCA applied to only prior support
% N.B. averaging across ts is done after estimating instantaneous speed
% (not before)

% 2018/2/7
% - simple distance = f(time after ready)
% - size/angle ready transient

% 2018/1/19
% estimating instantaneous speed during ts

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
initRSG2prior;
nTs=nTspp;
idDebug=  0;

dim= inf; % 3; % 5; % 4; % 2; % 3; %%%%%
useOptimD= 1; % 0; % if 1, use optimD
maxDim=16; % 11; % H/G, ER/EL/HR/HL: 8 7 10 9 7 7 9 6 for priorIncSL, 8 7 9 8 10 9 11 7

idStitch=0; % 1;
idIncShortInLong=1; % 0; % 1; if 2, use whole [ready set]
idExpFig=0; % 1;

idUseLowDimAngle=0; % 1; % 1; % 0; % 1; % 
dAngle= 2; % 3; % used only if idUseLowDimAngle==1

nPointDer=3; % 5; % 7; % num point for stencil derivative estimation (for tPeak; c.f. nPstencil=3 in estCurvature.m)
idSepPC=0; % 1; % 0;
idPlotPVAF=1; % if 0, use abs. amplitude from fitting for x axis (vs. peak times)
idPVAF=0; % 1; % 0; % if 0, common across S/L/SL; if 1, estimated separately
isUseMs=0; % 1; % 0; % if 1 use [ms], if 0 use [s]
idUseQuadFit=0; % 1;
idMed=1; % -1; % 0; % 1; % median for curvauture if 1; if -1, use max
idDebugCurv=0; % showing individual(PC&pr) quadFit & curvature
% idAngleReadout=1;
idAvgAtt=0;

idTrunBuffer=1; % estimating PVAF; same as runPSTH2traj_prior_conditionSpecific.m

idDataFig=[1 1 2]; % H EL

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

% figure's PaperUnits: 'inches'
optsExpFig.Height=6/2.54; % 2.1/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
optsExpFig.Width=6/2.54; % 3.15/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
optsExpFig.FontSize='8';
optsExpFig.FontMode='fixed';
optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
optsExpFig.LockAxes=1;
optsExpFig.LineMode='scaled'; % 'fixed'; % 'scaled';
optsExpFig.LineWidthMin=1;
optsExpFig.LineWidthMin=2;
% optsExpFig.LineWidth=.5; %1;
optsExpFig.Renderer='painters';

msize=6; % 3; % 5; % 10; % 4; % 6; % 2;
msize2=12; % 9;
lw=1.2; % 1; % 3; % 1.5;
lw2=0.5; % 1; % .5; % for each data set
% H/G, E/H, R/L
cmapHG=[0 0 0; .5 .5 .5]; % animal & E/H
cmapCond={[0 0 0; 1 1 1];
    [.5 .5 .5; 1 1 1]};
markers={'o','^','o','^'};

cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
cmapPr2=[rgb('IndianRed'); rgb('DodgerBlue'); rgb('ForestGreen')]; % for each data set

%% initial
binSize=20;
if idIncShortInLong==1
    if idSepPC
        nPr=nPr+1; % 3
        % short/long/shortInLong
        d=dir('trajKS_prior_S_*.mat'); % *:ER_G
        d2=dir('trajKS_prior_L_*.mat'); % *:ER_G
        d3=dir('trajKS_prior_SiL_*.mat'); % *:ER_G
        d=cat(1,d,d2,d3);
        T{end+1}=repmat(T{1}(end),nTspp,1); % SIL
        nBins=round([T{1}(:); T{2}(:); T{3}(:)]/binSize); % 24 28 ... 60
        tmpCmap{end+1,1}=repmat(rgb('DarkGreen'),nTspp,1); % SIL
        prNm{3}='Short in long';
    else %%%%% now used
        % use stitched trajKS
        if idStitch
            d=dir('trajKS_prior_*_stitchSiL.mat'); % ts_*bin20_smth40.mat');
        else
            nPr=nPr+1; % 3
            d=dir('trajKS_prior_*_bin20_smth40.mat'); % ts_*bin20_smth40.mat');
            T{end+1}=repmat(T{1}(end),nTspp,1); % SIL
            nBins=round([T{1}(:); T{2}(:); T{3}(:)]/binSize); % 24 28 ... 60
            tmpCmap{end+1,1}=repmat(rgb('DarkGreen'),nTspp,1); % SIL
            prNm{3}='Short in long';
        end
    end
elseif idIncShortInLong==2 % [RS]
    nPr=nPr+1; % 3
    d=dir('trajKS_ts_*_bin20_smth40.mat'); % ts_*bin20_smth40.mat');
    T{end+1}=repmat(T{1}(end),nTspp,1); % SIL
    nBins=round([T{1}(:); T{2}(:); T{3}(:)]/binSize); % 24 28 ... 60
    tmpCmap{end+1,1}=repmat(rgb('DarkGreen'),nTspp,1); % SIL
else
    d=dir('trajKS_prior_*_supportOnly.mat'); % ts_*bin20_smth40.mat');
    nBins=round([T{1}(:); T{2}(:)]/binSize); % 24 28 ... 60
end

try
    cd(psthDir);
catch
    psthDir='/Users/seonminahn/Dropbox (MIT)/psthDataHigh';
    cd(psthDir);
end

% speed & arc length (whole meas/overlap)
idExcB4pr=0; % 1; % 0; % 1; % speed weird before shortest ts
nExcB4pr=4; % for both prior, 4 data points
nT=(sum(diff(T{2}))+sum(T{2}(2)-T{2}(1)))/binSize; % 400/20=20
v=nan(nAnimal,nEH,nTarg,nPr,nTs,nT); % speed [2HG x 2EH x 2RL x 2SL x 5ts x 20timePoints]
vv=nan(nAnimal,nEH,nTarg,nPr,nTs,nT,maxDim); % velocity vector [2HG x 2EH x 2RL x 2SL x 5ts x 20timePoints x dim]
avv=nan(nAnimal,nEH,nTarg,nPr,nT,nT); % angle b/t velocity vector [2HG x 2EH x 2RL x 2SL x 20timePoints x 20timePoints]
avrv=nan(nAnimal,nEH,nTarg,nPr,nT); % angle b/t readout vector and velocity vector  [2HG x 2EH x 2RL x 2SL x 20timePoints]
v2=nan(nAnimal,nEH,nTarg,nPr,nTs,nT); % ||2nd derivative|| [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
c=nan(nAnimal,nEH,nTarg,nPr,nTs,nT); % curvature [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
nEpoch=3; % wholeMeas/overlap, short/long, SIL
x=nan(nAnimal,nEH,nTarg,nPr,nTs,nEpoch); % arc length [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]

% - simple distance = f(time after ready)
xPr=nan(nAnimal,nEH,nTarg,nPr,nTs,nT); % arc length [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
% - size/angle ready transient
asv=nan(nAnimal,nEH,nTarg,nPr*nTs,nPr*nTs);
msv=nan(nAnimal,nEH,nTarg,nPr,nTs);

% individual PC's pvaf & curvature
PVAF=nan(nAnimal,nEH,nTarg,nPr,maxDim); % common across S/L/SL
curv=nan(nAnimal,nEH,nTarg,nPr,nT,maxDim); % [222 x prior x time x dim]

pNm={'amplitude','preferred time','tuning width','baseline'};
nParam=length(pNm); % 4; % gaussEqn='a*exp(-((x-b)/c)^2)+d';
pFit=nan(nAnimal,nEH,nTarg,nPr,maxDim,nParam);

tPeak=nan(nAnimal,nEH,nTarg,nPr,maxDim); % -1,maxDim);
curvPeak=nan(nAnimal,nEH,nTarg,nPr,maxDim);

for iFile=1:length(d)
    %% for all conditions
    nm=d(iFile).name;disp(['===== ' nm ' =====']);
    load(nm); % D keep_neurons binSize smthWidth optimD use_sqrt proj_matrix eigenvalues
        
    % remove before shortest ts & build ShortInLong
    if idIncShortInLong==2 % [RS]
%         for iD=1:(nTspp*2)
%             % SL
%             if iD>=(nTspp+1)
%                 D(iD+nTspp)=D(iD);
%                 D(iD+nTspp).condition=['Short' D(iD+nTspp).condition];
%                 nBinTmp=T{1}(1)/binSize; % 480/20=24 bins
%                 iBin=nBinTmp-1; % from 23th(440-460) bin (include -1 to make 480trials useful for vv)
%                 jBin=T{1}(end)/binSize; % 800/20=40
%                 D(iD+nTspp).data=D(iD+nTspp).data(:,iBin:jBin);
%             end
%             iPr=1+floor((iD-1)/nTspp); % 1111122222
%             iBin=T{iPr}(1)/binSize; % 480/20=24 bins or 800/20=40 bins
%             D(iD).data=D(iD).data(:,iBin:end); 
%         end
    end
    
    if useOptimD % isinf(dim)
        load(nm,'optimD'); dim=optimD;
    end
    
    % finding animal, condition > use that projection
    if ~isempty(strfind(d(iFile).name,'_H_')) | ~isempty(strfind(d(iFile).name,'_H.mat'))
        idAnimal=1; iAnimalNm='H';
    else
        idAnimal=2; iAnimalNm='G';
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
    
    %%
    if idSepPC
        %% find out prior
        if ~isempty(strfind(d(iFile).name,'_S_'))
            iPr=1;
        elseif ~isempty(strfind(d(iFile).name,'_L_'))
            iPr=2;
        elseif ~isempty(strfind(d(iFile).name,'_SiL_'))
            iPr=3;
        end
        
        %% PVAF
        %     PVAF(idAnimal,iEH,iTarg,1:optimD)=eigenvalues(1:optimD)./sum(eigenvalues); %   same across S/L/SL
        disp('PVAF separate across S/L/SL:');
        disp(num2str(eigenvalues(1:optimD)'./sum(eigenvalues)));
        PVAF(idAnimal,iEH,iTarg,iPr,1:optimD)=eigenvalues(1:optimD)./sum(eigenvalues);
        
        
        %% avgAtt
        tmpDS=fillNan(D(1:nTs));
        tmpDcat=cat(3,tmpDS.data); % [nPC time ts]
        tmpD{1}=tmpDcat;
        tmpDcurv{1}=nanmean(tmpDcat,3); % [nPC time]
        
        
        if idAnimal==idDataFig(1) & idDataFig(2)==iEH & iTarg==idDataFig(3)
            if iPr==1
                dFig{1}=cat(3,tmpDS.data); % [nPC time ts]
            elseif iPr==2
                dFig{2}=cat(3,tmpDS.data);
            else
                dFig{3}=tmpDcat;
            end
            %          dFig=tmpDcurv;
        end
        
        %% peak time as a measure of curvature
        % (TBD) after separate PCA
        %     gaussEqn='a*exp(-((x-b)/c)^2)+d';
        % pFit=nan(nAnimal,nEH,nTarg,nPr,maxDim,nParam);
        
        %     for iPr=1:nPr
        % time vector common across priors
        t0=T{2-abs(iPr-2)}(1)-binSize/2; % 123> 2-101=121 % 470/790ms; shrotest ts 480>24th bin with bin size of 20ms
        t1=t0+binSize*(size(tmpD{1},2)-1); % 810/1201 for S/L
        tmpTmpT=t0:binSize:t1;
        tmpT=repmat(tmpTmpT,[size(tmpD{1},1) 1 nTspp]); % to match tmpD: [nPC time ts]
        for iDim=1:optimD
            tmpDD=tmpD{1}(iDim,:);
            tmpTT=tmpT(iDim,:);
            
            % remove nan with attrition
            tmpTT=tmpTT(~isnan(tmpDD));
            tmpDD=tmpDD(~isnan(tmpDD));
            
            % curve: coeffvalues (coeffnames), confint, feval, formula,
            % gof: sse, rsquare, dfe, adjrsquare, rmse
            % output: numobs, numparam, residuals
            [curve,gof,output]=fitGauss(tmpTT(:),tmpDD(:));
            
            if idDebug & iDim==1 % & idAnimal==2 & iEH==2 & iDim==1
                figure;
                plot(curve); ha; plot(tmpTT,tmpDD,'o','markersize',2);
                title([animalNm{idAnimal} ' ' ehNm{iEH} ' ' targNm{iTarg} ': PC' num2str(iDim) ', ' prNm{iPr}]);
                setFigPos(2-abs(iPr-2),(idAnimal-1)*nEH*nTarg+(iEH-1)*nTarg+iTarg);
                %                 setFigPos(1,2);
                text(mean(tmpTT),mean(tmpDD),num2str([curve.a curve.b curve.c curve.d]));
                xlabel('time after Ready (ms)');  ylabel('PC');
                applytofig(gcf,optsExpFig);
                %                 waitforbuttonpress; close;
            end
            
            pFit(idAnimal,iEH,iTarg,iPr,iDim,:)=[curve.a curve.b curve.c curve.d];
            
        end % for iDim=1:optimD
        %     end % for iPr=1:nPr
        
        
        %% curvature of each PC after avgAtt
        if idDebugCurv,figure(999); set(gcf,'position',[55 8 1057 797]); figure(1000); setFigPos(2,4); end
        for iDim=1:optimD
            %          for iPr=1:nPr
            if isUseMs
                dataCurv=[binSize*(0:(size(tmpDcurv{1},2)-1)); ... % making it into /s;         ;...%
                    tmpDcurv{1}(iDim,:)];
            else % use second
                dataCurv=[binSize*(0:(size(tmpDcurv{1},2)-1))/1000; ... % making it into /s;         ;...%
                    tmpDcurv{1}(iDim,:)];
            end
            % quadratic fit
            pQuad=polyfit(dataCurv(1,:),dataCurv(2,:),2); % 2nd order (1st element for quadratic coeff.)
            ax2_b=2*pQuad(1)*dataCurv(1,:)+pQuad(2);
            cTmp2=abs(2*pQuad(1)./((1+(ax2_b).^2).^(3/2))); % abs(pQuad(1));
            cQuad=estCurvature([dataCurv(1,:);polyval(pQuad,dataCurv(1,:))]);
            % main
            cTmp=estCurvature(dataCurv);
            if idUseQuadFit
                curv(idAnimal,iEH,iTarg,iPr,1:length(cTmp2),iDim)=cTmp2; % [222 x prior x time x dim]
            else
                curv(idAnimal,iEH,iTarg,iPr,1:length(cTmp),iDim)=cTmp; % [222 x prior x time x dim]
            end
            if idDebugCurv,
                figure(999); subplot(optimD,nPr,(iDim-1)*nPr+iPr); ha;
                yyaxis left;
                plot(dataCurv(1,:),dataCurv(2,:));
                plot(dataCurv(1,:),polyval(pQuad,dataCurv(1,:)),'--'); % quad fit
                yyaxis right;
                plot(dataCurv(1,1+[1:length(cTmp)]),cTmp(:)');
                plot(dataCurv(1,:),cTmp2,'--'); % curvature from quad fit
                
                title(['curv: ' num2str(nanmedian(cTmp(:)),2) '   quad: ' num2str(pQuad(1),2)]);
                figure(1000); plot(nanmedian(cTmp),nanmedian(cTmp2),'ko'); ha; xlabel('curvature from raw data'); ylabel('curvature from quadratic fit');
            end
            %          end % for iPr=1:nPr
        end % for iDim=1:optimD
        
        
        
        
        
        
%% if idSepPC %%%%% now used
    else
        
        
    %% PVAF
    %     PVAF(idAnimal,iEH,iTarg,1:optimD)=eigenvalues(1:optimD)./sum(eigenvalues); %   same across S/L/SL
    if idPVAF==0
        disp('PVAF common across S/L/SL:');
        disp(num2str(eigenvalues(1:optimD)'./sum(eigenvalues)));
        for iPr=1:nPr
            PVAF(idAnimal,iEH,iTarg,iPr,1:optimD)=eigenvalues(1:optimD)./sum(eigenvalues);
        end
    else % separtely across S/L/SL
        dPSTH=load(['PSTH_prior_' animalNm{idAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D nTrialExc tBuffer
        for iPr=1:nPr
            idD0=((iEH-1)*nEH+iTarg)+(iPr-1)*(nEH*nTarg*nTs);% 1234 + 20*(iPr-1)
            idD1=iPr*(nEH*nTarg*nTs); % 20 40 60
            
            % centering
            proj_matrixTmp=cat(3,proj_matrix,[meanPSTH(:) nan(size(proj_matrix,1),size(proj_matrix,2)-1)]); % [#cell x dim x 2(projMat/meanPSTH)]
            
            [Dtmp,proj_matrixTmp,keep_neuronsTmp,eigenvalues]=PSTH2traj(dPSTH.D(idD0:(nEH*nTarg):idD1),...
                proj_matrixTmp,keep_neurons,binSize,smthWidth,[],use_sqrt); %  proj_matrix,keep_neurons optimD
            PVAF(idAnimal,iEH,iTarg,iPr,1:optimD)=eigenvalues(1:optimD)./sum(eigenvalues);
            %         if sum(PVAF(idAnimal,iEH,iTarg,iPr,1:optimD)>0.5)
            %             disp('too high PVAF?');
            %         end
            disp('PVAF separate across S/L/SL:');
            disp(num2str(eigenvalues(1:optimD)'./sum(eigenvalues)));
        end
    end % if idPVAF==0
    
    
    %% avgAtt
    tmpDS=fillNan(D(1:nTs));
    tmpDcat=cat(3,tmpDS.data); % [nPC time ts]
    tmpD{1}=tmpDcat;
    tmpDcurv{1}=nanmean(tmpDcat,3); % [nPC time]
    
    tmpDL=fillNan(D((nTs+1):(2*nTs)));
    tmpDcat=cat(3,tmpDL.data); % [nPC time ts]
    tmpD{2}=tmpDcat;
    tmpDcurv{2}=nanmean(tmpDcat,3); % [nPC time]
    
    if ~idStitch
        tmpDcat=cat(3,D((2*nTs+1):end).data); % [nPC time ts]
        tmpD{3}=tmpDcat;
        tmpDcurv{3}=nanmean(tmpDcat,3); % [nPC time]
    end
    
    if idAnimal==idDataFig(1) & idDataFig(2)==iEH & iTarg==idDataFig(3)
        dFig{1}=cat(3,tmpDS.data); % [nPC time ts]
        dFig{2}=cat(3,tmpDL.data);
        if ~idStitch
            dFig{3}=tmpDcat;
        end
        %          dFig=tmpDcurv;
    end
    
    %% peak time w/ diff after avgAtt
    % simple max(abs)
    % peak across [short in long]+[long], rather than [long],
    % remove 1st data point of long as it's identical to last data point of [short in long]
    % tPeak=nan(nAnimal,nEH,nTarg,nPr-1,maxDim);
    
    % time vectors: short
    iPr=1;
    t0=T{2-abs(iPr-2)}(1)-binSize/2; % 123> 2-101=121 % 470/790ms; shrotest ts 480>24th bin with bin size of 20ms
    t1=t0+binSize*(size(tmpD{iPr},2)-1); % 810/1210 for S/L
    tVect{iPr}=t0:binSize:t1; % 470:20:810
    iPr=2; % long
    if ~idStitch
        t0=T{2-abs(iPr-2)}(1)-binSize/2;
    end
    t1=t0+binSize*(size(tmpD{iPr},2)-1); % +size(tmpD{iPr+1},2)-2-1); % remove 1st data point of long as it's identical to last data point of [short in long]
    tVect{iPr}=t0:binSize:t1; % 470:20:1210
    if ~idStitch
        iPr=3;
        t0=T{2-abs(iPr-2)}(1)-binSize/2; % 123> 2-101=121 % 470/790ms; shrotest ts 480>24th bin with bin size of 20ms
        t1=t0+binSize*(size(tmpD{iPr},2)-1); % 810/1210 for S/L
        tVect{iPr}=t0:binSize:t1; % 470:20:810
    end
        
    % tPeak=nan(nAnimal,nEH,nTarg,nPr-1,maxDim);
    for iPr=1:nPr % (nPr-1)
%         if iPr==nPr % (nPr-1)% stich SiL & L
%             tmpTmpDcurv=[tmpDcurv{2}(:,1:(end-1)) tmpDcurv{3}(:,2:end)];
%         else
%             tmpTmpDcurv=tmpDcurv{iPr}; % [nPC time]
%         end

        % no initial/final points
        if nPointDer>0
            tVect{iPr}=tVect{iPr}((1+(nPointDer-1)/2):(end-(nPointDer-1)/2));
        end

        for iDim=1:optimD
            % find point with dF(x)/dx=0 after smoothing
            tmpDD=tmpDcurv{iPr}(iDim,:); % still with initial/final points
            dFdx=deriv(tmpDD,binSize,1,nPointDer); %diff(tmpDD); % deriv(tmpDcurv{iPr}(iDim,:),binSize,1,7); % order,numPoint % no initial/final points
            tStat=find(diff(dFdx>0)~=0); %deriv(dFdx>0,1,1,nPointDer)~=0); %diff(dFdx>0)~=0); % no initial/final points
            if ~isempty(tStat) % if there exist deflection points
                % deal with multiple of them: pick up max amplitude from mean
                [maxTmp,iMax]=max(abs(tmpDD(tStat+(nPointDer-1)/2)-mean(tmpDD)));
                tPeak(idAnimal,iEH,iTarg,iPr,iDim)=tVect{iPr}(tStat(iMax)+1); % [nPC 1] % add 1 for use 1 time point later
            end

            if idDebug % & idAnimal==2 & iEH==2 & iDim==1
                figure;  ha; setFigPos(1,2);  % setFigPos(2-abs(iPr-2),(idAnimal-1)*nEH*nTarg+(iEH-1)*nTarg+iTarg);
                plot(tVect{iPr},tmpDD((1+(nPointDer-1)/2):(end-(nPointDer-1)/2)),'o'); 
                plotHorizon(gca,mean(tmpDD),[.5 .5 .5]);
                if ~isempty(tStat) % if there exist deflection points
                    plotVertical(gca,tVect{iPr}(tStat+1),[.5 .5 .5]); % all deflection points
                    plotVertical(gca,tVect{iPr}(tStat(iMax)+1),[1 0 0]); % deflection point with max amplitude
                end
                title([animalNm{idAnimal} ' ' ehNm{iEH} ' ' targNm{iTarg} ': PC' num2str(iDim) ', ' prNm{iPr}]);
                setFigPos(1,2); 
                waitforbuttonpress; close;
            end
        end
%         [maxTmp,iMax]=max(abs(tmpDcurv{iPr}),[],2); % tmpTmpDcurv,[],2);
%         tPeak(idAnimal,iEH,iTarg,iPr,1:length(iMax))=tVect{iPr}(iMax); % [nPC 1]
    end

    %% peak time as a measure of curvature 
    % @2018/7/9: stich SiL & L
    % (TBD) after separate PCA
    %     gaussEqn='a*exp(-((x-b)/c)^2)+d';
    % pFit=nan(nAnimal,nEH,nTarg,nPr,maxDim,nParam);

%     for iPr=1:nPr % (nPr-1)
%         % time vector common across priors
%         t0=T{1}(1)-binSize/2; % 123> 2-101=121 % 470/790ms; shrotest ts 480>24th bin with bin size of 20ms
%         if iPr==1, t1=t0+binSize*(size(tmpD{iPr},2)-1); % 810/1201 for S/L
%         else t1=t0+binSize*(size(tmpD{iPr},2)-1); end % +size(tmpD{iPr+1},2)-2-1);  end
%         tmpTmpT=t0:binSize:t1;
%         tmpT=repmat(tmpTmpT,[size(tmpD{iPr},1) 1 nTspp]); % to match tmpD: [nPC time ts]
%         for iDim=1:optimD
%              tmpDD=tmpD{iPr}(iDim,:);
% %             if iPr==1, tmpDD=tmpD{iPr}(iDim,:);
% %             else tmpDD=cat(2,tmpD{iPr}(iDim,1:(end-1),:),tmpD{iPr+2}(iDim,2:end,:)); tmpDD=tmpDD(:); end
%             tmpTT=tmpT(iDim,:);
%             
%             % remove nan with attrition
%             tmpTT=tmpTT(~isnan(tmpDD));
%             tmpDD=tmpDD(~isnan(tmpDD));
%             
%             % curve: coeffvalues (coeffnames), confint, feval, formula,
%             % gof: sse, rsquare, dfe, adjrsquare, rmse
%             % output: numobs, numparam, residuals
%             [curve,gof,output]=fitGauss(tmpTT(:),tmpDD(:));
%             
%             if idDebug & idAnimal==2 & iEH==2 & iDim==1
%                 figure;  ha; setFigPos(2-abs(iPr-2),(idAnimal-1)*nEH*nTarg+(iEH-1)*nTarg+iTarg);
%                 plot(tmpTT,tmpDD,'o'); plot(curve);  
%                 title([animalNm{idAnimal} ' ' ehNm{iEH} ' ' targNm{iTarg} ': PC' num2str(iDim) ', ' prNm{iPr}]);
%                 setFigPos(1,2); text(mean(tmpTT),mean(tmpDD),num2str([curve.a curve.b curve.c curve.d]));
%                 waitforbuttonpress; close;
%             end
%             
%             pFit(idAnimal,iEH,iTarg,iPr,iDim,:)=[curve.a curve.b curve.c curve.d];
%             
%         end % for iDim=1:optimD
%     end % for iPr=1:nPr
    
    % back up before stitch
%     for iPr=1:nPr
%         % time vector common across priors
%         t0=T{2-abs(iPr-2)}(1)-binSize/2; % 123> 2-101=121 % 470/790ms; shrotest ts 480>24th bin with bin size of 20ms
%         t1=t0+binSize*(size(tmpD{iPr},2)-1); % 810/1201 for S/L
%         tmpTmpT=t0:binSize:t1;
%         tmpT=repmat(tmpTmpT,[size(tmpD{iPr},1) 1 nTspp]); % to match tmpD: [nPC time ts]
%         for iDim=1:optimD
%             tmpDD=tmpD{iPr}(iDim,:);
%             tmpTT=tmpT(iDim,:);
%             
%             % remove nan with attrition
%             tmpTT=tmpTT(~isnan(tmpDD));
%             tmpDD=tmpDD(~isnan(tmpDD));
%             
%             % curve: coeffvalues (coeffnames), confint, feval, formula, 
%             % gof: sse, rsquare, dfe, adjrsquare, rmse
%             % output: numobs, numparam, residuals
%             [curve,gof,output]=fitGauss(tmpTT(:),tmpDD(:));
%             
%             if idDebug & idAnimal==2 & iEH==2 & iDim==1
%                 figure;  ha; setFigPos(2-abs(iPr-2),(idAnimal-1)*nEH*nTarg+(iEH-1)*nTarg+iTarg);
%                 plot(curve);  plot(tmpTT,tmpDD,'o'); 
%                 title([animalNm{idAnimal} ' ' ehNm{iEH} ' ' targNm{iTarg} ': PC' num2str(iDim) ', ' prNm{iPr}]);
%                 setFigPos(1,2); text(mean(tmpTT),mean(tmpDD),num2str([curve.a curve.b curve.c curve.d]));
%                 waitforbuttonpress; close;
%             end
%             
%             pFit(idAnimal,iEH,iTarg,iPr,iDim,:)=[curve.a curve.b curve.c curve.d];
%             
%         end % for iDim=1:optimD
%     end % for iPr=1:nPr
    
    
    %% curvature of each PC after avgAtt
     if idDebugCurv,figure(999); set(gcf,'position',[55 8 1057 797]); figure(1000); setFigPos(2,4); end
     for iDim=1:optimD
         for iPr=1:nPr
             if isUseMs
                 dataCurv=[binSize*(0:(size(tmpDcurv{iPr},2)-1)); ... % making it into /s;         ;...%
                     tmpDcurv{iPr}(iDim,:)];
             else % use second
                 dataCurv=[binSize*(0:(size(tmpDcurv{iPr},2)-1))/1000; ... % making it into /s;         ;...%
                     tmpDcurv{iPr}(iDim,:)];
             end
             % quadratic fit
             pQuad=polyfit(dataCurv(1,:),dataCurv(2,:),2); % 2nd order (1st element for quadratic coeff.)
             ax2_b=2*pQuad(1)*dataCurv(1,:)+pQuad(2);
             cTmp2=abs(2*pQuad(1)./((1+(ax2_b).^2).^(3/2))); % abs(pQuad(1));
             cQuad=estCurvature([dataCurv(1,:);polyval(pQuad,dataCurv(1,:))]);
             % main
             cTmp=estCurvature(dataCurv); % [1 size(dataCurv,2)-2]
             if idUseQuadFit
                 curv(idAnimal,iEH,iTarg,iPr,1:length(cTmp2),iDim)=cTmp2; % [222 x prior x time x dim]
             else
                 curv(idAnimal,iEH,iTarg,iPr,1:length(cTmp),iDim)=cTmp; % [222 x prior x time x dim]
             end
             
             % get curvuture at tPeak
             tmpTmpT=dataCurv(1,:)*1000+T{2-abs(iPr-2)}(1)-binSize/2; % [470:810]
             tmpTmpT=tmpTmpT(2:(end-1)); % nPstencil=3; in estCurvature
             tmpId=tPeak(idAnimal,iEH,iTarg,iPr,iDim)==tmpTmpT;
             if sum(tmpId)>0
                 curvPeak(idAnimal,iEH,iTarg,iPr,iDim)=cTmp(tPeak(idAnimal,iEH,iTarg,iPr,iDim)==tmpTmpT);
             else % no peak>median
                 curvPeak(idAnimal,iEH,iTarg,iPr,iDim)=squeeze(nanmedian(curv(idAnimal,iEH,iTarg,iPr,:,iDim),5)); % avg time
             end
             
             if idDebugCurv, 
                 figure(999); subplot(optimD,nPr,(iDim-1)*nPr+iPr); ha; 
                 yyaxis left;
                 plot(dataCurv(1,:),dataCurv(2,:)); 
%                  plot(dataCurv(1,:),polyval(pQuad,dataCurv(1,:)),'--'); % quad fit
                 yyaxis right;
                 plot(dataCurv(1,1+[1:length(cTmp)]),cTmp(:)'); 
%                  plot(dataCurv(1,:),cTmp2,'--'); % curvature from quad fit
                
                title(['curv: ' num2str(nanmedian(cTmp(:)),2) '   quad: ' num2str(pQuad(1),2)]);
                figure(1000); plot(nanmedian(cTmp),nanmedian(cTmp2),'ko'); ha; xlabel('curvature from raw data'); ylabel('curvature from quadratic fit');
             end
         end % for iPr=1:nPr
     end % for iDim=1:optimD
     if idDebugCurv, waitforbuttonpress; 
         close all; 
     end
    
%     if idDebug,
%         disp(squeeze(asv(idAnimal,iEH,iTarg,:,:)));
%     end
    
    
    end % if idSepPC
end % for i=1:length(d)

if idSepPC
    save('summaryGaussFitCurvPVAF_sepPCA.mat','PVAF','curv','pNm','nParam','pFit',...
        'useOptimD','idIncShortInLong','idSepPC','idPVAF','isUseMs','idMed','idTrunBuffer');
else
    save('summaryGaussFitCurvPVAF_stitchSiL.mat','PVAF','curv','pNm','nParam','pFit',...
        'useOptimD','idIncShortInLong','idSepPC','idPVAF','isUseMs','idMed','idTrunBuffer','tPeak','curvPeak');
%     save('summaryGaussFitCurvPVAF.mat','PVAF','curv','pNm','nParam','pFit',...
%         'useOptimD','idIncShortInLong','idSepPC','idPVAF','isUseMs','idMed','idTrunBuffer','tPeak');
end

%% statistical test b/t short vs long curvature
% if idMed==0
%     mCurv=nanmean(curv,5); % mean across time
% elseif idMed==1
%     mCurv=nanmedian(curv,5); % median across time
% elseif idMed==-1
%     mCurv=nanmax(curv,5); % max across time
% end
% mCshort=squeeze(mCurv(:,:,:,1,:,:));mCshort=mCshort(~isnan(mCshort));
% mClong=squeeze(mCurv(:,:,:,2,:,:));mClong=mClong(~isnan(mClong));
% mCSIL=squeeze(mCurv(:,:,:,3,:,:));mCSIL=mCSIL(~isnan(mCSIL));
% disp(['p(signrank(mCshort,mClong))= ' num2str(signrank(mCshort,mClong),3)]);
% % figure;scatterhist(mCshort,mClong); plotIdentity(gca); setFigPos(1,1); 
% % figure; histogram(mCshort-mClong,15,'DisplayStyle','stairs'); setFigPos(1,1); plotVertical(gca,0,[]); ha; plotVertical(gca,mean(mCshort-mClong),[1 0 0]); xlabel('curvature(short)-curvature(long)')
% disp(['p(signrank(mCshort,mCSIL))= ' num2str(signrank(mCshort,mCSIL),3)]);
% % figure;scatterhist(mCshort,mCSIL); plotIdentity(gca); setFigPos(2,1); 
% % figure; histogram(mCshort-mCSIL,15,'DisplayStyle','stairs'); setFigPos(2,1);  plotVertical(gca,0,[]); ha; plotVertical(gca,mean(mCshort-mCSIL),[1 0 0]); xlabel('curvature(short)-curvature(short in long)')
% disp(['p(signrank(mCSIL,mClong))= ' num2str(signrank(mCSIL,mClong),3)]);
% % figure;scatterhist(mCSIL,mClong); plotIdentity(gca); setFigPos(1,2); 
% % figure; histogram(mCSIL-mClong,15,'DisplayStyle','stairs'); setFigPos(1,2);  plotVertical(gca,0,[]); ha; plotVertical(gca,mean(mCSIL-mClong),[1 0 0]); xlabel('curvature(short in long)-curvature(long)')
% 
% PVAFTmp=squeeze(PVAF(:,:,:,1,:,:));
% pvafNoNan=PVAFTmp(~isnan(PVAFTmp));
% 
% figure; setFigPos(2,3);
% % simple scatter plot: long vs short/shortInLong
% gainPVAF=msize;
% scatter(mClong,mCshort,...%     normMinMax(pvafNoNan)*(gainPVAF-1)+1
%     msize,cmapPr(1,:),'linewidth',lw,'markerfacecolor','w'); ha;
% scatter(mClong,mCSIL,...%     normMinMax(pvafNoNan)*(gainPVAF-1)+1,
%     msize,cmapPr(3,:),'linewidth',lw,'markerfacecolor','w'); ha;
% axis square;
% plotIdentity(gca);
% set(gca,'ticklength',[0.01 0.01],'tickdir','out');
% xlabel('Curvature (Long)'); ylabel('Curvature');

%             figure;scatterhist(mCshort,mClong); plotIdentity(gca); setFigPos(1,1);
%             figure; histogram(mCshort-mClong,15,'DisplayStyle','stairs'); setFigPos(1,1); plotVertical(gca,0,[]); ha; plotVertical(gca,mean(mCshort-mClong),[1 0 0]); xlabel('curvature(short)-curvature(long)')
%             
%             figure;scatterhist(mCshort,mCSIL); plotIdentity(gca); setFigPos(2,1);
%             figure; histogram(mCshort-mCSIL,15,'DisplayStyle','stairs'); setFigPos(2,1);  plotVertical(gca,0,[]); ha; plotVertical(gca,mean(mCshort-mCSIL),[1 0 0]); xlabel('curvature(short)-curvature(short in long)')

% % ternary plot: not used as scaling issues
% 
% % Plot the ternary axis system
% [h,hg,htick]=terplot; % handles for patch(fill), grid, tick labels
% % Plot the data
% % First set the colormap (can't be done afterwards)
% colormap(flipud(bone));
% [hd,hcb]=ternaryc(mCshort,mCSIL,mClong,pvafNoNan,'o');
% % Add the labels
% hlabels=terlabel('Curvature (Short)','Curvature (Short in long)','Curvature (Long)');
% 
% set(hd,'linewidth',lw,'markersize',msize,'markerfacecolor','w');
% set(hlabels(3),'color',cmapPr(2,:));
% set(hlabels(2),'color',cmapPr(3,:)); % Short in long
% set(hlabels(1),'color',cmapPr(1,:));
% % %--  Modify the tick labels
% % set(htick(:,1),'color','y','linewidth',3)
% % set(htick(:,2),'color','c','linewidth',3)
% % set(htick(:,3),'color','m','linewidth',3)
% % %--  Change the colorbar
% % set(hcb,'xcolor','w','ycolor','w')

%% plot
% v=nan(nAnimal,nEH,nTarg,nPr,nTs,nT); % speed [2HG x 2EH x 2RL x 2SL x 5ts x 20timePoints]
% vv=nan(nAnimal,nEH,nTarg,nPr,nTs,nT,maxDim); % velocity vector [2HG x 2EH x 2RL x 2SL x 5ts x 20timePoints x dim]
% avv=nan(nAnimal,nEH,nTarg,nPr,nT,nT); % angle b/t velocity vector [2HG x 2EH x 2RL x 2SL x 20timePoints x 20timePoints]

% hFig1=figure(1); ha;setFigPos(1,1); % speed vs time (avg. across ts) two lines for short vs long
% hFig2=figure(2); ha;setFigPos(1,2); % avv
% hFig3=figure(3); ha;setFigPos(1,3); % avrv
% % figure(3); setFigPos(1,3); % avv colormap for each data set: short
% % figure(4); setFigPos(2,3); % smae for long
% hFig4=figure(4); ha;setFigPos(2,1); % curvature vs PVAF
% hFig5=figure(5); ha;setFigPos(2,2); % check curvature distribution across time

condNm={'ER','EL','HR','HL'};
nCmap=length(animalNm)*length(condNm);
cmap3=parula(nCmap); % hsv(nCmap);
icmap=1;

tmpMS=[];tmpML=[];tmpMSL=[];
AS=[]; AL=[];ASL=[];
AS2=[]; AL2=[];ASL2=[];

curvPVAFc=cell(nPr,1); % for linear regression [PVAF curv]

for i=1:nAnimal
    for j=1:nEH
        for k=1:nTarg
            iCond= (j-1)*nTarg+k;
            
            
            %% % fig3b1: only pc1
%             figure(3);ha;setFigPos(1,1);
%             
%             for iPr=1:nPr
%                 if idDataFig(1)==i & idDataFig(2)==j & idDataFig(3)==k
%                     pcSc=squeeze(dFig{iPr}(1,:,:)); % {S/L/SL}[nPC x time x ts]
%                     if iPr==2
%                         tmpX=T{iPr}(1):binSize:(T{iPr}(1)+binSize*(size(pcSc,1)-1));
%                     else
%                         tmpX=T{1}(1):binSize:(T{1}(1)+binSize*(size(pcSc,1)-1));
%                     end
%                     shadedErrorBar(tmpX(:)',nanmean(pcSc,2)',sem(pcSc,2,'omitnan'),...
%                         {'-','color',cmapPr(iPr,:),'linewidth',lw},1); drawnow; % 2
%                     
% %                     pcSc=dFig{iPr}(1,:) % {S/L/SL}[nPC x time] when dFig is tmpDcurv
% %                     if iPr==2
% %                         tmpX=T{iPr}(1):binSize:(T{iPr}(1)+binSize*(length(pcSc)-1));
% %                     else
% %                         tmpX=T{1}(1):binSize:(T{1}(1)+binSize*(length(pcSc)-1));
% %                     end
% %                     plot(tmpX,pcSc,'-','color',cmapPr(iPr,:),'linewidth',lw);
% 
% %                     disp(squeeze(curv(i,j,k,iPr,:,1)));
%                     
%                 end % if idDataFig(1)==i & idDataFig(2)==j & idDataFig(3)==k
%             end % iPr
%             axis tight;
            
            
            %% % fig3b2
            %% PVAF plot
%             % PVAF=nan(nAnimal,nEH,nTarg,nPr,maxDim); % common across S/L/SL
%             figure(7); ha; setFigPos(2,2);
%             for iPr=1:nPr
%                 plot(squeeze(PVAF(i,j,k,iPr,:)),'o-','markerfacecolor','w','color',cmapPr(iPr,:));
%             end % for iPr=1:nPr
%             xlabel('PC'); ylabel('% Var. explained');
%             
%             minPVAF=nanmin(PVAF(:));
%             maxPVAF=nanmax(PVAF(:));
            
            %% peak time distribution: PVAF/amplitude vs peak time (marker proportional to PVAF)
            figure(8); ha; setFigPos(1,5);
            for iPr=1:nPr
                for iDim=1:size(pFit,5)
                    if idPlotPVAF
                        tmpX=abs(squeeze(PVAF(i,j,k,iPr,iDim))); % pFit(i,j,k,iPr,iDim,1))); % a
                        tmpMsize=msize2; % (squeeze(PVAF(i,j,k,iPr,iDim))-minPVAF)/(maxPVAF-minPVAF)*(msize2-msize)+msize; % msize for min, msize2 for max
                    else
                        tmpX=abs(squeeze(pFit(i,j,k,iPr,iDim,1))); % a
                        tmpMsize=(squeeze(PVAF(i,j,k,iPr,iDim))-minPVAF)/(maxPVAF-minPVAF)*(msize2-msize)+msize; % msize for min, msize2 for max
                    end
%                     tmpY=squeeze(pFit(i,j,k,iPr,iDim,2)); % b
                    tmpY=squeeze(tPeak(i,j,k,iPr,iDim)); % b
                    if tmpX==0 & tmpY==0
                        disp([i j k iPr iDim]);
                    end
                    if ~isnan(tmpX) & ~isnan(tmpY) & ~isnan(tmpMsize) & tmpX~=0
                        plot(tmpX,tmpY,'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',tmpMsize);
                    end
                end
            end % for iPr=1:nPr
            %             xlabel('% Var. explained'); % amplitude');
            if idPlotPVAF
                xlabel('% Var. explained');
            else
                xlabel('amplitude');
            end
% %             ylabel('preferred time');
            
%             figure(9); ha; setFigPos(1,6); % tuning width vs baseline
%             for iPr=1:nPr
%                 for iDim=1:size(pFit,5)
%                     tmpX=squeeze(pFit(i,j,k,iPr,iDim,3)); % c
%                     tmpY=squeeze(pFit(i,j,k,iPr,iDim,4)); % c
%                     tmpMsize=(squeeze(PVAF(i,j,k,iPr,iDim))-minPVAF)/(maxPVAF-minPVAF)*(msize2-msize)+msize; % msize for min, msize2 for max
%                     if ~isnan(tmpX) & ~isnan(tmpY) & ~isnan(tmpMsize)
%                         plot(tmpX,tmpY,'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',tmpMsize);
%                     end
%                 end
%             end % for iPr=1:nPr
%             xlabel('tuning width'); ylabel('baseline');
            
            %% old curv=f(PVAF)
            hFig4=figure(4);       ha;setFigPos(2,1);
            %             % individual PC's pvaf & curvature
            % PVAF=nan(nAnimal,nEH,nTarg,maxDim); % common across S/L/SL
            % curv=nan(nAnimal,nEH,nTarg,nPr,nT,maxDim); % [222 x prior x time x dim]
            for iPr=1:nPr
                nDim=nnz(~isnan(PVAF(i,j,k,iPr,:)));
                tmpX=PVAF(i,j,k,iPr,1:nDim);
                if idMed==1
                    mCurv=squeeze(nanmedian(curv(i,j,k,iPr,:,1:nDim),5)); % avg time
                elseif idMed==0
                    mCurv=squeeze(nanmean(curv(i,j,k,iPr,:,1:nDim),5)); % avg time
                elseif idMed==-1
                    mCurv=squeeze(nanmax(curv(i,j,k,iPr,:,1:nDim),[],5)); % avg time
                end
                % remove outlier
%                 thOut=9*10^(-5);
%                 if sum(mCurv>thOut),
%                     disp(['outlier:' animalNm{i} ' ' ehNm{j} ' ' targNm{j} ' PC' num2str(find(mCurv>thOut))]); % outlier:G Hand Left PC4
%                     mCurv(mCurv>thOut)=nan;
%                 end
                tmpY=squeeze(curvPeak(i,j,k,iPr,1:nDim)); % mCurv;

                disp(length(mCurv)); % 8 7 10 9 7 7 9 6
                plot(tmpX(:),tmpY,'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize,'markerfacecolor','w');
                 if idDataFig(1)==i & idDataFig(2)==j & idDataFig(3)==k % examples in fig3b1
                     plot(tmpX(1),tmpY(1),'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize,'markerfacecolor',cmapPr(iPr,:));
                 end
%                 % plot PC number
%                 for iT=1:length(mCurv)
%                     text(tmpX(iT),mCurv(iT),[animalNm{i} ehNm{j}(1) targNm{k}(1) num2str(iT)],'color',cmapPr(iPr,:),'fontsize',8);
%                 end
                curvPVAFc{iPr}=[curvPVAFc{iPr}; tmpX(:),mCurv(:)];
            end
            if icmap==nCmap
                axis tight;
%                 for iPr=1:nPr
%                     s=plotReg(curvPVAFc{iPr}(:,1),curvPVAFc{iPr}(:,2),hFig4,cmapPr(iPr,:),0,1); % idRobust, idRedModel (PC1 if 1)
% %                     disp(['beta/t/pval/r2: ' num2str([s.beta(:)' s.tstat.t(:)' s.tstat.pval(:)' s.rsquare],2)]);
%                 end
                xlabel('% var. explained'); 
                ylabel('curvature@peak');
            end
            ylim([0 max(ylim)]);
            icmap=icmap+1;
            
            %% distribution of curvature
%             figure(5); ha; setFigPos(1,2);
%             for iPr=1:nPr
%                 for iDim=1:size(curv,ndims(curv))
%                     tmpCurv2=squeeze(curv(i,j,k,iPr,:,iDim)); % time x 1
%                     if sum(~isnan(tmpCurv2))
%                         [Ntmp,edgeTmp]=histcounts(tmpCurv2(~isnan(tmpCurv2)),'BinMethod','fd');
%                         plot(movavg(edgeTmp,2),Ntmp,'-','color',cmapPr(iPr,:),'linewidth',2);
% %                         histogram(tmpCurv2(~isnan(tmpCurv2)),'DisplayStyle','stairs','edgecolor',cmapPr(iPr,:)); ha;
%                     end
%                 end
%             end
%             xlabel('curvature'); ylabel('# time points');
%             
%             %% bar graphs: PC vs curv
%             figure(6); ha; setFigPos(2,2);
%             for iPr=1:nPr
%                 if idMed
%                     mCurv=squeeze(nanmedian(curv(i,j,k,iPr,:,1:nDim),5)); % avg time
%                 else
%                     mCurv=squeeze(nanmean(curv(i,j,k,iPr,:,1:nDim),5)); % avg time
%                 end
%                 plot(mCurv,'-','color',cmapPr(iPr,:)); ha;
%             end
%             xlabel('PC'); ylabel('median curvature');
            
        end % targ
    end % EH
end % animal

% figure(3); set(gca,'ticklength',[0.01 0.01],'tickdir','out','xtick',[T{1}(1) T{1}(3) T{2}(1) T{2}(3) T{2}(end)],'ytick',-1:1:1);
% figure(4); set(gca,'ticklength',[0.01 0.01],'tickdir','out');


%% plot me(di)an & error
figure(8);
for iPr=1:nPr
    if idPlotPVAF
        tmpX=abs(squeeze(PVAF(:,:,:,iPr,:))); % pFit(:,:,:,iPr,:,1))); % a
    else
        tmpX=abs(squeeze(pFit(:,:,:,iPr,:,1))); % a
    end
%     tmpY=squeeze(pFit(:,:,:,iPr,:,2)); % b
    tmpY=squeeze(tPeak(:,:,:,iPr,:)); % b % tPeak(:,:,:,iPr,:)
    plotErrorbar(nanmedian(tmpX(:)),sem(tmpX(:),1,'omitnan'),nanmedian(tmpY(:)),sem(tmpY(:),1,'omitnan'),cmapPr(iPr,:),2);
    plot(nanmedian(tmpX(:)),nanmedian(tmpY(:)),'o','color',cmapPr(iPr,:),'linewidth',2,'markersize',msize2*1.5,'markerfacecolor','w'); % msize for min, msize2 for max
end % for iPr=1:nPr
if idPlotPVAF
    xlabel('% Var. explained');
else
    xlabel('abs. amplitude');
end
ylabel('preferred time');
% axis([0 3 350 1150]);
plotHorizon(gca,[T{1}(1) T{1}(end) T{2}(end)],[]);

% figure(9); % tuning width vs baseline
% for iPr=1:nPr
%     tmpX=squeeze(pFit(:,:,:,iPr,:,3)); % c width
%     tmpY=squeeze(pFit(:,:,:,iPr,:,4)); % d baseline
%     plotErrorbar(nanmedian(tmpX(:)),sem(tmpX(:),1,'omitnan'),nanmedian(tmpY(:)),sem(tmpY(:),1,'omitnan'),cmapPr(iPr,:),2);
%     plot(nanmedian(tmpX(:)),nanmedian(tmpY(:)),'o','color',cmapPr(iPr,:),'linewidth',2,'markersize',msize2*1.5,'markerfacecolor','w'); % msize for min, msize2 for max
% end % for iPr=1:nPr
% xlabel('tuning width'); ylabel('baseline');

%             pFit(idAnimal,iEH,iTarg,iPr,iDim,:)=[curve.a curve.b curve.c curve.d]; 'a*exp(-((x-b)/c)^2)+d';
figure(10); ha; setFigPos(2,6); % peak distribution
for iParam=1:nParam
    subplot(nParam,1,iParam);ha;
    for iPr=1:nPr
        tmpX=squeeze(pFit(:,:,:,iPr,:,iParam)); % 
        histogram(tmpX(:),round(length(tmpX(:))/3),'DisplayStyle','stairs','edgecolor',cmapPr(iPr,:),'linewidth',2);
        plot(nanmean(tmpX(:)),0,'o','color',cmapPr(iPr,:),'linewidth',2,'markersize',msize2*1.5,'markerfacecolor','w'); % msize for min, msize2 for max
    end % for iPr=1:nPr
    xlabel(pNm{iParam});
end % for iParam=1:nParam

figure(11); ha; setFigPos(1,1); % histogram only for tPeak/prefferred time % tPeak=nan(nAnimal,nEH,nTarg,nPr-1,maxDim);
iParam=2;
bound=[0 1500];
for iPr=1:nPr
%     tmpX=squeeze(pFit(:,:,:,iPr,:,iParam)); %
    tmpX=squeeze(tPeak(:,:,:,iPr,:)); %
    tmpX=tmpX(~isnan(tmpX)); % remove nan
    disp(['% outlier: ' num2str(100-100*mean(tmpX>bound(1) & tmpX<bound(2)))]);
    tmpX=tmpX(tmpX>bound(1) & tmpX<bound(2));
    histogram(tmpX(:),round(length(tmpX(:))/3),'DisplayStyle','stairs','edgecolor',cmapPr(iPr,:),'linewidth',2);
    plot(nanmean(tmpX(:)),0,'o','color',cmapPr(iPr,:),'linewidth',2,'markersize',msize2*1.5,'markerfacecolor','w'); % msize for min, msize2 for max
end % for iPr=1:nPr
xlabel(pNm{iParam});
ylabel('# data sets');
