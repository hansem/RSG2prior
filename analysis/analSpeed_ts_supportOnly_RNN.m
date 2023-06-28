function analSpeed_ts_supportOnly_RNN

% 2018/5/17
% RNN: ls traj_prior_Weber*.mat (traj_prior_Weber_2Prior_Type2_005_v5)

% 2018/4/13
% curvature =f(varExplained), for each PC, colored by S/L/SL
% - %VarExplained calculated for each data set & S/L/SL
% - curvature: fit quadratic function

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
idDebug=0;

dim=inf; % 3; % 5; % 4; % 2; % 3; %%%%%
useOptimD=1; % 0; % if 1, use optimD
maxDim=11; % H/G, ER/EL/HR/HL: 8 7 10 9 7 7 9 6 for priorIncSL, 8 7 9 8 10 9 11 7

idIncShortInLong=1; % 0; % 1; if 2, use whole [ready set]
idExpFig=0; % 1;

idUseLowDimAngle=0; % 1; % 1; % 0; % 1; % 
dAngle= 2; % 3; % used only if idUseLowDimAngle==1

idUseQuadFit=0; % 1;
idMed=1; % median for curvauture if 1
idDebugCurv=0; % 1; % 0; % 1; % showing individual(PC&pr) quadFit & curvature
% idAngleReadout=1;

idTrunBuffer=1; % estimating PVAF; same as runPSTH2traj_prior_conditionSpecific.m


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
optsExpFig.Height=3/2.54; % 2.1/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
optsExpFig.Width=3/2.54; % 3.15/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
optsExpFig.FontSize='6';
optsExpFig.FontMode='fixed';
optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
optsExpFig.LockAxes=1;
optsExpFig.LineMode='scaled'; % 'fixed'; % 'scaled';
optsExpFig.LineWidthMin=0.5;
optsExpFig.LineWidthMin=1.2;
% optsExpFig.LineWidth=.5; %1;
optsExpFig.Renderer='painters';

msize=5; % 10; % 4; % 6; % 2;
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
initRSG2prior;
nTs=nTspp;
binSize=20;
if idIncShortInLong==1
    nPr=nPr+1; % 3
    d=dir('trajKS_prior_*_bin20_smth40.mat'); % ts_*bin20_smth40.mat');
    T{end+1}=repmat(T{1}(end),nTspp,1); % SIL
    nBins=round([T{1}(:); T{2}(:); T{3}(:)]/binSize); % 24 28 ... 60
    tmpCmap{end+1,1}=repmat(rgb('DarkGreen'),nTspp,1); % SIL
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

for iFile=1:length(d)
    %% for all conditions
    nm=d(iFile).name;disp(['===== ' nm ' =====']);
    load(nm); % D keep_neurons binSize smthWidth optimD use_sqrt proj_matrix eigenvalues
        
    % remove before shortest ts & build ShortInLong
    if idIncShortInLong==2 % [RS]
        for iD=1:(nTspp*2)
            % SL
            if iD>=(nTspp+1)
                D(iD+nTspp)=D(iD);
                D(iD+nTspp).condition=['Short' D(iD+nTspp).condition];
                nBinTmp=T{1}(1)/binSize; % 480/20=24 bins
                iBin=nBinTmp-1; % from 23th(440-460) bin (include -1 to make 480trials useful for vv)
                jBin=T{1}(end)/binSize; % 800/20=40
                D(iD+nTspp).data=D(iD+nTspp).data(:,iBin:jBin);
            end
            iPr=1+floor((iD-1)/nTspp); % 1111122222
            iBin=T{iPr}(1)/binSize; % 480/20=24 bins or 800/20=40 bins
            D(iD).data=D(iD).data(:,iBin:end); 
        end
    end
    
    if useOptimD % isinf(dim)
        load(nm,'optimD'); dim=optimD;
    end
    
    % finding animal, condition > use that projection
    if ~isempty(strfind(d(iFile).name,'_H_'))
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
    
    % PVAF
    %     PVAF(idAnimal,iEH,iTarg,1:optimD)=eigenvalues(1:optimD)./sum(eigenvalues); %   same across S/L/SL
    disp('PVAF common across S/L/SL:');
    disp(num2str(eigenvalues(1:optimD)'./sum(eigenvalues)));
    dPSTH=load(['PSTH_prior_' animalNm{idAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D nTrialExc tBuffer
    for iPr=1:nPr
        idD0=((iEH-1)*nEH+iTarg)+(iPr-1)*(nEH*nTarg*nTs);% 1234 + 20*(iPr-1)
        idD1=iPr*(nEH*nTarg*nTs); % 20 40 60
        
        [Dtmp,proj_matrixTmp,keep_neuronsTmp,eigenvalues]=PSTH2traj(dPSTH.D(idD0:(nEH*nTarg):idD1),...
            proj_matrix,keep_neurons,binSize,smthWidth,[],use_sqrt); %  proj_matrix,keep_neurons optimD
        PVAF(idAnimal,iEH,iTarg,iPr,1:optimD)=eigenvalues(1:optimD)./sum(eigenvalues);
%         if sum(PVAF(idAnimal,iEH,iTarg,iPr,1:optimD)>0.5)
%             disp('too high PVAF?');
%         end
        disp('PVAF separate across S/L/SL:');
        disp(num2str(eigenvalues(1:optimD)'./sum(eigenvalues)));
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
%         curv=cell(nPr,1);
    elseif strfind(nm,'ts')
        idProd=0; % 1; % plotting together
        idMeas=1; % 0; % 1;
        idAtt=0;        
%         curv=cell(nPr+1,1);
    elseif strfind(nm,'prior')
        idProd=0; % 1; % plotting together
        idMeas=1; % 0; % 1;
        idAtt=0; % 1; % 0; % 1; % 0; % 1;
    end    
    idCondSpecific=true; % false;
    conditionSpecificPCANm={'condSpecPCA'};
    
    %% curvature of each PC after avgAtt
    tmpDS=fillNan(D(1:nTs));
    tmpDcat=cat(3,tmpDS.data); % [nPC time ts]
    tmpDcurv{1}=nanmean(tmpDcat,3); % [nPC time]
    
    tmpDL=fillNan(D((nTs+1):(2*nTs)));
    tmpDcat=cat(3,tmpDL.data); % [nPC time ts]
    tmpDcurv{2}=nanmean(tmpDcat,3); % [nPC time]
    
     tmpDcat=cat(3,D((2*nTs+1):end).data); % [nPC time ts]
     tmpDcurv{3}=nanmean(tmpDcat,3); % [nPC time]
     
     if idDebugCurv,figure(999); set(gcf,'position',[55 8 1057 797]); end
     for iDim=1:optimD
         for iPr=1:nPr
             dataCurv=[binSize*(0:(size(tmpDcurv{iPr},2)-1)); tmpDcurv{iPr}(iDim,:)];
                 pQuad=polyfit(dataCurv(1,:),dataCurv(2,:),2); % 2nd order
                 cTmp2=abs(pQuad(1));
                 cTmp=estCurvature(dataCurv);
             if idUseQuadFit
                 curv(idAnimal,iEH,iTarg,iPr,1:length(cTmp),iDim)=cTmp2; % [222 x prior x time x dim]
             else
                 curv(idAnimal,iEH,iTarg,iPr,1:length(cTmp),iDim)=cTmp; % [222 x prior x time x dim]
             end
             if idDebugCurv, subplot(optimD,nPr,(iDim-1)*nPr+iPr); ha; 
                 plotyy(dataCurv(1,:),dataCurv(2,:),dataCurv(1,1+[1:length(cTmp)]),cTmp(:)'); 
                plot(dataCurv(1,:),polyval(pQuad,dataCurv(1,:)),'--');
                title([num2str(nanmedian(cTmp(:)),2) num2str(pQuad(1),2)]);
             end
         end % for iPr=1:nPr
     end % for iDim=1:optimD
     if idDebugCurv, waitforbuttonpress; 
         close; end
    
    %% avg attrition
    if idAtt
        %         if strfind(nm,'prior')
        %             % filling NaN at the end
        %             D=fillNan(D);
        %             tmpDcat=cat(3,D.data); % [nPC time nCond]
        %             meanTmpD=nanmean(tmpDcat,3); % [nPC time]
        %             for i=1:length(D) % 5 for prior, 10 for the others
        %                 D(i).data=meanTmpD(:,~isnan(D(i).data(1,:)));
        %             end % for i=1:length(D) % 5 for prior, 10 for the others
        %
        %         else % separately for two priors
        % short
        % filling NaN at the end
        tmpDS=fillNan(D(1:nTs));
        tmpDcat=cat(3,tmpDS.data); % [nPC time nCond]
        meanTmpD=nanmean(tmpDcat,3); % [nPC time]
        for i=1:nTs % 5 for prior, 10 for the others
            D(i).data=meanTmpD(:,~isnan(D(i).data(1,:)));
        end % for i=1:length(D) % 5 for prior, 10 for the others
        
        % long
        % filling NaN at the end
        tmpDL=fillNan(D((nTs+1):end));
        tmpDcat=cat(3,tmpDL.data); % [nPC time nCond]
        meanTmpD=nanmean(tmpDcat,3); % [nPC time]
        for i=(nTs+1):length(D) % 5 for prior, 10 for the others
            D(i).data=meanTmpD(:,~isnan(D(i).data(1,:)));
        end % for i=1:length(D) % 5 for prior, 10 for the others
        
        
        %         end
        
    end
    
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
    
    if idMeas
        %% measurement only
        %         hFig=figure;setFigPos(1,1);% set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+1)]));hold all;
        %         hFig2=figure;setFigPos(1,2);
        %         hFig3=figure;setFigPos(1,3);
        
        if idDebug
            figure; grid on; % debug
        end
        
        for i=1:length(tmpD)
            if length(tmpD(i).epochStarts)>1
                iSet=tmpD(i).epochStarts(2);
            else
                iSet=size(tmpD(i).data,2);
            end
            idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2 3 3 3 3 3
            iTs=rem(i-1,5)+1;
            
            if idAtt % with attrition
                if i~=1 & i~=(nTs+1) % ...
                    %                         & i~=(nTs*nPr+1) & i~=nTs*nPr+nTs+1 % with                        attrition to connect % for split tp
                    endPre=size(tmpD(i-1).data,2);
                    
                    tmpX=[tmpD(i-1).data(1,endPre) tmpD(i).data(1,(endPre+1):iSet)];
                    tmpY=[tmpD(i-1).data(2,endPre) tmpD(i).data(2,(endPre+1):iSet)];
                    tmpZ=[tmpD(i-1).data(3,endPre) tmpD(i).data(3,(endPre+1):iSet)];
                else
                    endPre=1;
                    tmpX=tmpD(i).data(1,1:iSet); % PC1
                    tmpY=tmpD(i).data(2,1:iSet);
                    tmpZ=tmpD(i).data(3,1:iSet);
                end
            else
                endPre=1;
                tmpX=tmpD(i).data(1,1:iSet);
                tmpY=tmpD(i).data(2,1:iSet);
                tmpZ=tmpD(i).data(3,1:iSet);
                % PC456
                %                 tmpX2=tmpD(i).data(4,1:iSet);
                %                 tmpY2=tmpD(i).data(5,1:iSet);
                %                 tmpZ2=tmpD(i).data(6,1:iSet);
                
                %                 plot3(tmpX,tmpY,tmpZ);ha;
                
            end % with attrition
            %             if i<=length(tmpD)/2 % bias+ % for split tp
            
            %% speed
            tmp=tmpD(i).data(1:dim,1:iSet); % [dim x time]
            if idDebug,            ha; plot3(tmp(1,:),tmp(2,:),tmp(3,:),'o-','color',tmpCmap{idPr}(iTs,:)); end % debug
            
            tmpV=sqrt(sum(diff(tmp,1,2).^2,1));       % [1 x (time-1)]
            nV=length(tmpV);
            
            v(idAnimal,iEH,iTh,idPr,iTs,1:nV)=tmpV; % [2HG x 2EH x 2RL x 2SL x 5ts x 60timePoints]
            %             disp(squeeze(v(idAnimal,iEH,iTh,idPr,iTs,:)))
            vv(idAnimal,iEH,iTh,idPr,iTs,1:nV,1:dim)=diff(tmp,1,2)';
            
            
            
        end % length(tmpD)
        
        
        
    end %  if idMeas
    
    %     % 2) check angle/magnitude b/t readout vectors across ts, prior: 2D cmap, if <90, parellel-not state dependent
    %     for iSV=1:size(sv,1) % # ts
    %         for jSV=1:size(sv,1)
    %             if iSV==jSV
    %                 asv(idAnimal,iEH,iTh,iSV,jSV)=0;
    %             else
    %                 asv(idAnimal,iEH,iTh,iSV,jSV)=acosd(dot(sv(iSV,:)',sv(jSV,:)'));
    %             end
    %             if ~isreal(asv(idAnimal,iEH,iTh,iSV,jSV)) % acosd(dot(sv(iSV,:)',sv(jSV,:)')))
    %                 disp(asv(idAnimal,iEH,iTh,iSV,jSV)); % acosd(dot(sv(iSV,:)',sv(jSV,:)')));
    %             end
    %         end
    %     end
    if idDebug,
        disp(squeeze(asv(idAnimal,iEH,iTh,:,:)));
    end
end % for i=1:length(d)

% if idAngleReadout==1
    % avrv=nan(nAnimal,nEH,nTarg,nPr,nT); % angle b/t readout vector and velocity vector  [2HG x 2EH x 2RL x 2SL x 20timePoints]
    for i=1:nAnimal
        for j=1:nEH
            for k=1:nTarg
                for iPr=1:nPr
                    vvTmp=squeeze(nanmean(vv(i,j,k,iPr,:,:,:),5)); % [20timePoints x dim]
                    vvTmp=vvTmp(~isnan(vvTmp(:,1)),~isnan(vvTmp(1,:)));
                    if idUseLowDimAngle
                        vvTmp=vvTmp(:,1:dAngle);
                    end
                    
                    % - define readout(480to800) vector
%                     rv=sum(vvTmp,1); % [1 x dim]    % only shortest to longest
                    tmpIdVV=round(linspace(1,nnz(~isnan(vvTmp(:,1))),nTspp)); % 1 5 9 13 17/1 6 11 16 21
                    rv=mean([sum(vvTmp(tmpIdVV(1):tmpIdVV(end),:),1);... % shortest to longest
                        sum(vvTmp(tmpIdVV(1+1):tmpIdVV(end-1),:),1)],1); % 2th to 4th

                    for iT=1:size(vvTmp,1) % across time
                        avrv(i,j,k,iPr,iT)=angleVectors(rv(:),vvTmp(iT,:)');
                    end % iT
                end % prior
            end % animal
        end % EH
    end % target


    
% else % idAngleReadout==1
    % angle b/t successive vectors
    %     avv=nan(nAnimal,nEH,nTarg,nPr,nT,nT);
    %  vv(idAnimal,iEH,iTh,idPr,iTs,1:nV,1:dim)=diff(tmp,1,2)';
    for i=1:nAnimal
        for j=1:nEH
            for k=1:nTarg
                iCond= (j-1)*nTarg+k;
                for iPr=1:nPr
                    % averaging vv across ts
                    if idExcB4pr
                        vvTmp=squeeze(nanmean(vv(i,j,k,iPr,:,nExcB4pr:end,:),5)); % [20timePoints x dim]
                    else
                        vvTmp=squeeze(nanmean(vv(i,j,k,iPr,:,:,:),5)); % [20timePoints x dim]
                    end
                    vvTmp=vvTmp(~isnan(vvTmp(:,1)),~isnan(vvTmp(1,:)));
                    if idUseLowDimAngle
                        vvTmp=vvTmp(:,1:dAngle);
                    end
                    for iT=1:size(vvTmp,1)
                        for jT=1:size(vvTmp,1)
                            if iT~=jT
                                avv(i,j,k,iPr,iT,jT)=angleVectors(vvTmp(iT,:)',vvTmp(jT,:)');
                            else
                                avv(i,j,k,iPr,iT,jT)=0;
                            end
                        end % jT
                    end % iT
                end % prior
            end % animal
        end % EH
    end % target
    
% end

%         condId=condId+1;
%         if idCondSpecific, break; end;
%         end %   if ~idMeasOnly & idMeasProd
%         if idCondSpecific, break; end;
%     end % targ
%     if idCondSpecific, break; end;
% end % EH

%% plot
% v=nan(nAnimal,nEH,nTarg,nPr,nTs,nT); % speed [2HG x 2EH x 2RL x 2SL x 5ts x 20timePoints]
% vv=nan(nAnimal,nEH,nTarg,nPr,nTs,nT,maxDim); % velocity vector [2HG x 2EH x 2RL x 2SL x 5ts x 20timePoints x dim]
% avv=nan(nAnimal,nEH,nTarg,nPr,nT,nT); % angle b/t velocity vector [2HG x 2EH x 2RL x 2SL x 20timePoints x 20timePoints]

hFig1=figure(1); ha;setFigPos(1,1); % speed vs time (avg. across ts) two lines for short vs long
hFig2=figure(2); ha;setFigPos(1,2); % avv
hFig3=figure(3); ha;setFigPos(1,3); % avrv
% figure(3); setFigPos(1,3); % avv colormap for each data set: short
% figure(4); setFigPos(2,3); % smae for long
hFig4=figure(4); ha;setFigPos(2,1); % curvature vs PVAF
% hFig5=figure(5); ha;setFigPos(2,2); % check curvature distribution across time

condNm={'ER','EL','HR','HL'};
nCmap=length(animalNm)*length(condNm);
cmap3=parula(nCmap); % hsv(nCmap);
icmap=1;

tmpMS=[];tmpML=[];tmpMSL=[];
AS=[]; AL=[];ASL=[];
AS2=[]; AL2=[];ASL2=[];

curvPVAFc=cell(nPr,1); % for linear regression

for i=1:nAnimal
    for j=1:nEH
        for k=1:nTarg
            iCond= (j-1)*nTarg+k;
            
            % speed
            figure(1);
            %             subplot(nEH*nTarg,nAnimal,iCond*nAnimal-(nAnimal-i)); ha;
            %             title([animalNm{i} ':' ehNm{j} targNm{k}]);
            for l=1:nPr
                if idExcB4pr
                    tmpM=nanmean(squeeze(v(i,j,k,l,:,nExcB4pr:end)),1); tmpM=tmpM(~isnan(tmpM));% avg across ts
                    tmpS=sem(squeeze(v(i,j,k,l,:,nExcB4pr:end)),1,'omitnan'); tmpS=tmpS(~isnan(tmpS));% nansem(squeeze(v(i,j,k,l,:,:)),1);
                    tMin=T{l}(1)-binSize*(nExcB4pr-1); % (size(tmpD(1+nTspp*(l-1)).data,2)-1); % 420 for short, 740 for long
                else
                    tmpM=nanmean(squeeze(v(i,j,k,l,:,:)),1); tmpM=tmpM(~isnan(tmpM));% avg across ts
                    tmpS=sem(squeeze(v(i,j,k,l,:,:)),1,'omitnan'); tmpS=tmpS(~isnan(tmpS));% nansem(squeeze(v(i,j,k,l,:,:)),1);
                    if l==3 % SIL
                        tMin=T{1}(1); 
                    else
                        tMin=T{l}(1);
                    end
                    if idIncShortInLong==2, tMin=tMin-1; end % if using whole ts, start from bin23[440-460]
                end
                tmpT=linspace(tMin,T{l}(end),length(tmpM)); % speed=diff (2 for shortest ts
                    % linspace(tMin,T{l}(end)-binSize,length(tmpM)); % ([1:length(tmpM)]-1)*binSize;
                plot(tmpT,tmpM(:)','-','color',cmapPr2(l,:),'linewidth',lw2); % for each data set; cmapHG(i,:)
                %                 errorbar(tmpT,tmpM(:)',tmpS(:)','o-','color',cmap{l},'linewidth',lw,'markerfacecolor',cmap{l},'markersize',msize);
                if l==1 % short
                    tmpMS=[tmpMS; tmpM(:)'];
                    %                     set(gca,'xticklabel',[],'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'ticklength',[0.03 0.03],'xlim',[T{1}(1) T{2}(end)]);
                elseif l==2 % long
                    tmpML=[tmpML; tmpM(:)'];
                    %                     set(gca,'xtick',[T{1}(1:2:end) T{2}(3:2:end)],'ticklength',[0.03 0.03],'xlim',[T{1}(1) T{2}(end)]);
                else % short in long
                    tmpMSL=[tmpMSL; tmpM(:)'];
                end
            end
            axis tight;
            if icmap==nCmap
                % avg across data sets
                for l=1:nPr
                    if l==1,tmpM=nanmean(tmpMS,1);tmpS=sem(tmpMS,1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    elseif l==2, tmpM=nanmean(tmpML,1); tmpS=sem(tmpML,1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);end;
                    else tmpM=nanmean(tmpMSL,1); tmpS=sem(tmpMSL,1,'omitnan'); end;
                    if idExcB4pr
                        tMin=T{l}(1)-binSize*(nExcB4pr-1); %(size(tmpD(1+nTspp*(l-1)).data,2)-1); %420 for short, 740 for long
                    else
                        if l==3 % SIL
                            tMin=T{1}(1);
                        else
                            tMin=T{l}(1);
                        end
                        if idIncShortInLong==2, tMin=tMin-1; end % if using whole ts, start from bin23[440-460]
                    end
                    tmpT=linspace(tMin,T{l}(end),length(tmpM)); % ([1:length(tmpM)]-1)*binSize;
                    shadedErrorBar(tmpT,tmpM(:)',tmpS(:)',{'-','color',cmapPr(l,:),'linewidth',lw},1); % tmpCmap{l,1}(1,:)
                    for iT=1:length(T{l})
                        if l==3
                            idTs=tmpT==((length(T{1})~=iT)*T{1}(iT)+(length(T{1})==iT)*(T{1}(iT))); % -binSize
                        else
                            idTs=tmpT==((length(T{l})~=iT)*T{l}(iT)+(length(T{l})==iT)*(T{l}(iT))); % -binSize
                        end
                        plot(tmpT(idTs),tmpM(idTs),'o',...
                            'markerfacecolor',tmpCmap{l,1}(iT,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2);
                    end
                end
                set(gca,'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xlim',[T{1}(1)-binSize T{2}(end)+binSize],'xticklabel',{T{1}(1);[];T{1}(end);[];T{2}(end)},...
                    'ylim',[0 0.2],'ytick',0:0.1:0.2);
            end %  if icmap==nCmap
            if idExpFig & i==nAnimal & j==nEH & k==nTarg
                drawnow;
                if idIncShortInLong==1
                    savefig(hFig1,fullfile(figDir,'4link',['priorSpeed.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                    xlabel([]);ylabel([]);
                    exportfig(hFig1,fullfile(figDir,'4link',['priorSpeed.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
                elseif idIncShortInLong==2
                     savefig(hFig1,fullfile(figDir,'4link',['priorSpeed_ts.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                    xlabel([]);ylabel([]);
                    exportfig(hFig1,fullfile(figDir,'4link',['priorSpeed_ts.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
                else
                    savefig(hFig1,fullfile(figDir,'4link',['priorSpeed_supportOnly.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                    xlabel([]);ylabel([]);
                    exportfig(hFig1,fullfile(figDir,'4link',['priorSpeed_supportOnly.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
                end
            end
            
            %%
            % avv=nan(nAnimal,nEH,nTarg,nPr,nT,nT); % angle b/t velocity vector [2HG x 2EH x 2RL x 2SL x 20timePoints x 20timePoints]
            figure(2);
            
            %             subplot(nEH*nTarg,nAnimal,iCond*nAnimal-(nAnimal-i)); ha;
            %             title([animalNm{i} ':' ehNm{j} targNm{k}]);
            for l=1:nPr
                tmpData=squeeze(avv(i,j,k,l,:,:)); % [nTs x nTs]
                vB=0:(sum(~isnan(tmpData(1,:)))-1); % distance difference
                distM=toeplitz(vB);
                
                tmpData=tmpData(~isnan(tmpData(:,1)),~isnan(tmpData(1,:)));
                
                A=[]; % sorted by angle
                % mean across pair9
                for iDist=1:length(vB)
                    A=[A nanmean(tmpData(distM==vB(iDist)))];
                end
                plot(vB*binSize,A,'-','color',cmapPr2(l,:),'linewidth',lw2); % cmapHG(i,:)
                if l==1
                    AS=[AS; A(:)'];
                elseif l==2
                    AL=[AL; A(:)'];
                else % SIL
                    ASL=[ASL; A(:)'];
                end
            end % pr
            axis tight;
            if icmap==nCmap
                % avg across data sets
                for l=1:nPr
                    if l==1,tmpM=nanmean(AS,1);tmpS=sem(AS,1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    elseif l==2, tmpM=nanmean(AL,1); tmpS=sem(AL,1,'omitnan');  % nansem(squeeze(v(i,j,k,l,:,:)),1);end
                    else tmpM=nanmean(ASL,1); tmpS=sem(ASL,1,'omitnan'); end
                    tmpT=binSize*(0:(size(tmpM,2)-1));
                    
                    shadedErrorBar(tmpT,tmpM(:)',tmpS(:)',{'-','color',cmapPr(l,:),'linewidth',lw},1); % tmpCmap{l,1}(1,:)
                end
                dtTmp=T{1}(2)-T{1}(1);
                set(gca,'tickdir','out','xtick',0:dtTmp:sum(diff(T{2})),'xticklabel',{'0','','160','','320'},'ytick',[0 45 90 130]);
            end %  if icmap==nCmap
            if idExpFig & i==nAnimal & j==nEH & k==nTarg
                drawnow;
                if idIncShortInLong==1
                savefig(hFig2,fullfile(figDir,'4link',['priorAngle.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                xlabel([]);ylabel([]);
                exportfig(hFig2,fullfile(figDir,'4link',['priorAngle.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
                elseif idIncShortInLong==2
                savefig(hFig2,fullfile(figDir,'4link',['priorAngle_ts.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                xlabel([]);ylabel([]);
                exportfig(hFig2,fullfile(figDir,'4link',['priorAngle_ts.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);                    
                else
                                    savefig(hFig2,fullfile(figDir,'4link',['priorAngle_supportOnly.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                xlabel([]);ylabel([]);
                exportfig(hFig2,fullfile(figDir,'4link',['priorAngle_supportOnly.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);

                end
            end
            
            %             for l=1:nPr
            %                 figure(3+(l-1)); % % avv=nan(nAnimal,nEH,nTarg,nPr,nT,nT); % angle b/t velocity vector [2HG x 2EH x 2RL x 2SL x 20timePoints x 20timePoints]
            %                 subplot(nEH*nTarg,nAnimal,iCond*nAnimal-(nAnimal-i)); ha;
            %                 title([animalNm{i} ':' ehNm{j} targNm{k}]);
            %                 tmpData=squeeze(avv(i,j,k,l,:,:));
            %                 imagesc(tmpData(1:(sum(~isnan(tmpData(1,:)))),1:(sum(~isnan(tmpData(1,:)))))); colorbar;
            %                 axis tight;
            %             end
            %
            %
            % %             if iCond==3 & i==1
            % %                 figure(1);                ylabel('speed (a.u.)');
            % %                 figure(4);                ylabel('norm of 2nd derivative (a.u.)');
            % %                 figure(5);                ylabel('curvature (a.u.)');
            % %                 figure(2);                ylabel('travelling distance from ready (a.u.)');
            % %                 figure(3);                ylabel('travelling distance during prior (a.u.)');
            % %
            % %                 figure(9); ylabel('distance across priors');
            % %                 figure(10); ylabel(['angle of time vectors [' num2str(t0sv*binSize) ',' num2str(t1sv*binSize)  ']']);
            % %                 figure(11); ylabel(['|time vectors [' num2str(t0sv*binSize) ',' num2str(t1sv*binSize)  ']|']);
            % %
            % %
            % %             end
            
                        %%
            % avrv=nan(nAnimal,nEH,nTarg,nPr,nT); % angle b/t readout vector and velocity vector  [2HG x 2EH x 2RL x 2SL x 20timePoints]
            figure(3);
            
            for l=1:nPr
                tmpData=squeeze(avrv(i,j,k,l,:)); % [nTs x 1]
                
                plot([0:(length(tmpData)-1)]*binSize,tmpData,'-','color',cmapPr2(l,:),'linewidth',lw2); % cmapHG(i,:)
                if l==1
                    AS2=[AS2; tmpData(:)'];
                elseif l==2
                    AL2=[AL2; tmpData(:)'];
                else % SIL
                    ASL2=[ASL2; tmpData(:)'];
                end
            end % pr
            axis tight;
            if icmap==nCmap
                % avg across data sets
                for l=1:nPr
                    if l==1,tmpM=nanmean(AS2,1);tmpS=sem(AS2,1,'omitnan'); % nansem(squeeze(v(i,j,k,l,:,:)),1);
                    elseif l==2, tmpM=nanmean(AL2,1); tmpS=sem(AL2,1,'omitnan');  % nansem(squeeze(v(i,j,k,l,:,:)),1);end
                    else tmpM=nanmean(ASL2,1); tmpS=sem(ASL2,1,'omitnan'); end
                    tmpT=binSize*(0:(size(tmpM,2)-1));
                    
                    shadedErrorBar(tmpT,tmpM(:)',tmpS(:)',{'-','color',cmapPr(l,:),'linewidth',lw},1); % tmpCmap{l,1}(1,:)
                end
                dtTmp=T{1}(2)-T{1}(1);
                set(gca,'tickdir','out','xtick',0:dtTmp:sum(diff(T{2})),'xticklabel',{'0','','160','','320',''},'ytick',[0 45 90 130]);
            end %  if icmap==nCmap
            if idExpFig & i==nAnimal & j==nEH & k==nTarg
                drawnow;
                if idIncShortInLong==1
                savefig(hFig2,fullfile(figDir,'4link',['priorAngle2.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                xlabel([]);ylabel([]);
                exportfig(hFig2,fullfile(figDir,'4link',['priorAngle2.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
                elseif idIncShortInLong==2
                savefig(hFig2,fullfile(figDir,'4link',['priorAngle2_ts.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                xlabel([]);ylabel([]);
                exportfig(hFig2,fullfile(figDir,'4link',['priorAngle2_ts.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);                    
                else
                                    savefig(hFig2,fullfile(figDir,'4link',['priorAngle2_supportOnly.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                xlabel([]);ylabel([]);
                exportfig(hFig2,fullfile(figDir,'4link',['priorAngle2_supportOnly.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);

                end
            end
            
            
            %% 
%             figure(5); % check distribution of curvature across time
%             for iPr=1:nPr
%                 subplot(nPr,1,iPr);ha;
%                 for iDim=1:nnz(~isnan(curv(i,j,k,iPr,1,:)))
%                     histogram(curv(i,j,k,iPr,:,iDim),'DisplayStyle','stairs');
%                 end
%             end
%             waitforbuttonpress; clf;
            
            figure(4);            
%             % individual PC's pvaf & curvature
% PVAF=nan(nAnimal,nEH,nTarg,maxDim); % common across S/L/SL
% curv=nan(nAnimal,nEH,nTarg,nPr,nT,maxDim); % [222 x prior x time x dim]
            
            for iPr=1:nPr
                nDim=nnz(~isnan(PVAF(i,j,k,iPr,:)));
                tmpX=PVAF(i,j,k,iPr,1:nDim);
                if idMed
                    mCurv=squeeze(nanmedian(curv(i,j,k,iPr,:,1:nDim),5)); % avg time
                else
                    mCurv=squeeze(nanmean(curv(i,j,k,iPr,:,1:nDim),5)); % avg time
                end
                if sum(mCurv>9*10^(-5)),
                    disp('outlier');
                end
                plot(tmpX(:),mCurv,'o','color',cmapPr(iPr,:),'linewidth',lw,'markersize',msize,'markerfacecolor','w');
                curvPVAFc{iPr}=[curvPVAFc{iPr}; tmpX(:),mCurv(:)];
            end
                
            if icmap==nCmap
                axis tight;
                for iPr=1:nPr
                     s=plotReg(curvPVAFc{iPr}(:,1),curvPVAFc{iPr}(:,2),hFig4,cmapPr(iPr,:));
                     disp(num2str([s.beta(:)' s.tstat.t(:)' s.tstat.pval(:)' s.rsquare],2));
                end
                xlabel('% var. explained'); ylabel('curvature');
            end


            icmap=icmap+1;
            
        end
    end
end


% % arc length comparison
% % longest of short vs long: similar
% % overlap: longer for short
%
%             condNm={'ER','EL','HR','HL'};
%             nCmap=length(animalNm)*length(condNm);
%             cmap=parula(nCmap); % hsv(nCmap);
%              icmap=1;
%
%              figure(1000);ha; setFigPos(2,4); xlabel('trajectory length (t_s=800,short)');ylabel('trajectory length (t_s=1200,long)');
%              figure(1001);ha; setFigPos(2,5);xlabel('trajectory length (t_s=800,short)');ylabel('trajectory length (t_s=800,long)');
% for i=1:nAnimal
%     for j=1:nEH
%         for k=1:nTarg
%             figure(1000);
%             % longest of short vs long: similar
%             plot(squeeze(x(i,j,k,1,end,1)),squeeze(x(i,j,k,2,end,1)),'o','color',cmap(icmap,:),'linewidth',lw,'markerfacecolor',cmap(icmap,:),'markersize',msize);
%
%             figure(1001);
%             % overlap: longer for short
%             plot(squeeze(x(i,j,k,1,end,1)),squeeze(x(i,j,k,2,1,1)),'o','color',cmap(icmap,:),'linewidth',lw,'markerfacecolor',cmap(icmap,:),'markersize',msize);
%              icmap= icmap+1;
%             end % for i=1:nAnimal
%         end % for j=1:nEH
%     end % for k=1:nTarg
%      figure(1000);axis tight; plotIdentity(gca);
% figure(1001);axis tight; plotIdentity(gca);