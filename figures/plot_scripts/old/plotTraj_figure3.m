function varargout=plotTraj_figure3(fname,varargin)

% 2018/5/8
% modified to deal with D(DataHigh format) as varargin

% RNN
% cc; plotTraj_figure3('traj_ts_Weber_2Prior_Type2_000_v3.mat');
% cc; plotTraj_figure3('traj_ts_Weber_2Prior_Type2_000_v4.mat');
% cc; plotTraj_figure3('traj_ts_Weber_2Prior_Type2_005_v2.mat');
% cc; plotTraj_figure3('traj_ts_Weber_2Prior_Type2_005_v5.mat');

% for COSYNE
% plotTraj_figure3('trajKS_ts_EL_H_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_tp_EL_H_bin20_smth40.mat'); 

% PCA done separately for meas/prod
% plotTraj_figure3('trajKS_ts_EL_H_bin20_smth40.mat'); cc; % chosen for representative
% plotTraj_figure3('trajKS_ts_ER_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_ts_HL_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_ts_HR_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_ts_EL_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_ts_ER_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_ts_HL_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_ts_HR_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_tp_EL_H_bin20_smth40.mat'); cc; % chosen for representative
% plotTraj_figure3('trajKS_tp_ER_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_tp_HL_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_tp_HR_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_tp_EL_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_tp_ER_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_tp_HL_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_tp_HR_G_bin20_smth40.mat'); cc;

% plotTraj_figure3('trajKS_periSet_EL_H_bin20_smth40.mat'); cc; % chosen for representative
% plotTraj_figure3('trajKS_periSet_ER_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_periSet_HL_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_periSet_HR_H_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_periSet_EL_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_periSet_ER_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_periSet_HL_G_bin20_smth40.mat'); cc;
% plotTraj_figure3('trajKS_periSet_HR_G_bin20_smth40.mat'); cc;

% not including short in long
% plotTraj_figure3('trajKS_prior_EL_H_supportOnly.mat'); 
% plotTraj_figure3('trajKS_prior_ER_H_supportOnly.mat'); 
% plotTraj_figure3('trajKS_prior_HL_H_supportOnly.mat'); 
% plotTraj_figure3('trajKS_prior_HR_H_supportOnly.mat'); 
% plotTraj_figure3('trajKS_prior_EL_G_supportOnly.mat'); 
% plotTraj_figure3('trajKS_prior_ER_G_supportOnly.mat'); 
% plotTraj_figure3('trajKS_prior_HL_G_supportOnly.mat'); 
% plotTraj_figure3('trajKS_prior_HR_G_supportOnly.mat'); 

% including short in long
% plotTraj_figure3('trajKS_prior_EL_H_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_prior_ER_H_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_prior_HL_H_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_prior_HR_H_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_prior_EL_G_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_prior_ER_G_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_prior_HL_G_bin20_smth40.mat'); 
% plotTraj_figure3('trajKS_prior_HR_G_bin20_smth40.mat'); 

% old
% plotTraj_figure3('trajKS_periSet_E2_H_bin20_smth40.mat');
% plotTraj_figure3('trajKS_periSet_H2_H_bin20_smth40.mat');

% difference from plotTraj
%   noSet

if isempty(fname)
    fname='trajKS_periSet_E2_H_bin20_smth40.mat';
end

% plotTraj(fname)

% plot state space trajectory for each condition (of EH and targets) & two epoches
% data in DataHigh format, 40 conditions (2 pr x 5 ts x 2 EH x 2 targets)
% input: filename (assuming in /Users/hansem/Dropbox (MIT)/psthDataHigh/)
%   e.g. plotTraj('traj_periSet_G_attrition.mat'); plotTraj('traj_periSet_H_attrition.mat')
% if fname contains '_attrition', attrition
% assuming two epoches exist (epochStart(2)) as in periset
% save both in png and fig formats
% e.g. plotTraj('traj_periSet_ER_H.mat')

%% update data after kilosort @9/1/2017
% cc;load PSTH_periSet_H_SUMU.mat % runExportDataHigh_CCN17_woOutlier.m
% DataHigh(D,'DimReduce'); % 20ms bin; 1 sp threshold; PCA, 1:25 dimension; smthKernel: 25ms; optimD:
%           plotTraj('trajKS_periSet_H.mat');
%           plotTraj('trajKS_periSet_H_SUkilosort.mat')
%           plotTraj('trajKS_periSet_G.mat');                 % N.B. kilosort not finished

%% initial
initRSG2prior;
idExpFig=0; % 1; % 0; % 1;
idKilosort=0; % 1;

idCosyne=0; % 1;

cd(psthDir);
if isempty(varargin)
    try
        load(fname); % dimReduce_periSet4_splitTp_attrition.mat'); % load('dimReduce_periSet4_attrition.mat'); % _eyeRightAttrition _attrition
    catch
        load(fullfile(singleTDir,fname));
        %     load(['./sess/' fname]);
    end
else
    D=varargin{1};
end

if strfind(fname,'periSet') % trajKS_periSet_ER_H
    idProd=1; % 0; % 1;
    idMeas=1; % 0;
    idAtt=0; % 1; % 1; % 0;
elseif strfind(fname,'ts')
    idProd=0; % 1; % plotting together
    idMeas=1; % 0; % 1;
    idAtt= 0; % 1;
elseif strfind(fname,'tp')
    idProd=1; % plotting together
    idMeas=0; % 1;
    idAtt=0;
elseif ~isempty(strfind(fname,'prior')) | ~isempty(strfind(fname,'_pr_'))
    idProd=0; % 1; % plotting together
    idMeas=1; % 0; % 1;
    idAtt=1;
%     %     if strfind(fname,'priorSL')
%     %         idAtt=0;
%     %     else
%     %     idAtt=1;
%     %     end
else
     idProd=1; % 1; % plotting together
    idMeas=0; % 0; % 1;
    idAtt= 0; % 1;
end

if idKilosort
    fnmEnd='_Kilosort';
else
    fnmEnd=[];
end
idFlipPC3=0; % 1;
idCondSpecific=true; % false;
conditionSpecificPCANm={'condSpecPCA'};

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

if idCosyne
    % figure's PaperUnits: 'inches'
    optsExpFig.Height=22/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
    optsExpFig.Width=22/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
    optsExpFig.FontSize='24';
    optsExpFig.FontMode='fixed';
    optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
    optsExpFig.LockAxes=1;
    optsExpFig.LineMode='scaled';
    optsExpFig.LineWidthMin=0.5;
    optsExpFig.LineWidthMax=1.5;
    optsExpFig.Renderer='painters';
    
    lw=2; % trajectory
    lw2=1; % marker
    msize=14; % 5; % 10; % 4;
    msize2=6; % for speed markers
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
    
    lw=1; % 1.5;
    lw2=1;
    msize=5; % 10; % 4;
    msize2=1;
end

nPr=2;    prNm={'Short','Long'};
nEH=2;ehNm={'Eye','Hand'};
nTarg=2; targNm={'Right','Left'}; % [0 180]
nTs=5;
nSplit=2; splitNm={'bias+','bias-'}; % more/less biased; for ts@prior mean, early/late

% identifying attrition & epochs & animal
% idAtt=true;
% % if ~isempty(strfind(fname,'_attrition'))
% %     idAtt=true;
% % else
% %     idAtt=false;
% % end
epNm={'fixation';'targetOn';'periReady';'periSet';'production';'ts';'tp';'_priorS_';'_priorL_';'_priorSL_';'t_s';'t_p';'prior'};
epNm2={'preFix','postFix';'preTOn','postTOn';'preReady','postReady';'','';'prod','postProd';'','';'','';'','';'','';'','';'','';'','';'',''};
for i=1:length(epNm)
    if ~isempty(strfind(fname,epNm{i})), break; end;
end
idEp=i;

animalNm={'H','G'}; % ,'t_s','t_p'};
for i=1:length(animalNm)
    if ~isempty(strfind(fname,['_' animalNm{i} '_'])), break; end;
end
idAnimal=i;

condNm={'ER','EL';'HR','HL';'E2','H2'};
for iEH=1:size(condNm,1)
    for iTh=1:nTarg
        if ~isempty(strfind(fname,condNm{iEH,iTh})),
            idCondSpecific=true; tmpId=[iEH; iTh];
            break;
        end;
    end
end
if ~exist('tmpId','var') % for RNN
    iEH=1;iTh=1;tmpId=[iEH; iTh];
end
% individual sessions' data
if ~isempty(strfind(fname,'_prior_')) % & ~isempty(strfind(fname,'avgDir'))
%     if ~isempty(strfind(D(1).condition,'Eye'))
%         idCondSpecific=true; tmpId=[1; 1];
%     else % Hand
%         idCondSpecific=true; tmpId=[2; 1];
%     end
    tmpCmap{end+1,1}=[0 0 1;0 0 1;0 0 1;0 0 1;0 0 1]; % tmpCmap{1,1};
    linestyle={'-','-','--'};
else
    linestyle={'-','-'};
end % if ~isempty(strfind(fname,'prior')) & ~isempty(strfind(fname,'avgDir'))

%% avg attrition
if idAtt
    if ~isempty(strfind(fname,'_priorS_')) | ~isempty(strfind(fname,'_priorL_')) | ~ isempty(strfind(fname,'_priorSL_'))
        % filling NaN at the end
        D=fillNan(D);
        tmpDcat=cat(3,D.data); % [nPC time nCond]
        meanTmpD=nanmean(tmpDcat,3); % [nPC time]
        for i=1:length(D) % 5 for prior, 10 for the others
            D(i).data=meanTmpD(:,~isnan(D(i).data(1,:)));
        end % for i=1:length(D) % 5 for prior, 10 for the others
        
    else % separately for two priors
        if length(D(1).epochStarts)>1 & ~isempty(strfind(fname,'periSet'))
            % avg w/ attrition only before set
            tmpDS=fillNan(D(1:nTs));
            tmpDcat=cat(3,tmpDS.data); % [nPC time nCond]
            for i=1:nTs % short
                iSet=D(i).epochStarts(2);
                D(i).data(:,1:iSet)=mean(tmpDcat(:,1:iSet,i:end),3); % mean for ts>current
            end
            
            tmpDL=fillNan(D((nTs+1):nTs*2));
            tmpDcat=cat(3,tmpDL.data); % [nPC time nCond]
            for i=(nTs+1):(nTs*2) % long
                iSet=D(i).epochStarts(2);
                D(i).data(:,1:iSet)=mean(tmpDcat(:,1:iSet,(i-nTs):end),3); % mean for ts>current
            end
            
        else % for ts-only
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
            tmpDL=fillNan(D((nTs+1):nTs*2));
            tmpDcat=cat(3,tmpDL.data); % [nPC time nCond]
            meanTmpD=nanmean(tmpDcat,3); % [nPC time]
            for i=(nTs+1):(nTs*2) % 5 for prior, 10 for the others
                D(i).data=meanTmpD(:,~isnan(D(i).data(1,:)));
            end % for i=1:length(D) % 5 for prior, 10 for the others
        end
        
        if strfind(fname,'_prior_')
            if isempty(strfind(fname,'_supportOnly')) & isempty(strfind(fname,'Weber'))
                % S,L,SL together>for SL
                tmpDSL=fillNan(D((nTs*2+1):end));
                tmpDcat=cat(3,tmpDSL.data); % [nPC time nCond]
                meanTmpD=nanmean(tmpDcat,3); % [nPC time]
                for i=(nTs*2+1):length(D) % 5 for prior, 10 for the others
                    D(i).data=meanTmpD(:,~isnan(D(i).data(1,:)));
                end % for i=1:length(D) % 5 for prior, 10 for the others
            end
        end % if strfind(fname,'_prior_')
    end
    
end

% no set

if ~isempty(strfind(fname,'noSet'))
    idNoSet=1;
    linestyle{end+1}='-';
else
    idNoSet=0;
    
end

%% main
condId=1; % eyeRight
for iEH=1:nEH
    for iTh=1:nTarg
        if idCondSpecific
            tmpD=D;iEH=tmpId(1);iTh=tmpId(2);
            if iEH>2, iEH=2; end; % E2, H2 avgDir
        else
            tmpD=D(condId:4:length(D));
        end
        
        tmpMin=[]; % for setting measure/produce in the same range
        tmpMax=[];
        
        diffV=[];
        sSet=[]; % for curvature
        sSetAll=[]; % all data (binSize20) except before 1st set
        
        if idMeas 
%             ha;
            %% measurement only
            hFig=figure;set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+1)]));hold all;
            
            % H
            if idAnimal==1
                if strfind(fname,'periSet')
                    set(gca,'View',[29 60]); % [49 24]);
                    grid on; % off;
                else % ts/tp
                    if strfind(fname,'_prior_')
                        if strfind(fname,'supportOnly')
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[-24    12]); %12 -56]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[-18.0000  -42.4000]); % -10 -60]); % 3 -60]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[   51.0000  -44.4000]); %19 -82]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[ -80.8000   37.2000]); %30 -70]);
                            end
                        else
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[-41 -54]); %12 -56]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[-110    67.6]); % -10 -60]); % 3 -60]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[47 36]); %19 -82]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[34 -39]); %30 -70]);
                            end
                        end
                    elseif strfind(fname,'_ts_')
                        if strfind(fname,'depth') %%%%%%
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[-17.2000   79.2000]); %12 -56]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[-15 -64]); % -10 -60]); % 3 -60]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[39.0000  -75.6000]); %19 -82]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[37.2000  -89.2000]); %30 -70]);
                            end
%                             if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
%                                 set(gca,'View',[138 -78]); %12 -56]);
%                             elseif iEH==1 & iTh==2
%                                 set(gca,'View',[53 80]); % -10 -60]); % 3 -60]);
%                             elseif iEH==2& iTh==1
%                                 set(gca,'View',[17 -84]); %19 -82]);
%                             elseif iEH==2 & iTh==2
%                                 set(gca,'View',[153 -76]); %30 -70]);
%                             end
                        else
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[-17.2000   79.2000]); %12 -56]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[-15 -64]); % -10 -60]); % 3 -60]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[39.0000  -75.6000]); %19 -82]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[37.2000  -89.2000]); %30 -70]);
                            end
                        end
                    end
                end
                
            else % G
                if strfind(fname,'periSet')
                    set(gca,'View',[57 4]);
                    grid on; % off;
                else
                    if strfind(fname,'_prior_')
                        if strfind(fname,'supportOnly')
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[-20.6000    6.8000]); %-9   -34]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[ 16.4000   13.6000]); %-22   -48]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[  -66.6000   27.6000]); %5    74]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[34.8000   15.2000]); %22   -52]);
                            end
                        else
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[33 77]); %12 -56]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[32 76]); % -10 -60]); % 3 -60]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[-10 31]); %19 -82]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[172 19]); %30 -70]);
                            end
                        end
                    elseif strfind(fname,'_ts_')
                        if strfind(fname,'depth') %%%%%%
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[ 7.8000   54.8000]); %12 -56]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[ -203.6000  -48.8000]); % -10 -60]); % 3 -60]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[-32.2000  -62.8000]); %19 -82]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[-69.6000  -76.0000]); %30 -70]);
                            end
%                             if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
%                                 set(gca,'View',[ -211 78]); %12 -56]);
%                             elseif iEH==1 & iTh==2
%                                 set(gca,'View',[ -168 75]); % -10 -60]); % 3 -60]);
%                             elseif iEH==2& iTh==1
%                                 set(gca,'View',[166 78]); %19 -82]);
%                             elseif iEH==2 & iTh==2
%                                 set(gca,'View',[-176 -60]); %30 -70]);
%                             end
                        else
                            if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                                set(gca,'View',[ 7.8000   54.8000]); %12 -56]);
                            elseif iEH==1 & iTh==2
                                set(gca,'View',[ -203.6000  -48.8000]); % -10 -60]); % 3 -60]);
                            elseif iEH==2& iTh==1
                                set(gca,'View',[-32.2000  -62.8000]); %19 -82]);
                            elseif iEH==2 & iTh==2
                                set(gca,'View',[-69.6000  -76.0000]); %30 -70]);
                            end
                        end
                    end
                end % periSet
            end % animal
            
            for i=1:length(tmpD)
                if length(tmpD(i).epochStarts)>1
                    iSet=tmpD(i).epochStarts(2);
                else
                    iSet=size(tmpD(i).data,2);
                end
                idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2 3 3 3 3 3 [for SL]
                
                if idAtt % with attrition
                    if i~=1 & i~=(nTs+1) & ...
                            i~=(2*nTs+1) % ... % SL
                        %                         & i~=(nTs*nPr+1) & i~=nTs*nPr+nTs+1 % with                        attrition to connect % for split tp
                        if length(tmpD(i-1).epochStarts)>1
                            endPre=tmpD(i-1).epochStarts(end); % size(tmpD(i-1).data,2);
                        else
                            endPre=size(tmpD(i-1).data,2);
                        end
                        
                        tmpX=[tmpD(i-1).data(1,endPre) tmpD(i).data(1,(endPre+1):iSet)];
                        tmpY=[tmpD(i-1).data(2,endPre) tmpD(i).data(2,(endPre+1):iSet)];
                        tmpZ=[tmpD(i-1).data(3,endPre) tmpD(i).data(3,(endPre+1):iSet)];
                    elseif i>(2*nTs+1) % no plot for after avg for SL
                        tmpX=[];tmpY=[];tmpZ=[];
                    else
                        tmpX=tmpD(i).data(1,1:iSet); % PC1
                        tmpY=tmpD(i).data(2,1:iSet);
                        tmpZ=tmpD(i).data(3,1:iSet);
                    end
                    sSet=[sSet tmpD(i).data(:,iSet)]; %  1:3,iSet)];
                    
                    if rem(i,nTs)==0 % longest ts for each prior
                        tShortestTs=size(tmpD((idPr-1)*nTs+1).data,2);
                        sSetAll{idPr}=tmpD(i).data(:,tShortestTs:end);
                    end
                    
                else
                    tmpX=tmpD(i).data(1,1:iSet);
                    tmpY=tmpD(i).data(2,1:iSet);
                    tmpZ=tmpD(i).data(3,1:iSet);
                end % with attrition
                %             if i<=length(tmpD)/2 % bias+ % for split tp
                
                plot3(tmpX,tmpY,tmpZ,'-',... % ,linestyle{idPr},...
                    'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw); % lines tmpCmap(i,:)       iEH
                
                % for setting measure/produce in the same range
                tmpMin=min([tmpMin min(tmpD(i).data(1:3,:),[],2)],[],2); % [3 x1]
                tmpMax=max([tmpMax max(tmpD(i).data(1:3,:),[],2)],[],2);
                
            end % for i=1:length(tmpD)
            
            
            % marker
            for i=1:length(tmpD)
                if length(tmpD(i).epochStarts)>1
                    iSet=tmpD(i).epochStarts(2);
                else
                    iSet=size(tmpD(i).data,2);
                end
                idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2 3 3 3 3 3 [for SL]
                
                if idAtt % with attrition
                    if i~=1 & i~=(nTs+1) & ...
                            i~=(2*nTs+1) % ... % SL
                        %                         & i~=(nTs*nPr+1) & i~=nTs*nPr+nTs+1 % with                        attrition to connect % for split tp
                        if length(tmpD(i-1).epochStarts)>1
                            endPre=tmpD(i-1).epochStarts(end); % size(tmpD(i-1).data,2);
                        else
                            endPre=size(tmpD(i-1).data,2);
                        end
                        
                        tmpX=[tmpD(i-1).data(1,endPre) tmpD(i).data(1,(endPre+1):iSet)];
                        tmpY=[tmpD(i-1).data(2,endPre) tmpD(i).data(2,(endPre+1):iSet)];
                        tmpZ=[tmpD(i-1).data(3,endPre) tmpD(i).data(3,(endPre+1):iSet)];
                        
                        
                        
                    elseif i>(2*nTs+1) % no plot for after avg for SL
                        tmpX=[];tmpY=[];tmpZ=[];
                    else
                        tmpX=tmpD(i).data(1,1:iSet); % PC1
                        tmpY=tmpD(i).data(2,1:iSet);
                        tmpZ=tmpD(i).data(3,1:iSet);
                    end
                    if i~=1
                        diffV=[diffV tmpD(i).data(:,iSet)-tmpD(i-1).data(:,endPre)]; % [dim x 1]
                    end
                else
                    tmpX=tmpD(i).data(1,1:iSet);
                    tmpY=tmpD(i).data(2,1:iSet);
                    tmpZ=tmpD(i).data(3,1:iSet);
                end % with attrition
                %             if i<=length(tmpD)/2 % bias+ % for split tp
                
                % marker
                
                plot3(tmpD(i).data(1,iSet),tmpD(i).data(2,iSet),tmpD(i).data(3,iSet),'o',...
                    'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2);
                %                     'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw2); % marker tmpCmap(i,:)          iEH
                if isempty(strfind(fname,'prior'))
                    plot3(tmpD(i).data(1,1),tmpD(i).data(2,1),tmpD(i).data(3,1),'^',...
                        'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2); % marker tmpCmap(i,:)          iEH
                    %                     'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw2); % marker tmpCmap(i,:)          iEH
                end
                

            end % for i=1:length(tmpD)
            
            
            tmp=[tmpMin tmpMax]';
            axis tight; %(tmp(:)'); % tight; %
            if max(abs(tmpMin))<10 & max(abs(tmpMax))<10
                set(gca,'tickdir','out','TickLength', [.02 .02],...
                    'xtick',unique(round(tmpMin(1):1:tmpMax(1))),'ytick',unique(round(tmpMin(2):1:tmpMax(2))),'ztick',unique(round(tmpMin(3):1:tmpMax(3))));
            end
            xlabel('PC1');ylabel('PC2'); zlabel('PC3');
            grid on; 
            

            drawnow;
            
%             if ~idMeasProd
                if idExpFig
                    drawnow;
                    savefig(hFig,fullfile(figDir,'3PC',['PC123_' epNm{idEp} '_readySet_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                    xlabel([]);ylabel([]);zlabel([]);
                    exportfig(hFig,fullfile(figDir,'3PC',['PC123_' epNm{idEp} '_readySet_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.eps']),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
                    
                end
%             end
            
varargout{1}=diffV;
varargout{2}=sSet;
varargout{3}=sSetAll;


        end %  if idMeas
        
        if idProd 
            %% production only
            %          if ~idMeasProd
                    hFig2=figure;set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+2)]));
                    title([condNm{iEH,iTh} ' ' animalNm{idAnimal}]);
            %          end
            % H
            if idAnimal==1
                if strfind(fname,'periSet')
                    set(gca,'View',[20 26]); % [26 38]);
                    grid on; % off;
                else
                    if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                        set(gca,'View',[50    44]);
                    elseif iEH==1 & iTh==2
                        set(gca,'View',[80.4000   50.2000]); % 78 59]); % 70    48]);
                    elseif iEH==2& iTh==1
                        set(gca,'View',[-162    -54]);
                    elseif iEH==2 & iTh==2
                        set(gca,'View',[15    -82]);
                    end
                end
            else
                if strfind(fname,'periSet')
                    set(gca,'View',[12 30]);
                    grid on; % off;
                else
                    if iEH==1 & iTh==1 % {'EyeLeft','EyeRight','HandLeft','HandRight'};
                        set(gca,'View',[12    24]);
                    elseif iEH==1 & iTh==2
                        set(gca,'View',[170 -80]);
                    elseif iEH==2& iTh==1
                        set(gca,'View',[9 24]);
                    elseif iEH==2 & iTh==2
                        set(gca,'View',[-15 74]);
                    end
                end
            end
            
            hold all;
            
            sSet=[];
            sIC=[]; % [length(tmpD) x 3D]
            if ~exist('binSize')
            binSize=20;
            end
            dtIC=200;
            
            tmpMin=[];
            tmpMax=[];
            
            for i=1:length(tmpD)
                if ~idNoSet | i<=nPr*nTs
                idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
                
                if length(tmpD(i).epochStarts)>1
                    iSet=tmpD(i).epochStarts(2);
                else
                    iSet=1; % tmpD(i).epochStarts(2);
                end
%                 sSet=[sSet; tmpD(i).data(1:3,iSet)'];
                
                iIC=iSet+round(dtIC/binSize);
%                 sIC=[sIC; tmpD(i).data(1:3,iIC)'];
               
                    plot3(tmpD(i).data(1,iSet:end),tmpD(i).data(2,iSet:end),tmpD(i).data(3,iSet:end),linestyle{idPr},...
                        'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw); % lines tmpCmap(i,:)       iEH
                    plot3(tmpD(i).data(1,iIC:end),tmpD(i).data(2,iIC:end),tmpD(i).data(3,iIC:end),'o',...
                        'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw2,'markersize',msize2); % lines tmpCmap(i,:)       iEH
                    
%                 end
                tmpMin=min([tmpMin min(tmpD(i).data(1:3,iSet:end),[],2)],[],2); % [3 x1]
                tmpMax=max([tmpMax max(tmpD(i).data(1:3,iSet:end),[],2)],[],2);

                end
            end
            % marker
           
            for i=1:length(tmpD)
                if ~idNoSet | i<=nPr*nTs
                idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
                
                if length(tmpD(i).epochStarts)>1
                    iSet=tmpD(i).epochStarts(2);
                else
                    iSet=1; % tmpD(i).epochStarts(2);
                end
                sSet=[sSet; tmpD(i).data(1:3,iSet)'];
                
                iIC=iSet+round(dtIC/binSize);
                 if ~strfind(fname,'set2IC')
                sIC=[sIC; tmpD(i).data(1:3,iIC)'];
                end % if ~strfind(fname,'set2IC')
               
                    plot3(tmpD(i).data(1,iSet),tmpD(i).data(2,iSet),tmpD(i).data(3,iSet),'o',...
                        'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2);
%                         'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw2); % marker tmpCmap(i,:)      iEH
                    plot3(tmpD(i).data(1,end),tmpD(i).data(2,end),tmpD(i).data(3,end),'s',...
                        'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2);
%                         'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw2); % marker tmpCmap(i,:)      iEH
%                 end
                tmpMin=min([tmpMin min(tmpD(i).data(1:3,iSet:end),[],2)],[],2); % [3 x1]
                tmpMax=max([tmpMax max(tmpD(i).data(1:3,iSet:end),[],2)],[],2);

                end
            end
            
            
            tmp=[tmpMin tmpMax]';
            
             if ~strfind(fname,'set2IC')
            % connecting line for sSet or sIC
            for iPr=1:nPr
                tsId=[1:nTs]+(iPr-1)*nTs;
%                 plot3(sSet(tsId,1),sSet(tsId,2),sSet(tsId,3),':','color',pplot.cmap{nPr*iPr-1},'linewidth',lw);
                plot3(sIC(tsId,1),sIC(tsId,2),sIC(tsId,3),':','color',pplot.cmap{nPr*iPr-1},'linewidth',lw);
            end
            end % if ~strfind(fname,'set2IC')
            
            axis tight; %(tmp(:)'); % tight;
            set(gca,'tickdir','out','TickLength', [.02 .02]); % ,...
%                 'xtick',unique(round(tmpMin(1):1:tmpMax(1))),'ytick',unique(round(tmpMin(2):1:tmpMax(2))),'ztick',unique(round(tmpMin(3):1:tmpMax(3))));
            xlabel('PC1');ylabel('PC2'); zlabel('PC3');
            grid on;
            
            
            
            if idExpFig
                drawnow;
                savefig(hFig2,fullfile(figDir,'3PC',['PC123_' epNm{idEp} '_setGo_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_prodOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                xlabel([]);ylabel([]);zlabel([]);
                exportfig(hFig2,fullfile(figDir,'3PC',['PC123_' epNm{idEp} '_setGo_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.eps']),optsExpFig);%             'periSet4_prodOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);

            end
            
            condId=condId+1;
            if idCondSpecific, break; end;
        end %   if idProd
        if idCondSpecific, break; end;
    end % targ
    if idCondSpecific, break; end;
end % EH

% if exist('hFig')
% varargout{1}=hFig;
% end

%% back up 20170702
%             tmp1=spring(round(nTspp*1.5));
%             tmp2=summer(round(nTspp*1.5));
%             tmpCmap={flipud(autumn(nTspp)), flipud(tmp1(1:nTspp,:));... % shortEye, shortHand
%                 winter(nTspp), tmp2(1:nTspp,:)}; % longEye, longHand
%         tmpCmap=[flipud(autumn(nTspp));flipud(winter(nTspp))];

% cmap for split tp :80 conditions

%         if condId==1||condId==2 % eye
%             tmpCmap=[    1.0000    1.0000         0;...
%                 1.0000    0.7500         0;...
%                 1.0000    0.5000         0;...
%                 1.0000    0.2500         0;...
%                 1.0000         0         0;...
%                 0         0    1.0000;...
%                 0    0.2500    0.8750;...
%                 0    0.5000    0.7500;...
%                 0    0.7500    0.6250;...
%                 0    1.0000    0.5000;...
%                  1.0000    1.0000         0;...
%                 1.0000    0.7500         0;...
%                 1.0000    0.5000         0;...
%                 1.0000    0.2500         0;...
%                 1.0000         0         0;...
%                 0         0    1.0000;...
%                 0    0.2500    0.8750;...
%                 0    0.5000    0.7500;...
%                 0    0.7500    0.6250;...
%                 0    1.0000    0.5000];
%         else
%             tmpCmap=[    1.0000    0.5714    0.4286;...
%                 1.0000    0.4286    0.5714;...
%                 1.0000    0.2857    0.7143;...
%                 1.0000    0.1429    0.8571;...
%                 1.0000         0    1.0000;...
%                 0    0.5000    0.4000;...
%                 0.1429    0.5714    0.4000;...
%                 0.2857    0.6429    0.4000;...
%                 0.4286    0.7143    0.4000;...
%                 0.5714    0.7857    0.4000;...
%                 1.0000    0.5714    0.4286;...
%                 1.0000    0.4286    0.5714;...
%                 1.0000    0.2857    0.7143;...
%                 1.0000    0.1429    0.8571;...
%                 1.0000         0    1.0000;...
%                 0    0.5000    0.4000;...
%                 0.1429    0.5714    0.4000;...
%                 0.2857    0.6429    0.4000;...
%                 0.4286    0.7143    0.4000;...
%                 0.5714    0.7857    0.4000];
%         end
        
    %% modification 20170429
%     % mean across ts with attribution for measurement
%     nTs=5;
%     nPr=2;
%     tmpD=handles.D(handles.selected_conds);
%     tmpD0=tmpD;
%     % short
%     for j=1:(nTs-1)
%         if j==1
%             endPre=1;
%         else
%             endPre=tmpD0(j-1).epochStarts(2)+1;
%         end
%         iSet=tmpD0(j).epochStarts(2);
%         tmp=tmpD0(j).data(:,endPre:iSet);
%         
%         for k=(j+1):nTs
%             tmp=cat(3,tmp,tmpD0(k).data(:,endPre:iSet));
%         end
%         tmpD(j).data=[mean(tmp,3) tmpD0(j).data(:,iSet+1:end)];
%         tmpD(j).epochStarts(2)=iSet-endPre+1;        
%     end
%     j=nTs;
%     endPre=tmpD0(j-1).epochStarts(2)+1;
%     iSet=tmpD0(j).epochStarts(2);
%     tmpD(j).data=tmpD0(j).data(:,endPre:end);
%     tmpD(j).epochStarts(2)=iSet-endPre+1;
%     
%     % long
%     for j=(nTs+1):(nTs*nPr-1)
%         if j==(nTs+1)
%             endPre=1;
%         else
%             endPre=tmpD0(j-1).epochStarts(2)+1;
%         end
%         iSet=tmpD0(j).epochStarts(2);
%         tmp=tmpD0(j).data(:,endPre:iSet);
%         
%         for k=(j+1):nTs
%             tmp=cat(3,tmp,tmpD0(k).data(:,endPre:iSet));
%         end
%         tmpD(j).data=[mean(tmp,3) tmpD0(j).data(:,iSet+1)];
%         tmpD(j).epochStarts(2)=iSet-endPre+1;        
%     end
%         j=nTs*nPr;
%     endPre=tmpD0(j-1).epochStarts(2)+1;
%     iSet=tmpD0(j).epochStarts(2);
%     tmpD(j).data=tmpD0(j).data(:,endPre:end);
%     tmpD(j).epochStarts(2)=iSet-endPre+1;
%     
%     handles.D(handles.selected_conds)=tmpD;