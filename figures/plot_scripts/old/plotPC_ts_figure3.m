function plotPC_ts_figure3(fname)

% 2018/7/12: avgAttB4PCA

% plotPC_ts_figure3('trajKS_ts_EL_H_avgAttB4PCA.mat'); cc; % chosen for representative
% plotPC_ts_figure3('trajKS_ts_ER_H_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HL_H_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HR_H_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_EL_G_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_ER_G_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HL_G_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HR_G_avgAttB4PCA.mat'); cc;

% plotPC_ts_figure3('trajKS_prior_EL_H_avgAttB4PCA.mat'); cc; % chosen for representative
% plotPC_ts_figure3('trajKS_prior_ER_H_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_prior_HL_H_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_prior_HR_H_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_prior_EL_G_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_prior_ER_G_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_prior_HL_G_avgAttB4PCA.mat'); cc;
% plotPC_ts_figure3('trajKS_prior_HR_G_avgAttB4PCA.mat'); cc;

% plotPC_ts_figure3('trajKS_ts_EL_H_bin20_smth40.mat'); cc; % chosen for representative
% plotPC_ts_figure3('trajKS_ts_ER_H_bin20_smth40.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HL_H_bin20_smth40.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HR_H_bin20_smth40.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_EL_G_bin20_smth40.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_ER_G_bin20_smth40.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HL_G_bin20_smth40.mat'); cc;
% plotPC_ts_figure3('trajKS_ts_HR_G_bin20_smth40.mat'); cc;

% plotPC_ts_figure3('trajKS_prior_EL_H_bin20_smth40.mat'); waitforbuttonpress; cc; % chosen for representative
% plotPC_ts_figure3('trajKS_prior_ER_H_bin20_smth40.mat'); waitforbuttonpress; cc;
% plotPC_ts_figure3('trajKS_prior_HL_H_bin20_smth40.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_HR_H_bin20_smth40.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_EL_G_bin20_smth40.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_ER_G_bin20_smth40.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_HL_G_bin20_smth40.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_HR_G_bin20_smth40.mat'); waitforbuttonpress; cc;

% plotPC_ts_figure3('trajKS_prior_SEL_H_bin10_smth20.mat'); waitforbuttonpress; cc; % chosen for representative
% plotPC_ts_figure3('trajKS_prior_LEL_H_bin10_smth20.mat'); waitforbuttonpress; cc; % chosen for representative
% plotPC_ts_figure3('trajKS_prior_SER_H_bin10_smth20.mat'); waitforbuttonpress; cc;
% plotPC_ts_figure3('trajKS_prior_LER_H_bin10_smth20.mat'); waitforbuttonpress; cc;
% plotPC_ts_figure3('trajKS_prior_SHL_H_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_LHL_H_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_SHR_H_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_LHR_H_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_SEL_G_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_LEL_G_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_SER_G_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_LER_G_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_SHL_G_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_LHL_G_bin10_smth20.mat');  waitforbuttonpress;cc;
% plotPC_ts_figure3('trajKS_prior_SHR_G_bin10_smth20.mat'); waitforbuttonpress; cc;
% plotPC_ts_figure3('trajKS_prior_LHR_G_bin10_smth20.mat'); waitforbuttonpress; cc;

% old
% plotPC_SFN17('trajKS_periSet_E2_H_bin20_smth40.mat');
% plotPC_SFN17('trajKS_periSet_H2_H_bin20_smth40.mat');
% plotPC_SFN17('trajKS_periSet_E2_G_bin20_smth40.mat');
% plotPC_SFN17('trajKS_periSet_H2_G_bin20_smth40.mat');

if isempty(fname)
    fname='trajKS_periSet_EL_H_bin20_smth40.mat'; % trajKS_periSet_E2_H_bin20_smth40.mat';
end

% plot state space trajectory for each condition (of EH and targets) & two epoches
% data in DataHigh format, 40 conditions (2 pr x 5 ts x 2 EH x 2 targets)
% input: filename (assuming in /Users/hansem/Dropbox (MIT)/psthDataHigh/)
%   e.g. plotTraj('traj_periSet_G_attrition.mat'); plotTraj('traj_periSet_H_attrition.mat')
% if fname contains '_attrition', attrition
% assuming two epoches exist (epochStart(2)) as in ts
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
idExpFig=1; % 0; % 1;
idKilosort=0; % 1;

idCosyne=0; % 1;

cd(psthDir);
try
    load(fname); % dimReduce_periSet4_splitTp_attrition.mat'); % load('dimReduce_periSet4_attrition.mat'); % _eyeRightAttrition _attrition
catch
    load(['./sess/' fname]);
end

if strfind(fname,'periSet') % trajKS_periSet_ER_H
    idMeasProd=1;
    idMeasOnly=0;
    idAtt=0;
elseif strfind(fname,'ts') 
        idMeasProd=0; % 1; % plotting together
    idMeasOnly=1; % 0; % 1;
    idAtt=0; % 1;
% elseif strfind(fname,'tp') 
%         idMeasProd=0; % 1; % plotting together
%     idMeasOnly=1; % 0; % 1;
%     idAtt=0;
elseif strfind(fname,'prior')
    idMeasProd=0; % 1; % plotting together
    idMeasOnly=1; % 0; % 1;
        idAtt=0; % 1;
% %     if strfind(fname,'priorSL')
% %         idAtt=0;
% %     else
% %     idAtt=1;
% %     end
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
    lw=2; % trajectory
    lw2=1; % marker
    msize=10;
else
% figure's PaperUnits: 'inches'
optsExpFig.Height=2.1/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
optsExpFig.Width=3.15/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
% Width2=2.34/2.54; % 1.7*2*1.2; %
% Width3=2.1/2.54; % 1.3*2*1.2;
% optsExpFig.FontSize='6';
% optsExpFig.FontMode='fixed';
%                 optsExpFig.Format='eps'; % 'eps'; % 'eps'; % 'tiff'; % 'pdf'; % 'png';
% optsExpFig.LockAxes=1;
% optsExpFig.LineMode='fixed'; % 'scaled';
% optsExpFig.LineWidth=.5; % 1;
% % optsExpFig.LineWidthMin=0.5;
% % optsExpFig.LineWidthMax=1.2;
% 
% optsExpFig.Renderer='painters';

lw=.75; % 1;
msize=2; % 4;
lw2=.5; % 1; % marker
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
epNm={'fixation';'targetOn';'periReady';'periSet';'production';'ts';'tp';'_priorS_';'_priorL_';'_priorSL_';'t_s';'t_p';'_prior_'};
epNm2={'preFix','postFix';'preTOn','postTOn';'preReady','postReady';'','';'prod','postProd';'','';'','';'','';'','';'','';'','';'',''};
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
% individual sessions' data
if ~isempty(strfind(fname,'_prior_')) % & ~isempty(strfind(fname,'avgDir'))
    if ~isempty(strfind(D(1).condition,'Eye'))
        idCondSpecific=true; %tmpId=[1; 1];
    else % Hand
        idCondSpecific=true;% tmpId=[2; 1];
    end
    tmpCmap{end+1,1}=[0 0.6 0;0 0.6 0;0 0.6 0;0 0.6 0;0 0.6 0;]; % 1;0 0 1;0 0 1;0 0 1;0 0 1]; % [0 0 1;0 0 1;0 0 1;0 0 1;0 0 1]; % tmpCmap{1,1};
    linestyle={'-','-','--'};
else
    linestyle={'-','-'};
end % if ~isempty(strfind(fname,'prior')) & ~isempty(strfind(fname,'avgDir'))

%% avg attrition
if idAtt
    if strfind(fname,'prior')
        % filling NaN at the end
        D=fillNan(D);
        tmpDcat=cat(3,D.data); % [nPC time nCond]
        meanTmpD=nanmean(tmpDcat,3); % [nPC time]
        for i=1:length(D) % 5 for prior, 10 for the others
            D(i).data=meanTmpD(:,~isnan(D(i).data(1,:)));
        end % for i=1:length(D) % 5 for prior, 10 for the others
        
    else % separately for two priors
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
        
        
    end
    
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
        
        if idMeasOnly || idMeasProd
        
        %% measurement only
        hFig=figure;setFigPos(1,1);% set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+1)]));hold all;
        hFig2=figure;setFigPos(1,2);
        hFig3=figure;setFigPos(1,3);
        for i=1:length(tmpD)
            if length(tmpD(i).epochStarts)>1
                iSet=tmpD(i).epochStarts(2);
            else
                iSet=size(tmpD(i).data,2);
            end
            idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
            
            
            
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
                if size(tmpCmap,1)==3 & idPr==2 % long prior offset time
                    endPre=size(tmpD(nTs).data,2);
                else
                    endPre=1;
            end
%                 endPre=1;
                tmpX=tmpD(i).data(1,1:iSet);
                tmpY=tmpD(i).data(2,1:iSet);
                tmpZ=tmpD(i).data(3,1:iSet);
                % PC456
%                 tmpX2=tmpD(i).data(4,1:iSet);
%                 tmpY2=tmpD(i).data(5,1:iSet);
%                 tmpZ2=tmpD(i).data(6,1:iSet);
            end % with attrition
%             if i<=length(tmpD)/2 % bias+ % for split tp

%%%%% note PC3 filpped for visualization
if idFlipPC3
%                 plot3(tmpX,tmpY,-tmpZ,'-',...
%                     'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw); % lines tmpCmap(i,:)      iEH
%                  plot3(tmpD(i).data(1,iSet),tmpD(i).data(2,iSet),-tmpD(i).data(3,iSet),'o',...
%                     'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize); % marker tmpCmap(i,:)      iEH
               
else
    figure(hFig);
%     subplot(3,1,1),
    plot(endPre+[0:(length(tmpX)-1)],tmpX,'-',...
        'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
%     if strfind(fname,'priorSL')
%     else
%         plot(iSet,tmpD(i).data(1,iSet),'o',...
%             'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
%     end
    
    figure(hFig2);
%     subplot(3,1,2),
    plot(endPre+[0:(length(tmpX)-1)],tmpY,'-',...
        'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
%     if strfind(fname,'priorSL')
%     else
%         plot(iSet,tmpD(i).data(2,iSet),'o',...
%             'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
%     end
    
    figure(hFig3);
%     subplot(3,1,3),
    plot(endPre+[0:(length(tmpX)-1)],tmpZ,'-',...
        'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
%     if strfind(fname,'priorSL')
%     else
%         plot(iSet,tmpD(i).data(3,iSet),'o',...
%             'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
%     end
    % PC456
%     figure(hFig2);
%     subplot(3,1,1),plot(endPre+[0:(length(tmpX2)-1)],tmpX2,'-',...
%         'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
%     if strfind(fname,'priorSL')
%     else
%         plot(iSet,tmpD(i).data(4,iSet),'o',...
%             'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
%     end
%     ylabel(['PC' num2str(4)]);
%     subplot(3,1,2),plot(endPre+[0:(length(tmpY2)-1)],tmpY2,'-',...
%         'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
%     if strfind(fname,'priorSL')
%     else
%         plot(iSet,tmpD(i).data(5,iSet),'o',...
%             'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
%     end
%     ylabel(['PC' num2str(5)]);
%     subplot(3,1,3),plot(endPre+[0:(length(tmpZ2)-1)],tmpZ2,'-',...
%         'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
%     if strfind(fname,'priorSL')
%     else
%         plot(iSet,tmpD(i).data(6,iSet),'o',...
%             'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
%     end
%     ylabel(['PC' num2str(6)]);
% %                  if strfind(fname,'tp') 
% %                 plot3(tmpD(i).data(1,1),tmpD(i).data(2,1),tmpD(i).data(3,1),'o',...
% %                     'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize); % marker tmpCmap(i,:)          iEH
% %                  elseif strfind(fname,'priorSL')
% %                  else
% %                plot3(tmpD(i).data(1,iSet),tmpD(i).data(2,iSet),tmpD(i).data(3,iSet),'o',...
% %                     'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize); % marker tmpCmap(i,:)          iEH
% %                  end
                
                
%            
end
%             else
%                 plot3(tmpX,tmpY,tmpZ,':',...
%                     'color',tmpCmap(i,:),'linewidth',2); % lines
%                 plot3(tmpD(i).data(1,iSet),tmpD(i).data(2,iSet),tmpD(i).data(3,iSet),'^',...
%                     'markerfacecolor',tmpCmap(i,:),'markeredgecolor','k','markersize',10); % marker
%             end
            % for setting measure/produce in the same range
            tmpMin=min([tmpMin min(tmpD(i).data(1:3,:),[],2)],[],2); % [3 x1]
            tmpMax=max([tmpMax max(tmpD(i).data(1:3,:),[],2)],[],2);
            
        end
        
        % plot dot
        for i=1:length(tmpD)
            if length(tmpD(i).epochStarts)>1
                iSet=tmpD(i).epochStarts(2);
            else
                iSet=size(tmpD(i).data,2);
            end
            idPr=floor((i-1)/nTs)+1; %  1     1     1
            
            if size(tmpCmap,1)==3 & idPr==2 % long prior offset time
                tmpX2=iSet+size(tmpD(nTs).data,2);
            else
                tmpX2=iSet;
            end
            
            tmpX=tmpD(i).data(1,1:iSet);
            tmpY=tmpD(i).data(2,1:iSet);
            tmpZ=tmpD(i).data(3,1:iSet);
            
            if idPr~=3 % not for short in long
                figure(hFig);
                plot(tmpX2,tmpD(i).data(1,iSet),'o',...
                    'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2);
                
                figure(hFig2);
                plot(tmpX2,tmpD(i).data(2,iSet),'o',...
                    'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2);
                
                figure(hFig3);
                plot(tmpX2,tmpD(i).data(3,iSet),'o',...
                    'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2);
                
                if ~isempty(strfind(fname,'prior')) & ~idCosyne
                    figure(hFig);
                    plot(0,tmpD(i).data(1,1),'^',...
                        'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw2);
                    
                    figure(hFig2);
                    plot(0,tmpD(i).data(2,1),'^',...
                        'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw2);
                    
                    figure(hFig3);
                    plot(0,tmpD(i).data(3,1),'^',...
                        'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw2);
                end
            end
        end
                
                
        tmp=[tmpMin tmpMax]'; 
%         axis tight; %(tmp(:)'); % tight; %
%         set(gca,'xtick',floor(tmpMin(1)):1:ceil(tmpMax(1)),'ytick',floor(tmpMin(2)):1:ceil(tmpMax(2)),'ztick',floor(tmpMin(3)):1:ceil(tmpMax(3)));
%         xlabel('PC1');ylabel('PC2'); zlabel('PC3');
%         grid on;

        
        end %  if idMeasOnly || idMeasProd
        
        if ~idMeasOnly & idMeasProd
        %% production only
        periSetWindow=200; nSetWin=round(periSetWindow/binSize);
%          if ~idMeasProd
%         figure;set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+2)]));
        hFig1=figure;setFigPos(2,1);% set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+1)]));hold all;
        hFig12=figure;setFigPos(2,2);
        hFig13=figure;setFigPos(2,3);
%          end
        hold all;
        for i=1:length(tmpD)
            iSet=tmpD(i).epochStarts(2);
            idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
            %%%%% note PC3 filpped for visualization
            if idFlipPC3
                %             plot3(tmpD(i).data(1,iSet:end),tmpD(i).data(2,iSet:end),-tmpD(i).data(3,iSet:end),'-',...
                %                 'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw); % lines tmpCmap(i,:)       iEH
                %             plot3(tmpD(i).data(1,iSet),tmpD(i).data(2,iSet),-tmpD(i).data(3,iSet),'o',...
                %                 'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize); % marker tmpCmap(i,:)      iEH
            else
                endPre=-nSetWin; % 1;
                tmpX=tmpD(i).data(1,(iSet+endPre):end);
                tmpY=tmpD(i).data(2,(iSet+endPre):end);
                tmpZ=tmpD(i).data(3,(iSet+endPre):end);
                %      subplot(3,1,1),
                figure(hFig1);
                plot(endPre+[0:(length(tmpX)-1)],tmpX,'-',...
                    'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
                %                 plot(endPre,tmpX(endPre),'o',...
                %                     'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
                %     subplot(3,1,2),
                figure(hFig12);
                plot(endPre+[0:(length(tmpY)-1)],tmpY,'-',...
                    'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
                %                                 plot(endPre,tmpY(endPre),'o',...
                %                     'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
                
                %     subplot(3,1,3),
                figure(hFig13);
                plot(endPre+[0:(length(tmpZ)-1)],tmpZ,'-',...
                    'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;
                
                %                 plot(endPre,tmpZ(endPre),'o',...
                %                     'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize);
                
                %     plot3(tmpD(i).data(1,iSet:end),tmpD(i).data(2,iSet:end),tmpD(i).data(3,iSet:end),'-',...
                %                 'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw); % lines tmpCmap(i,:)       iEH
                %             plot3(tmpD(i).data(1,iSet),tmpD(i).data(2,iSet),tmpD(i).data(3,iSet),'o',...
                %                 'markerfacecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markeredgecolor','k','markersize',msize); % marker tmpCmap(i,:)      iEH
            end
            tmpMin=min([tmpMin min(tmpD(i).data(1:3,(iSet+endPre):end),[],2)],[],2); % [3 x1]
            tmpMax=max([tmpMax max(tmpD(i).data(1:3,(iSet+endPre):end),[],2)],[],2);
        end
        
        % plot dot
        for i=1:length(tmpD)
            iSet=tmpD(i).epochStarts(2);
            idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
            
            tmpX=tmpD(i).data(1,iSet);
            tmpY=tmpD(i).data(2,iSet);
            tmpZ=tmpD(i).data(3,iSet);
            
            figure(hFig1);
            plot(0,tmpX,'o',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
            
            figure(hFig12);
            plot(0,tmpY,'o',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
            
            figure(hFig13);
            plot(0,tmpZ,'o',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
            
                        tmpX=tmpD(i).data(1,end);
            tmpY=tmpD(i).data(2,end);
            tmpZ=tmpD(i).data(3,end);
            
            figure(hFig1);
            plot(endPre+(length(tmpD(i).data(1,(iSet+endPre):end))-1),tmpX,'s',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
            
            figure(hFig12);
            plot(endPre+(length(tmpD(i).data(1,(iSet+endPre):end))-1),tmpY,'s',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
            
            figure(hFig13);
            plot(endPre+(length(tmpD(i).data(1,(iSet+endPre):end))-1),tmpZ,'s',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
        end
        
%         for iSubplot=1:3
% subplot(3,1,iSubplot),
%         axis tight; %(tmp(:)'); % tight; %
%         set(gca,'xticklabel',[]);
%         ylabel(['PC' num2str(iSubplot)]);
% %         grid on;
%         end

        
        %% production time-locked to go
        
        hFig4=figure;setFigPos(2,4);% set(gcf,'position',pplot.(['rect' num2str(iEH) '_' num2str(2*(iTh-1)+1)]));hold all;
        hFig42=figure;setFigPos(2,5);
        hFig43=figure;setFigPos(2,6);
        %          end
        hold all;
        
        % find longest time
        iSet=tmpD(end).epochStarts(2);
        durProd=size(tmpD(end).data(1,iSet:end),2);
        
        for i=1:length(tmpD)
            iSet=tmpD(i).epochStarts(2);
            idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
            
            durProdCurrent=size(tmpD(i).data(1,iSet:end),2);
            endPre=durProd-durProdCurrent;
            tmpX=tmpD(i).data(1,iSet:end);
            tmpY=tmpD(i).data(2,iSet:end);
            tmpZ=tmpD(i).data(3,iSet:end);
            
            figure(hFig4);
            plot(endPre+[0:(length(tmpX)-1)],tmpX,'-',...
                'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;           
            
            figure(hFig42);
            plot(endPre+[0:(length(tmpY)-1)],tmpY,'-',...
                'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;            
            
            figure(hFig43);
            plot(endPre+[0:(length(tmpZ)-1)],tmpZ,'-',...
                'color',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'linewidth',lw);ha;           
                        
            tmpMin=min([tmpMin min(tmpD(i).data(1:3,(iSet):end),[],2)],[],2); % [3 x1]
            tmpMax=max([tmpMax max(tmpD(i).data(1:3,(iSet):end),[],2)],[],2);
        end
        
        % plot dot
        for i=1:length(tmpD)
            iSet=tmpD(i).epochStarts(2);
            idPr=floor((i-1)/nTs)+1; %  1     1     1     1     1     2     2     2     2     2
            
             durProdCurrent=size(tmpD(i).data(1,iSet:end),2);
            endPre=durProd-durProdCurrent;
            
            tmpX=tmpD(i).data(1,end);
            tmpY=tmpD(i).data(2,end);
            tmpZ=tmpD(i).data(3,end);
            
            figure(hFig4);
            plot(endPre+(length(tmpD(i).data(1,iSet:end))-1),tmpX(1),'s',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
            
            figure(hFig42);
            plot(endPre+(length(tmpD(i).data(2,iSet:end))-1),tmpY(1),'s',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);
            
            figure(hFig43);            
            plot(endPre+(length(tmpD(i).data(3,iSet:end))-1),tmpZ(1),'s',...
                'markeredgecolor',tmpCmap{idPr,1}(i-(idPr-1)*nTs,:),'markerfacecolor','w','markersize',msize,'linewidth',lw);

        end
        
        condId=condId+1;
        if idCondSpecific, break; end;
        end %   if ~idMeasOnly & idMeasProd
        if idCondSpecific, break; end;
    end % targ
    if idCondSpecific, break; end;
end % EH

%% final touch & export
if idMeasOnly || idMeasProd
% ready
figure(hFig);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(1) tmpMax(1)]);
plotVertical(gca,0,[]);
set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',unique([0 min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})])/binSize,...
    'xticklabel',{0,min(T{1})/1000,[],max(T{1})/1000,[],max(T{2})/1000},'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig,fullfile(figDir,'3','PC',['PC1_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    applytofig(hFig,optsExpFig);
    print(hFig,fullfile(figDir,'3','PC',['PC1_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.eps']),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
%     exportfig(hFig,fullfile(figDir,'3','PC',['PC1_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end

figure(hFig2);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(2) tmpMax(2)]);
plotVertical(gca,0,[]);
set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',unique([0 min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})])/binSize,...
    'xticklabel',{0,min(T{1})/1000,[],max(T{1})/1000,[],max(T{2})/1000},'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig2,fullfile(figDir,'3','PC',['PC2_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    applytofig(hFig2,optsExpFig);
    print(hFig2,fullfile(figDir,'3','PC',['PC2_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.eps']),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
%     exportfig(hFig2,fullfile(figDir,'3','PC',['PC2_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end

figure(hFig3);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(3) tmpMax(3)]);
plotVertical(gca,0,[]);
set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',unique([0 min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})])/binSize,...
    'xticklabel',{0,min(T{1})/1000,[],max(T{1})/1000,[],max(T{2})/1000},'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig3,fullfile(figDir,'3','PC',['PC3_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    applytofig(hFig3,optsExpFig);
    print(hFig3,fullfile(figDir,'3','PC',['PC3_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.eps']),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
%     exportfig(hFig3,fullfile(figDir,'3','PC',['PC3_' epNm{idEp} '_ready_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end

end % if idMeasOnly || idMeasProd

if ~idMeasOnly & idMeasProd
% set
 optsExpFig.Width=Width2; 
figure(hFig1);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(1) tmpMax(1)]);
plotVertical(gca,0,[]);
set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',[0 max(T{1})/2-binSize max(T{1})-binSize]/binSize, ...% unique([0 -min(T{1}) -mean(T{1}) -max(T{1}) -mean(T{2}) -max(T{2}) min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
    'xticklabel',[0 max(T{1})/2/1000 max(T{1})/1000],'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig1,fullfile(figDir,'3PC',['PC1_' epNm{idEp} '_set_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    exportfig(hFig1,fullfile(figDir,'3PC',['PC1_' epNm{idEp} '_set_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end
            
figure(hFig12);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(2) tmpMax(2)]);
plotVertical(gca,0,[]);
set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',[0 max(T{1})/2-binSize max(T{1})-binSize]/binSize, ...% unique([0 -min(T{1}) -mean(T{1}) -max(T{1}) -mean(T{2}) -max(T{2}) min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
    'xticklabel',[0 max(T{1})/2/1000 max(T{1})/1000],'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig12,fullfile(figDir,'3PC',['PC2_' epNm{idEp} '_set_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    exportfig(hFig12,fullfile(figDir,'3PC',['PC2_' epNm{idEp} '_set_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end

figure(hFig13);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(3) tmpMax(3)]);
plotVertical(gca,0,[]);
set(gca,'tickdir','out','TickLength', [.02 .02],'xtick',[0 max(T{1})/2-binSize max(T{1})-binSize]/binSize, ...% unique([0 -min(T{1}) -mean(T{1}) -max(T{1}) -mean(T{2}) -max(T{2}) min(T{1}) mean(T{1}) max(T{1}) mean(T{2}) max(T{2})]),...
    'xticklabel',[0 max(T{1})/2/1000 max(T{1})/1000],'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig13,fullfile(figDir,'3PC',['PC3_' epNm{idEp} '_set_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    exportfig(hFig13,fullfile(figDir,'3PC',['PC3_' epNm{idEp} '_set_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end


% go
optsExpFig.Width=Width3;
figure(hFig4);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(1) tmpMax(1)]);
plotVertical(gca,max(get(gca,'xlim')),[]);
set(gca,'xtick',max(get(gca,'xlim'))+([-max(T{1})+binSize -(max(T{1})/2-binSize) 0])/binSize,'tickdir','out','TickLength', [.02 .02],...
    'xticklabel',[-max(T{1})/1000 -max(T{1})/2/1000 0],'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig4,fullfile(figDir,'3PC',['PC1_' epNm{idEp} '_go_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    exportfig(hFig4,fullfile(figDir,'3PC',['PC1_' epNm{idEp} '_go_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end

figure(hFig42);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(2) tmpMax(2)]);
plotVertical(gca,max(get(gca,'xlim')),[]);
set(gca,'xtick',max(get(gca,'xlim'))+([-max(T{1})+binSize -max(T{1})/2+binSize 0])/binSize,'tickdir','out','TickLength', [.02 .02],...
    'xticklabel',[-max(T{1})/1000 -max(T{1})/2/1000 0],'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig42,fullfile(figDir,'3PC',['PC2_' epNm{idEp} '_go_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    exportfig(hFig42,fullfile(figDir,'3PC',['PC2_' epNm{idEp} '_go_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end

figure(hFig43);box off;
axis tight; %(tmp(:)'); % tight; %
ylim([tmpMin(3) tmpMax(3)]);
plotVertical(gca,max(get(gca,'xlim')),[]);
set(gca,'xtick',max(get(gca,'xlim'))+([-max(T{1})+binSize -max(T{1})/2+binSize 0])/binSize,'tickdir','out','TickLength', [.02 .02],...
    'xticklabel',[-max(T{1})/1000 -max(T{1})/2/1000 0],'ytick',unique(round(get(gca,'ytick'))),'yticklabel',{unique(round(get(gca,'ytick')))});
if idExpFig
    drawnow;
    savefig(hFig43,fullfile(figDir,'3PC',['PC3_' epNm{idEp} '_go_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.fig'])); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    exportfig(hFig43,fullfile(figDir,'3PC',['PC3_' epNm{idEp} '_go_' ehNm{iEH}(1) targNm{iTh}(1) '_' animalNm{idAnimal} '.' optsExpFig.Format]),optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);
end
end

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