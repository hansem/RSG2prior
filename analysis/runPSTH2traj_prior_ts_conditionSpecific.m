function runPSTH2traj_prior_ts_conditionSpecific

% 2019/1/20
% re-run after manual update of behavioral idOut
% run for both idIncShortInLong= 0 or 1
% debug for idIncShortInLong (D>oldD)
% smthWidth 20 %%%%%
% no buffer after Set before PCA

% 2018/7/11
% avg across ts w/ attrition before PCA
% after getting PC coeff(projMatrix), project original data w/o avgAtt
% for ts, prior (short,long,shortInLong)
% filename: avgAttB4PCA
% buffer removed for SiL end & L start
% caveat: assuming # trials/ts is identical (TBD: new PSTH)

% 2018/5/9
% save offset of meanPSTH for centering in PCA

% 2018/3/27: correct wrong tInfo (12/7, 8/22, 8/23)
% tBuffer(80ms) also added after end of prior: cut off 60ms after smoothing
% e.g. for 480ms, PSTH=161ms > trajKS=8-3*2(before/after)=2
% assuming visual response latency to set is longer than tBuffer
% now shortInLong included: idIncShortInLong=1

% 2018/2/18
% PCA to prior support (not including shortInLong)
% to test importance of set-induced responses to read out rotating trajectory into proper initial conditions: what if randomly projected?
% use PSTH_prior_H/G_full_SUMU.mat (which include shortInLong)

% 2017/1/16
%   PCA separately for ts and tp

% 2017/12/8
%   after kilosort(+redo)
%   condition-specific for both H/G
%   after PCA, getting PVAF(#PC=3) by projecting data across conditions' PC
%   binSize=20;smthWidth=40;
% save PVAF across conditions

% 2017/9/25
% debug: project_matrix, keep_neurons should be empty
%       runPSTH2traj_periSet_conditionSpecific_H(5,40);cc;
%       runPSTH2traj_periSet_conditionSpecific_H(10,40);cc;
%       runPSTH2traj_periSet_conditionSpecific_H(5,25);cc;

% project PSTH into trajectory
% psth_prior_H_splitTp > traj_prior_ER_H, traj_prior_EL_H,traj_prior_HR_H,traj_prior_HL_H

%%
idScree=1; % 0;
idPlotTraj=1; % 0;

idIncShortInLong= 1; % 0; % 1; 
idTrunBuffer=1;
idSave=1;

initRSG2prior; % psthDir animalNm
cd(psthDir);
% optsExpFig.Width=3.5;
% optsExpFig.Height=3.5;

if idIncShortInLong, nPr=3; end % 0;

binSize=20;
smthWidth=20; % 40; %%%%%

condNm={'ER','EL','HR','HL'};
condNm2={'Eye Left','Eye Right','Hand Left','Hand Right'};
nCond=numel(condNm);

% parameters for dimreduce
optimD=[]; %%%%%
use_sqrt=0; %%%%%
proj_matrix=[];
keep_neurons=[];


%% prior
for iAnimal=1:length(animalNm)    
    load(['PSTH_prior_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D nTrialExc tBuffer
    oldD=D; % data for after-PCA
    
    if idIncShortInLong==0
        D=D(1:(nEH*nTarg*nPr*nTspp)); % not including shrot in long (D(41:60))
        oldD=D;
        
        % no buffer after Set
        for iD=1:(nEH*nTarg*nTspp*nPr) % 1:40
            D(iD).data=D(iD).data(:,1:(end-tBuffer));
        end
        
    else
        % remove overlapped buffer b/t SiL & L
        % no buffer after Set > no data for shortest ts for long > not remove overlap then
        for iD=1:(nEH*nTarg*nTspp*(nPr-1)) % 1:40
            D(iD).data=D(iD).data(:,1:(end-tBuffer));
        end
%         for iD=(1+nEH*nTarg*nTspp):(nEH*nTarg*nTspp*(nPr-1)) % 21:40 L 
%             D(iD).data=D(iD).data(:,(tBuffer+1):end);
%         end
        for iD=(1+nEH*nTarg*nTspp*(nPr-1)):(nEH*nTarg*nTspp*nPr) % 41:60 SiL
            D(iD).data=D(iD).data(:,1:(end-tBuffer));
        end
    end
    
    % avgAttB4PCA
%     Davg=struct([]);
    for iPr=1:nPr
        for iHE=1:nEH
            for iTarg=1:nTarg
                d0=(iPr-1)*nEH*nTarg*nTspp+... % 0 20 40
                    ((iHE-1)*nEH+iTarg); % 1 2 3 4
                d1=iPr*nEH*nTarg*nTspp; % 20 40 60
                tmp=struct2mat(D(d0:nCond:d1),'data'); % [neuron x time x 5ts]
                
                Davg((iPr-1)*nEH*nTarg+(iHE-1)*nEH+iTarg)=D(d0); % 12 struct arrays
                Davg((iPr-1)*nEH*nTarg+(iHE-1)*nEH+iTarg).data=squeeze(nanmean(tmp,3));
            end
        end
    end
    
    % main
    for iHE=1:nEH
        for iTarg=1:nTarg
            % PCA with avgAtt
            [Dtmp,proj_matrix,keep_neurons,eigenvaluesAvgAtt,meanPSTH]=PSTH2traj(Davg(((iHE-1)*nEH+iTarg):nCond:length(Davg)),...
                [],[],binSize,smthWidth,[],use_sqrt); %  proj_matrix,keep_neurons optimD
            meanPSTH=meanPSTH(:);
            if idIncShortInLong==0 
                fnTmp=['trajKS_prior_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_supportOnly_avgAttB4PCA_smth20.mat']; %  '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
            else %%%%%
                fnTmp=['trajKS_prior_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_avgAttB4PCA_smth20.mat']; % '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
            end
            optimD=size(Dtmp(1).data,1);
            
            % projection with full data (separately ts)
            proj_matrixTmp=cat(3,proj_matrix,[meanPSTH(:) nan(size(proj_matrix,1),size(proj_matrix,2)-1)]); % [#cell x dim x 2(projMat/meanPSTH)]
            [D,proj_matrix,keep_neurons,eigenvalues,totalVar]=PSTH2traj(oldD(((iHE-1)*nEH+iTarg):nCond:length(oldD)),...
                proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt); %  proj_matrix,keep_neurons optimD
            
            % no buffer after Set before PCA
            if idTrunBuffer % tBuffer(80ms) also added after end of prior: cut off after smoothing
                nBinBuffer=round(tBuffer/binSize)-1; % 4-1; if 4, no data from shortest ts
                if idIncShortInLong==1 % overlapped buffer b/t SiL & L removed      
                    for iD=1:nTspp % short
                        D(iD).data=D(iD).data(:,(nBinBuffer+1):(end-nBinBuffer)); % )); % -nBinBuffer));
                    end
                    for iD=(nTspp+1):(nTspp+nTspp) % long
                        D(iD).data=D(iD).data(:,(nBinBuffer+1):(end-nBinBuffer)); % )); % 1:(end)); % -nBinBuffer));
                    end
                    for iD=(2*nTspp+1):(3*nTspp) % SiL
                        D(iD).data=D(iD).data(:,(nBinBuffer+1):(end-nBinBuffer)); %end);
                    end
                else          
                    for iD=1:length(D)
                        D(iD).data=D(iD).data(:,(nBinBuffer+1):(end-nBinBuffer)); % )); % -nBinBuffer));
                    end
                end
            end
            
            if idSave %  ~exist(fnTmp,'file')
                save(fnTmp,...
                    'D','proj_matrix','keep_neurons','binSize','smthWidth','optimD','use_sqrt','eigenvalues','meanPSTH','totalVar','eigenvaluesAvgAtt');
            end
            
            if idPlotTraj
                % plot
                plotTraj_figure3(fnTmp); ha; setFigPos(iAnimal,(iHE-1)*nEH+iTarg); 
%                 title(['prior ' animalNm{iAnimal} ' ' ehNm{iHE}(1) targNm{iTarg}(1)]);
                applytofig(gcf,optsExpFig); % plotTraj_SFN17(fnTmp);
                if idIncShortInLong==0
                    savefig(gcf,fullfile(figDir,'3','PC',['PC123_prior_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_supportOnly']));
                    remAxLabel;
                    print(gcf,fullfile(figDir,'3','PC',['PC123_prior_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_supportOnly.eps']),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);                    
                else
                    savefig(gcf,fullfile(figDir,'3','PC',['PC123_prior_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal}]));
                    remAxLabel;
                    print(gcf,fullfile(figDir,'3','PC',['PC123_prior_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '.eps']),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
                end
            end
        end % iTarg
    end % iHE
    clear Davg;
end % for iAnimal=1:length(animalNm)


%% ts
for iAnimal=1:length(animalNm)    
    load(['PSTH_t_s_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D
    oldD=D;
 
    nPr=2;
    
    % avgAttB4PCA
%     Davg=struct([]);
    for iPr=1:nPr % (nPr-1) % not short in long
        for iHE=1:nEH
            for iTarg=1:nTarg
                d0=(iPr-1)*nEH*nTarg*nTspp+... % 0 20
                    ((iHE-1)*nEH+iTarg); % 1 2 3 4
                d1=iPr*nEH*nTarg*nTspp; % 20 40
                tmp=struct2mat(oldD(d0:nCond:d1),'data'); % [neuron x time x 5ts]
                
                Davg((iPr-1)*nEH*nTarg+(iHE-1)*nEH+iTarg)=oldD(d0); % 8 struct arrays
                Davg((iPr-1)*nEH*nTarg+(iHE-1)*nEH+iTarg).data=squeeze(nanmean(tmp,3));
            end
        end
    end
    
    % main
    for iHE=1:nEH
        for iTarg=1:nTarg
            % PCA with avgAtt
            [Dtmp,proj_matrix,keep_neurons,eigenvaluesAvgAtt,meanPSTH]=PSTH2traj(Davg(((iHE-1)*nEH+iTarg):nCond:length(Davg)),...
                [],[],binSize,smthWidth,[],use_sqrt); %  proj_matrix,keep_neurons optimD
            meanPSTH=meanPSTH(:);
            fnTmp=['trajKS_ts_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_avgAttB4PCA.mat']; %'_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
            optimD=size(Dtmp(1).data,1);
            
            % projection with full data (separately ts)
            proj_matrixTmp=cat(3,proj_matrix,[meanPSTH(:) nan(size(proj_matrix,1),size(proj_matrix,2)-1)]); % [#cell x dim x 2(projMat/meanPSTH)]
            [D,proj_matrix,keep_neurons,eigenvalues,totalVar]=PSTH2traj(oldD(((iHE-1)*nEH+iTarg):nCond:length(oldD)),...
                proj_matrixTmp,keep_neurons,binSize,smthWidth,optimD,use_sqrt); %  proj_matrix,keep_neurons optimD
            
            if idSave % ~exist(fnTmp,'file')
                save(fnTmp,...
                    'D','proj_matrix','keep_neurons','binSize','smthWidth','optimD','use_sqrt','eigenvalues','meanPSTH','totalVar','eigenvaluesAvgAtt');
            end
            
            if idPlotTraj
                % plot
                plotTraj_figure3(fnTmp);  ha; setFigPos(iAnimal,(iHE-1)*nEH+iTarg); 
%                 title(['ts ' animalNm{iAnimal} ' ' ehNm{iHE}(1) targNm{iTarg}(1)]);
                applytofig(gcf,optsExpFig); % % plotTraj_SFN17(fnTmp);
                savefig(gcf,fullfile(figDir,'3','PC',['PC123_ts_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal}]));
                remAxLabel;
                print(gcf,fullfile(figDir,'3','PC',['PC123_ts_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '.eps']),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);

                
            end
        end % iTarg
    end % iHE
    clear Davg;
end % for iAnimal=1:length(animalNm)

%%   scree plot
if idScree
    legendMat=cell(size(condNm));
    load pplot.mat;
    msize=2;
    lw=1;
    
    nPCmax=10;
    thVar=.75;
    
    nCmap=length(animalNm)*length(condNm);
    cmap=parula(nCmap); % hsv(nCmap);
    
    % prior
    icmap=1;
    fnTmp1='trajKS_prior_';
    figure;ha;setFigPos(1,6);

    for iAnimal=1:length(animalNm)
        if idIncShortInLong==0
            fnTmp2=['_' animalNm{iAnimal}  '_supportOnly_avgAttB4PCA.mat']; % '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
        else
            fnTmp2=['_' animalNm{iAnimal} '_avgAttB4PCA.mat']; % _bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
        end
        for i=1:length(condNm)
            load([fnTmp1 condNm{i} fnTmp2]); tmpD=load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues
            disp([fnTmp1 condNm{i} fnTmp2]);
            if isfield(tmpD,'eigenvaluesAvgAtt') % use eigenvalues from avgAtt before PCA
                plot(1:nPCmax,cumsum(eigenvaluesAvgAtt(1:nPCmax))./sum(eigenvaluesAvgAtt)*100,'o-','color',cmap(icmap,:),'markerfacecolor',cmap(icmap,:),'markeredgecolor',cmap(icmap,:),...
                    'markersize',msize,'linewidth',lw);
                optD=find(cumsum(eigenvaluesAvgAtt(1:nPCmax))./sum(eigenvaluesAvgAtt)>thVar,1,'first');
            else % if isfield(tmp,'totalVar')
                plot(1:nPCmax,cumsum(eigenvalues(1:nPCmax))./totalVar*100,'o-','color',cmap(icmap,:),'markerfacecolor',cmap(icmap,:),'markeredgecolor',cmap(icmap,:),...
                    'markersize',msize,'linewidth',lw);
                optD=find(cumsum(eigenvalues(1:nPCmax))./totalVar>thVar,1,'first');
            end
            disp(optD);
            legendMat{icmap}=[animalNm{iAnimal} ' ' condNm2{i}]; % ':' num2str(optD) ];
            icmap=icmap+1;
        end
        xlabel('Number of PCs'); ylabel('Cum. var. explained');
        set(gca,'ytick',0:25:100,'yticklabel',{0;[];50;[];100},'ylim',[0 100],'xlim',[0 nPCmax],...
            'xtick',[0:3:9],'tickdir','out','ticklength',[.02 .02]);
        %         plot([0 3],ones(1,2)*100*sum(eigenvalues(1:3))./sum(eigenvalues),'k:');
        %
        %         plotHorizon(gca,thVar,[]);
        %         plotVertical(gca,3,[]);
%         legend(legendMat,'location','best'); legend boxoff;
    end % for iAnimal=1:length(animalNm)
    applytofig(gcf,optsExpFig); % plotTraj_SFN17(fnTmp);
    if idIncShortInLong==0
        savefig(gcf,fullfile(figDir,'3','PC','screePlot_prior_supportonly_avgAttB4PCA'));
        remAxLabel;
        print(gcf,fullfile(figDir,'3','PC','screePlot_prior_supportonly_avgAttB4PCA.eps'),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    else
        savefig(gcf,fullfile(figDir,'3','PC','screePlot_prior_incShortInLong_avgAttB4PCA'));
        remAxLabel;
        print(gcf,fullfile(figDir,'3','PC','screePlot_prior_incShortInLong_avgAttB4PCA.eps'),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    end
    
    % ts
    icmap=1;
    fnTmp1='trajKS_ts_';
    figure;ha;setFigPos(2,6);

    for iAnimal=1:length(animalNm)
        fnTmp2=['_' animalNm{iAnimal} '_avgAttB4PCA.mat']; %_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
        for i=1:length(condNm)
            load([fnTmp1 condNm{i} fnTmp2]); tmpD=load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues
            disp([fnTmp1 condNm{i} fnTmp2]);
            if isfield(tmpD,'eigenvaluesAvgAtt') % use eigenvalues from avgAtt before PCA
                plot(1:nPCmax,cumsum(eigenvaluesAvgAtt(1:nPCmax))./sum(eigenvaluesAvgAtt)*100,'o-','color',cmap(icmap,:),'markerfacecolor',cmap(icmap,:),'markeredgecolor',cmap(icmap,:),...
                    'markersize',msize,'linewidth',lw);
                optD=find(cumsum(eigenvaluesAvgAtt(1:nPCmax))./sum(eigenvaluesAvgAtt)>thVar,1,'first');
            else % if isfield(tmp,'totalVar')
                plot(1:nPCmax,cumsum(eigenvalues(1:nPCmax))./totalVar*100,'o-','color',cmap(icmap,:),'markerfacecolor',cmap(icmap,:),'markeredgecolor',cmap(icmap,:),...
                    'markersize',msize,'linewidth',lw);
                optD=find(cumsum(eigenvalues(1:nPCmax))./totalVar>thVar,1,'first');
            end
            disp(optD);
            legendMat{icmap}=[animalNm{iAnimal} ' ' condNm2{i}]; % ':' num2str(optD) ];
            icmap=icmap+1;
        end
        xlabel('Number of PCs'); ylabel('Cum. var. explained');
        set(gca,'ytick',0:25:100,'yticklabel',{0;[];50;[];100},'ylim',[0 100],'xlim',[0 nPCmax],...
            'xtick',[0:3:9],'tickdir','out','ticklength',[.02 .02]);
        %         plot([0 3],ones(1,2)*100*sum(eigenvalues(1:3))./sum(eigenvalues),'k:');
        %
        %         plotHorizon(gca,thVar,[]);
        %         plotVertical(gca,3,[]);
%         legend(legendMat,'location','best'); legend boxoff;
    end % for iAnimal=1:length(animalNm)
                    applytofig(gcf,optsExpFig); % plotTraj_SFN17(fnTmp);
                savefig(gcf,fullfile(figDir,'3','PC','screePlot_ts_avgAttB4PCA'));
                remAxLabel;
                print(gcf,fullfile(figDir,'3','PC','screePlot_ts_avgAttB4PCA.eps'),'-depsc'); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);

    
%     for iAnimal=1:length(animalNm)
%         fnTmp2=['_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
%         
%         figure;ha;setFigPos(iAnimal,6);
%         for i=1:length(condNm)
%             load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues
%             disp([fnTmp1 condNm{i} fnTmp2]);
%             plot(cumsum(eigenvalues(1:nPCmax))./sum(eigenvalues),'.-','color',pplot.cmap{i});
%             optD=find(cumsum(eigenvalues(1:20))./sum(eigenvalues)>thVar,1,'first');
%             disp(optD);
%             legendMat{i}=[condNm{i} ':' num2str(optD)];
%         end
%         xlabel('# PC'); ylabel('% Var. Explained');
%         set(gca,'ytick',0:.25:1,'ylim',[0 1],'xlim',[1 nPCmax]);
%         plotHorizon(gca,thVar,[]);
%         plotVertical(gca,3,[]);
%         legend(legendMat,'location','best'); legend boxoff;
%     end % for iAnimal=1:length(animalNm)

   

end


%%   after PCA, getting PVAF(#PC=3) by projecting data across conditions' PC
% PVAF=cell(nEH*nTarg,nEH*nTarg,length(animalNm));
% PVAF3=nan(nEH*nTarg,nEH*nTarg,length(animalNm));
% nPC=3;
% for iAnimal=1:length(animalNm)    
% 
%     load(['PSTH_periSet_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D
%     
%     % main
%     oldD=D;
%     for iHE=1:nEH
%         for iTarg=1:nTarg
%             ii=(iHE-1)*nEH+iTarg;            
%             
%             fnTmp=['trajKS_periSet_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
%             load(fnTmp); % D  proj_matrix,keep_neurons optimD
%             P=proj_matrix;
%             K=keep_neurons;
%             O=optimD;
%             
%             for jHE=1:nEH
%                 for jTarg=1:nTarg
%                     jj=(jHE-1)*nEH+jTarg;
%                     [D,proj_matrix,keep_neurons,eigenvalues,totalVar]=PSTH2traj(oldD(((jHE-1)*nEH+jTarg):nCond:length(oldD)),...
%                         P,K,binSize,smthWidth,O,use_sqrt); %  proj_matrix,keep_neurons optimD
%                     
%                     PVAF{ii,jj,iAnimal}=100*eigenvalues./totalVar; % (sum(eigenvalues));
%                     PVAF3(ii,jj,iAnimal)=sum(PVAF{ii,jj,iAnimal}(1:nPC));
%                 end
%             end
%             
%             D=oldD;
%             
%         end % iTarg
%     end % iHE
%     
%     figure; setFigPos(iAnimal,1);
%     imagesc(PVAF3(:,:,iAnimal)); colorbar;
%     set(gca,'xtick',[],'ytick',[]);
% end % for iAnimal=1:length(animalNm)
