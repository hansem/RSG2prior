function plotScree_ts_tp

% 2018/2/21
% for COSYNE

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

initRSG2prior; % psthDir animalNm
cd(psthDir);

binSize=20;
smthWidth=40;

condNm={'ER','EL','HR','HL'};
condNm2={'Eye Left','Eye Right','Hand Left','Hand Right'};
nCond=numel(condNm);

% parameters for dimreduce
optimD=[]; %%%%%
use_sqrt=0; %%%%%
proj_matrix=[];
keep_neurons=[];

% % ts
% for iAnimal=1:length(animalNm)    
%     load(['PSTH_t_s_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D
% 
%     % main
%     oldD=D;
%     for iHE=1:nEH
%         for iTarg=1:nTarg
%             [D,proj_matrix,keep_neurons,eigenvalues]=PSTH2traj(D(((iHE-1)*nEH+iTarg):nCond:length(D)),...
%                 [],[],binSize,smthWidth,[],use_sqrt); %  proj_matrix,keep_neurons optimD
%             fnTmp=['trajKS_ts_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
%             optimD=size(D(1).data,1);
%             if ~exist(fnTmp,'file')
%                 save(fnTmp,...
%                     'D','proj_matrix','keep_neurons','binSize','smthWidth','optimD','use_sqrt','eigenvalues');
%             end
%             D=oldD;
%             
%             if idPlotTraj
%                 % plot
%                 plotTraj_figure3(fnTmp); % plotTraj_SFN17(fnTmp);
%             end
%         end % iTarg
%     end % iHE
%     
% end % for iAnimal=1:length(animalNm)
% 
% for iAnimal=1:length(animalNm)    
%     load(['PSTH_t_p_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D
% 
%     % main
%     oldD=D;
%     for iHE=1:nEH
%         for iTarg=1:nTarg
%             [D,proj_matrix,keep_neurons,eigenvalues]=PSTH2traj(D(((iHE-1)*nEH+iTarg):nCond:length(D)),...
%                 [],[],binSize,smthWidth,[],use_sqrt); %  proj_matrix,keep_neurons optimD
%             fnTmp=['trajKS_tp_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
%             optimD=size(D(1).data,1);
%             if ~exist(fnTmp,'file')
%                 save(fnTmp,...
%                     'D','proj_matrix','keep_neurons','binSize','smthWidth','optimD','use_sqrt','eigenvalues');
%             end
%             D=oldD;
%             
%             if idPlotTraj
%                 % plot
%                 plotTraj_figure3(fnTmp); % plotTraj_SFN17(fnTmp);
%             end
%         end % iTarg
%     end % iHE
%     
% end % for iAnimal=1:length(animalNm)


%%   scree plot
if idScree
    legendMat=cell(size(condNm));
    load pplot.mat;
    msize=4; % 6; % 2;
    lw=1.5;
    
    nPCmax=12;
    thVar=.75;
    
    nCmap=length(animalNm)*length(condNm);
    cmap=[0 0 0; .5 .5 .5]; % animal & E/H
    cmap2={[0 0 0; 1 1 1];
        [.5 .5 .5; 1 1 1]};
    markers={'o','^','o','^'}; % condNm={'ER','EL','HR','HL'};
%     cmap=parula(nCmap); % hsv(nCmap);
    
    % ts
    icmap=1;
    fnTmp1='trajKS_ts_';
    figure;ha;setFigPos(1,6);

    for iAnimal=1:length(animalNm)
        fnTmp2=['_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
        for i=1:length(condNm)
            load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues '-'
            disp([fnTmp1 condNm{i} fnTmp2]);
            plot(1:nPCmax,cumsum(eigenvalues(1:nPCmax))./sum(eigenvalues)*100,[markers{i} ],'color',cmap(iAnimal,:),'markerfacecolor',cmap2{iAnimal}(floor((i-1)/2)+1,:),'markeredgecolor',cmap(iAnimal,:),...
                'markersize',msize,'linewidth',lw);
            optD=find(cumsum(eigenvalues(1:20))./sum(eigenvalues)>thVar,1,'first');
            disp(optD);
            legendMat{icmap}=[animalNm{iAnimal} ' ' condNm2{i}]; % ':' num2str(optD) ];
            icmap=icmap+1;
        end
        xlabel('Number of PCs'); ylabel('Cum. var. explained');
        set(gca,'ytick',0:25:100,'yticklabel',{0;[];50;[];100},'ylim',[0 100],'xlim',[0 nPCmax],...
            'xtick',[0:3:nPCmax],'tickdir','out','ticklength',[.02 .02]);
        %         plot([0 3],ones(1,2)*100*sum(eigenvalues(1:3))./sum(eigenvalues),'k:');
        %
        %         plotHorizon(gca,thVar,[]);
        %         plotVertical(gca,3,[]);
%         legend(legendMat,'location','best'); legend boxoff;
    end % for iAnimal=1:length(animalNm)
    
    % tp
    legendMat=cell(size(condNm));
    icmap=1;
    fnTmp1='trajKS_tp_';
    figure;ha;setFigPos(2,6);

    for iAnimal=1:length(animalNm)
        fnTmp2=['_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
        for i=1:length(condNm)
            load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues '-'
            disp([fnTmp1 condNm{i} fnTmp2]);
            plot(1:nPCmax,cumsum(eigenvalues(1:nPCmax))./sum(eigenvalues)*100,[markers{i}],'color',cmap(iAnimal,:),'markerfacecolor',cmap2{iAnimal}(floor((i-1)/2)+1,:),'markeredgecolor',cmap(iAnimal,:),...
                'markersize',msize,'linewidth',lw);
            optD=find(cumsum(eigenvalues(1:20))./sum(eigenvalues)>thVar,1,'first');
            disp(optD);
            legendMat{icmap}=[animalNm{iAnimal} ' ' condNm2{i}]; % ':' num2str(optD) ];
            icmap=icmap+1;
        end
        xlabel('Number of PCs'); ylabel('Cum. var. explained');
        set(gca,'ytick',0:25:100,'yticklabel',{0;[];50;[];100},'ylim',[0 100],'xlim',[0 nPCmax],...
            'xtick',[0:3:nPCmax],'tickdir','out','ticklength',[.02 .02]);
        %         plot([0 3],ones(1,2)*100*sum(eigenvalues(1:3))./sum(eigenvalues),'k:');
        %
        %         plotHorizon(gca,thVar,[]);
        %         plotVertical(gca,3,[]);
%         legend(legendMat,'location','best'); legend boxoff;
    end % for iAnimal=1:length(animalNm)
    
    
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
