function runPSTH2traj_ts_tp_conditionSpecific_bootstrap

% 2018/6/19
% filename: xxx_boostrap > xxx_bootstrap
% only for tp, for now
%
% below are inherited from runPSTH2traj_prior_conditionSpecific_bootstrap..m
%   apply projecitonWeight/offsetMeanPSTH from original data, rather than doing PCA on bootstrapped data
%   correct wrong tInfo (12/7, 8/22, 8/23)

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

%% random projections ????
% % gaussian, unit length, orthogonalized
% nD=2;
% nRand=1000;


%%
idScree=0; % 1; % 0;
idPlotTraj=0; % 1; % 0;
idSave=1;

initRSG2prior; % psthDir animalNm
% cd(psthDir);
bootDir='/Users/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'; % '/Volumes/hansem/Desktop/bootstrapPSTH_RSG2prior_conditionSpecific/'; % macmini
cd(bootDir);
d=dir('PSTH_periSet*.mat'); % 2000 [1000 boot x 2 animals] 'PSTH_periSet_G_1000_full_SUMU'

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
    
    
for i=1:length(d) % 2000
    % extract info from file name
    if strfind(d(i).name,'_H_')
        iAnimal=1;
    else
        iAnimal=2;
    end
    iBoot=str2num(d(i).name(regexp(d(i).name,'\d'))); % 1... 1000
    
    load(d(i).name); % D[1x40].data([neuron x time]) nTrialExc
    
    % truncate for tp from periSet
    for iD=1:length(D)
        t0=D(iD).epochStarts(2)+1;
%         t1=tmpTs+tBuffer;
        D(iD).data=D(iD).data(:,t0:end); % t1);
        D(iD).epochStarts=1;
        D(iD).epochColors=D(iD).epochColors(1,:);
    end
    
    % main
    oldD=D;
    for iHE=1:nEH
        for iTarg=1:nTarg
            
            % load projectionWeight/offset from trajKS_tp_EL_G_bin20_smth40.mat
            traj=load(fullfile(psthDir,['trajKS_tp_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_bin20_smth40.mat'])); % keep_neurons, meanPSTH, proj_matrix
            
            % match cell ID
            % if not exist in original but exist in bootstrap, put zero
            nN=size(oldD(1).id,1);
            P0=zeros(nN,traj.optimD);
            M0=zeros(nN,1);
            K=false(nN,1);
            for iN=1:nN
                idN=oldD(1).id(iN,1)==traj.D(1).id(:,1) & oldD(1).id(iN,2)==traj.D(1).id(:,2); % [555 x 1]
                if sum(idN)>0,K(iN)=traj.keep_neurons(idN);end
                if sum(idN(traj.keep_neurons))>0 % [519 x 1]
                    P0(iN,:)=traj.proj_matrix(idN(traj.keep_neurons),:); % as size(proj_matrix)==size(keep_neurons)
                    M0(iN)=traj.meanPSTH(idN(traj.keep_neurons));
                end
            end        
            P0=P0(K,:); M0=M0(K);
            P=cat(3,P0,[M0 nan(size(P0,1),size(P0,2)-1)]);
            
            [D,proj_matrix,keep_neurons,eigenvalues]=PSTH2traj(D(((iHE-1)*nEH+iTarg):nCond:length(D)),...
                P,K,binSize,smthWidth,traj.optimD,use_sqrt); %  proj_matrix,keep_neurons optimD
            fnTmp=['trajKS_tp_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '_' num2str(iBoot) '.mat'];
            optimD=size(D(1).data,1);
                        
            if idSave %  ~exist(fnTmp,'file')
                save(fnTmp,...
                    'D','proj_matrix','keep_neurons','binSize','smthWidth','optimD','use_sqrt','eigenvalues','iBoot');
            end
            D=oldD;
            
            if idPlotTraj
                % plot
                plotTraj_figure3(fnTmp); % plotTraj_SFN17(fnTmp);
            end
        end % iTarg
    end % iHE
    
end % d=dir
    
    
    
%     % ts
%     for iAnimal=1:length(animalNm)
%         load(['PSTH_t_s_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D
%         
%         % main
%         oldD=D;
%         for iHE=1:nEH
%             for iTarg=1:nTarg
%                 
%                 fnTmp=['trajKS_ts_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
%                 load(fnTmp,'keep_neurons');
%                
%                 nNeuron=nnz(keep_neurons);
%                 
%                 for iRand=1:nRand
%                     % gaussian, unit length, orthogonalized
%                     randProj=randn(nNeuron,nD);
%                     randProjOrth=orth(randProj);
%                     proj_matrix=normalize(randProjOrth,1); % normCol(proj_matrix)
%                     
%                     [D,proj_matrix,keep_neurons,eigenvalues]=PSTH2traj(D(((iHE-1)*nEH+iTarg):nCond:length(D)),...
%                         proj_matrix,keep_neurons,binSize,smthWidth,optimD,use_sqrt); %  proj_matrix,keep_neurons optimD
%                     
%                     optimD=size(D(1).data,1);
%                     
%                                     fnTmp2=['trajKS_ts_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '_boot' num2str(iRand) '.mat'];
% 
% %                     if ~exist(fnTmp2,'file')
%                         save(fullfile(pwd,'randProj_trajKS_ts',fnTmp2),...
%                             'D','proj_matrix','keep_neurons','binSize','smthWidth','optimD','use_sqrt','eigenvalues');
% %                     end
%                     D=oldD;
%                     
%                     if idPlotTraj
%                         % plot
%                         plotTraj_figure3(fnTmp); % plotTraj_SFN17(fnTmp);
%                     end
%                 
%                 
%                 end % rand
%                 
%                 
%             end % iTarg
%         end % iHE
%         
%     end % for iAnimal=1:length(animalNm)

    
%     % check data
%     fnTmp='trajKS_ts_HL_G_bin20_smth40_boot1000.mat';
%     load(fnTmp);    
%     figure; ha; 
%     for i=1:length(D)
%         if i==length(D)
%             lw=3;
%         else
%             lw=1;
%         end
%         plot(D(i).data(1,:),D(i).data(2,:),'.-','linewidth',lw);
%     end
%     axis tight; xlabel('PC1'); ylabel('PC2');    
%     c=estCurvTraj(nm,dCurv); % [1 x S/L/SL] cell
    

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
    msize=2;
    lw=1;
    
    nPCmax=10;
    thVar=.75;
    
    nCmap=length(animalNm)*length(condNm);
    cmap=parula(nCmap); % hsv(nCmap);
    
    % ts
    icmap=1;
    fnTmp1='trajKS_ts_';
    figure;ha;setFigPos(1,6);

    for iAnimal=1:length(animalNm)
        fnTmp2=['_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
        for i=1:length(condNm)
            load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues
            disp([fnTmp1 condNm{i} fnTmp2]);
            plot(1:nPCmax,cumsum(eigenvalues(1:nPCmax))./sum(eigenvalues)*100,'o-','color',cmap(icmap,:),'markerfacecolor',cmap(icmap,:),'markeredgecolor',cmap(icmap,:),...
                'markersize',msize,'linewidth',lw);
            optD=find(cumsum(eigenvalues(1:20))./sum(eigenvalues)>thVar,1,'first');
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
        legend(legendMat,'location','best'); legend boxoff;
    end % for iAnimal=1:length(animalNm)
    
    % tp
    legendMat=cell(size(condNm));
    icmap=1;
    fnTmp1='trajKS_tp_';
    figure;ha;setFigPos(2,6);

    for iAnimal=1:length(animalNm)
        fnTmp2=['_' animalNm{iAnimal} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat'];
        for i=1:length(condNm)
            load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues
            disp([fnTmp1 condNm{i} fnTmp2]);
            plot(1:nPCmax,cumsum(eigenvalues(1:nPCmax))./sum(eigenvalues)*100,'o-','color',cmap(icmap,:),'markerfacecolor',cmap(icmap,:),'markeredgecolor',cmap(icmap,:),...
                'markersize',msize,'linewidth',lw);
            optD=find(cumsum(eigenvalues(1:20))./sum(eigenvalues)>thVar,1,'first');
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
        legend(legendMat,'location','best'); legend boxoff;
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
