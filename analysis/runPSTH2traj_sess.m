function runPSTH2traj_sess % (binSize,smthWidth)

% 1. project sessions' trial-averaged PSTH to traj
% 2. plot trajectory
% 3. scree plot

% only for H
% epoch: prior (S,L,SL) but common across S/L/SL (tBuffer=80);
%       Short480Eye Short480Hand Short560Eye Short560Hand Short640Eye Short640Hand Short720Eye Short720Hand Short800Eye Short800Hand Long800Eye Long800Hand Long900Eye Long900Hand Long1000Eye Long1000Hand Long1100Eye Long1100Hand Long1200Eye Long1200Hand ShortLong800Eye ShortLong800Hand ShortLong900Eye ShortLong900Hand ShortLong1000Eye ShortLong1000Hand ShortLong1100Eye ShortLong1100Hand ShortLong1200Eye ShortLong1200Hand 
% no attrition before PCA (if done, noisier at longer ts not considered)
% trial-averaged PSTH, not individual trials
% EH-specific (direction avg. when getting PSTH)
% curvature estimation separately

% Q: use projection matrix from all-session?

%% init
iAnimal=1;
initRSG2prior; % psthDir animalNm sessDir

cd(sessDir);
fnTmp=dir(sessDir);

%% parameters for dimreduce
binSize=10; % 5; % 20; %%%%%
smthWidth=40; % 25; % 40; % 40; %%%%%
optimD=[]; %%%%%
use_sqrt=0; %%%%%
proj_matrix=[];
keep_neurons=[]; 

nNeuron=[];

%% main
for i=1:length(fnTmp)
    if  strfind(fnTmp(i).name,'PSTH_') & strfind(fnTmp(i).name,'_prior_') & strfind(fnTmp(i).name,'.mat')
        tmp=fnTmp(i).name;
        load(tmp); % D nTrialExc tBuffer
        oldD=D;
        disp(['===== ' tmp ' =====']);
        nNeuron=[nNeuron;size(D(1).data,1)]
        
        for iHE=1:nEH
            disp(['--- ' ehNm{iHE} ' ---']);
%             try
            [D,proj_matrix,keep_neurons,eigenvalues]=PSTH2traj(D((iHE):2:length(D)),... % (iHE-1)*nEH+iTarg
                [],[],binSize,smthWidth,[],use_sqrt);
%             catch
%                 disp('');
%             end
            if ~isempty(D)            
                optimD=size(D(1).data,1);
                filenameTmp=['trajKS' tmp(5:(end-4)) '_' ehNm{iHE} '.mat'];
                save(filenameTmp,...
                    'D','proj_matrix','keep_neurons','binSize','smthWidth','optimD','use_sqrt','eigenvalues');

                % plot            
                try
                plotTraj(filenameTmp);
                catch
                    disp('');
                end

                savefig(gcf,['PC123_' filenameTmp(1:(end-4)) '.fig']); % ['periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3.fig']);
    %             exportfig(hFig,['PC123_' fnTmp(1:(end-4)) '.' optsExpFig.Format],optsExpFig);%             'periSet4_measOnly_' ehNm{iEH} targNm{iTh} '_PC1PC2PC3_splitTp.png'],optsExpFig);

                
            end
            D=oldD;
        end %         for iHE=1:nEH
%         waitforbuttonpress;         
        close all;
    end % if strfind(fnTmp(i).name,'prior')
end % for i=1:fnTmp





%% check PVAF from already done data set
% condNm={'E2','H2'}; % {'ER','EL','HR','HL'};
% legendMat=cell(size(condNm));
% load pplot.mat;
% fnTmp1='trajKS_periSet_';
% fnTmp2=['_' animalNm{i} '_bin' num2str(binSize) '_smth' num2str(smthWidth) '.mat']; % '_H.mat';
% nPCmax=25;
% thVar=.70;
% figure;ha;
% for i=1:length(condNm)
%     load([fnTmp1 condNm{i} fnTmp2]); % eigenvalues
%     disp([fnTmp1 condNm{i} fnTmp2]);
%     plot(cumsum(eigenvalues(1:nPCmax))./sum(eigenvalues),'.-','color',pplot.cmap{i});  
%     optD=find(cumsum(eigenvalues(1:20))./sum(eigenvalues)>thVar,1,'first');
%     disp(optD);
%     legendMat{i}=[condNm{i} ':' num2str(optD)];
% end
% xlabel('# PC'); ylabel('% Var. Explained'); 
% set(gca,'ytick',0:.25:1,'ylim',[0 1],'xlim',[1 nPCmax]);
% plotHorizon(gca,thVar,[]);
% legend(legendMat,'location','best'); legend boxoff;

