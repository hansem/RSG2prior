function runExportDataHigh_singleSession_woOutlier_kilosort

% 2018/3/26
% re-run for sessions with incorrect tInfo: 12/7, 8/22, 8/23

% 2017/10/17
%  extract single-trial PSTH for all sessions & save  psthDataHigh/singleTrial

% 2017/9/25
%   effectorSpecific
%   only H for now
%   outlier removed
%   (TBD) trial-level

% trial wise (session specific)
% ex1_spikecounts: D(7*30).Data[61 400]
% 	61 neurons, 400ms spike trains, 1ms resolution, 7 conditions, 30 trials/condition
% D
%     data
%     condition
%     epochStarts
%     epochColors

% note
% for H only, for now
%  extract PSTH for all sessions (peri-set)
%  effectorSpecific (Avg. across modalities)

% (TBD) ts-modulated sessions
sessList=[170822;170823]; % 161207; %%%%%
% sess=[161204;161211;161218;121221];
sessId=[2; 8; 12; 15]; % sessions with # ts-mod cells>10

idSUonly=0; % 1; % 0; % 1; % 0; % 1;
nAnimal=1; % length(animalNm) % only H for now

%% exclusion crteria in exportDataHigh
nTrialExc=0; % 5; % if # trials per condition<nTrialExc, cell is removed from analysis

    tBuffer=80; % include 80ms before shortest ts of each prior
    nCondTmp=2; % EH

%%
initRSG2prior;

dirName=neuDir2; %     neuDir2='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/neuralData/kilosort';
    cd([dirName]);


    prNm=[prNm(:); 'ShortLong']; T{end+1}=T{end};tmpCmap{end+1,1}=tmpCmap{1,1};

    nEpoch=5+2+1; nEpoch2=5; epochId=[1:5 3 4 3];
    epochNm={'fixation','targetOn','periReady','periSet','production','t_s','t_p','prior'};
idEpoch=4;

%% 
load pplot.mat;
for iAnimal=1:nAnimal % length(animalNm):(-1):1
    % session
    for i=1:length(fn{iAnimal}) % 1:length(iSess) % sessId:sessId % 1:length(fname)
        fnm=fn{iAnimal}{i}; % H_RSGprior_20161203
        fid=fnm(end-5:end); % '161203'
        
        if sum(str2num(fid)==sessList)
        disp(['===== ' fnm ' =====']);
        load(fid); % kilosort, eg. 161203.mat: sp, idClust, idSU[# cluster x (id,idSU)], tInfo
                        
        % extract behavior
        beh=load([behDir animalNm{iAnimal} '_RSGprior_DMFC.mat']);
        idSess=str2num(fid)==beh.sessId; iSess=find(idSess); movingBallTrials=unique(tInfo(:,1)); % size same as idNoise
        if strcmp(fid,'161207'),idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
        nTmax=nnz(idSess);
        idOutSess=beh.idOut(idSess);
%         idOut2=false(size(tInfo,1),1);% remove outlier
%         movingBallTrials=unique(tInfo(:,1)); % size same as idNoise
%         iOutSess=find(beh.idOut(idSess));
%         for iOut=1:length(iOutSess)
%             idOut2(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
%         end
        
% consistent with plotSDFexport
        NEV=kilosort2NEV(sp,idClust,idSU);
        
        % loop through each unit
        nUnit=size(idSU,1); % both SU & MU
        
        D=[];
        iD=1;
        for iT=1:nTmax
            % remove outlier
            if ~idOutSess(iT)
                
                idPrior=2-beh.idShortTrial(iSess(iT));
                idHE=2-beh.idHandEye(iSess(iT));
                idT=find(T{idPrior}==beh.T(iSess(iT)));
                idTheta=(beh.theta(iSess(iT))==180)+1;
                
                D(iD).condition=[prNm{idPrior} ...
                    num2str(beh.T(iSess(iT))) ...
                    ehNm{idHE} ...
                    targNm{idTheta}];
                D(iD).epochColors=[tmpCmap{idPrior,1}(idT,:);...
                    tmpCmap{idPrior,1}(idT,:)]; % darker for 2nd phase of each epoch
                D(iD).epochStarts=[1 beh.T(iSess(iT))]; % ts; 1ms resolution
                D(iD).data=[];D(iD).id=[];
                D(iD).tp=beh.t(iSess(iT));
                D(iD).idShortTrial=beh.idShortTrial(iSess(iT));
                D(iD).ts=beh.T(iSess(iT));
                D(iD).idHandEye=beh.idHandEye(iSess(iT));
                D(iD).theta=beh.theta(iSess(iT));
                
                
                for iE=1:nUnit
                    if idSUonly & idSU(j,2)==0 % MU when idSUonly==1
                        
                    else
                        tmpId=[idSU(iE,1) idSU(iE,1) 1 length(movingBallTrials)]; % id(electrode,Unit,start trial,end trial)
                        try % exportDataHigh_conditionspecific
                            tmpD=exportDataHigh_conditionspecific(NEV,tmpId,tInfo,idEpoch,iT); % [1 x time]
%                             tmpD=exportDataHigh_effectorSpecific(NEV,tmpId,tInfo,idEpoch,iT); % [1 x time]
                        catch
                            disp('');
                        end
                        if size(tmpD,2)~=size(D(iD).data,2) && size(D(iD).data,2)>0 % debug
                            disp([num2str(size(tmpD)) ', ' num2str(size(D(iD).data))]);
                        end
                        D(iD).data=[D(iD).data;tmpD(:)']; % [neuron x time]
                        D(iD).id=[D(iD).id;str2num(fid) idSU(iE,:) ]; % [neuron x time]
                    end
                    
                end % for iE=1:nUnit
                
                disp(['trial#' num2str(iD) ': ts=' num2str(D(iD).ts) ', tp=' num2str(D(iD).tp) ',#spikes=' num2str(nnz(D(iD).data))]);
                iD=iD+1;
            end % if beh.t(iT)>0
        end % for iT=1:nTmax
    
    save([psthDir '/singleTrial/PSTH_' fid '_' epochNm{idEpoch} '_' animalNm{iAnimal} '_avgDir_SUMU.mat'],'D','-v7.3');
    
        end % if sum(fid==sessList)
    
    end % for i=1:length(fname)
    
end % iAnimal

% % check results
% figure;
% for i=1:length(D)
%     imagesc(D(i).data); title(['trial#' num2str(i)]);
%     drawnow;waitforbuttonpress;
% end

% % check results: looking at one cells' PSTH
% load('/Users/hansem/Dropbox (MIT)/psthDataHigh/singleTrial/PSTH_161221_periSet_H_avgDir_SUMU.mat')
% tmpId=find(D(1).id(:,2)==2155) % 17
% % load('/Users/hansem/Dropbox (MIT)/psthDataHigh/singleTrial/PSTH_161218_periSet_H_avgDir_SUMU.mat')
% % tmpId=find(D(1).id(:,2)==2036) % 17
% % load('/Users/hansem/Dropbox (MIT)/psthDataHigh/singleTrial/PSTH_161211_periSet_H_avgDir_SUMU.mat')
% % find(D(1).id(:,2)==1042) % 17
% tmp=[];j=1;for i=1:length(D),if D(i).ts==720,tmp{j}=D(i).data(tmpId,:);j=j+1;end;end;
% tmp2=fillNanCell(tmp);
% figure; imagesc(tmp2)
% tmp22=smoother(tmp2,40,1);
% figure; imagesc(tmp22);
% figure; plot(mean(tmp22));

% % debug: correcting idShortTrial
% initRSG2prior; % psthDir animalNm sessDir
% cd(singleTDir);
% fnTmp=dir(singleTDir);
% for i=1:length(fnTmp)
%     fn2=fnTmp(i).name;
%     if  strfind(fn2,'PSTH_') & strfind(fn2,'.mat')
%          load(fn2); % D (# trials).data [neurons x (ts+tp)] .id tp, idShortTrial ts idHandEye
%         disp(['===== ' fn2 ' =====']);
%         for i=1:length(D)
%             if strfind(D(i).condition,'Long')
%                 D(i).idShortTrial=0;
%             else
%                 D(i).idShortTrial=1;
%             end
%         end % for length (D)
%         save(fn2,'D');
%         clear D;
%     end
% end % for i=1:length(fnTmp)
        
% function D=initD
% iD=1;
% for i=1:nPr
%     for j=1:length(T{i})
%         for k=1:nEH
%             for m=1:nTarg
%                 D(iD).condition=[prNm{i} num2str(T{i}(j)) '/' ehNm{k} '/' targNm{m}];
%                 D(iD).data=[];
%                 iD=iD+1;
%             end % m
%         end % k
%     end % j
% end % i
% 
% 
% function chooseBestSession
% 
% %% best sessions for single-trial analysis
% %% pick up common trials across cells
% listSess=[161221;...
%       161220;...
%       161222;...
%       161223;...
%       161211];
%   load pplot.mat;
%   for i=1:length(listSess)
%       idSess=fidMat(:,1)==listSess(i);
%       nC=nnz(idSess);
%       nT=max(cnameMat(fidMat(:,1)==listSess(i),end));
%       tmp=zeros(nC,nT);
%       for j=1:nC
%           startT=cnameValid{strcmp(fname,['H_RSGprior_20' num2str(listSess(i))])}(j,end-1);
%           endT=cnameValid{strcmp(fname,['H_RSGprior_20' num2str(listSess(i))])}(j,end);
%           tmp(j,startT:endT)=1;
%       end
%       figure; set(gcf,'position',pplot.(['rect1_' num2str(i)]));
%       imagesc(tmp); colormap('gray');
%       xlabel('trials'); ylabel('cell');
%       
%       edges=[1:nT]-0.5;
%       nStartT=histc(fidMat(idSess,end-1),edges);
%       nTlength=histc(fidMat(idSess,end)-fidMat(idSess,end-1),edges);
%       nEndT=histc(fidMat(idSess,end),edges);
%       figure; set(gcf,'position',pplot.(['rect2_' num2str(i)]));
%       plot(1:nT,cumsum(nStartT),'r',1:nT,cumsum(nTlength),'g',1:nT,cumsum(nEndT),'b'); 
%       axis tight; xlabel('trials'); ylabel('cumulative # cells');
% %       figure; plot([1:nT]'.*(max(cumsum(nTlength))-cumsum(nTlength)));
%   end

