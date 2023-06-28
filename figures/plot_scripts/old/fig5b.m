function fig5b

% 2019/4/11
% save overlap distance into '/Users/hansem/Dropbox (MIT)/fileFromNHPrig/neuralData/kilosort_wf/overlapDistSp.mat'

% 2019/3/23
% stat: 1) no modulation during prior support, 2) ready transient, 3) different speed between short vs long
% figure 3d: discretize distance b/t two priors (from ready 240, 240, 320-prior support, 240, 240, 240)
% supplementary figure for speed during set-go

% 2019/3/12
% fig5b: speed (combined across animals)
% sFig XX: overlap distance

% 2019/1/29
% measure speed in full state space
% binningSmoothing PSTH > speed > avgAcrossTs with attrition
% idAvgBeforeSp=1: avgAcrossTs with attrition > binningSmoothing PSTH > speed
% idUseStencil, nP: # of points in the finite difference method

% 2019/1/24
% measuring distance between two overlap ts across prior
% in full state space (noisy): PSTH_periSet_H/G_full_SUMU.mat
% modified from runPSTH2traj_ts_tp_conditionSpecific.m
% alternative: do PCA over ts&tp together-  trajKS_periSet_XX_X.mat

%%

binSize=20;
smthWidth=40; 

overlapT=800;

initRSG2prior; % psthDir animalNm
cd(psthDir);

condNm={'ER','EL','HR','HL'};
condNm2={'Eye Left','Eye Right','Hand Left','Hand Right'};
nCond=numel(condNm);

% key data: distance
dist=cell(length(animalNm)*nEH*nTarg,1);
distPC=cell(length(animalNm)*nEH*nTarg,1);

win=[160 160 160 160 160 160 160 160 160 160]; % [480 320 480]; % [160 160 160 160 160 160 160 160 160 160]; % [160 160 160 320 160 160 160 160 160]; %  [240 240 320 240 240 240]; % bin size for figure 3d
iStart=cumsum(win/binSize)-win/binSize+1; % 1         241         481         801        1041        1281
iEnd=cumsum(win/binSize);          %   240         480         800        1040        1280        1520
dist2=nan(length(animalNm)*nEH*nTarg,length(win));
distPC2=nan(length(animalNm)*nEH*nTarg,length(win));
idMean=1; % if 1, Root-mean-sqaured-difference; if 0, Root-sum-squared-difference

% speed
idAvgBeforeSp=0; % 0; % 1;
idUseStencil=0; % 1; % 0;
nP=5; % 3; % 5; % 3;
if idUseStencil==0, nP=1; end;
idSanityCheck=0; % 1; % 0;
sp=cell(length(animalNm)*nEH*nTarg,2); % last for ts/tp separately; aligned to ready/set
smthWidthSp= 20; % 80; % 40; % 20; % if 0, no smoothing

% save for laminar project (keep condition-specificity)
saveDir='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/neuralData/kilosort_wf/';
distCell=cell(length(animalNm)*nEH*nTarg,1);
spCell=cell(length(animalNm)*nEH*nTarg,nPr*nTspp);
idDist=cell(length(animalNm)*nEH*nTarg,1); % [neuron x session/neuronId/idSU]
idSp=cell(length(animalNm)*nEH*nTarg,1); % [neuron x session/neuronId/idSU]

% main
for iAnimal=1:length(animalNm)    
    
    %% in full state space
    load(['PSTH_periSet_' animalNm{iAnimal} '_full_SUMU.mat']); % PSTH_t_s_H_SUMU.mat'); % D(40) nTrialExc tBuffer
    disp(animalNm{iAnimal});
    
    % find id for overlap among 40 conditions
    tmpId=nEH*nTarg*(nTs-1)+[1:8]; % 17:24  short800ER>EL>HR>HL>Long800...
    tmpD0=D(tmpId);
    
    % smoothing & binning
    tmpD=smoothBin(tmpD0,binSize,smthWidth);
    
    for iHE=1:nEH
        for iTarg=1:nTarg
            disp(tmpD(1).condition);
            iEHTarg=(iHE-1)*nTarg+iTarg; % 1:4
            iTmp=(iAnimal-1)*nEH*nTarg+iEHTarg; % 1:8
            
            Ds=tmpD(iEHTarg); % [neuron x
            Dl=tmpD(iEHTarg+nEH*nTarg);
            
            % find minimum truncated by shorter tp
            minLength=min([size(Ds.data,2) size(Dl.data,2)]);
            
            % dist
            dist{iTmp}=sqrt(sum((Ds.data(:,1:minLength)-Dl.data(:,1:minLength)).^2,1))'; % vector
            
            distCell{iTmp}=Ds.data(:,1:minLength)-Dl.data(:,1:minLength); % short-long [neuron x time]
            idDist{iTmp}=Ds.id;
        end
    end
    
    % figure 3d: discretize distance b/t two priors (from ready 240, 240, 320-prior support, 240, 240, 240)
    % no smoothing, binSize win=[240 240 320 240 240 240];
    % for the last bin, use minimum 
    tmpD=smoothBin2(tmpD0,win);
    for iHE=1:nEH
        for iTarg=1:nTarg
            iEHTarg=(iHE-1)*nTarg+iTarg; % 1:4
            iTmp=(iAnimal-1)*nEH*nTarg+iEHTarg; % 1:8
            
            Ds=tmpD(iEHTarg); % [neuron x nWin]
            Dl=tmpD(iEHTarg+nEH*nTarg);
            
            % find minimum truncated by shorter tp
            minLength=min([size(Ds.data,2) size(Dl.data,2)]);
            
            % dist
            if idMean
                dist2(iTmp,1:minLength)=sqrt(mean((Ds.data(:,1:minLength)-Dl.data(:,1:minLength)).^2,1))'; % vector
            else
                dist2(iTmp,1:minLength)=sqrt(sum((Ds.data(:,1:minLength)-Dl.data(:,1:minLength)).^2,1))'; % vector
            end
            
        end
    end
    
    %% SPEED
    disp('===== speed =====');
    if idAvgBeforeSp %  avgAcrossTs with attrition > binningSmoothing PSTH > speed
        
        for iHE=1:nEH
            for iTarg=1:nTarg
                iEHTarg=(iHE-1)*nTarg+iTarg; % 1:4
                iTmp=(iAnimal-1)*nEH*nTarg+iEHTarg; % 1:8
                
                % separating ts/tp
                iTmpD=1;
                for iD=iEHTarg:nCond:length(D) % for each pr/ts
                    
                    Dts(iTmpD)=D(iD); Dtp(iTmpD)=D(iD);
                    tmpTs=D(iD).epochStarts(2); % ts/binSize
                    disp([D(iD).condition ', ' num2str(tmpTs)]);
                    Dts(iTmpD).data=Dts(iTmpD).data(:,1:tmpTs);
                    Dtp(iTmpD).data=Dtp(iTmpD).data(:,tmpTs:end);
                    Dts(iTmpD).epochStarts=1; Dtp(iTmpD).epochStarts=1;
                    iTmpD=iTmpD+1;
                 end
                
                % avgAcrossTs with attrition
                tmpD=fillNan(Dts(1:nTspp)); % % for each pr/ts
                tmpD2=fillNan(Dts(nTspp+[1:nTspp])); % % for each pr/ts
                for iD=1:nTspp
                    Dts(iD).data=squeeze(nanmean(cat(3,tmpD.data),3)); % [neuron x time x 5tsShort] >  [neuron x time]
                    Dts(nTspp+iD).data=squeeze(nanmean(cat(3,tmpD2.data),3));  % [neuron x time x 5tsLong] >  [neuron x time]
                end
                
                % binningSmoothing
                Dts=smoothBin(Dts,binSize,smthWidthSp);
                Dtp=smoothBin(Dtp,binSize,smthWidthSp);
                
                % sanity check
                if idSanityCheck
                    % find(Dts(1).id(:,1)==161207 & Dts(1).id(:,2)==2055 & Dts(1).id(:,3)==1) 164
                    idNeuron=164; disp(Dts(1).id(idNeuron,:));
                    figure;  cmap=[tmpCmap{1,1}; tmpCmap{2,1}];
                    for iD=1:(nPr*nTspp) % for each pr/ts
                        plot(Dts(iD).data(idNeuron,:),'color',cmap(iD,:),'linewidth',2);ha;
                    end
                    figure; % sanity check
                    for iD=1:(nPr*nTspp) % for each pr/ts
                        if idUseStencil
                            plot(sqrt(sum(deriv(Dts(iD).data(idNeuron,:)',binSize/1000,1,nP).^2,1)),'color',cmap(iD,:),'linewidth',2);ha; % sanity check
                        else
                            plot(sqrt(sum(diff(Dts(iD).data(idNeuron,:),1,2).^2,1)),'color',cmap(iD,:),'linewidth',2);ha; % sanity check
                        end
                    end
                end
                
                % speed
                sp{iTmp,1}=cell(1,1); % ts
                sp{iTmp,2}=cell(1,1); % tp
                
                for iD=1:(nPr*nTspp) % for each pr/ts
                    if idUseStencil
                        tmpSp1=sqrt(sum(deriv(Dts(iD).data',binSize/1000,1,nP).^2,1)); % [1 x time-1]
                        tmpSp2=sqrt(sum(deriv(Dtp(iD).data',binSize/1000,1,nP).^2,1)); % [1 x time-1]
                        
                    else
                        tmpSp1=sqrt(sum(diff(Dts(iD).data,1,2).^2,1)); % [1 x time-1]
                        tmpSp2=sqrt(sum(diff(Dtp(iD).data,1,2).^2,1)); % [1 x time-1]
                    end
                    %                 if iD==1 % isempty(sp{iTmp,1})
                    %                     sp{iTmp,1}=tmpSp1(:)';
                    %                 else
                    sp{iTmp,1}=[sp{iTmp,1}; tmpSp1(:)'];
                    %                 end
                    %                 if iD==1 % isempty(sp{iTmp,2})
                    %                     sp{iTmp,2}=tmpSp2(:)';
                    %                 else
                    sp{iTmp,2}=[sp{iTmp,2}; tmpSp2(:)'];
                    %                 end
                end % iD
                sp{iTmp,1}=sp{iTmp,1}(2:end);
                sp{iTmp,2}=sp{iTmp,2}(2:end);
                
                % sanity check
                if idSanityCheck
                    if idUseStencil
                        figure;setFigPos(1,4); imagesc(deriv(Dts(iD).data',binSize/1000,1,nP).^2); colorbar;applytofig4keynote
                        figure;setFigPos(1,5); plot(transpose(deriv(Dts(iD).data',binSize/1000,1,nP).^2)); applytofig4keynote
                         tmp=deriv(Dts(iD).data',binSize/1000,1,nP).^2;
                        figure;setFigPos(1,6); imagesc(normalize(tmp,2)); colorbar;
                        figure;setFigPos(2,6); plot(sqrt(sum(normalize(tmp,2).^2,1))); colorbar;
                    else
                        figure; imagesc(diff(Dts(iD).data,1,2).^2); colorbar;applytofig4keynote
                        figure; plot(transpose(diff(Dts(iD).data,1,2).^2)); applytofig4keynote
                    end
%                     [x,y]=find(diff(Dts(iD).data,1,2).^2>0.01);
% %                     x(y>50)
%                     Dts(iD).id(x(y>50),:)
%
                    % 161218        1172           0
                    % 161221        1165           0
                    % 161223        1229           1
                    % 161206        1228           0
                    % 161219        1151           1
                    %
                    % 161206        1228           0
                    % 161214        2231           0
                    % 161218        1114           0
                    % 161218        1167           1  0
                    
                    idNeuron=500; % 499; % find(Dts(iD).id(:,1)==161223 & Dts(iD).id(:,2)==1229)
%                                         idNeuron=335; % find(Dts(iD).id(:,1)==161219 & Dts(iD).id(:,2)==1151)
                    %                     idNeuron=284; % find(Dts(iD).id(:,1)==161218 & Dts(iD).id(:,2)==1114)
                    %                     idNeuron=250; % find(Dts(iD).id(:,1)==161214 & Dts(iD).id(:,2)==2159)
                    figure;cmap=[tmpCmap{1,1}; tmpCmap{2,1}];
                    tmpSpM=cell(nPr,1);
                    for iD=1:nTspp:(nTspp+1) % eye right 1,6
                        iDtmp=(iD-1)/nTspp+1; % 1,2
                        if idUseStencil
                            tmpSp=sqrt(sum(deriv(Dts(iD).data(idNeuron,:)',binSize/1000,1,nP).^2,1));tmpSpM{iDtmp}=tmpSp(:)';
                        else
                            tmpSp=sqrt(sum(diff(Dts(iD).data(idNeuron,:),1,2).^2,1));tmpSpM{iDtmp}=tmpSp(:)';
                        end
                        setFigPos(1,2); % yyaxis left;
                        subplot(2,1,1); plot(transpose(tmpSp),'color',cmap(iD,:)); ha; % FR
                        % yyaxis right;
                        subplot(2,1,2); plot(transpose((Dts(iD).data(idNeuron,:))),'color',cmap(iD,:)); ha;
                    end
                    applytofig4keynote
                    
                    
                end
                
            end % iTarg
        end % iEH
        
        
     %%   USED
    else % idAvgBeforeSp==0: binningSmoothing PSTH > speed > avgAcrossTs with attrition
        % smoothing & binning
        Dsp=smoothBin(D,binSize,smthWidthSp);
        for iHE=1:nEH
            for iTarg=1:nTarg
                iEHTarg=(iHE-1)*nTarg+iTarg; % 1:4
                iTmp=(iAnimal-1)*nEH*nTarg+iEHTarg; % 1:8
                sp{iTmp,1}=cell(1,1); % ts
                sp{iTmp,2}=cell(1,1); % tp
                
                iiD=1;
                for iD=iEHTarg:nCond:length(Dsp) % for each pr/ts
                    tmpTs=Dsp(iD).epochStarts(2); % ts/binSize
                    if idUseStencil
                        tmpSp1=sqrt(sum(deriv(Dsp(iD).data(:,1:tmpTs)',binSize/1000,1,nP).^2,1)); % [1 x time-1]
                        tmpSp2=sqrt(sum(deriv(Dsp(iD).data(:,tmpTs:end)',binSize/1000,1,nP).^2,1)); % [1 x time-1]
                    else %%%%% used
                        tmpSp1=sqrt(sum(diff(Dsp(iD).data(:,1:tmpTs),1,2).^2,1)); % [1 x time-1]
                        tmpSp2=sqrt(sum(diff(Dsp(iD).data(:,tmpTs:end),1,2).^2,1)); % [1 x time-1]
                    end
                    %                 if iD==1 % isempty(sp{iTmp,1})
                    %                     sp{iTmp,1}=tmpSp1(:)';
                    %                 else
                    sp{iTmp,1}=[sp{iTmp,1}; tmpSp1(:)'];
                    %                 end
                    %                 if iD==1 % isempty(sp{iTmp,2})
                    %                     sp{iTmp,2}=tmpSp2(:)';
                    %                 else
                    sp{iTmp,2}=[sp{iTmp,2}; tmpSp2(:)'];
                    %                 end
                    
                    % saving for laminar project
                    spCell{iTmp,iiD}=sqrt(diff(Dsp(iD).data,1,2).^2); % [neuron x time] for both ts/tp
                    idSp{iTmp}=Dsp(iD).id;
                    iiD=iiD+1;
                end % iD
                
                sp{iTmp,1}=sp{iTmp,1}(2:end);
                sp{iTmp,2}=sp{iTmp,2}(2:end);
                
                % sanity check
                if idSanityCheck
                    if idUseStencil
                        figure;setFigPos(1,4); imagesc(deriv(Dsp(iD).data(:,1:tmpTs)',binSize/1000,1,nP).^2); colorbar;applytofig4keynote
                        figure;setFigPos(1,5); plot(transpose(deriv(Dsp(iD).data(:,1:tmpTs)',binSize/1000,1,nP).^2)); applytofig4keynote
                         tmp=deriv(Dsp(iD).data(:,1:tmpTs)',binSize/1000,1,nP).^2;
                        figure;setFigPos(1,6); imagesc(normalize(tmp,2)); colorbar;
                        figure;setFigPos(2,6); plot(sqrt(sum(normalize(tmp,2).^2,1)));
                    else
                        figure; imagesc(diff(Dsp(iD).data,1,2).^2); colorbar;applytofig4keynote
                        figure; plot(transpose(diff(Dsp(iD).data,1,2).^2)); applytofig4keynote
                    end
                                        idNeuron=506; % find(Dsp(iD).id(:,1)==161223 & Dsp(iD).id(:,2)==1229)
%                                         idNeuron=340; % find(Dsp(iD).id(:,1)==161219 & Dsp(iD).id(:,2)==1151)
                    %                     idNeuron=287; % find(Dsp(iD).id(:,1)==161218 & Dsp(iD).id(:,2)==1114)
                    %                     idNeuron=253; % find(Dsp(iD).id(:,1)==161214 & Dsp(iD).id(:,2)==2159)
                    figure;cmap=[tmpCmap{1,1}; tmpCmap{2,1}];
                    tmpSpM=cell(nPr*nTspp,1);
                    for iD=1:nCond:length(Dsp) % eye right
                        iDtmp=(iD-1)/nCond+1;
                        tmpTs=Dsp(iD).epochStarts(2);
                        if idUseStencil
                            tmpSp=sqrt(deriv(Dsp(iD).data(idNeuron,1:tmpTs)',binSize/1000,1,nP).^2);tmpSpM{iDtmp}=tmpSp(:)';
                        else
                            tmpSp=sqrt(diff(Dsp(iD).data(idNeuron,1:tmpTs),1,2).^2);tmpSpM{iDtmp}=tmpSp(:)';
                        end
                        setFigPos(1,2); % yyaxis left;
%                         subplot(2,1,1); plot(transpose(tmpSp),'color',cmap(iDtmp,:)); ha; % FR
                        % yyaxis right;
                        subplot(2,1,2); plot(transpose((Dsp(iD).data(idNeuron,1:tmpTs))),'color',cmap(iDtmp,:)); ha;
                    end
                    tmpSp=fillNanCell(tmpSpM);
                    subplot(2,1,1);ha; shadedErrorBar([],nanmean(tmpSp(1:nTspp,:),1),sem(tmpSp(1:nTspp,:),1,'omitnan'),'r',1);
                    shadedErrorBar([],nanmean(tmpSp(nTspp+[1:nTspp],:),1),sem(tmpSp(1:nTspp,:),1,'omitnan'),{'b','linewidth',3},1)
                    applytofig4keynote
                    
                    % 161218        1172           0
                    % 161221        1165           0
                    % 161223        1229           1
                    % 161206        1228           0
                    % 161219        1151           1
                    %
                    % 161206        1228           0
                    % 161214        2231           0
                    % 161218        1114           0
                    % 161218        1167           1  0
                    
%                     figure; setFigPos(1,1); imagesc(diff(Dsp(iD).data(:,1:tmpTs),1,2).^2); colorbar;applytofig4keynote
%                     [x,y]=find(diff(Dsp(iD).data,1,2).^2>0.01)
%                     x(y>50)
                end % sanity check
                

            end % iTarg
        end % iEH
    end % if idAvgBeforeSp
    
    %% PCA
    for iHE=1:nEH
        for iTarg=1:nTarg
            fnTmp=['trajKS_periSet_' ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal} '.mat']; %  '_bin' num2str(binSize) '_smth' num2str(smthWidth)
            disp([ehNm{iHE}(1) targNm{iTarg}(1) '_' animalNm{iAnimal}]);
            
            load(fnTmp); % D optimD
            disp(['optimD: ' num2str(optimD)]);
            
            iCond=1+(nEH*nTarg); % short800
            iCond2=1+(nEH*nTarg)+1; % long800
            disp([D(iCond).condition ' vs ' D(iCond2).condition]);
            
             iEHTarg=(iHE-1)*nTarg+iTarg; % 1:4
            iTmp=(iAnimal-1)*nEH*nTarg+iEHTarg; % 1:8
            
            % find minimum truncated by shorter tp
            minLength=min([size(D(iCond).data,2) size(D(iCond2).data,2)]);
            
            tmpDist=sqrt(sum((D(iCond).data(:,1:minLength)-D(iCond2).data(:,1:minLength)).^2,1))';
            distPC{iTmp}=tmpDist; % vector
            
            for iWin=1:length(win)
                if iStart(iWin)<=length(tmpDist) % safe
                    iEndTmp=min([iEnd(iWin) length(tmpDist)]);
                    distPC2(iTmp,iWin)=mean(tmpDist(iStart(iWin):iEndTmp));
                     
                end
            end
        end % iTarg
    end % iHE
    
end % for iAnimal=1:length(animalNm)

%% save for laminar project
save(fullfile(saveDir,'overlapDistSp.mat'),'idSp','idDist','distCell','spCell','binSize','smthWidth');

%% plot speed
% sp=cell(length(animalNm)*nEH*nTarg,2); % last for ts/tp separately; aligned to ready/set

% ts (short vs long)
figure; setFigPos(2,1); ha; % speed(whole ts): short vs long
mV{1}=cell(1,1); mV{2}=cell(1,1);
for i=1:size(sp,1) % 1:size(sp,1)/2 %(1+size(sp,1)/2):size(sp,1) %  1:size(sp,1)
    vTs=fillNanCell(sp{i,1}); % [10ts x time]
    vTss=nanmean(vTs(1:nTspp,:),1); % [1 x time]
    vTsl=nanmean(vTs(nTspp+[1:nTspp],:),1); % [1 x time]
    mV{1}=[mV{1}; vTss];mV{2}=[mV{2}; vTsl];
    if idUseStencil
        tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(length(vTss)-1))]; % (binSize/2)+[0:binSize:(binSize*(length(vTss)-1))];
        plot(tmpT,vTss,'r-');
        tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(length(vTsl)-1))]; % (binSize/2)+[0:binSize:(binSize*(length(vTsl)-1))];
        plot(tmpT,vTsl,'b-');
    else
        tmpT=(binSize/2)+[0:binSize:(binSize*(length(vTss)-1))];
        plot(tmpT,vTss,'r-');
        tmpT=(binSize/2)+[0:binSize:(binSize*(length(vTsl)-1))];
        plot(tmpT,vTsl,'b-');
    end
end
% plot avg
mmv1=fillNanCell(mV{1}); % [8data x time]
if idUseStencil
    tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(size(mmv1,2)-1))]; % %(binSize/2)+[0:binSize:(binSize*(size(mmv1,2)-1))]; % 10:20:....
else
    tmpT=(binSize/2)+[0:binSize:(binSize*(size(mmv1,2)-1))]; % 10:20:....
end
shadedErrorBar(tmpT(:)',nanmean(mmv1,1),sem(mmv1,1,'omitnan'),{'-','color','r','linewidth',2},1); drawnow;
mmv2=fillNanCell(mV{2}); % [8data x time]
if idUseStencil
    tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(size(mmv2,2)-1))]; %(binSize/2)+[0:binSize:(binSize*(size(mmv2,2)-1))]; % 10:20:....
else
    tmpT=(binSize/2)+[0:binSize:(binSize*(size(mmv2,2)-1))]; % 10:20:....
end
shadedErrorBar(tmpT(:)',nanmean(mmv2,1),sem(mmv2,1,'omitnan'),{'-','color','b','linewidth',2},1); drawnow;

plotVertical(gca,[T{1}(1) overlapT T{2}(end)],[]);
set(gca,'xtick',[0 T{1}(1:2:3) T{2}(1:2:end)],'tickdir','out','ticklength',[0.015 0.015]); %  overlapT overlapT+200]);
xlabel('Time from ready (ms)');
ylabel('Speed (a.u.)');

% stat: 1) no modulation during prior support, 2) ready transient, 3) different speed between short vs long
nBin=[320 400]/binSize; % 16 20
B=nan(nAnimal*nEH*nTarg,nPr); % stat: 1) no modulation during prior support
nBinReady=200/binSize; % 10
mV_ReadySupport=nan(nAnimal*nEH*nTarg*nPr,2); % ready vs priorSupport
for iPr=1:nPr
    tmp=cell2mat(mV{iPr}); % [8 x 59]
    tmp=tmp(:,~isnan(tmp(1,:))); % remove nan
    % 2) ready transient
    mV_ReadySupport(nAnimal*nEH*nTarg*(iPr-1)+[1:(nAnimal*nEH*nTarg)],1)=mean(tmp(:,1:nBinReady),2);
    tmp=tmp(:,(end-nBin(iPr)+1):end); % only durong prior support
    mV_ReadySupport(nAnimal*nEH*nTarg*(iPr-1)+[1:(nAnimal*nEH*nTarg)],2)=mean(tmp,2);
    for iTmp=1:size(tmp,1)
        tmpB=regress(tmp(iTmp,:)',[[1:size(tmp,2)]' ones(size(tmp,2),1)]);
        B(iTmp,iPr)=tmpB(1);
    end
end
% stat: 1) no modulation during prior support
signrank(B(:,1)) % short
signrank(B(:,2)) % long
% 2) ready transient
signrank(mV_ReadySupport(:,1),mV_ReadySupport(:,2))
% 3) different speed b/t short vs long
signrank(mV_ReadySupport(1:nAnimal*nEH*nTarg,2),mV_ReadySupport(nAnimal*nEH*nTarg+[1:nAnimal*nEH*nTarg],2))

% ts (short vs long) H
figure; setFigPos(1,5); ha; % speed(whole ts): short vs long
mV{1}=cell(1,1); mV{2}=cell(1,1);
for i=1:size(sp,1)/2 %(1+size(sp,1)/2):size(sp,1) %  1:size(sp,1)
    vTs=fillNanCell(sp{i,1}); % [10ts x time]
    vTss=nanmean(vTs(1:nTspp,:),1); % [1 x time]
    vTsl=nanmean(vTs(nTspp+[1:nTspp],:),1); % [1 x time]
    mV{1}=[mV{1}; vTss];mV{2}=[mV{2}; vTsl];
    
    tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(length(vTss)-1))]; % (binSize/2)+[0:binSize:(binSize*(length(vTss)-1))];
    plot(tmpT,vTss,'r-');
    tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(length(vTsl)-1))]; % (binSize/2)+[0:binSize:(binSize*(length(vTsl)-1))];
    plot(tmpT,vTsl,'b-');
end
% plot avg
mmv1=fillNanCell(mV{1}); % [8data x time]
tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(size(mmv1,2)-1))]; % %(binSize/2)+[0:binSize:(binSize*(size(mmv1,2)-1))]; % 10:20:....
shadedErrorBar(tmpT(:)',nanmean(mmv1,1),sem(mmv1,1,'omitnan'),{'-','color','r','linewidth',2},1); drawnow;
mmv2=fillNanCell(mV{2}); % [8data x time]
tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(size(mmv2,2)-1))]; %(binSize/2)+[0:binSize:(binSize*(size(mmv2,2)-1))]; % 10:20:....
shadedErrorBar(tmpT(:)',nanmean(mmv2,1),sem(mmv2,1,'omitnan'),{'-','color','b','linewidth',2},1); drawnow;

plotVertical(gca,[T{1}(1) overlapT],[]);
set(gca,'xtick',[T{1}(1:3) T{2}(1:3:end)]); %  overlapT overlapT+200]);
xlabel('time from ready (ms)');
ylabel('speed');

% ts (short vs long) G
figure; setFigPos(1,6); ha; % speed(whole ts): short vs long
mV{1}=cell(1,1); mV{2}=cell(1,1);
for i=(1+size(sp,1)/2):size(sp,1) %  1:size(sp,1)
    vTs=fillNanCell(sp{i,1}); % [10ts x time]
    vTss=nanmean(vTs(1:nTspp,:),1); % [1 x time]
    vTsl=nanmean(vTs(nTspp+[1:nTspp],:),1); % [1 x time]
    mV{1}=[mV{1}; vTss];mV{2}=[mV{2}; vTsl];
    
    tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(length(vTss)-1))]; % (binSize/2)+[0:binSize:(binSize*(length(vTss)-1))];
    plot(tmpT,vTss,'r-');
    tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(length(vTsl)-1))]; % (binSize/2)+[0:binSize:(binSize*(length(vTsl)-1))];
    plot(tmpT,vTsl,'b-');
end
% plot avg
mmv1=fillNanCell(mV{1}); % [8data x time]
tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(size(mmv1,2)-1))]; % %(binSize/2)+[0:binSize:(binSize*(size(mmv1,2)-1))]; % 10:20:....
shadedErrorBar(tmpT(:)',nanmean(mmv1,1),sem(mmv1,1,'omitnan'),{'-','color','r','linewidth',2},1); drawnow;
mmv2=fillNanCell(mV{2}); % [8data x time]
tmpT=(binSize*round(nP/2)-binSize/2)+[0:binSize:(binSize*(size(mmv2,2)-1))]; %(binSize/2)+[0:binSize:(binSize*(size(mmv2,2)-1))]; % 10:20:....
shadedErrorBar(tmpT(:)',nanmean(mmv2,1),sem(mmv2,1,'omitnan'),{'-','color','b','linewidth',2},1); drawnow;

plotVertical(gca,[T{1}(1) overlapT],[]);
set(gca,'xtick',[T{1}(1:3) T{2}(1:3:end)]); %  overlapT overlapT+200]);
xlabel('time from ready (ms)');
ylabel('speed');

% ts (for each ts)
figure; setFigPos(2,5); ha; % speed(whole tp): for each ts
cmap=[tmpCmap{1,1}; tmpCmap{2,1}];
subplot(2,1,1); ha; % H
for i=1:(size(sp,1)/2) % 1:size(sp,1) % 8data sets
    vTp=fillNanCell(sp{i,1}); % [10ts x time]
    for j=1:size(vTp,1) % 10ts
        tmpT=(binSize/2)+[0:binSize:(binSize*(length(vTp(j,:))-1))];
        plot(tmpT,vTp(j,:),'-','linewidth',1,'color',cmap(j,:));
    end
end
% plot avg
ylabel('speed');
subplot(2,1,2); ha; % G
for i=(1+size(sp,1)/2):size(sp,1) % 1:size(sp,1) % 8data sets
    vTp=fillNanCell(sp{i,1}); % [10ts x time]
    for j=1:size(vTp,1) % 10ts
        tmpT=(binSize/2)+[0:binSize:(binSize*(length(vTp(j,:))-1))];
        plot(tmpT,vTp(j,:),'-','linewidth',1,'color',cmap(j,:));
    end
end
ylabel('speed');
xlabel('time from ready (ms)');

% tp
vTpMat=cell(nAnimal,nTs*nPr,(size(sp,1)/2));
figure; setFigPos(2,6); ha; % speed(whole tp): for each ts
cmap=[tmpCmap{1,1}; tmpCmap{2,1}];
subplot(2,1,1); ha;
for i=1:(size(sp,1)/2) % 1:size(sp,1) % 8data sets
    vTp=fillNanCell(sp{i,2}); % [10ts x time]
    for j=1:size(vTp,1) % 10ts
        tmpT=(binSize/2)+[0:binSize:(binSize*(length(vTp(j,:))-1))];
        plot(tmpT,vTp(j,:),'-','linewidth',1,'color',cmap(j,:));
        vTpMat{1,j,i}=vTp(j,:);
    end
end
% plot avg
ylabel('speed');
subplot(2,1,2); ha;
for i=(1+size(sp,1)/2):size(sp,1) % 1:size(sp,1) % 8data sets
    vTp=fillNanCell(sp{i,2}); % [10ts x time]
    for j=1:size(vTp,1) % 10ts
        tmpT=(binSize/2)+[0:binSize:(binSize*(length(vTp(j,:))-1))];
        plot(tmpT,vTp(j,:),'-','linewidth',1,'color',cmap(j,:));
        vTpMat{2,j,i-size(sp,1)/2}=vTp(j,:);
    end
end
ylabel('speed');
xlabel('time from set (ms)');

% sFig 13: averaged across conditions, smoothing
for iAnimal=1:nAnimal
%     figure; setFigPos(2,iAnimal); ha; box off;
    for iTs=1:(nTspp*nPr)
        if iTs==1 || iTs==(nTspp+1)
            figure; setFigPos(2,iAnimal); ha; box off;
        end
        mVts=nanmean(fillNanCell(vTpMat(iAnimal,iTs,:)),1); % 1 x timePoint
        smVts=smoother(mVts,binSize,binSize);
        tmpT=(binSize/2)+[0:binSize:(binSize*(length(smVts)-1))];
        plot(tmpT,smVts,'-','linewidth',1,'color',cmap(iTs,:));
        set(gca,'xtick',0:200:1000,'ytick',20:10:50,'tickdir','out','ticklength',[0.015 0.015]);
    end
%     set(gca,'xtick',0:200:1000,'ytick',20:10:50,'tickdir','out','ticklength',[0.015 0.015]);
end


%% plot distance
% full state space
figure; setFigPos(1,1); ha;
for i=1:length(dist) % plot individual condition/data set
    tmpT=(binSize/2)+[0:binSize:(binSize*(length(dist{i}(:))-1))]; % 10:20:....
    plot(tmpT,dist{i}(:),'k-');
end
% plot avg
D=fillNanCell(dist,0); % [#dataSet/conditions x time] [8 70]
tmpT=(binSize/2)+[0:binSize:(binSize*(size(D,2)-1))]; % 10:20:....
shadedErrorBar(tmpT(:)',nanmean(D,1),sem(D,1,'omitnan'),{'-','color','k','linewidth',2},1); drawnow;
plotVertical(gca,[T{1}(1) overlapT],[]);
set(gca,'xtick',[0 T{1}(1) overlapT overlapT+200],'tickdir','out','ticklength',[0.015 0.015]);
xlabel('Time from Ready (ms)');
ylabel('Distance between two priors (a.u.)');
% applytofig4paper;

figure; setFigPos(1,2); ha; % discrete
tmpT=cumsum(win)-(win/2); % 10:20:....
for i=1:size(dist2,1) % plot individual condition/data set
    plot(tmpT,dist2(i,:),'k-');
end
% plot avg
shadedErrorBar(tmpT(:)',nanmean(dist2,1),sem(dist2,1,'omitnan'),{'-','color','k','linewidth',2},1); drawnow;
plotVertical(gca,[T{1}(1) overlapT],[]);
set(gca,'xtick',cumsum(win),'tickdir','out','ticklength',[0.015 0.015]);
if ~idMean, set(gca,'ylim',[45 110]); else set(gca,'ylim',[1.9 4.8]); end
xlabel('Time from Ready (ms)');
ylabel('Distance between two priors (a.u.)');

% PC space
figure; setFigPos(2,1); ha;
for i=1:length(distPC) % plot individual condition/data set
    tmpT=(binSize/2)+[0:binSize:(binSize*(length(dist{i}(:))-1))]; % 10:20:....
    plot(tmpT,distPC{i}(:),'k-');
end
% plot avg
D=fillNanCell(distPC,0); % [#dataSet/conditions x time] [8 70]
tmpT=(binSize/2)+[0:binSize:(binSize*(size(D,2)-1))]; % 10:20:....
shadedErrorBar(tmpT(:)',nanmean(D,1),sem(D,1,'omitnan'),{'-','color','k','linewidth',2},1); drawnow;
plotVertical(gca,[T{1}(1) overlapT],[]);
set(gca,'xtick',[T{1}(1) overlapT T{2}(end)],'tickdir','out','ticklength',[0.015 0.015]);
xlabel('time from ready (ms)');
ylabel('distance between two priors');
% applytofig4paper;

figure; setFigPos(2,2); ha; % discrete
tmpT=cumsum(win)-(win/2); % 10:20:.
for i=1:size(distPC2,1) % plot individual condition/data set
    plot(tmpT,distPC2(i,:),'k-');
end
% plot avg
shadedErrorBar(tmpT(:)',nanmean(distPC2,1),sem(distPC2,1,'omitnan'),{'-','color','k','linewidth',2},1); drawnow;
plotVertical(gca,[T{1}(1) overlapT],[]);
set(gca,'xtick',cumsum(win),'tickdir','out','ticklength',[0.015 0.015]);
xlabel('Time from Ready (ms)');
ylabel('Distance between two priors (a.u.)');


function D=smoothBin(D,binSize,smthWidth)
handles.binWidth=binSize;
handles.use_sqrt=0;
    
handles.kern=smthWidth;

m = mean([D.data],2) * 1000;
mean_thresh=1; % remove neurons with FR<1 sp/s
keep_neurons = m>=mean_thresh; %-Inf; %m >= mean_thresh; % -Inf; % mean_thresh; %%%%%
handles.keep_neurons=keep_neurons;
handles.mean_thresh=mean_thresh;
% disp(nnz(keep_neurons));

% Remove low firing rate neurons
for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(handles.keep_neurons,:);
    D(itrial).id = D(itrial).id(handles.keep_neurons,:);
end
disp(['# neurons: ' num2str(size(D(1).data,1))]);

% Bin the spikes and use square root transform (use gpfa's pack)
if (handles.binWidth ~= 1 || handles.use_sqrt)
    [d(1:length(D)).spikes] = deal(D.data);
    for itrial = 1:length(D)
        d(itrial).trialId = itrial;
    end
    s = getSeq(d, handles.binWidth, 'useSqrt', handles.use_sqrt); % spike counts
    [D.data] = deal(s.y);
end

% Smooth data if necessary (automatically zero if GPFA selected)
for i = 1:length(D)
    D(i).data = smoother(D(i).data, handles.kern, handles.binWidth);
end

% convert into sp/s
for i = 1:length(D)
    D(i).data = D(i).data./handles.binWidth*1000;
end

% binning epochStarts
if length(D(i).epochStarts)>1
    for i = 1:length(D)
        D(i).epochStarts(2) = round(D(i).epochStarts(2)/handles.binWidth);
    end
end

%%
function D2=smoothBin2(D,win)
% no smoothing
% use window size specified by win=[240 240 320 240 240 240] last window
% doesn't have to be full
handles.use_sqrt=0;
    
m = mean([D.data],2) * 1000;
mean_thresh=1; % remove neurons with FR<1 sp/s
keep_neurons = m>=mean_thresh; %-Inf; %m >= mean_thresh; % -Inf; % mean_thresh; %%%%%
handles.keep_neurons=keep_neurons;
handles.mean_thresh=mean_thresh;
% disp(nnz(keep_neurons));

% Remove low firing rate neurons
for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(handles.keep_neurons,:);
    D(itrial).id = D(itrial).id(handles.keep_neurons,:);
end
disp(['# neurons: ' num2str(size(D(1).data,1))]);

iStart=cumsum(win)-win+1; % 1         241         481         801        1041        1281
iEnd=cumsum(win);          %   240         480         800        1040        1280        1520

D2=D;

for iD=1:length(D)
    D2(iD).data=[];
    D2(iD).epochStart=[1 find(iEnd==800)];
    for iBin=1:length(win)
        % deal with variable ending
        nT=size(D(iD).data,2);
        if iStart(iBin)<=nT % safe
            
            iEndTmp=min([iEnd(iBin) nT]);
            nSpike=sum(D(iD).data(:,iStart(iBin):iEndTmp),2);
            D2(iD).data=[D2(iD).data nSpike/length(iStart(iBin):iEndTmp)*1000]; % sp/s
            
        
        end % if iEnd(iBin)>nT
    end % for iD=1:length(D)
end % for iBin=1:length(win)