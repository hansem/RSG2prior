function analSetTransient

% 2019/3/12
% fig5d: distance during set transient (combined across data sets+shadedError)

% magnitude: distance between initial & end points (or mean length?)
% 2019/1/28

% measure angle between prior (with across-ts & within-prior as reference)
% 2019/1/24

% analyzing set transient using simple distance/nearest distance(kinet)
% trajKS_set2IC_XX_X.mat
% tIC=200
% 2019/1/22

%% init
initRSG2prior;
cd(psthDir);
load pplot.mat; % pplot.cmap for separate figure

idUseLowDim=1; % for angle

% plot
cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
cmapRot=[rgb('DeepPink'); rgb('Aqua')];
msize=4; % 10; % 4; % 6; % 2;
msize2=1;
lw=.75;% 1.2; % 1; % 3; % 1.5;
lw2=2;
cmapMat=[tmpCmap{1,1}; tmpCmap{2,1}];

%% main
d=dir('trajKS_set2IC_*.mat'); % prior*_supportOnly.mat'); % IC specific % _bin20_smth40
nDS=length(d);

load(d(1).name,'binSize','durIC');
nT=round(durIC/binSize);

% distance matrix
x=nan(nDS,nTspp+nTspp,nT); % ref: short800
y=nan(nDS,nPr,nTspp,nT); % ref: mean of each prior

m=nan(nDS,nPr,nTspp); % magnitude
a=nan(nDS,nTspp+nTspp); % angle bewteen two priors - ref: short800

for iDS=1:nDS
    
    disp(['===== ' d(iDS).name ' =====']);
    load(d(iDS).name); % binSize smthWidth optimD use_sqrt proj_matrix keep_neurons D eigenvalues meanPSTH
    disp(['optimD: ' num2str(optimD)]);
    
    dsName=d(iDS).name(15:18);
    iAnimalNm=d(iDS).name(18);
    if strcmp(iAnimalNm,animalNm{1}), iAnimal=1; else iAnimal=2; end
    if strcmp(dsName(1:2),'ER'), iCond=1; elseif strcmp(dsName(1:2),'EL'), iCond=2; elseif strcmp(dsName(1:2),'HR'), iCond=3; else iCond=4; end
    
    nPC=size(D(1).data,1);
    
    if idUseLowDim
        nPCangle=3;
    else
        nPCangle=nPC;
    end
    
    % ref: short800
    ref=D(nTspp).data(1:nPC,:); % [nPC x #time]
    tmp=cat(3,D.data); % [nPC x #time x #trajectories]
    x(iDS,:,:)=squeeze(sqrt(sum((tmp-ref).^2,1)))';  % [#time x #trajectories]'
    
    % ref: mean of each prior
    for iPr=1:nPr
        ref=D((iPr-1)*nTspp+round(nTspp/2)).data(1:nPC,:); % [nPC x #time]
        tmp=cat(3,D((iPr-1)*nTspp+[1:nTspp]).data); % [nPC x #time x #trajectories]
        y(iDS,iPr,:,:)=squeeze(sqrt(sum((tmp-ref).^2,1)))';  % [#time x #trajectories]'
    end
    
    % magnitude
    for iPr=1:nPr
        tmp=cat(3,D((iPr-1)*nTspp+[1:nTspp]).data); % [nPC x #time x #trajectories]
%         m(iDS,iPr,:)=sqrt(sum((tmp(:,end,:)-tmp(:,1,:)).^2,1)); % distance between initial & end points
        m(iDS,iPr,:)=sqrt(sum(mean(diff(tmp,1,2).^2,2),1)); % mean length?
    end
    
    % angle % ref: short800
    ref=D(nTspp).data(1:nPCangle,end)-D(nTspp).data(1:nPCangle,1); % [nPC x 1/#time]
    tmp=cat(3,D.data); % [nPC x #time x #trajectories]
    tmp=tmp(1:nPCangle,end,:)-tmp(1:nPCangle,1,:);
    for iTmp=1:size(tmp,3)
        a(iDS,iTmp)=estAngle(ref(:),squeeze(tmp(:,1,iTmp)));
    end
    
    % plot
    figure; setFigPos(iAnimal,iCond); ha; title([dsName(end) ' ' dsName(1:2)]);
    for iT=1:(2*nTspp)
        if iT>nTspp % flip sign for long
            plot(1:nT,-squeeze(x(iDS,iT,:)),'linewidth',2,'color',cmapMat(iT,:)); 
        else
            plot(1:nT,squeeze(x(iDS,iT,:)),'linewidth',2,'color',cmapMat(iT,:));
        end
    end
    xlabel('time after Set'); 
    ylabel('distance from short800');
    set(gca,'xtick',[0 nT/2 nT],'xticklabel',[0 durIC/2 durIC]);
    applytofig4keynote;
    
    figure; setFigPos(iAnimal,iCond); ha; title([dsName(end) ' ' dsName(1:2)]);
    for iPr=1:nPr
        for iT=1:nTspp
            if iT>nTspp/2 % flip sign for long
                plot(1:nT,-squeeze(y(iDS,iPr,iT,:)),'linewidth',2,'color',tmpCmap{iPr,1}(iT,:));
            else
                plot(1:nT,squeeze(y(iDS,iPr,iT,:)),'linewidth',2,'color',tmpCmap{iPr,1}(iT,:));
            end
        end
    end
    xlabel('time after Set'); 
    ylabel('distance from prior mean');
    set(gca,'xtick',[0 nT/2 nT],'xticklabel',[0 durIC/2 durIC]);
    applytofig4keynote;
    
end % for iDS=1:nDS

%% avg across data sets
mx=squeeze(mean(x,1)); % [ts x time]
my=squeeze(mean(y,1)); % [pr x ts x time]

sx=squeeze(sem(x,1)); % [ts x time]
sy=squeeze(sem(y,1)); % [pr x ts x time]

% plot
figure; setFigPos(1,1); ha; % distance from short800
for iT=1:(2*nTspp)
    if iT>nTspp % flip sign for long
        plot(1:nT,-squeeze(mx(iT,:)),'linewidth',2,'color',cmapMat(iT,:));
    else
        plot(1:nT,squeeze(mx(iT,:)),'linewidth',2,'color',cmapMat(iT,:));
    end
end
xlabel('time after Set');
ylabel('distance from short800');
set(gca,'xtick',[0 nT/2 nT],'xticklabel',[0 durIC/2 durIC]);
applytofig4keynote;

% fig5d %%%%%
figure; setFigPos(1,1);ha;box off;  % distance from prior mean
tmpX=binSize/2+binSize*([1:nT]-1);
for iPr=1:nPr
    for iT=1:nTspp
        if iT>nTspp/2 % flip sign for long
%             plot(1:nT,-squeeze(my(iPr,iT,:)),'linewidth',2,'color',tmpCmap{iPr,1}(iT,:));
            shadedErrorBar(tmpX,-squeeze(my(iPr,iT,:)),squeeze(sy(iPr,iT,:)),{'linewidth',lw,'color',tmpCmap{iPr,1}(iT,:)},1);
        else
%             plot(1:nT,squeeze(my(iPr,iT,:)),'linewidth',2,'color',tmpCmap{iPr,1}(iT,:));
            shadedErrorBar(tmpX,squeeze(my(iPr,iT,:)),squeeze(sy(iPr,iT,:)),{'linewidth',lw,'color',tmpCmap{iPr,1}(iT,:)},1);
        end
    end
end
xlabel('Time after Set');
ylabel('Distance from prior mean');
set(gca,'xtick',[0 durIC/2 durIC],'xticklabel',[0 durIC/2 durIC],'tickDir','out','tickLength',[0.015 0.015],'ytick',[-1:1:1]);
xlim([0 durIC]);
% applytofig4keynote;

return;
% close all;

% magnitude
figure; setFigPos(2,1); ha;
for iDS=1:nDS
    plot(T{1}(:),squeeze(m(iDS,1,1:nTspp)),'o-','color','r'); 
    plot(T{2}(:),squeeze(m(iDS,2,1:nTspp)),'o-','color','b'); 
end
% plot mean
shadedErrorBar(T{1}(:),squeeze(nanmean(m(:,1,1:nTspp),1)),squeeze(sem(m(:,1,1:nTspp),1,'omitnan')),{'o-','color','r'},1);
shadedErrorBar(T{2}(:),squeeze(nanmean(m(:,2,1:nTspp),1)),squeeze(sem(m(:,2,1:nTspp),1,'omitnan')),{'o-','color','b'},1);
xlabel('t_s'); ylabel('magnitude');
applytofig4keynote;

% angle: scatterplot across prior vs within prior
figure; setFigPos(2,2); ha; % pattern across TS
for iDS=1:nDS
    plot(T{1}(:),a(iDS,1:nTspp),'o-','color','r'); 
    plot(T{2}(:),a(iDS,nTspp+[1:nTspp]),'o-','color','b'); 
end
xlabel('t_s'); ylabel('angle from short800');
applytofig4keynote;

figure; setFigPos(2,3); ha; % scatter after avg across ts withi prior
for iDS=1:nDS
    mX=nanmean(a(iDS,1:(nTspp-1)));
    sX=sem(a(iDS,1:(nTspp-1)),2,'omitnan');
    
    mY=nanmean(a(iDS,nTspp+[1:nTspp]));
    sY=sem(a(iDS,nTspp+[1:nTspp]),2,'omitnan');
    
    plotErrorbar(mX,sX,mY,sY,'k');
end
axisEqual;
plotIdentity(gca);
xlabel('within-prior angle');
ylabel('across-prior angle');
applytofig4keynote;

