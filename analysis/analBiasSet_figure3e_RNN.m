function analBiasSet_figure3e_RNN

% 2018/5/16
% RNN: ls traj_prior_Weber*.mat (traj_prior_Weber_2Prior_Type2_005_v5)
% change into form of tp-ts
%
% 2018/5/9
% add bootstrap data as error bar
% figure 3e
% trim a lot of codes from analSet_randProj.m

% 2018/4/12
% generate BLS from projection to readout vector
% from analSet_randProj

%% initial
idTpMTs=1;  % plot te-ts vs ts
idPlot=0;
idDebug=0; %1;

dim=inf; % 3; % 5; % 4; % 2; % 3; %%%%%
useOptimD=1; % 0; % if 1, use optimD
maxDim=11; % 8; %  trajKS_prior_(CondNm)_(H/G)_supportOnly.mat
% H/G, ER/EL/HR/HL: 8 7 10 9 7 7 9 6 for priorIncSL, 8 7 9 8 10 9 11 7

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
cmap3=[tmpCmap{1,1};tmpCmap{2,1}];

tmpMarker={'^','v','<','>','o','s','d','h'};

% figure's PaperUnits: 'inches'
optsExpFig.Height=3/2.54; %2.1/2.54; % 2*1.2; % 0.8; % '1.2'; % '2'; % 7;
optsExpFig.Width=3/2.54; %3.15/2.54; % 3.6; % 1.5*1.2*2; % '2.4'; % '4';
Width2=2.34/2.54; % 1.7*2*1.2; %
Width3=2.1/2.54; % 1.3*2*1.2;
optsExpFig.FontSize='6';
optsExpFig.FontMode='fixed';
optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
optsExpFig.LockAxes=1;
optsExpFig.LineMode='scaled';
optsExpFig.LineWidthMin=0.5;
optsExpFig.LineWidthMin=1.2;
optsExpFig.Renderer='painters';

lw=2; % 1;
lw2=1;
msize=4;
msize2= 2;

idExpFig=0; % 0; % 1;

cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
cmapPr2=[rgb('IndianRed'); rgb('DodgerBlue'); rgb('ForestGreen')]; % for each data set


%% initial
initRSG2prior;nTs=nTspp;
try
    cd(psthDir);
catch
    psthDir='/Users/seonminahn/Dropbox (MIT)/psthDataHigh';
    cd(psthDir);
end

binSize=20;
smthWidth=40;

% RNN specific
d=dir('traj_prior_Weber*.mat');
nRNN=length(d);

for idAnimal=1:nRNN
    iAnimalNm=['Weber' d(idAnimal).name((end-9):(end-4))]; % Weber000_v3
    disp(iAnimalNm);
    
    hFig2=figure; ha; setFigPos(1,idAnimal); % % te by projection to readout vector vs BLS =f(ts)
    teRv=cell(nPr,1); % 1/speed fit to BLS to avg
    
    tsNm=d(idAnimal).name; % ['trajKS_prior_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_supportOnly.mat']; % '_bin20_smth40.mat']; % ['trajKS_ts_' ehNm{iEH}(1) targNm{iTarg}(1) '_' iAnimalNm '_bin20_smth40.mat'];
    [teRv,teBLS]=estTeRv(teRv,tsNm,nPr,T);
    
    % output teRv{nPr}[ 1 x timePoints]
    
    %% plot
    
    for iPr=1:nPr
        % BLS
        tmpT=min(T{iPr}):max(T{iPr});
        if ~idTpMTs
            plot(tmpT,teBLS{iPr},'-','linewidth',lw2,'color',tmpCmap{iPr,1}(1,:));
            if iPr==nPr, plotIdentity(gca); end; %  axis tight;
        else
            plot(tmpT,teBLS{iPr}-tmpT,'-','linewidth',lw2,'color',tmpCmap{iPr,1}(1,:));
            if iPr==nPr, plotHorizon(gca,0,[]); end; %  axis tight;
        end
        % data: now no std
        mTeRv=squeeze(mean(teRv{iPr},1));
        sTeRv=squeeze(std(teRv{iPr},0,1)); % sTeRv=squeeze(sem(teRv{iPr},1));
        tmpT=linspace(min(T{iPr}),max(T{iPr}),length(mTeRv));
        if ~idTpMTs
            plot(tmpT(:),mTeRv(:),'o','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize,'markerfacecolor','w'); %
            c=corr(mTeRv(:),BLS(tmpT(:),0.05,[min(tmpT) max(tmpT)],'uniform'));
            disp(['R^2(BLS): ' num2str(c.^2,2)]);
%             h=shadedErrorBar(tmpT,mTeRv,sTeRv,{'.','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize},1); drawnow; % 2
        else
            plot(tmpT(:),mTeRv(:)-tmpT(:),'o','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize,'markerfacecolor','w'); %
            c=corr(mTeRv(:)-tmpT(:),BLS(tmpT(:),0.05,[min(tmpT) max(tmpT)],'uniform')-tmpT(:));
            disp(['R^2(BLS-ts): ' num2str(c.^2,2)]);
%             h=shadedErrorBar(tmpT,mTeRv-tmpT,sTeRv,{'.','color',tmpCmap{iPr,1}(1,:),'linewidth',lw2,'markersize',msize},1); drawnow; % 2
        end
    end
    if idTpMTs
        set(gca,'ticklength',[0.02 0.02],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000}); % ,...
        %                         'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
    else
        set(gca,'ticklength',[0.02 0.02],'tickdir','out','xtick',[T{1}(1:2:end) T{2}(3:2:end)],'xticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000},...
            'ytick',[T{1}(1:2:end) T{2}(3:2:end)],'yticklabel',{T{1}(1)/1000;[];T{1}(end)/1000;[];T{2}(end)/1000});
    end
    xlabel('t_s (s)'); ylabel('transformed states projected to readout (s)');
    axis tight; title(iAnimalNm);
%     legend(iAnimalNm,'location','best'); legend boxoff;
    
end % animal



function [teRv,teBLS]=estTeRv(teRv,tsNm,nPr,T)

% output teRv{nPr}[ 1 x timePoints]

[diffV,sSet,sSetAll]=extStateSet(tsNm); % plotTraj_figure3(tsNm); ha; % diffV:[dim x 9ts] sSet:[optimD x 10ts] for curvature; sSetAll:[optimD x 17(short)/21(long)]

% 2) check R2 for actual data is high
% more constraints by middle points between sampled ts?
isUseMid=1; %%%%%
% BLS
wm=0.05; %%%%% 0.048 for G, 0.05 for H
idUseSp=0; % 1; %%%%% use speed, not BLS
if idUseSp
    if isUseMid
        tmpT=linspace(min(T{1}),max(T{1}),length(sSetAll{1}));tmpTs{1}=tmpT;
        te{1}=1./BLS(tmpT,wm,[min(T{1}) max(T{1})],'uniform');
        tmpT=linspace(min(T{2}),max(T{2}),length(sSetAll{2}));tmpTs{2}=tmpT;
        te{2}=1./BLS(tmpT,wm,[min(T{2}) max(T{2})],'uniform');
    else
        te{1}=1./BLS(T{1},wm,[min(T{1}) max(T{1})],'uniform');tmpTs{1}=T{1};
        te{2}=1./BLS(T{2},wm,[min(T{2}) max(T{2})],'uniform');tmpTs{2}=T{2};
    end
else
    if isUseMid
        tmpT=linspace(min(T{1}),max(T{1}),length(sSetAll{1}));tmpTs{1}=tmpT;
        te{1}=BLS(tmpT,wm,[min(T{1}) max(T{1})],'uniform');
        tmpT=linspace(min(T{2}),max(T{2}),length(sSetAll{2}));tmpTs{2}=tmpT;
        te{2}=BLS(tmpT,wm,[min(T{2}) max(T{2})],'uniform');
    else
        te{1}=BLS(T{1},wm,[min(T{1}) max(T{1})],'uniform');tmpTs{1}=T{1};
        te{2}=BLS(T{2},wm,[min(T{2}) max(T{2})],'uniform');tmpTs{2}=T{2};
    end
end
% - BLS (1ms resolution): wm defined above
teBLS{1}=BLS(min(T{1}):max(T{1}),wm,[min(T{1}) max(T{1})],'uniform');
teBLS{2}=BLS(min(T{2}):max(T{2}),wm,[min(T{2}) max(T{2})],'uniform');


% only need sSet, sSetAll from above (plotTraj_figure3)
% use midPoints (isUseMid) & not speed (idUseSp)
sPrjRV=cell(nPr,1); sPrjRVFit=cell(nPr,1);
% - define readout(480to800) vector
% sSet:[optimD x 10ts] for curvature; sSetAll:[optimD x 17(short)/21(long)]
vR{1}=mean([sSet(:,5)-sSet(:,1) sSet(:,4)-sSet(:,2)],2); % [optimD x 1]
vR{1}=normalize(vR{1},1); % [optimD x 1]
vR{2}=mean([sSet(:,10)-sSet(:,6) sSet(:,9)-sSet(:,7)],2); % [optimD x 1]
vR{2}=normalize(vR{2},1); % [optimD x 1]
% - projection
for iPr=1:nPr
    sPrjRV{iPr}=vR{iPr}'*sSetAll{iPr}; % [1 x optimD] x sSetAll:[optimD x 17(short)/21(long)]
    stats=regstats(te{iPr}(:),sPrjRV{iPr}(:),'linear',{'beta','rsquare','yhat'}); % rsquare
    sPrjRVFit{iPr}=stats.yhat(:);
%     disp(stats.rsquare);
    teRv{iPr}=cat(1,teRv{iPr},sPrjRVFit{iPr}(:)'); % [8dataSet x timePoints] to Avg
end

