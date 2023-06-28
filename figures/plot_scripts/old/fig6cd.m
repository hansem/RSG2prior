function fig6cd

% 19/3/19
% combine across sessions & animals
% 

% 19/3/8
% fig4bc (mean Xu, varXu): H 161218
% sfig4bc (mean Xu, varXu): G 17050  op[
% re-run for G's GPFA 170506: corr b/t Xu & v didn't work

% 19/3/6
% putting behaviora data into GPFA data set & sSet (Xu), xro(u)

% 19/3/5
% - check variance across ts (zscoring for each conditino/prior across ts)
% - correlation with tp
% (TBD) d-prime > if works, example trajectory > crossVal_sess_anal_randPerturb_fig3e.m
% no IC/production
% data: H161218 G170506
% note: no stitching across sessions, all neurons need to be simulataneously recorded

% 18/9/27
% G170506
% idAllMeas: also just to check trajectory for whole measurement period

% 2018/9/26
% variance decrease by changing uts connecting only ts1-ts2, ts4-ts5? idUtsNew
% GPFA

% dependence: nData.m; estReadout.m

%% init
close all;
initRSG2prior; 
% other sessions: (criteria:At least, 40 cells with at least 2000 trials)
%         H: 161218 161222 (161211 >55 cells with >1700 trials)
%         G: 170817 170822 170818 170823 170511 170821 170507 170506 (more trials<>more neurons)
%               H/G: sess#12 16 (8) / 8 11 9 12 5 10 2 1 for varargin
% crossVal_sess(1,12); crossVal_sess(2,1); crossVal_sess(1,16); crossVal_sess(2,8); crossVal_sess(1,8);

iAnimal=2; % 2; % 1; % animalNm
if iAnimal==1
sessId=161218; % 161223; % 170506; % 161218 170506 161222 170817 170507 170818
    maxT0=7; % 18; % 5; % 7; % maximum of start trials to maximize # neuron&# trials
    minT1=1603; % 1839; % 1402; % 1603; % minimum of end trials
else
    sessId=170511; % 170823; % 170822; % 170821; % 170818; %170817; % 170511; % 170508; % 170507; % 170506; % 161218 170506 161222 170817 170507 170818
    maxT0=32; % 27; % 3; % 79; % 6; % 16; % 32; % 29; % 2; % 5; % 7; % maximum of start trials to maximize # neuron&# trials
    minT1=2828; % 2863; % 3586; % 2717; % 3573; % 3181; % 2828; % 1567; % 1732; % 1402; % 1603; % minimum of end trials
end

% session with # neurons>50
    % done
% 161206: 303        1556  % done
% 161211: 13        1477 %
    % 161218: 7 1603 %%%%%%
% 161221:17        1514 % 
    % 161222: 37 1810 
% 161223:18 1839

    % 170506: 5  1402 
    % 170507: 2 1732
    % 170508: 29 1567
    % 170511: 32 2828 %%%%%%
    % 170817: 16 3181
    % 170818: 6 3573
    % 170821:79 2717
    % 170822:3 3586 
% 170823: 27        2863 %

% load PSTH_170823_periSet_G_avgDir_SUMU.mat;
% disp([D(1).id t0 t1]); figure; setFigPos(1,1); hist(t0,50); figure; setFigPos(1,2); hist(t1,50);

idAllMeas=0; % 1; % if 0, only prior support

tBuffer=80; % to save trials with shortest ts

nH=20;
optsExpFig.Width=3;
optsExpFig.Height=3;
optsExpFig.FontSize='6';
optsExpFig.FontMode='fixed';
optsExpFig.Format='eps'; % 'tiff'; % 'pdf'; % 'png';
optsExpFig.LockAxes=1;
optsExpFig.LineMode='fixed'; % 'scaled';
%     optsExpFig.LineWidthMin=0.5;
%     optsExpFig.LineWidthMax=1.5;
optsExpFig.LineWidth=.5; %1;
optsExpFig.Renderer='painters';

lw=1.5;
lw2=0.25;
msize=4; % 5; % 10; % 4;
msize2=2;

% H/G, E/H, R/L
cmapHG=[0 0 0; .5 .5 .5]; % animal & E/H
%     cmapCond={[0 0 0; 1 1 1];
%         [.5 .5 .5; 1 1 1]};
markers={'o','^','o','^'};

nCmap=8; % nAnimal*nEH*nTarg;
cmapCond=parula(nCmap); % hsv(nCmap);
icmap=1;

cmapPr=[rgb('FireBrick'); rgb('RoyalBlue'); rgb('DarkGreen')]; %
cmapPr2=[rgb('IndianRed'); rgb('DodgerBlue'); rgb('ForestGreen')]; % for each data set

%% main

pPoly=nan(nPr,nEH,nTarg);

% only simulataneously recorded
D=PSTH2onlyFullRecord(['PSTH_' num2str(sessId) '_periSet_' animalNm{iAnimal} '_avgDir_SUMU.mat'],maxT0,minT1);

% task Ps
ts=[D.ts];
tp=[D.tp];
idPr=2-[D.idShortTrial]; % 1 for short, 2 for long
idEH=2-[D.idHandEye]; % 1 for eye, 2 for hand
theta=2-(0==[D.theta]); % 1 for right, 2 for left

% separate prior support & tp
if idAllMeas
    trajFn=['trajKS_' num2str(sessId) '_ts_' animalNm{iAnimal} '_condSpec_GPFA_CV_cellFullOnly.mat'];
else
    trajFn=['trajKS_' num2str(sessId) '_pr_' animalNm{iAnimal} '_condSpec_GPFA_CV_cellFullOnly.mat'];
end
disp(trajFn);
if exist(trajFn)
    D1=load(trajFn); % D proj_matrix keep_neurons eigenvalues binSize smthWidth thFR optimD
    
    % est uts for each prior (using mean states)
    u=nan(nPr,nTspp+1,D1.optimD); % second: difference vectors for all ts + uts
    Xu=cell(nEH,nTarg,nPr,nTspp,2); % last for uts0 vs uts
    dA=nan(nEH,nTarg,nPr,nTspp); % angle b/t uts- vs uts
    
    %     hF=figure; setFigPos(1,5); % % var(uts0 vs uts)
    
    % check variance across ts
    Vset=nan(nPr,nEH,nTarg,nTspp);
    nVset=nan(nPr,nEH,nTarg,nTspp);
    Mset=nan(nPr,nEH,nTarg,nTspp);
    
    h=figure; setFigPos(1,1);% Set
    h2=figure; setFigPos(1,2); % IC>mu
    h3=figure; setFigPos(1,3); % normalized variance
    h4=figure; setFigPos(1,4); % tp corr
    
    for iEH=1:nEH % condition-specific
        for iTarg=1:nTarg
            
            for iPr=1:nPr
                id=(idPr==iPr) & (idEH==iEH) & (theta==iTarg);
                
                % original uts
                u(iPr,end,:)=estReadout(D1.D(id));
                
                % ts-specific u
                Davg=avgD(D1.D(id)); % D(trials) > D(5 ts)
                [u(iPr,1,:),u(iPr,2,:),u(iPr,3,:),u(iPr,4,:),u(iPr,5,:)]=estReadout(Davg,1); % last for idUtsNew
                
                %                 figure; setFigPos(iPr,(iEH-1)*nTarg+iTarg); % histogram of Xu
                
                % projection
                for iTs=1:nTspp
                    tmpId=id&ts==T{iPr}(iTs);
                    Dtmp=cat(3,D1.D(tmpId).data); % [#LV x time x trials]
                    X=squeeze(Dtmp(:,end,:)); % neural states [#LV x trials]
                    
                    % original uts
                    uts0=squeeze(u(iPr,end,:)); % [#LV x 1]
                    Xu{iEH,iTarg,iPr,iTs,1}=X'*uts0; % [trials x 1]
                    
                    if ~isfield(D1,'sSet')
                        D1.sSet=nan(length(D1.D),1);
                        D1.xro=nan(length(D1.D),D1.optimD);
                        D1.tp=nan(length(D1.D),1);
                        D1.ts=nan(length(D1.D),1);
                        D1.idPr=nan(length(D1.D),1);
                        D1.idHandEye=nan(length(D1.D),1);
                        D1.theta=nan(length(D1.D),1);
                    end
                    D1.sSet(tmpId)=Xu{iEH,iTarg,iPr,iTs,1}; % [trials x 1]
                    D1.xro(tmpId,:)=repmat(uts0',nnz(tmpId),1); % [trials x #LV]
                    D1.tp(tmpId)=tp(tmpId);
                    D1.ts(tmpId)=ts(tmpId);
                    D1.idPr(tmpId)=iPr; % 1 for short, 2 for long
                    D1.idHandEye(tmpId)=iEH; % 1 for eye, 2 for hand
                    D1.theta(tmpId)=iTarg; % 1 for right, 2 for left
                    
                end % ts
                
                % z scoring for each cond (across ts)
                [tmpM,tmpSd]=meanSD(cell2mat(squeeze(Xu(iEH,iTarg,iPr,:,1))));
                for iTs=1:nTspp
                    % check variance across ts
                    Mset(iPr,iEH,iTarg,iTs)=nanmean((Xu{iEH,iTarg,iPr,iTs,1}-tmpM)./tmpSd);
                    Vset(iPr,iEH,iTarg,iTs)=nanvar((Xu{iEH,iTarg,iPr,iTs,1}-tmpM)./tmpSd);
                    
                    %                     % diference vector
                    %                     uts=squeeze(u(iPr,iTs,:)); % [#LV x 1]
                    %                     Xu{iEH,iTarg,iPr,iTs,2}=X'*uts; % [trials x 1]
                    %
                    %                     dA(iEH,iTarg,iPr,iTs)=angleVectors(uts0,uts); % angle
                    
                    %                     % plot
                    %                     subplot(nTspp,1,iTs); ha;
                    %                     histogram(X'*uts0,nH,'displaystyle','stairs','edgecolor',tmpCmap{iPr,1}(iTs,:)); ha;
                    %                     histogram(X'*uts,nH,'displaystyle','stairs','edgecolor',tmpCmap{iPr,1}(iTs,:),'LineStyle','-.'); drawnow;
                    %
                    %                     if iTs==nTspp
                    %                         xlabel('Xu'); ylabel('# trials');
                    %                     end
                    
                end % for iTs=1:nTspp
                
                % for normalized var
                tmpGrad=gradient(squeeze(Mset(iPr,iEH,iTarg,:))); % [1 x 5ts]
                for iTs=1:nTspp
                    nVset(iPr,iEH,iTarg,iTs)=Vset(iPr,iEH,iTarg,iTs)./tmpGrad(iTs);
                end
                figure(h);
                %                 plot(T{iPr}(:),squeeze(Mset(iPr,iEH,iTarg,:)),['-'],'color',pplot.cmap{(iEH-1)*nTarg+iTarg},'markerfacecolor','w');ha; drawnow; % pplot.marker8_2{im}  pplot.cmap{(i-1)*2+1}
                plot(T{iPr}(:),squeeze(Mset(iPr,iEH,iTarg,:)),['-'],'color',pplot.cmap{(iPr-.5)*2},'markerfacecolor','w','linewidth',lw2);ha; drawnow; % pplot.marker8_2{im}  pplot.cmap{(i-1)*2+1}
                figure(h2);
                %                 plot(T{iPr}(:),squeeze(Vset(iPr,iEH,iTarg,:)),['-'],'color',pplot.cmap{(iEH-1)*nTarg+iTarg},'markerfacecolor','w');ha; drawnow; % pplot.marker8_2{im}
                plot(T{iPr}(:),squeeze(Vset(iPr,iEH,iTarg,:)),['-'],'color',pplot.cmap{(iPr-.5)*2},'markerfacecolor','w','linewidth',lw2);ha; drawnow; % pplot.marker8_2{im}
                
                % polynomial fit
                [polyP,polyS]=polyfit(T{iPr}(:),squeeze(Vset(iPr,iEH,iTarg,:)),2);
                plot(T{iPr}(:),polyval(polyP,T{iPr}(:)),[':'],'color',pplot.cmap{(iEH-1)*nTarg+iTarg},'markerfacecolor','w','linewidth',lw2);ha; drawnow; % pplot.marker8_2{im}
                pPoly(iPr,iEH,iTarg)=polyP(1);
                
                figure(h3);
                plot(T{iPr}(:),squeeze(nVset(iPr,iEH,iTarg,:)),['-'],'color',pplot.cmap{(iEH-1)*nTarg+iTarg},'markerfacecolor','w','linewidth',lw2);ha; drawnow; % pplot.marker8_2{im}
                
                %                 applytofig(gcf,optsExpFig);
            end % for iPr=1:nPr
        end % . for iEH=1:nEH % condition-specific
    end %    for iTarg=1:nTarg
    
    % stat on quadratic coefficient
signrank(pPoly(:),0,'tail','left')
    %     signrank(pPoly(:))
    
    % save
%      v2struct(D1);
%      save(trajFn,'sSet','Xu','tp','ts','idPr','idHandEye','theta','-append');
    
    % plot grand averages
    Mset2=squeeze(mean(mean(Mset,2),3)); % nan(nPr,nTspp);
    Mset3=squeeze(sem(reshape(shiftdim(Mset,1),nEH*nTarg,nTspp,nPr),1));  % nEH,nTarg,nTspp,nPr > nEH*nTarg,nTspp,nPr > nTspp,nPr
    Vset2=squeeze(mean(mean(Vset,2),3)); % nan(nPr,nTspp);
    Vset3=squeeze(sem(reshape(shiftdim(Vset,1),nEH*nTarg,nTspp,nPr),1));  % nEH,nTarg,nTspp,nPr > nEH*nTarg,nTspp,nPr > nTspp,nPr
    nVset2=squeeze(mean(mean(nVset,2),3)); % nan(nPr,nTspp);
    nVset3=squeeze(sem(reshape(shiftdim(nVset,1),nEH*nTarg,nTspp,nPr),1));  % nEH,nTarg,nTspp,nPr > nEH*nTarg,nTspp,nPr > nTspp,nPr
    for i=1:nPr % short long
        %     for l=1:nTspp
        %         id=d.idPr==i & d.ts==T{i}(l);
        %
        %         Vset2(i,l)=var(sSet(id));
        %         Mset2(i,l)=mean(sSet(id));
        %     end
        figure(h);
        %         plot(T{i}(:),squeeze(Mset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}  pplot.cmap{(i-1)*2+1}
        shadedErrorBar(T{i}(:),squeeze(Mset2(i,:)),squeeze(Mset3(:,i)),{['-'],'color',pplot.cmap{(i-.5)*2},'markerfacecolor','w','linewidth',1,'markersize',msize},1);ha; drawnow; % pplot.marker8_2{im}
        plotVertical(gca,median(T{i}),[]);
        for iTs=1:nTspp
            plot(T{i}(iTs),squeeze(Mset2(i,iTs)),['o'],'color',tmpCmap{i,1}(iTs,:),'markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
        end
        figure(h2);
        %         plot(T{i}(:),squeeze(Vset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
        shadedErrorBar(T{i}(:),squeeze(Vset2(i,:)),squeeze(Vset3(:,i)),{['-'],'color',pplot.cmap{(i-.5)*2},'markerfacecolor','w','linewidth',1,'markersize',msize},1);ha; drawnow; % pplot.marker8_2{im}
        plotVertical(gca,median(T{i}),[]);
        for iTs=1:nTspp
            plot(T{i}(iTs),squeeze(Vset2(i,iTs)),['o'],'color',tmpCmap{i,1}(iTs,:),'markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
        end
        figure(h3);
        %     plot(T{i}(:),squeeze(nVset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
        shadedErrorBar(T{i}(:),squeeze(nVset2(i,:)),squeeze(nVset3(:,i)),{['-'],'color','k','markerfacecolor','w','linewidth',2},1);ha; drawnow; % pplot.marker8_2{im}
        
    end % prior
    figure(h); box off;
    set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out','ytick',-2:.5:2);
    xlabel('t_s (ms)'); ylabel('mean across trials'); axis tight;
    applytofig4paper; %applytofig4keynote;
    figure(h2); box off;
    set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out','ytick',-2:.5:2); % ,'ytick',-3:0.5:3
    xlabel('t_s (ms)'); ylabel('var. across trials');  axis tight; applytofig4paper; %applytofig4keynote;
    figure(h3);
    xlabel('t_s (ms)'); ylabel('norm. var. across trials');  axis tight; applytofig4keynote;
    
    %     optsExpFig.FontSize=12;
    %
    %     % var(uts0 vs uts): after combining across conditions
    %     dtXu=cellfun(@(x)(x-mean(x)),Xu,'UniformOutput',false);
    %     figure; setFigPos(1,5);
    %     for iPr=1:nPr
    %         for iTs=1:nTspp
    %             %             x=squeeze(cellfun(@var,Xu(:,:,iPr,iTs,1))); % EH x RL
    %             %             y=squeeze(cellfun(@var,Xu(:,:,iPr,iTs,2))); % EH x RL
    %             x=dtXu(:,:,iPr,iTs,1); x=cell2mat(x(:)); x=x(:);
    %             x=var(x); % EH x RL
    %             y=dtXu(:,:,iPr,iTs,2); y=cell2mat(y(:)); y=y(:);
    %             y=var(y); % EH x RL
    %
    %             plot(x(:),y(:),'o','color',tmpCmap{iPr,1}(iTs,:)); ha;
    %
    %         end % for iTs=1:nTspp
    %     end % pr
    %     plotIdentity(gca); drawnow;
    %     xlabel('var(u_t_s)'); ylabel('var(\Deltats)');
    %     applytofig(gcf,optsExpFig);
    %
    %     % angle uts0 vs uts
    %     figure; setFigPos(2,5);
    %     imagesc(reshape(dA,nEH*nTarg*nPr,nTspp)); colorbar;
    %     xlabel('ts'); ylabel('conditions (EH, RL, SL)');
    %     applytofig(gcf,optsExpFig);
    
else
    D1=D; % prior support
    for i=1:length(D1)
        D1(i).epochStarts=1;
        D1(i).epochColors=D1(i).epochColors(1,:);
        idPr=2-D1(i).idShortTrial; % 1 for short, 2 for long
        if idAllMeas
            D1(i).data=D1(i).data(:,1:D1(i).ts);
        else
            D1(i).data=D1(i).data(:,(T{idPr}(1)-tBuffer+1):D1(i).ts);
        end
    end % for i=1:length(D1)
    
    % delete extra fields in data to prevent error in cvreducedims
    % : only data, conditions, epochStarts, epochColors
    D12=rmfield(D1,{'id','tp','idShortTrial','ts','idHandEye','theta'});
    DataHigh(D12,'DimReduce');
    
    % production
    trajFn2=['trajKS_' num2str(sessId) '_tp_' animalNm{iAnimal} '_condSpec_GPFA_CV_cellFullOnly.mat'];
    disp(trajFn2);
    if exist(trajFn2)
        D2=load(trajFn2); % D proj_matrix keep_neurons eigenvalues binSize smthWidth thFR optimD
    else
        D2=D; % tp
        for i=1:length(D2)
            D2(i).epochStarts=1;
            D2(i).epochColors=D2(i).epochColors(1,:);
            D2(i).data=D2(i).data(:,D2(i).ts:end);
        end % for i=1:length(D2)
        
        % delete extra fields in data to prevent error in cvreducedims :
        % only data, conditions, epochStarts, epochColors
        D22=rmfield(D2,{'id','tp','idShortTrial','ts','idHandEye','theta'});
        DataHigh(D22,'DimReduce');
    end


end



function Davg=avgD(D)

% averaging across ts in D(trials)
nT=nData(D,2); % [trials x 1]
nTm=unique(nT); % [5ts x 1]

for i=1:length(nTm)
    id=nT==nTm(i);
    Davg(i).data=squeeze(mean(cat(3,D(id).data),3)); % [LV x time x trials] > [LV x time]
end % for i=1:nTm


function D2=PSTH2onlyFullRecord(fname,T0,T1)

% convert PSTH to include only neurons simultaneously recorded from t0 to t1
% fn: e.g. PSTH_161218_periSet_H_avgDir_SUMU.mat: D(trials) in DataHigh format, t0 (start trial for neurons), t1 (end trial)
% t0: maximum ID of start trial
% t1: minimum ID of end trial
%
% D=PSTH2onlyFullRecord('PSTH_161218_periSet_H_avgDir_SUMU.mat',7,1603);

%% 
initRSG2prior;
cd(singleTDir);

load(fname); % D t0 t1

id=t0<=T0 & t1>=T1; % [neuron x 1]

% exclude trials
D2=D(T0:T1);

% exclude neurons
for i=1:length(D2)
    D2(i).data=D2(i).data(id,:);
end

% %% 2019/3/6
% % putting behavioral parameters to GPFA data set
% d1=load('trajKS_161218_pr_H_condSpec_GPFA_CV_cellFullOnly.mat'); % 1597 trials
% d2=load('trajKS_161218_pr_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly_newUts_sqrtSingleT.mat'); % 1610 trials
% 
% iT1=1;
% for iT=1:length(d2.ts)
%     strcmp(d2.D(iT).condition,d1.D(iT1).condition)
%     
% end

