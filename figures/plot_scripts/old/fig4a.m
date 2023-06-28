function fig4a

% 19/3/10
% plot example trajectory
% crossVal_sess_plotExTraj_fig3e.m was only for PCA
% Trial175/Trial1610 was best > 175-(maxT0-1)=169

% 19/3/8
% fig4bc (mean Xu, varXu): H 161218
% sfig4bc (mean Xu, varXu): G 170506

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
initRSG2prior; 
% other sessions: (criteria:At least, 40 cells with at least 2000 trials)
%         H: 161218 161222 (161211 >55 cells with >1700 trials)
%         G: 170817 170822 170818 170823 170511 170821 170507 170506 (more trials<>more neurons)
%               H/G: sess#12 16 (8) / 8 11 9 12 5 10 2 1 for varargin
% crossVal_sess(1,12); crossVal_sess(2,1); crossVal_sess(1,16); crossVal_sess(2,8); crossVal_sess(1,8);

iAnimal=1; % 1; % animalNm
if iAnimal==1
sessId=161218; % 170506; % 161218 170506 161222 170817 170507 170818
    maxT0=7; % 5; % 7; % maximum of start trials to maximize # neuron&# trials
    minT1=1603; % 1402; % 1603; % minimum of end trials
else
    sessId=170506; % 170506; % 161218 170506 161222 170817 170507 170818
    maxT0=5; % 5; % 7; % maximum of start trials to maximize # neuron&# trials
    minT1=1402; % 1402; % 1603; % minimum of end trials
end
% 161218: 7 1603
% 170506: 5  1402
% 161222:
% 170817:

idTrial=169; % 152; % 175; % short hand - best for IC

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
% if exist(trajFn)
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
    
    % 13 15 16 
%     for idTrial=7:length(D1.D)
    idTrial=17;

    % particualar trial of interest: idTrial
    iPr=idPr(idTrial); % 1 for short, 2 for long
    tmpIdPr=iPr;
    iEH=idEH(idTrial); % 1 for eye, 2 for hand
    iTarg=theta(idTrial); % 1 for right, 2 for left
    
    disp([prNm{iPr} ' ' ehNm{iPr} ' ' targNm{iTarg} ' trial#' num2str(idTrial)]);
    
%     for iEH=1:nEH % condition-specific
%         for iTarg=1:nTarg
%             
%             for iPr=1:nPr
                id=(idPr==iPr) & (idEH==iEH) & (theta==iTarg);
                
                % original uts
                u(iPr,end,:)=estReadout(D1.D(id));
                uts0=squeeze(u(iPr,end,:)); % [#LV x 1]
                
                % ts-specific u
                Davg=avgD(D1.D(id)); % D(trials) > D(5 ts).data [#LV x time]
                [u(iPr,1,:),u(iPr,2,:),u(iPr,3,:),u(iPr,4,:),u(iPr,5,:)]=estReadout(Davg,1); % last for idUtsNew
                
                %                 figure; setFigPos(iPr,(iEH-1)*nTarg+iTarg); % histogram of Xu
                
                % projection
%                 for iTs=1:nTspp
%                     tmpId=id&ts==T{iPr}(iTs);
                    Dtmp=cat(3,D1.D(idTrial).data); % [#LV x time x trials]
                    X=squeeze(Dtmp(:,end,:)); % neural states [#LV x trials]
                    
                    % original uts
                    Xu{iEH,iTarg,iPr,iTs,1}=X'*uts0; % [trials x 1]
                    
%                 end % ts
                
                %% actual plot
                % rotating trajectory+readout vector in state space
                
                % other prior
                if iPr==1
                    Davg2=avgD(D1.D((idPr==(3-iPr)) & (idEH==iEH) & (theta==iTarg))); % D(trials) > D(5 ts)
                    M1=struct2mat(Davg,'data'); % short
                    M2=struct2mat(Davg2,'data'); % long
                    
                    Dtmp=cat(3,cat(2,M1,nan(size(M1,1),size(M2,2)-size(M1,2),size(M1,3))),M2); % [9dim x timePoints x 10ts]
                else
                    Davg2=avgD(D1.D((idPr==(3-iPr)) & (idEH==iEH) & (theta==iTarg))); % D(trials) > D(5 ts)
                    M1=struct2mat(Davg,'data'); % long
                    M2=struct2mat(Davg2,'data'); % short
                    
                    Dtmp=cat(3,cat(2,M2,nan(size(M2,1),size(M1,2)-size(M2,2),size(M2,3))),M1); % [9dim x timePoints x 10ts]
                end
                
                % making 3D
                Dtmp=Dtmp(1:3,:,:); % [3dim x 24timePoints x 10ts]

%                 if idCheckTraj

                    Dtmp=permute(Dtmp,[3 2 1]); % [10cond x 25time x 3]
                    figure; setFigPos(1,1); ha; grid on;
                    
%                     if k==1 % avg-ts traj
                        DtmpMean=cat(1,nanmean(Dtmp(1:nTspp,:,:),1),nanmean(Dtmp(nTspp+[1:nTspp],:,:),1)); % [2pr x 25time x 3]
                        
%                         if idSmoothAfterAvgAtt
                            binWidth=1; % 20 ms bin in PC
                            kernSD=2; % 40ms smoothing was used in PC
                            for iDtmpMean=1:size(DtmpMean,1) % pr
                                tmpTmpTmp=squeeze(DtmpMean(iDtmpMean,:,:))'; % [3pc 24time]
                                isNanTmp=~isnan(tmpTmpTmp(1,:));
                                DtmpMean(iDtmpMean,isNanTmp,:)=smoother(tmpTmpTmp(:,isNanTmp),kernSD,binWidth)';
                            end
%                         end
                                
                        plot3Dmultiple(DtmpMean,[tmpCmap{1,1}(1,:);tmpCmap{2,1}(1,:)],[],lw); % avg-ts traj
%                         plot3(DtmpMean(1,1,1),DtmpMean(1,1,2),DtmpMean(1,1,3),'^','color',tmpCmap{1,1}(1,:)); % start
%                         plot3(DtmpMean(2,1,1),DtmpMean(2,1,2),DtmpMean(2,1,3),'^','color',tmpCmap{2,1}(1,:));
%                     end
                    
                    % states for each ts
                    cmapAll=cell2mat(tmpCmap(:,1)); % [10ts x 3]
                    for iTsTmp=1:size(Dtmp,1)
                        idPrTmp=floor((iTsTmp-1)./nTspp)+1; % 1111122222
%                         if k==1 % rotation
                            tTmp=nnz(~isnan(Dtmp(iTsTmp,:,1)));
                            plot3(DtmpMean(idPrTmp,tTmp,1),DtmpMean(idPrTmp,tTmp,2),DtmpMean(idPrTmp,tTmp,3),...
                                'o','markerfacecolor',cmapAll(iTsTmp,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2); % states for each ts
%                         else
%                             tTmp=iIC;
% %                             plotTraj_figure3('trajKS_tp_EL_H_bin20_smth40.mat',tmpDtraj);
%                             if rem(iTsTmp,nTspp)~=0 % connecting lines
%                                 plot3(Dtmp(iTsTmp+[0:1],tTmp,1),Dtmp(iTsTmp+[0:1],tTmp,2),Dtmp(iTsTmp+[0:1],tTmp,3),...
%                                     '-','color',tmpCmap{idPrTmp}(1,:),'linewidth',lw); % states for each ts
%                             end
%                             plot3(Dtmp(iTsTmp,tTmp,1),Dtmp(iTsTmp,tTmp,2),Dtmp(iTsTmp,tTmp,3),...
%                                 'o','markerfacecolor',cmapAll(iTsTmp,:),'markeredgecolor','k','markersize',msize,'linewidth',lw2); % states for each ts
%                         end
                    end
                    
                    % readout
%                     if k==1
                        tmpX=nanmean(DtmpMean(tmpIdPr,[1 nnz(~isnan(DtmpMean(tmpIdPr,:,1)))],1)); % plot readout around center of rotation
                        tmpY=nanmean(DtmpMean(tmpIdPr,[1 nnz(~isnan(DtmpMean(tmpIdPr,:,2)))],2));
                        tmpZ=nanmean(DtmpMean(tmpIdPr,[1 nnz(~isnan(DtmpMean(tmpIdPr,:,3)))],3));
%                     else
%                         tmpX=nanmean(Dtmp([1:nTspp]+(tmpIdPr-1)*nTspp,iIC,1)); % plot readout around center of rotation
%                         tmpY=nanmean(Dtmp([1:nTspp]+(tmpIdPr-1)*nTspp,iIC,2));
%                         tmpZ=nanmean(Dtmp([1:nTspp]+(tmpIdPr-1)*nTspp,iIC,3));
%                     end
                    plot3(tmpX+0.1*[-1;1]*D1.xro(idTrial,1),...
                        tmpY+0.1*[-1;1]*D1.xro(idTrial,2),...
                        tmpZ+0.1*[-1;1]*D1.xro(idTrial,3),'--','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % readout
%                     if idIncEarly
%                         plot3(DtmpMean(1,24,1),DtmpMean(1,24,2),DtmpMean(1,24,3),'x','color',tmpCmap{1,1}(1,:)); % 480
%                         plot3(DtmpMean(2,40,1),DtmpMean(2,40,2),DtmpMean(2,40,3),'x','color',tmpCmap{2,1}(1,:));
%                     end
%                     plot3Dmultiple(Dtmp,cell2mat(tmpCmap(:,1)),[]); grid on; % individual ts traj
                    
                    % this trial's state        
                    cmap3=tmpCmap{tmpIdPr,1}(D1.ts(idTrial)==T{tmpIdPr},:);
%                     plot3(Dtraj.data(1,:),Dtraj.data(2,:),Dtraj.data(3,:),'-','linewidth',2,'color',cmap3); % leaved trial's trajectory
%                     if k==1
                        plot3(D1.D(idTrial).data(1,end),D1.D(idTrial).data(2,end),D1.D(idTrial).data(3,end),'o','linewidth',lw2,'markeredgecolor',cmap3,'markerfacecolor',cmap3,'markersize',msize); % start leaved trial's trajectory
%                     else
%                         plot3(D1.D(idTrial).data(1,iIC),D1.D(idTrial).data(2,iIC),D1.D(idTrial).data(3,iIC),'o','linewidth',lw2,'markeredgecolor',cmap3,'markerfacecolor',cmap3,'markersize',msize); % start leaved trial's trajectory
%                     end
%                     
                    % trial state to readout
%                     % xb=x(x'x)^(-1)x'y
%                     x=D(iTrial).xro(1:3);
%                     if k==1,y=Dtraj.data(1:3,end); else y=Dtraj.data(1:3,iIC); end
%                     xb=x*((x'*x)^(-1))*x'*y;
                    % (a-p)-((a-p)'*n)*n
                    a=[tmpX;tmpY;tmpZ]; n=D1.xro(idTrial,1:3)';
%                     if k==1,
                        p=D1.D(idTrial).data(1:3,end); 
%                     else p=D1.D(idTrial).data(1:3,iIC); end
                    xb=(a-p)-((a-p)'*n)*n;
%                     if k==1
                        %                         plot3([Dtraj.data(1,end);xb(1)],...
                        %                             [Dtraj.data(2,end);xb(2)],...
                        %                             [Dtraj.data(3,end);xb(3)],':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
                        plot3(D1.D(idTrial).data(1,end)+[0;1]*xb(1),...
                            D1.D(idTrial).data(2,end)+[0;1]*xb(2),...
                            D1.D(idTrial).data(3,end)+[0;1]*xb(3),':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
%                     else
%                         %                         plot3([Dtraj.data(1,iIC);xb(1)],...
%                         %                             [Dtraj.data(2,iIC);xb(2)],...
%                         %                             [Dtraj.data(3,iIC);xb(3)],':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
%                         plot3(Dtraj.data(1,iIC)+[0;1]*xb(1),...
%                             Dtraj.data(2,iIC)+[0;1]*xb(2),...
%                             Dtraj.data(3,iIC)+[0;1]*xb(3),':','color',tmpCmap{tmpIdPr,1}(1,:),'linewidth',lw2);                    % projection
%                     end
                    
                    
                    axis tight;
                    xlabel('Latent Factor 1'); ylabel('Latent Factor 2'); zlabel('Latent Factor 3');
                    set(gca,'view',[-29.6000  -18.0000],... % -31.2000  -10.8000],... % 40 74],... % 66 48],...
                        'xtick',0:0.2:0.2,'ytick',-0.2:0.1:(-0.1),'ztick',-0.1:0.1:0.1);
                    
                    
%                     if k==1


% % %                         set(gca,'view',[-21 83],... % 40 74],... % 66 48],...
% % %                             'xtick',-1:1,'ytick',-1:1,'ztick',-3:1:3);
% % %                         xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
                        
                        
                        
                        
                        
%                     else
% %                         axis tight;
%                         set(gca,'view',[-18.8000  -56.4000],...  % [147.6000   71.6000] % [-108 34],...
%                             'xtick',-3:1:3,'ytick',-2:1:3,'ztick',-2:1:1);
%                         xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
%                     end
%                     waitforbuttonpress; close;
%                 end % idCheckTraj
                
                %%
                
%                 tW=0;
%                 while 1-tW
%                     tW=waitforbuttonpress; 
%                 end
%                 close all;
                
%     end % for idTrial=1:length(D1.D)
                

                                applytofig(gcf,optsExpFig);
%             end % for iPr=1:nPr
%         end % . for iEH=1:nEH % condition-specific
%     end %    for iTarg=1:nTarg
    
    % save
%      v2struct(D1);
%      save(trajFn,'sSet','Xu','tp','ts','idPr','idHandEye','theta','-append');
    
%     % plot grand averages
%     Mset2=squeeze(mean(mean(Mset,2),3)); % nan(nPr,nTspp);
%     Mset3=squeeze(sem(reshape(shiftdim(Mset,1),nEH*nTarg,nTspp,nPr),1));  % nEH,nTarg,nTspp,nPr > nEH*nTarg,nTspp,nPr > nTspp,nPr
%     Vset2=squeeze(mean(mean(Vset,2),3)); % nan(nPr,nTspp);
%     Vset3=squeeze(sem(reshape(shiftdim(Vset,1),nEH*nTarg,nTspp,nPr),1));  % nEH,nTarg,nTspp,nPr > nEH*nTarg,nTspp,nPr > nTspp,nPr
%     nVset2=squeeze(mean(mean(nVset,2),3)); % nan(nPr,nTspp);
%     nVset3=squeeze(sem(reshape(shiftdim(nVset,1),nEH*nTarg,nTspp,nPr),1));  % nEH,nTarg,nTspp,nPr > nEH*nTarg,nTspp,nPr > nTspp,nPr
%     for i=1:nPr % short long
%         %     for l=1:nTspp
%         %         id=d.idPr==i & d.ts==T{i}(l);
%         %
%         %         Vset2(i,l)=var(sSet(id));
%         %         Mset2(i,l)=mean(sSet(id));
%         %     end
%         figure(h);
%         %         plot(T{i}(:),squeeze(Mset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}  pplot.cmap{(i-1)*2+1}
%         shadedErrorBar(T{i}(:),squeeze(Mset2(i,:)),squeeze(Mset3(:,i)),{['-'],'color',pplot.cmap{(i-.5)*2},'markerfacecolor','w','linewidth',1,'markersize',msize},1);ha; drawnow; % pplot.marker8_2{im}
%         plotVertical(gca,median(T{i}),[]);
%         for iTs=1:nTspp
%             plot(T{i}(iTs),squeeze(Mset2(i,iTs)),['o'],'color',tmpCmap{i,1}(iTs,:),'markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
%         end
%         figure(h2);
%         %         plot(T{i}(:),squeeze(Vset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
%         shadedErrorBar(T{i}(:),squeeze(Vset2(i,:)),squeeze(Vset3(:,i)),{['-'],'color',pplot.cmap{(i-.5)*2},'markerfacecolor','w','linewidth',1,'markersize',msize},1);ha; drawnow; % pplot.marker8_2{im}
%         plotVertical(gca,median(T{i}),[]);
%         for iTs=1:nTspp
%             plot(T{i}(iTs),squeeze(Vset2(i,iTs)),['o'],'color',tmpCmap{i,1}(iTs,:),'markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
%         end
%         figure(h3);
%         %     plot(T{i}(:),squeeze(nVset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
%         shadedErrorBar(T{i}(:),squeeze(nVset2(i,:)),squeeze(nVset3(:,i)),{['-'],'color','k','markerfacecolor','w','linewidth',2},1);ha; drawnow; % pplot.marker8_2{im}
%         
%     end % prior
%     figure(h); box off;
%     set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out','ytick',-2:.5:2);
%     xlabel('t_s (ms)'); ylabel('mean across trials'); axis tight;
%     applytofig4paper; %applytofig4keynote;
%     figure(h2); box off;
%     set(gca,'xtick',unique([T{1}(1:2:end) T{2}(3:end)]),'ticklength',[0.01 0.01],'tickdir','out','ytick',-2:.5:2); % ,'ytick',-3:0.5:3
%     xlabel('t_s (ms)'); ylabel('var. across trials');  axis tight; applytofig4paper; %applytofig4keynote;
%     figure(h3);
%     xlabel('t_s (ms)'); ylabel('norm. var. across trials');  axis tight; applytofig4keynote;
%     
%     %     optsExpFig.FontSize=12;
%     %
%     %     % var(uts0 vs uts): after combining across conditions
%     %     dtXu=cellfun(@(x)(x-mean(x)),Xu,'UniformOutput',false);
%     %     figure; setFigPos(1,5);
%     %     for iPr=1:nPr
%     %         for iTs=1:nTspp
%     %             %             x=squeeze(cellfun(@var,Xu(:,:,iPr,iTs,1))); % EH x RL
%     %             %             y=squeeze(cellfun(@var,Xu(:,:,iPr,iTs,2))); % EH x RL
%     %             x=dtXu(:,:,iPr,iTs,1); x=cell2mat(x(:)); x=x(:);
%     %             x=var(x); % EH x RL
%     %             y=dtXu(:,:,iPr,iTs,2); y=cell2mat(y(:)); y=y(:);
%     %             y=var(y); % EH x RL
%     %
%     %             plot(x(:),y(:),'o','color',tmpCmap{iPr,1}(iTs,:)); ha;
%     %
%     %         end % for iTs=1:nTspp
%     %     end % pr
%     %     plotIdentity(gca); drawnow;
%     %     xlabel('var(u_t_s)'); ylabel('var(\Deltats)');
%     %     applytofig(gcf,optsExpFig);
%     %
%     %     % angle uts0 vs uts
%     %     figure; setFigPos(2,5);
%     %     imagesc(reshape(dA,nEH*nTarg*nPr,nTspp)); colorbar;
%     %     xlabel('ts'); ylabel('conditions (EH, RL, SL)');
%     %     applytofig(gcf,optsExpFig);
%     
% else
%     D1=D; % prior support
%     for i=1:length(D1)
%         D1(i).epochStarts=1;
%         D1(i).epochColors=D1(i).epochColors(1,:);
%         idPr=2-D1(i).idShortTrial; % 1 for short, 2 for long
%         if idAllMeas
%             D1(i).data=D1(i).data(:,1:D1(i).ts);
%         else
%             D1(i).data=D1(i).data(:,(T{idPr}(1)-tBuffer+1):D1(i).ts);
%         end
%     end % for i=1:length(D1)
%     
%     % delete extra fields in data to prevent error in cvreducedims
%     % : only data, conditions, epochStarts, epochColors
%     D12=rmfield(D1,{'id','tp','idShortTrial','ts','idHandEye','theta'});
%     DataHigh(D12,'DimReduce');
% end

%% production
% trajFn2=['trajKS_' num2str(sessId) '_tp_' animalNm{iAnimal} '_condSpec_GPFA_CV_cellFullOnly.mat'];
% disp(trajFn2);
% if exist(trajFn2)
%     D2=load(trajFn2); % D proj_matrix keep_neurons eigenvalues binSize smthWidth thFR optimD
% else
%     D2=D; % tp
%     for i=1:length(D2)
%         D2(i).epochStarts=1;
%         D2(i).epochColors=D2(i).epochColors(1,:);
%         D2(i).data=D2(i).data(:,D2(i).ts:end);
%     end % for i=1:length(D2)
%     
%     % delete extra fields in data to prevent error in cvreducedims
%     % : only data, conditions, epochStarts, epochColors
%     D22=rmfield(D2,{'id','tp','idShortTrial','ts','idHandEye','theta'});
%     DataHigh(D22,'DimReduce');
% end

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

