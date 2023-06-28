function analVarSingleT

% 2019/3/4
% check mean(sSet) & sigma(sSet)
% w/o vs w/ sqrt transformation
%         H: 161218 161222 (161211 for w/o >55 cells with >1700 trials)
%         G: 170817 170822 170818 170507 170506 (more trials<>more neurons)
%               H/G: sess#12 16 (8) / 8 11 9 2 1 for varargin
% ignore IC (return)
% normlized variance by gradient

% 18/9/26
% sqrt transform applied
% trajKS_161218_pr_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly_newUts_sqrtSingleT.mat
% trajKS_161218_t_p_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly_newUts_sqrtSingleT.mat
%
% checking noise structure: var(Xu) vs var(XPC1), var(XPC2)....

% 2018/9/24
% for all ts, uts1-uts5; sSet1-sSet5
% checking d' by randomly rotating vectors

% 2018/9/19
% variance decrease by changing uts connecting only ts1-ts2, ts4-ts5? idUtsNew
%    trajKS_161218_pr_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly_newUts
%    D.sSet1, D.sSet2 (D.uts1, D.uts2)

%%
anNm='H'; % 'H'; % 'G';
iAnimal=1; % 2; % 1;
sessNum=12; % 12 16 / 8 11 9 2 1

idSqrt=1; % 0; % 1;
if idSqrt
    fnEnd='_newUts_sqrtSingleT';
else
    fnEnd='';
end

initRSG2prior; % T
idPlotHist=0; % 1;
nH=20;
optsExpFig.Width=5;
optsExpFig.Height=5;
optsExpFig.FontSize=16;
optsExpFig.Linewidth=2;

disp(['trajKS_' num2str(sessDate{iAnimal}(sessNum)) '_pr_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly' fnEnd '.mat']);
d=load(['trajKS_' num2str(sessDate{iAnimal}(sessNum)) '_pr_' anNm '_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly' fnEnd '.mat']); % D.sSet1, D.sSet2
% d0=load('trajKS_161218_pr_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly'); % D.sSet
 
if idSqrt
    for iTs=1:nTspp
        sSetP{iTs}=[d.D.(['sSet' num2str(iTs)])]; % [1610 trials x 1]
    end
    
    xro=[d.D.xro]; % [1610 trials x 3]
    uts1=[d.D.uts1];  % [1610 trials x 3]
    uts2=[d.D.uts2]; % [1610 trials x 3]
    uts3=[d.D.uts3];  % [1610 trials x 3]
    uts4=[d.D.uts4]; % [1610 trials x 3]
    uts5=[d.D.uts5];  % [1610 trials x 3]
end

sSet=[d.D.sSet]; % same as d.D.sSet % [1610 trials x 1]

sSet3D=d.sSet; % [1610 trials x 3]

%% test new BLS prediction: variance profile
% 2019/3/4: mean & sigma of sSet
% sSet sIC

% d2=load('trajKS_161218_t_p_H_condSpec_poolSessNew_CV_avgAttB4PCA_cellFullOnly'); % D.sSet
% sIC=[d2.D.sIC]; % [1610 trials x 1]

Vset=nan(nPr,nEH,nTarg,nTspp);
nVset=nan(nPr,nEH,nTarg,nTspp);
Mset=nan(nPr,nEH,nTarg,nTspp);
% VIC=nan(nPr,nEH,nTarg,nTspp);

h=figure; setFigPos(1,1); im=1; % Set
h2=figure; setFigPos(1,2); % IC>mu
h3=figure; setFigPos(1,3); % normalized variance

for i=1:nPr % short long
    for j=0:(nEH-1) % eye hand (RG BO)
        for k=1:nTarg % right left
            for l=1:nTspp
                id=d.idPr==i & (1-d.idHandEye)==j & d.theta==180*(k-1) & d.ts==T{i}(l);
                
                Vset(i,j+1,k,l)=nanvar(sSet(id));
                Mset(i,j+1,k,l)=nanmean(sSet(id));
                %             VIC(i,j+1,k,l)=var(sIC(id));
                
            end % ts
            tmpGrad=gradient(squeeze(Mset(i,j+1,k,:))); % [1 x 5ts]
            for l=1:nTspp
                nVset(i,j+1,k,l)=Vset(i,j+1,k,l)./tmpGrad(l);
            end % ts
            figure(h);
            plot(T{i}(:),squeeze(Mset(i,j+1,k,:)),['-'],'color',pplot.cmap{j*nTarg+k},'markerfacecolor','w');ha; drawnow; % pplot.marker8_2{im}  pplot.cmap{(i-1)*2+1}
            figure(h2);
            plot(T{i}(:),squeeze(Vset(i,j+1,k,:)),['-'],'color',pplot.cmap{j*nTarg+k},'markerfacecolor','w');ha; drawnow; % pplot.marker8_2{im}
            figure(h3);
            plot(T{i}(:),squeeze(nVset(i,j+1,k,:)),['-'],'color',pplot.cmap{j*nTarg+k},'markerfacecolor','w');ha; drawnow; % pplot.marker8_2{im}
            
            %             figure(h2);
            %             plot(T{i}(:),squeeze(VIC(i,j+1,k,:)),[pplot.marker8_2{im} '-'],'color',pplot.cmap{(i-1)*2+1},'markerfacecolor','w');ha; drawnow;
            %             xlabel('t_s (ms)'); ylabel('var. across trials');
            im=im+1;
        end % targ
    end % EH
end % pr
% plot grand averages
Vset2=squeeze(mean(mean(Vset,2),3)); % nan(nPr,nTspp);
Mset2=squeeze(mean(mean(Mset,2),3)); % nan(nPr,nTspp);
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
    plot(T{i}(:),squeeze(Mset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}  pplot.cmap{(i-1)*2+1}
    figure(h2);
    plot(T{i}(:),squeeze(Vset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
    figure(h3);
%     plot(T{i}(:),squeeze(nVset2(i,:)),['-'],'color','k','markerfacecolor','w','linewidth',2);ha; drawnow; % pplot.marker8_2{im}
    shadedErrorBar(T{i}(:),squeeze(nVset2(i,:)),squeeze(nVset3(:,i)),{['-'],'color','k','markerfacecolor','w','linewidth',2},1);ha; drawnow; % pplot.marker8_2{im}
end
figure(h);
xlabel('t_s (ms)'); ylabel('mean across trials'); axis tight; applytofig4keynote;
figure(h2);
xlabel('t_s (ms)'); ylabel('var. across trials');  axis tight; applytofig4keynote;
figure(h3);
xlabel('t_s (ms)'); ylabel('norm. var. across trials');  axis tight; applytofig4keynote;
            
return;





%% compare var(u) and var(u') for all ts1-ts5
V=[]; % first column for new uts, 2nd for old
V2=nan(nPr,nEH,nTarg,nTspp,2); % last for two uts, last-1 for ts1/ts5
dA=nan(nPr,nEH,nTarg,nTspp); % last for two uts or ts1/ts5

% checking noise structure
nPC=size(sSet3D,2); V3=nan(nPr,nEH,nTarg,nTspp,nPC+1); %  last for original uts

h=figure; im=1;


for i=1:nPr % short long
    for j=0:(nEH-1) % hand eye
        for k=1:nTarg % right left
            
            h3=figure; setFigPos(i,j*nTarg+k); ha;
            
            id=d.idPr==i & d.idHandEye==j & d.theta==180*(k-1);
            if idPlotHist, h2=figure; setFigPos(i,j*nTarg+k); end % HR, HL, ER, EL
            
            for iTs=1:nTspp
                x=sSetP{iTs}(id&d.ts==T{i}(iTs));
                x2=sSet(id&d.ts==T{i}(iTs));
                
                V=[V; var(x) var(x2)]; % ts1
                
                V2(i,j+1,k,iTs,:)=[var(x),var(x2)];
                
                % checking noise structure
                V3(i,j+1,k,iTs,end)=var(x2);
                sSetTmp=sSet3D(id&d.ts==T{i}(iTs),:); % [trials x 3PC]
                [coeff score latent]=pca(sSetTmp);
                V3(i,j+1,k,iTs,1:nPC)=latent;
                
                % check distribution
                if idPlotHist
                    figure(h2);
                    subplot(nTspp,1,iTs); histogram(x,nH,'displaystyle','stairs','edgecolor','m'); ha;  histogram(x2,nH,'displaystyle','stairs','edgecolor','r'); ha;
                    drawnow;
                    applytofig(gcf,optsExpFig);drawnow;
                end
                
                figure(h);pause(0.1);
                if iTs==round(nTspp/2) % 3 prior mean
                    plot(V2(i,j+1,k,iTs,1),V2(i,j+1,k,iTs,2),pplot.marker8_2{im},'color',tmpCmap{i,1}(iTs,:),'linewidth',3,'markersize',15); ha; drawnow;
                else
                    plot(V2(i,j+1,k,iTs,1),V2(i,j+1,k,iTs,2),pplot.marker8_2{im},'color',tmpCmap{i,1}(iTs,:)); ha; drawnow;
                end
                
            end % for iTs=1:nTspp
            
            figure(h3); pause(0.1);
            for iPC=1:(nPC+1)
                plot(T{i}(:),squeeze(V3(i,j+1,k,:,iPC)),[pplot.marker8_2{im} '-'],'linewidth',2,'color',pplot.cmap{iPC}); ha; drawnow;
            end % for iPC=1:(nPC+1)
            xlabel('t_s'); ylabel('Var(Xu) or Var(X*PC)');
            applytofig(gcf,optsExpFig);drawnow;
            
            im=im+1;
        end
    end
end
figure(h);
plotIdentity(gca); xlabel('var(new uts)'); ylabel('var(old uts)');
% figure; plot(V(:,1),V(:,2),'ko'); axis tight; plotIdentity(gca);



%% (TBD) checking d' by randomly rotating vectors
% (old) checking sum variance over ts by randomly rotating vectors

nRand=1000; aRand=rand(nRand,1)*90;

Vall=nan(nPr,nEH,nTarg,nRand+1); 

im=1; ifig=1; % var sum over ts vs angle(u,u')
nWin=50; % window for moving average

for i=1:nPr % short long
    for j=0:(nEH-1) % hand eye
        for k=1:nTarg % right left
            % original uts
            tmpV=[];
            for l=1:nTspp
                id2=d.idPr==i & d.idHandEye==j & d.theta==180*(k-1) & d.ts==T{i}(l);
                tmpV=[tmpV; var(sSet(id2))];
            end % ts
            Vall(i,j+1,k,1)=sum(tmpV);
            
            % perturbing uts : common across ts
            id=d.idPr==i & d.idHandEye==j & d.theta==180*(k-1);
            uts=d.D(find(id,1,'first')).xro(:);
            
            for iRand=1:nRand
                u2=randVect(uts,aRand(iRand),1); % [3 x 1]
                
                sSetNew=sSet3D(id,:)*u2; % [trials(id) x 1]
            
                tmpV=[];
                for l=1:nTspp
                    id2=d.ts(id)==T{i}(l);
                    tmpV=[tmpV; var(sSetNew(id2))];
                end % ts
                Vall(i,j+1,k,1+iRand)=sum(tmpV);
                
            end % for iRand=1:nRand
            
            Vall2=squeeze(Vall(i,j+1,k,2:end));
            figure(ifig); setFigPos(i,j*nTarg+k); % HR, HL, ER, EL
            plot(aRand,Vall2,'.','color',pplot.cmap{(i-1)*2+1},'markersize',2);ha; 
            plot(0,squeeze(Vall(i,j+1,k,1)),pplot.marker8_2{im},'color',pplot.cmap{(i-1)*2+1},'markerfacecolor','w'); ha;
            
            [tmpXsort,idSort]=sort(aRand);
            x2=smooth(tmpXsort(:),nWin);y2=smooth(Vall2(idSort),nWin);
            plot(x2,y2,'-','color',pplot.cmap{(i-1)*2+1},'linewidth',2);
            
            drawnow;
            xlabel('angle between u_t_s vs rotated u_t_s'); ylabel('var(Xu) sum across t_s');
          applytofig(gcf,optsExpFig);
          
            im=im+1;ifig=ifig+1;
        end % targ
    end % EH
end % pr


return;


%% compare var(u) and var(u') only for ts1, ts5
V=[]; % first column for new uts, 2nd for old
V2=nan(nPr,nEH,nTarg,2,2); % last for two uts, last-1 for ts1/ts5
dA=nan(nPr,nEH,nTarg,2); % last for two uts or ts1/ts5

h=figure; im=1;

for i=1:nPr % short long
    for j=0:(nEH-1) % hand eye
        for k=1:nTarg % right left
            id=d.idPr==i & d.idHandEye==j & d.theta==180*(k-1);
            
            x=sSet1(id&d.ts==T{i}(1));
            y=sSet2(id&d.ts==T{i}(end));
            x2=sSet(id&d.ts==T{i}(1));
            y2=sSet(id&d.ts==T{i}(end));
            
            V=[V; var(x) var(x2)]; % ts1
            V=[V; var(y) var(y2)]; % ts5
            
            V2(i,j+1,k,1,:)=[var(x),var(x2)];
            V2(i,j+1,k,2,:)=[var(y),var(y2)];
            
            % check distribution 
            if idPlotHist
                figure; setFigPos(i,j*nTarg+k); % HR, HL, ER, EL
                subplot(2,1,1); histogram(x,nH,'displaystyle','stairs','edgecolor','m'); ha;  histogram(x2,nH,'displaystyle','stairs','edgecolor','r'); ha;
                subplot(2,1,2); histogram(y,nH,'displaystyle','stairs','edgecolor','c'); ha;  histogram(y2,nH,'displaystyle','stairs','edgecolor','b'); ha;
                drawnow;
                applytofig(gcf,optsExpFig);
            end
            
            figure(h);
            plot(V2(i,j+1,k,1,1),V2(i,j+1,k,1,2),pplot.marker8_2{im},'color','r'); ha;
            plot(V2(i,j+1,k,2,1),V2(i,j+1,k,2,2),pplot.marker8_2{im},'color','b'); ha;drawnow;
            im=im+1;
        end
    end
end

plotIdentity(gca); xlabel('var(new uts)'); ylabel('var(old uts)');
% figure; plot(V(:,1),V(:,2),'ko'); axis tight; plotIdentity(gca);

% check angle b/t encoding vectors

