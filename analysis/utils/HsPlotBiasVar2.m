function HsPlotBiasVar2 % (sessId,idHandEye,theta,ts,tp,idShortTrial,wFitSessCond,idOut)

% 2019/3/27
% separate figure for bias & var (model vs data)

% 2018/9/17
% not session-specific: wFitCond

% 2018/9/3 plot bias vs variance for each session
% now session-specific (line data-model, filled for model, circle for data)
% wFitSessCond: 17sess x 2EH x 2RL: w_m, w_p, offset1(long), offset2

% plot bias versus variance of each prior
% color code depends on Colors2P (default: Colors2P(3))
    
%% init
S=initRSG2prior;
v2struct(S);
% iS=unique(sessId); % 161203
% iEH=unique(idHandEye); % 1hand0eye
% iTarg=unique(theta); % 0Right180left
% idShortTrial: 1short0long

% simulation BLS nSim times and average them
nSim=500;
nRep=1; % 50 times more trials than data: now condition-specific, 1600 trials/40 cond: 40 trials/cond

% plot
lw=0.5;
ms=2;
tmpMarker={'s','d';'o','^'}; % eyeRight, eyeLeft; handRight, handLeft

idDebug=0; % 1;

%% main
iFig=1;
for iAnimal=1:nAnimal
    
    if iAnimal==1
        tmpMarker2='^';
    else
        tmpMarker2='o';
    end
    load([animalNm{iAnimal} '_RSGprior_DMFC.mat']);
    figure(iFig); setFigPos(iAnimal,1); ha; % bias var
    figure(iFig+10); setFigPos(iAnimal,2); ha; %
    figure(iFig+11); setFigPos(iAnimal,3); ha; %
    
    for j=1:nEH
        for k=1:nTarg
            disp([ehNm{j} targNm{k}]);
            
            wm=wFitCond(j,k).w_m;
            wp=wFitCond(j,k).w_p;
            bL=wFitCond(j,k).offset1;
            bS=wFitCond(j,k).offset2;
            
            %         for jS=1:length(iS) % session
            %             disp(['===== ' num2str(iS(jS)) ' =====']);
            %             wm=wFitSessCond(jS,j,k).w_m;
            %             wp=wFitSessCond(jS,j,k).w_p;
            %             bL=wFitSessCond(jS,j,k).offset1;
            %             bS=wFitSessCond(jS,j,k).offset2;
            
            %% simulation BLS nSim times and average them
            % do simulation & calculate mean/SD
            mutSim=nan(nPr,nTspp,nSim); % mean [2priors x 5ts/pr x #simulation]
            sdtSim=nan(nPr,nTspp,nSim); % SD  [2priors x 5ts/pr x #simulation]
            rmsBiasSim=nan(nPr,nSim); % mean Bias  [2priors x #simulation]
            rmsVSim=nan(nPr,nSim); % mean sqrtVar  [2priors x #simulation]
            for iSim=1:nSim
                for i=1:nPr % for each prior
                    tmpId= idHandEye==(2-j) & theta==180*(k-1) & idShortTrial==(2-i) & ~idOut; % sessId==iS(jS) &
                    tmpT=repmat(T(tmpId),nRep,1);
                    
                    % simulation
                    offset=bS*(2-i)+bL*(i-1);
                    tm=tmpT+wm*randn(size(tmpT)).*tmpT;
                    te=BayesEst(tm,wm,[min(tmpT) max(tmpT)],'uniform')+offset;
                    tpModel=te+wp*randn(size(te)).*te;
                    
                    % calculate mean/SD
                    tmpTlist=unique(tmpT);
                    for jT=1:nnz(tmpTlist)
                        id=(tmpT==tmpTlist(jT));
                        
                        mutSim(i,jT,iSim)=mean(tpModel(id));
                        sdtSim(i,jT,iSim)=std(tpModel(id));
                    end
                    
                    % prior-specific meanBias & sqrtVar
                    rmsBiasSim(i,iSim)=sqrt(mean((mutSim(i,:,iSim)-tmpTlist(:)').^2));
                    rmsVSim(i,iSim)=sqrt(mean((sdtSim(i,:,iSim)).^2));
                end % for i=1:nPr
            end % for iSim=1:nSim
            
            %% calculate mean and SD for each ts
            mut=nan(nPr,nTspp); % mean [2priors x 5ts/pr]
            sdt=nan(nPr,nTspp); % SD
            rmsBias=nan(nPr,1); % mean Bias [2priors x 1]
            rmsV=nan(nPr,1); % mean sqrtVar
            for i=1:nPr
                for iTs=1:length(Tmat{i})
                    
                    tmpId=idHandEye==(2-j) & theta==180*(k-1) & idShortTrial==(2-i) & T==Tmat{i}(iTs)  & ~idOut; % sessId==iS(jS) &
                    
                    mut(i,iTs)=mean(t(tmpId));
                    sdt(i,iTs)=std(t(tmpId));
                    
                end % for iTs=1:length(Tmat{i})
                
                % prior-specific meanBias & sqrtVar
                rmsBias(i)=sqrt(mean((mut(i,:)-Tmat{i}(:)').^2));
                rmsV(i)=sqrt(mean((sdt(i,:)).^2));
            end % i pr
            
            %% plot
            
            for i=1:nPr
                x1=rmsBias(i); % data
                y1=rmsV(i);
                % choose best
                tmpE=sqrt(((rmsBiasSim(i,:)-x1).^2)+((rmsVSim(i,:)-y1).^2));
                iBest=min(tmpE)==tmpE;
                x2=rmsBiasSim(i,iBest); % mean(rmsBiasSim(i,:)); % model
                y2=rmsVSim(i,iBest); % mean(rmsVSim(i,:));
                
                figure(iFig); % bias vs var
                plot([x1;x2],[y1;y2],'-','color',pplot.cmap{(i-1)*2+1},'linewidth',lw); % line %%%%%
                plot(x1,y1,tmpMarker{j,k},'color',pplot.cmap{(i-1)*2+1},'linewidth',lw,'markersize',ms,'markerfacecolor',pplot.cmap{(i-1)*2+1}); % data
                plot(x2,y2,tmpMarker{j,k},'color',pplot.cmap{(i-1)*2+1},'linewidth',lw,'markersize',ms,'markerfacecolor','w'); % model
                drawnow;
                
                figure(iFig+10); % data vs model
                plot(x1,x2,tmpMarker2,'color',pplot.cmap{(i-1)*2+1},'linewidth',lw,'markersize',ms,'markerfacecolor','w'); % bias
                drawnow;
                
                figure(iFig+11); % data vs model
                plot(y1,y2,tmpMarker2,'color',pplot.cmap{(i-1)*2+1},'linewidth',lw,'markersize',ms,'markerfacecolor','w'); % var
                drawnow;
                
                %% debug
                if idDebug
                    threshold=50;
                    if abs(x1-x2)>threshold | abs(y1-y2)>threshold
                        % raw data
                        tmpId= idHandEye==(2-j) & theta==180*(k-1) & idShortTrial==(2-i) & ~idOut; % sessId==iS(jS) &
                        hFigTmp=figure;
                        plotTpTs(T(tmpId),t(tmpId),idDebug,pplot.cmap{(i-1)*2+1},hFigTmp);ha;
                        errorbar(Tmat{i}(:),squeeze(mutSim(i,:,iBest)),squeeze(sdtSim(i,:,iBest)),'color',pplot.cmap{(i-1)*2+1}); % BLS
                        waitforbuttonpress; close;
                    end % if abs(x1-x2)>threshold | abs(y1-y2)>threshold
                end % idDebug
                
            end          % i pr
            %         end % for jS=1:length(iS) % session
            
            figure(iFig); % bias vs var
            % plot    circle
            axis tight;
            %         tmpMax=1.3*max([xlim ylim]);
            %         for i=linspace(0,round(tmpMax),5), plotCircleHalf(gca,i,0.5+[0 0 0]); end;
            %         if iS(1)<170000 % H
%             if iAnimal==1
                axis([0 150 0 150]);  % tight;
                for i=50:50:150, plotCircleHalf(gca,i,0.5+[0 0 0]); end;
%             else
%                 axis([0 200 0 200]);  % tight;
%                 for i=50:50:200, plotCircleHalf(gca,i,0.5+[0 0 0]); end;
%             end
            set(gca,'xtick',0:50:250,'ytick',0:50:250,'tickdir','out','ticklength',[0.02 0.02]); xlabel('Average Bias (ms)');ylabel('Average Var^1^/^2 (ms)'); % ylabel('$\sqrt{Var}$ (ms)','Interpreter','Latex');
            %         applytofig(gcf,optsExpFig);
            
            figure(iFig+10); % data vs model
            %         if iS(1)<170000 % H
%             if iAnimal==1
                axis([0 150 0 150]);  % tight;
%             else
%                 axis([0 200 0 200]);  % tight;
%             end
            plotIdentity(gca);
            set(gca,'xtick',0:50:250,'ytick',0:50:250,'tickdir','out','ticklength',[0.02 0.02]); xlabel('Data Average Bias (ms)');ylabel('Bayesian Model Average Bias (ms)');
            %         applytofig(gcf,optsExpFig);
            
            
            figure(iFig+11); % data vs model
            %         if iS(1)<170000 % H
%             if iAnimal==1
                axis([0 150 0 150]);  % tight;
%             else
%                 axis([0 200 0 200]);  % tight;
%             end
            plotIdentity(gca);
            set(gca,'xtick',0:50:250,'ytick',0:50:250,'tickdir','out','ticklength',[0.02 0.02]); xlabel('Data Var^1^/^2 (ms)');ylabel('Bayesian Model Var^1^/^2 (ms)');
            %         applytofig(gcf,optsExpFig);
            
            
        end % for k=1:nTarg
    end % for j=1:nEH
    
    iFig=iFig+100;
    
end % for iAnimal=1:nAnimal


%% subfunctions



function estimate = BayesEst(...
    measurements, ...
    weberMeasure, ...
    priorParams,  ...
    prior)
%% Using Matlab's integral() function

%measurements must be a column vector of the measurement values
%
%weberMeasure is the Weber fraction of the measurement. Pass in one value
%or a column vector with length equal to measurements.
%
%priorParams defines the distribution. For a uniform prior this is the just
%the support (start and end). For a Gaussian prior it's the mean and
%variance.
%
%prior is the prior type. Leave blank for uniform.

wm = weberMeasure;
if numel(wm) == 1
    wm = repmat(wm,size(measurements,1),1);
end

if ~exist('prior','var') || isempty(prior)
    prior = 'uniform';
end

switch prior
    case 'uniform'
        alpha = 1./integral(@(s)posteriorNum(s,measurements,wm), ...
            priorParams(1), priorParams(2),'ArrayValued',true);
        estimate = integral(@(s)alpha.*s.*posteriorNum(s,measurements,wm), ...
            priorParams(1), priorParams(2),'ArrayValued',true);

end


function numValue = posteriorNum(s,m,wm)
%% This is the simplest version of the calculation. 
%Does not allow the measurements to covary.

sigmaDet = ((s*wm).^2).^size(m,2);
numValue = 1./sqrt(2*pi*sigmaDet).*exp(-.5*sum((s-m).^2./(s*wm).^2,2));



function plotCircleHalf(h,radius,color)
% assuming x ranging from negative to positive
x=linspace(0,radius,100);
y=sqrt((radius^2)-(x.^2));
plot(x,y,':','color',color); hold on;

