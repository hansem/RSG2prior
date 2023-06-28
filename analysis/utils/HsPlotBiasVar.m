function HsPlotBiasVar(d,ctidx,varargin)

% plot bias versus variance of each prior
% color code depends on Colors2P (default: Colors2P(3))

% Input structure: d
% d.ts - Nx1 where N: no of trials - sample
% d.tp - Nx1 where N: no of trials - responses
% ctidx - for each trial in d.ts specifies whether it belongs to prior 1 or
%         2. Enteries 1 or 2. 1 for short

% optional input to show BLS predictions: dots by simulating BLS multiple times
    % dBLS.wm - Wm at which BLS was calculated
    % dBLS.wp - Wp at which BLS was calculated
    % dBLS.offset1 - movement offset for long (can be omitted)
    % dBLS.offset2 - movement offset for short (can be omitted)
    % dBLS.nSim - # of simulations/dots, 50 by default
    
%% main
figure; hold all; hF=[]; % handles for legends    

% initial
iPr=unique(ctidx); % 1 or 2
nPr=length(iPr);
nTsPr=length(unique(d.ts(ctidx==iPr(1)))); % # ts per prior: assuming identical across priors

cmap=Colors2P(3); % use overlap colors to maximize contrast
cmap=squeeze([cmap(1,end,:);cmap(2,1,:)]); % [2priors x 3]

%% whether to do BLS simulation

if ~isempty(varargin)
    dBLS=varargin{1}; % field: wm, wp, offset1, offset2
    
    % if no offset
    if ~isfield(dBLS,'offset1') & ~isfield(dBLS,'offset2')
        dBLS.offset1=0;        dBLS.offset2=0;
    end % if ~isfield(dBLS,'offset1') & ~isfield(dBLS,'offset2')
    dBLS.offset=[dBLS.offset2 dBLS.offset1]; % offset1 for long
    
    if ~isfield(dBLS,'nSim')
        dBLS.nSim=50;
    end
    
    % do simulation & calculate mean/SD
    mutSim=nan(nPr,nTsPr,dBLS.nSim); % mean [2priors x 5ts/pr x #simulation]
    sdtSim=nan(nPr,nTsPr,dBLS.nSim); % SD  [2priors x 5ts/pr x #simulation]
    rmsBiasSim=nan(nPr,dBLS.nSim); % mean Bias  [2priors x #simulation]
    rmsVSim=nan(nPr,dBLS.nSim); % mean sqrtVar  [2priors x #simulation]
    for iSim=1:dBLS.nSim
        for i=1:nPr % for each prior
            tmpT=d.ts(ctidx==iPr(i));
            
            % simulation
            tm=tmpT+dBLS.wm*randn(size(tmpT)).*tmpT;
            te=BayesEst(tm,dBLS.wm,[min(tmpT) max(tmpT)],'uniform')+dBLS.offset(i);
            tp=te+dBLS.wp*randn(size(te)).*te;
            
            % calculate mean/SD
            tmpTlist=unique(tmpT);
            for j=1:nnz(tmpTlist)
                id=(tmpT==tmpTlist(j));
                
                mutSim(i,j,iSim)=mean(tp(id));
                sdtSim(i,j,iSim)=std(tp(id));                
            end
            
            % prior-specific meanBias & sqrtVar
            rmsBiasSim(i,iSim)=sqrt(mean((mutSim(i,:,iSim)-tmpTlist(:)').^2));
            rmsVSim(i,iSim)=sqrt(mean((sdtSim(i,:,iSim)).^2));
            
            % plot
            plot(rmsBiasSim(i,iSim),rmsVSim(i,iSim),'.','color',cmap(i,:),'markersize',3);drawnow;
        end % for i=1:nPr
    end % for iSim=1:nSim
    
    % plot    
    axis tight;
    tmpMax=1.3*max([xlim ylim]);
    for i=linspace(0,round(tmpMax),5), plotCircleHalf(gca,i,0.5+[0 0 0]); end;
    
end % if ~isempty(varargin)

%% calculate mean and SD for each ts
mut=nan(nPr,nTsPr); % mean [2priors x 5ts/pr]
sdt=nan(nPr,nTsPr); % SD
rmsBias=nan(nPr,1); % mean Bias [2priors x 1]
rmsV=nan(nPr,1); % mean sqrtVar
for i=1:nPr
    tsSet=unique(d.ts(ctidx==iPr(i)));
    for j=1:numel(tsSet)
        id=(d.ts==tsSet(j) & ctidx==iPr(i));
        
        mut(i,j)=mean(d.tp(id));
        sdt(i,j)=std(d.tp(id));        
        
    end % for j=1:numel(tsSet)
    
    % prior-specific meanBias & sqrtVar
    rmsBias(i)=sqrt(mean((mut(i,:)-tsSet(:)').^2));
    rmsV(i)=sqrt(mean((sdt(i,:)).^2));
            
    % plot
    hTmp=plot(rmsBias(i),rmsV(i),'o','color',cmap(i,:),'markersize',6,'markerfacecolor','w');drawnow;hF=[hF; hTmp];
    
end % for i=1:nPr
        
% if no simulation, draw circle here
if isempty(varargin)
    axis tight;
    tmpMax=1.3*max([xlim ylim]);
    for i=linspace(0,round(tmpMax),5), plotCircleHalf(gca,i,0.5+[0 0 0]); end;    
end

% legend(hF,'short','long','location','southwest'); legend boxoff;
xlabel('mean Bias (ms)'); ylabel('$\sqrt{Var}$ (ms)','Interpreter','Latex');


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

