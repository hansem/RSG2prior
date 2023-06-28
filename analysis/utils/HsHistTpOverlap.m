function HsHistTpOverlap(d,ctidx,varargin)

% function HsHistTpOverlap(d,ctidx,varargin)

% plot histogram of tp for each prior @ overlap ts
% color code depends on Colors2P (default: Colors2P(3))

% Input structure: d
% d.ts - Nx1 where N: no of trials - sample
% d.tp - Nx1 where N: no of trials - responses
% ctidx - for each trial in d.ts specifies whether it belongs to prior 1 or
%         2. Enteries 1 or 2. 1 for short
% varargin for idStim

%% init
if ~isempty(varargin)
    idStim=varargin{1};
else
    idStim=0;
end

%% main
figure; hold all; hF=[]; % handles for legends    

% figure settings
nBin=20;
DisplayStyle='stairs';
barWidth=0.4;
markersize=10;

% initial
iPr=unique(ctidx); % 1 or 2
nPr=length(iPr);
nTsPr=length(unique(d.ts(ctidx==iPr(1)))); % # ts per prior: assuming identical across priors

cmap=Colors2P(1+4*idStim); % use overlap colors to maximize contrast
cmap=squeeze([cmap(1,end,:);cmap(2,1,:)]); % [2priors x 3]

% identifying overlap ts
tOverlap=max(unique(d.ts(ctidx==iPr(1))));
if tOverlap~=min(unique(d.ts(ctidx==iPr(2)))) % double check
    disp('more than 1 overlap?');
end

% getting bin locations
tpOveralp=d.tp(d.ts==tOverlap);
x=linspace(min(tpOveralp),max(tpOveralp),nBin);
dx=x(2)-x(1);

% plot
tpOverlap=cell(nPr,1);
for i=1:nPr
    id = d.ts==tOverlap & ctidx==iPr(i);    
    tpOverlap{i}=d.tp(id);
    
    h=histogram(tpOverlap{i},x); % nBin); % ,'DisplayStyle',DisplayStyle);
    
     Xval    = h.BinEdges + h.BinWidth*0.5 +h.BinWidth*(i-(nPr+1)/2)/2; % binWidth*0.25 for short, 0.75 for long
     Xval    = Xval(1:end-1);
     Yval    = h.Values;
     delete(h);
                
     hTmp      = bar(Xval,Yval,'BarWidth',barWidth,'FaceColor',cmap(i,:),'EdgeColor','none');
    hF=[hF; hTmp];
    
end
    
axis tight; 
maxYlim=max(ylim);
% plot mean
for i=1:nPr
    plot(mean(tpOverlap{i}),1.1*maxYlim,'color',cmap(i,:),'markerfacecolor',cmap(i,:),'marker','v','markersize',markersize);
end % for i=1:nPr
xlabel('t_p (ms)'); ylabel('# trials');
plotVertical(gca,tOverlap,[]);

%% subfunctions


function h=plotVertical(hAx,varargin)
% plot vertical lines into the plot
% input: hAx for gca of current plot
% varargin: first cell element for position (if empty, zero)
%               , 2nd cell for color (if empty, black)
% 2014/10/2
% hansem sohn

hold on;

if isempty(varargin)
    cmap=[0 0 0]; x=0;
else
    if isempty(varargin{1})
        x=0;
    else
        x=varargin{1};
    end
    if isempty(varargin{2})
        cmap=[0 0 0];
    else
        cmap=varargin{2};
    end
end

y=get(hAx,'ylim'); h=[];
for i=1:length(x)
    hTmp=plot([x(i); x(i)],[min(y); max(y)],':','color',cmap); hold on;
    h=[h; hTmp];
end


