function h=plotHorizon(hAx,varargin)
% plot vertical lines into the plot
% input: hAx for gca of current plot
% varargin: first cell element for position (if empty, zero)
%               , 2nd cell for color (if empty, black)
% 2014/10/2
% hansem sohn

hold on;

if isempty(varargin)
    cmap=[0 0 0]; y=0;
else
    if isempty(varargin{1})
        y=0;
    else
        y=varargin{1};
    end
    if isempty(varargin{2})
        cmap=[0 0 0];
    else
        cmap=varargin{2};
    end
end

x=get(hAx,'xlim'); h=[];
for i=1:length(y)
    hTmp=plot([min(x); max(x)],[y(i); y(i)],':','color',cmap); hold on;
    h=[h; hTmp];
end

%% debug
% plot(rand(10,1));
% plotHorizon(gca);
% clf; plotHorizon(gca,0,[1 0 0]);
% plotHorizon(gca,[0 1],[1 0 1]);
% plotHorizon(gca,[],[1 0 1]);
% plotHorizon(gca,[0 1],[]);
