function varargout=errorbarXY_meanSem(mx,my,sx,sy,varargin)
% varargout=errorbarXY_meanSem(mx,my,sx,sy,varargin)
% plot errorbars for both X and Y
% input: mean/SEM data for mx,my,sx,sy; cmap,
% marker,markersize,markerfacecolor,markeredgecolor
% hansem.sohn@gmail.com 2018/1/18

if isempty(varargin)
    marker='o';
    cmap=[0 0 0];    
    msize=4;
            markerfacecolor='w';
        markeredgecolor=cmap;
else
    if nargin-4==1 % cmap
        cmap=varargin{1};
        marker='o';
        msize=4;
                markerfacecolor='w';
        markeredgecolor=cmap;
    elseif nargin-4==2 % cmap & marker
        cmap=varargin{1};
        marker=varargin{2};
        msize=4;
                markerfacecolor='w';
        markeredgecolor=cmap;
    elseif nargin-4==3 % cmap & marker
        cmap=varargin{1};
        marker=varargin{2};
        msize=varargin{3};
                markerfacecolor='w';
        markeredgecolor=cmap;
    elseif nargin-4==5 % cmap & marker
        cmap=varargin{1};
        marker=varargin{2};
        msize=varargin{3};
        markerfacecolor=varargin{4};
        markeredgecolor=varargin{5};
    end % if nargin-2==1 % cmap
end % if isempty(varargin)


% plot
hold all;
plot(mx+sx*[-1 1],my*[1 1],'-','color',cmap,'linewidth',1);ha; % horizontal
plot(mx*[1 1],my+sy*[-1 1],'-','color',cmap,'linewidth',1);ha; % vertical
h=plot(mx,my,marker,'color',cmap,'markerfacecolor',markerfacecolor,'markeredgecolor',markeredgecolor,'linewidth',2,'markersize',msize);ha;
varargout(1)={h};


   