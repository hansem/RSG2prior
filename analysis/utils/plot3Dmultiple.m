function h=plot3Dmultiple(x,cmap,linespec,varargin)

% function h=plot3Dmultiple(X,cmap,linespec)
%   draw n (m-time) trajectories in 3D state space
% input 
%       X: n x m x 3 matrix 
%       cmap: (n or 1) x (3 or 1) matrix
%       linespec: n x 1 matrix
%       varargin: linewidth
% output
%       h: handle matrix [n x 1]

ha;
% if ~isempty(x)
    n=size(x,1);
    m=size(x,2);
%     if size(x,1)~=size(y,1)
%         x=repmat(x,size(y,1),1);
%     end

    if size(cmap,2)==1
        cmap=repmat(cmap,size(cmap,1),3);
    end
    if size(cmap,1)==1
        cmap=repmat(cmap,n,1);
    end
        
    if isempty(linespec)
        linespec=repmat({'-'},n,1);
    end
    
    if isempty(varargin)
        linewidth=1;
    else
        linewidth=varargin{1};
    end
    h=nan(n,1);
    
    %% main
    for i=1:n
        h(i)=plot3(x(i,:,1),x(i,:,2),x(i,:,3),linespec{i},'color',cmap(i,:),'linewidth',linewidth); % ,'markersize',3);
    end
% else
%     n=size(y,1);
%     h=nan(n,1);
%     for i=1:n
%         h(i)=plot(y(i,:),linespec,'color',cmap(i,:));
%     end
% end

