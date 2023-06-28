function [n,x,hStair]=histStairs(data,nBin,idPDF,h,varargin)
% plot histogram with stairs: histStairs(data,nBin,idPDF,h,varargin)
% varargin for color, linestyle
% Hansem Sohn 2016/7/12
if idPDF==1 % PDF as output
    [n,x]=histp(data,nBin);
else
    if exist('histcounts') & isempty(nBin)
        [n,x]=histcounts(data);
        n=[n(:);0]; x=x(:);
    else
        [n,x]=hist(data,nBin);
    end
end
figure(h);
if ~isempty(varargin)
    cmap=varargin{1};
    if length(varargin)==2% linestyle
        linestyle=varargin{2};
    else
        linestyle='-';
    end
else
    cmap='k';linestyle='-';
end
hStair=stairs(x,n,'color',cmap,'linestyle',linestyle); hold on;