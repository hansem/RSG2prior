function semOut=sem(data,varargin)

% function semOut=sem(data,varargin)
% input: data, dimension, 'omitnan', idMoment

% calculate standard error mean = std/sqrt(N)
if numel(varargin)==0 % assuming data is a (column) vector
    data=data(:);
    semOut=std(data,0,1)/sqrt(length(data));
elseif numel(varargin)==1
    dim=varargin{:};
    semOut=std(data,0,dim)/sqrt(size(data,dim));
elseif numel(varargin)==2 % including NaN
    dim=varargin{1};
    str=varargin{2}; %'omitnan'
    try
    semOut=std(data,0,dim,str)./sqrt(sum(~isnan(data),dim));
    catch
        semOut=nanstd(data,0,dim)./sqrt(sum(~isnan(data),dim));
    end
elseif numel(varargin)==3 % including NaN
    idMoment=varargin{1};
    dim=varargin{2};
    str=varargin{3}; %'omitnan'
    try
    semOut=std(data,idMoment,dim,str)./sqrt(sum(~isnan(data),dim));
    catch
        semOut=nanstd(data,idMoment,dim)./sqrt(sum(~isnan(data),dim));
    end
end
