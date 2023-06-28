function dataNanFill=fillNanCell(data,varargin)
% dataNanFill=fillNanCell(data,varargin)
% filling NaN in cell arrays that has different length in each cell
%
% input: data (vector of cells), each cell has vectors of different length
% varargin for % idLock indicate how data are aligned (0 from front or 1 from end; default:0)
%
% output: matrix [nCell x nDataPoint]

%% 
if ~isempty(varargin)
    idLock=varargin{1};
else
    idLock=0;
end


%%
data=data(:)'; %[1 x trials]
n=max(cellfun(@length,data));

if idLock==0 % from front
    dataNanFill=cellfun(@(x) [x(:); nan(n-length(x),1)],data,'UniformOutput',0);
    dataNanFill=cell2mat(dataNanFill)'; % [#timePoints x #trials]
else
    dataNanFill=cellfun(@(x) [nan(n-length(x),1); x(:)],data,'UniformOutput',0);
    dataNanFill=cell2mat(dataNanFill)'; % [#timePoints x #trials]
end