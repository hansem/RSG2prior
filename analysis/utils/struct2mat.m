function M=struct2mat(D,fieldname,varargin)

% function M=struct2mat(D,fieldname)
% convert D(dataHigh format)'s field to matrix
% if sizes are different, filling Nan to match size
% input: D(dataHigh format), fieldname(string), varargin for % idLock
% indicate how data are aligned (0 from front or 1 from end; default:0)
% only for columns (for row, alway aligned from top)
% output: M[max # row x max # columns]
% 2018/4/18 hansem@mit.edu

%% 
if ~isempty(varargin)
    idLock=varargin{1};
else
    idLock=0;
end

[X{1:length(D)}]=deal(D.(fieldname)); % convert to cells

nR=cellfun(@(x)size(x,1),X); % nData(D,1); % # rows to find max
nC=cellfun(@(x)size(x,2),X); % nData(D,2);

mnR=max(nR);
mnC=max(nC);

xFillR=cellfun(@(x)[x; nan(mnR-size(x,1),size(x,2))],X,'UniformOutput',0);
if idLock==0 % from front
    xFillRC=cellfun(@(x)[x nan(size(x,1),mnC-size(x,2))],xFillR,'UniformOutput',0);
else
    xFillRC=cellfun(@(x)[nan(size(x,1),mnC-size(x,2)) x],xFillR,'UniformOutput',0);
end

M=cat(3,xFillRC{:,:}); % cell2mat(xFillRC); 

%% debug
% D=[];
% D(1).data=rand(3,4);
% D(2).data=rand(2,5);
% D(3).data=rand(1,7);
% 
% M=struct2mat(D,'data');

