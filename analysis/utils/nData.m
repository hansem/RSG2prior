function N=nData(D,id)

% function nData(D,id)
% get # rows/columns of D.data 
% input: D(dataHigh format), id(1 for row, 2 for column)
% output: N [length(D) x 1]
% 20184/18 hansem@mit.edu

[X{1:length(D)}]=deal(D.data);
N=cellfun(@(x)size(x,id),X);