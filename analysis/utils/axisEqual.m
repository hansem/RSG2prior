function axisEqual

% function axisEqual
%
% making current figure's axis equal
% depends on xlim, ylim

m=min([xlim ylim]);
M=max([xlim ylim]);

axis([m M m M]);