function [chi2,pval]=chi2testIndep(x)
% chi^2 to test independence in data x
% input x is a matrix with each element for # of observed samples

nR=size(x,1);
nC=size(x,2);

marginalRow=repmat(sum(x,1),[nR 1]);
marginalCol=repmat(sum(x,2),[1 nC]);
expectedX=marginalRow.*marginalCol./sum(sum(x,1),2);

chi2mat=((x-expectedX).^2)./expectedX;
chi2=sum(sum(chi2mat,1),2);
dof=(nR-1)*(nC-1);
pval=1-chi2cdf(chi2,dof);
return;
% % debug
% x=[31 14 45; 2 5 53; 53 45 2];
% [chi2,pval]=chi2testIndep(x)
% chi2 should be 125.516
% http://math.hws.edu/javamath/ryan/ChiSquare.html

% x=[16 68-16;6 29-6]; [chi2,pval]=chi2testIndep(x)
% x=[16 68-16;11 40-11]; [chi2,pval]=chi2testIndep(x)