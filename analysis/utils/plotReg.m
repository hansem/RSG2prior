function s=plotReg(x,y,hFig,cmap,varargin)

% function s=plotReg(x,y,hFig,cmap,varargin)
% draw a linear regression line b/t x and y in figure hFig with color cmap
% input: x y hFig cmap varargin(id for robust fig; idRegModel: 0 for no error in variable x, 1 otherwise)
% output: s, structure from regstats 'beta','yhat','r','rsquare','tstat'

if ~isempty(varargin)
    idRobust=varargin{1};
    if length(varargin)>1
        idRegModel=varargin{2};
    else
        idRegModel=0;
    end
    if idRegModel==0
        if idRobust==1
            [B,s]=robustfit(x(:),y(:));
            s.beta=B;
        else
            s=regstats(y(:),x(:),'linear',{'beta','yhat','r','rsquare','tstat'});
        end
    else % 'model-2' regression: PC1 %%%%%
        [m,b,r,sm,sb]=lsqfitma(x,y);
        s.beta=[m;b]; B=s.beta;
        s.yhat=m*x+b;
%         % use SD
%         BINT=[m-sm m+sm;...
%             b-sb b+sb];
        piInv=tinv(0.975,length(x)-2); %/sqrt(length(x)) %  norminv(0.975,0,1); % 1.96
        BINT=[m-sm*piInv m+sm*piInv;...
            b-sb*piInv b+sb*piInv]; % assuming sm/sb is standard error
        
        X=[linspace(min(x),max(x),10)' ones(10,1)]; %[x(:) ones(length(x),1)];
        [xSort,idXsort]=sort(linspace(min(x),max(x),10)'); % x
        ha;
        h=shadedErrorBar(xSort(:)',(X(idXsort,:)*B)',[X(idXsort,:)*BINT(:,2)-X(idXsort,:)*B X(idXsort,:)*B-X(idXsort,:)*BINT(:,1)]',...
            {'-','color',cmap,'linewidth',2},1); drawnow; % 2
        
%         set(h.mainLine,'ZData',100000*ones(size(h.mainLine.XData))); % assuming surf <100000
%         set(h.patch,'ZData',100000*ones(size(h.patch.XData)));
    end
else
    s=regstats(y(:),x(:),'linear',{'beta','yhat','r','rsquare','tstat'});
    idRobust=0;
     idRegModel=0;
end

if idRegModel==0
    if ~isempty(hFig)
        figure(hFig);
        hold all;
    end
    
    if idRobust~=1
        X=[x(:) ones(length(x),1)];
        [xSort,idXsort]=sort(x);
        [B,BINT] = regress(y(:),X); % 1st row for slope, 2nd for intercept
        if ~isempty(hFig)
            shadedErrorBar(xSort(:)',(X(idXsort,:)*B)',[X(idXsort,:)*BINT(:,2)-X(idXsort,:)*B X(idXsort,:)*B-X(idXsort,:)*BINT(:,1)]',...
                {'-','color',cmap,'linewidth',2},1); drawnow; % 2
        end
    else
        if ~isempty(hFig)
            plot(x,s.beta(:)'*[ones(size(x(:)')); x(:)'],'-','color',cmap,'linewidth',2);
        end
    end
end
%% test
% n=30; x = linspace(1,10,n)+3*randn(1,n); % (1:10)';
% y = 10 - 2*x + 3*randn(1,n); y(10) = 0;
% hFig=figure;cmap=[1 0 0]; plot(x,y,'o','color',cmap'); idRobust=1;
% s=plotReg(x,y,hFig,cmap,idRobust)
% hFig=figure; plot(x,y,'o','color',cmap');s=plotReg(x,y,hFig,cmap,0)
% hFig=figure; plot(x,y,'o','color',cmap');s=plotReg(x,y,hFig,cmap,0,1)




%%
function [m,b,r,sm,sb]=lsqfitma(X,Y)
% lsqfitma.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2016 Mar 17.
% 
% M-file to calculate a "MODEL-2" least squares fit.
%
%     The line is fit by MINIMIZING the NORMAL deviates.
%
%     The equation of the line is:     y = mx + b.
%
%     This line is called the MAJOR AXIS.  All points are given EQUAL
%       weight.  The units and range for X and Y must be the same.
%     Equations are from York (1966) Canad. J. Phys. 44: 1079-1086;
%       re-written from Kermack & Haldane (1950) Biometrika 37: 30-41;
%       after a derivation by Pearson (1901) Phil. Mag. V2(6): 559-572.
%
%     Data are input and output as follows:
%
%	    [m,b,r,sm,sb] = lsqfitma(X,Y)
%
%             X    =    x data (vector)
%             Y    =    y data (vector)
%
%             m    =    slope
%             b    =    y-intercept
%             r    =    correlation coefficient
%             sm   =    standard deviation of the slope
%             sb   =    standard deviation of the y-intercept
%
%     Note that the equation passes through the centroid:  (x-mean, y-mean)

% Determine the size of the vector
 
n = length(X);
 
% Calculate sums and other re-used expressions
 
Sx = sum(X);
Sy = sum(Y);
xbar = Sx/n;
ybar = Sy/n;
U = X - xbar;
V = Y - ybar;
 
Suv = sum(U .* V);
Su2 = sum(U .^2);
Sv2 = sum(V .^2);
 
sigx = sqrt(Su2/(n-1));
sigy = sqrt(Sv2/(n-1));
 
% Calculate m, b, r, sm, and sb
 
m = (Sv2 - Su2 + sqrt(((Sv2 - Su2)^2) + (4 * Suv^2)))/(2 * Suv);
b = ybar - m * xbar;
r = Suv / sqrt(Su2 * Sv2);
 
sm = (m/r) * sqrt((1 - r^2)/n);
sb1 = (sigy - sigx * m)^2;
sb2 = (2 * sigx * sigy) + ((xbar^2 * m * (1 + r))/r^2);
sb = sqrt((sb1 + ((1 - r) * m * sb2))/n);

%%
function [m,b,r,sm,sb]=lsqfitgm(X,Y)
% lsqfitgm.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2016 Mar 17.
% 
% M-file to calculate a "MODEL-2" least squares fit.
%
%     The SLOPE of the line is determined by calculating the GEOMETRIC MEAN
%       of the slopes from the regression of Y-on-X and X-on-Y.
%
%     The equation of the line is:     y = mx + b.
%
%     This line is called the GEOMETRIC MEAN or the REDUCED MAJOR AXIS.
%
%     See Ricker (1973) Linear regressions in Fishery Research, J. Fish.
%       Res. Board Can. 30: 409-434, for the derivation of the geometric
%       mean regression.
%
%     Since no statistical treatment exists for the estimation of the
%       asymmetrical uncertainty limits for the geometric mean slope,
%       I have used the symmetrical limits for a model I regression
%       following Ricker's (1973) treatment.  For ease of computation,
%       equations from Bevington and Robinson (1992) "Data Reduction and
%       Error Analysis for the Physical Sciences, 2nd Ed."  pp: 104, and
%       108-109, were used to calculate the symmetrical limits: sm and sb.
%
%     Data are input and output as follows:
%
%	    [m,b,r,sm,sb] = lsqfitgm(X,Y)
%
%             X    =    x data (vector)
%             Y    =    y data (vector)
%
%             m    =    slope
%             b    =    y-intercept
%             r    =    correlation coefficient
%             sm   =    standard deviation of the slope
%             sb   =    standard deviation of the y-intercept
%
%     Note that the equation passes through the centroid:  (x-mean, y-mean)
%
%     WARNING:  Both lsqfitx.m and lsqfity.m must be present in a directory
%               that is in your MATLABPATH in order for this algorithm to
%               execute properly.

% Determine slope of Y-on-X regression

[my] = lsqfity(X,Y);

% Determine slope of X-on-Y regression

[mx] = lsqfitx(X,Y);

% Calculate geometric mean slope

m = sqrt(my * mx);

if (my < 0) && (mx < 0)
	m = -m;
end

% Determine the size of the vector
 
n = length(X);
 
% Calculate sums and means
 
Sx = sum(X);
Sy = sum(Y);
xbar = Sx/n;
ybar = Sy/n;

% Calculate geometric mean intercept

b = ybar - m * xbar;

% Calculate more sums

Sx2 = sum(X.^2);

% Calculate re-used expressions

den = n * Sx2 - Sx^2;

% Calculate r, sm, sb and s2

r = sqrt(my / mx);

if (my < 0) && (mx < 0)
	r = -r;
end

diff = Y - b - m .* X;

s2 = sum(diff .* diff) / (n-2);
sm = sqrt(n * s2 / den);
sb = sqrt(Sx2 * s2 / den);

