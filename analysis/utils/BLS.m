function te = BLS(tm,wm,support,prior,varargin)
% te = BLS(tm,wm,support,prior,varargin)
% 
% support:[min max] for 'uniform' prior
% if prior is gaussian mixture, support is vector for mu of mixtures,
% varargin for {[mixture proportion]},{[sigmas]}

% 2019/3/15 add option of uniformL1cost
% algorithm: not exact, monte carlo sampling from posterior > median of samples
% dependence: randpdf.m

%% Using Matlab's integral() function if exists
% if size(tm,1)==1,idRow=true; else idRow=false; end;
% tm=tm(:);
% wm = wm(:);
% if numel(wm) == 1
%     wm = repmat(wm,size(tm,1),1);
% end

% examples
% BLS(tmeas,wm,gmobj.mu,'gaussianmixture',gmobj.ComponentProportion,sqrt(squeeze(gmobj.Sigma)))...
% BLS(tmeas,wm,meanSD,'gaussian'

%%
if ~exist('prior','var') || isempty(prior)
    prior = 'uniform';
end

switch prior
    case 'gaussMotor'
        if exist('integral','file')
            alpha = 1./integral(@(s)posteriorGaussianMotor(s,tm,wm,support(1),support(2)), ...
                support(1)-4*support(2), support(1)+ 4*support(2), ...
                'ArrayValued',true);
            te = integral(@(s)alpha.*s.* ...
                posteriorGaussianMotor(s,tm,wm,support(1),support(2)), ...
                support(1)-4*support(2), support(1)+ 4*support(2), ...
                'ArrayValued',true);
        else
            alpha = 1./quadv(@(s)posteriorGaussianMotor(s,tm,wm,support(1),support(2)), ...
                support(1)-4*support(2), support(1)+ 4*support(2), ...
                'ArrayValued',true);
            te = quadv(@(s)alpha.*s.* ...
                posteriorGaussianMotor(s,tm,wm,support(1),support(2)), ...
                support(1)-4*support(2), support(1)+ 4*support(2), ...
                'ArrayValued',true);
        end
        
    case 'product'
        a = support(1,1);
        b = support(1,2);
        c = support(2,1);
        d = support(2,2);
        support = [a*c b*d];
        if exist('integral','file')
            alpha = 1./integral(@(s)ProdUniformPdf( s, a, b, c, d ).* ...
                likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
            te = integral(@(s)alpha.*s.*ProdUniformPdf( s, a, b, c, d ).* ...
                likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
        else
            alpha = 1./quadv(@(s)ProdUniformPdf( s, a, b, c, d ).* ...
                likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
            te = quadv(@(s)alpha.*s.*ProdUniformPdf( s, a, b, c, d ).* ...
                likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
        end
        
    case 'uniform'
        if exist('integral','file')
            alpha = 1./integral(@(s)likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
            te = integral(@(s)alpha.*s.*likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
        else
            alpha = 1./quadv(@(s)likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
            te = quadv(@(s)alpha.*s.*likelihood(s,tm,wm), support(1), support(2), ...
                'ArrayValued',true);
        end

    case 'uniformL1cost'
        nSample=10000; % per each tm
        te=nan(length(tm),1);
        
        for iTm=1:length(tm)
            alpha = 1./integral(@(s)likelihood(s,tm(iTm),wm), support(1), support(2), ...
                'ArrayValued',true); % [1 x 1]
            tmpPdf=@(s)alpha.*likelihood(s,tm(iTm),wm).*(s>support(1)&s<support(2));
            
            teTmp=randpdf(tmpPdf(tm),tm,[nSample 1]);
            
            te(iTm)=median(teTmp);
        
        end % for iTm=1:length(tm)
        
        return;
        
        figure; 
        tm=400:900; 
        te=BLS(tm,0.1,[480 800],'uniformL1cost');
        plot(tm,te,'.');
        te2=BLS(tm,0.1,[480 800],'uniform');
        ha; plot(tm,te2,'r.');
        plotIdentity(gca);
        legend('medianPosterior(L1 cost)','meanPosterior(L2 cost)','location','best'); legend boxoff;
        xlabel('t_m (ms)'); ylabel('t_e (ms)'); xlim([min(tm) max(tm)]); applytofig4keynote

    case 'gaussian'
        if exist('integral','file')
            alpha = 1./integral(@(s)posteriorGaussian(s,tm,wm,support(1),support(2)), ...
                max(1e-6,support(1)-4*support(2)), support(1)+ 4*support(2), ...
                'ArrayValued',true);
            te = integral(@(s)alpha.*s.* ...
                posteriorGaussian(s,tm,wm,support(1),support(2)), ...
                max(1e-6,support(1)-4*support(2)), support(1)+ 4*support(2), ...
                'ArrayValued',true);
        else
            alpha = 1./quadv(@(s)posteriorGaussian(s,tm,wm,support(1),support(2)), ...
                max(1e-6,support(1)-4*support(2)), support(1)+ 4*support(2), ...
                'ArrayValued',true);
            te = quadv(@(s)alpha.*s.* ...
                posteriorGaussian(s,tm,wm,support(1),support(2)), ...
                max(1e-6,support(1)-4*support(2)), support(1)+ 4*support(2), ...
                'ArrayValued',true);
        end
        
    case 'gaussianmixture'
        mu=support(:);
        mixP=varargin{1}(:)'; % disp(mixP);
        if length(mixP)<length(mu)
            mixP=[mixP 1-sum(mixP)]; % last one is determined
        end
        sigma=varargin{2}(:); sigma2=(reshape(sigma,[1 1 length(sigma)])).^2; % vector; only diagonal in covMat but need not to be the same
        
        gmobj=gmdistribution(mu,sigma2,mixP); % plot([mu(1)-4*sigma(1):1:mu(end)+ 4*sigma(end)]',pdf(gmobj,[mu(1)-4*sigma(1):1:mu(end)+ 4*sigma(end)]'));
        
        if exist('integral','file')
            alpha = 1./integral(@(s)pdf2d(gmobj,s).*likelihood(s,tm,wm), ...
                max(1e-6,mu(1)-4*sigma(1)), mu(end)+ 4*sigma(end), ...
                'ArrayValued',true);
            te = integral(@(s)alpha.*s.*pdf2d(gmobj,s).*likelihood(s,tm,wm), ...
                max(1e-6,mu(1)-4*sigma(1)), mu(end)+ 4*sigma(end), ...
                'ArrayValued',true);
        else
            alpha = 1./quadv(@(s)pdf2d(gmobj,s).*likelihood(s,tm,wm), ...
                max(1e-6,mu(1)-4*sigma(1)), mu(end)+ 4*sigma(end), ...
                'ArrayValued',true);
            te = quadv(@(s)alpha.*s.*pdf2d(gmobj,s).*likelihood(s,tm,wm), ...
                max(1e-6,mu(1)-4*sigma(1)), mu(end)+ 4*sigma(end), ...
                'ArrayValued',true);
        end
end
% if idRow, te=te(:)'; end;


end
%%
function numValue = likelihood(s,m,wm)
%% This is the simplest version of the calculation. 
%Does not allow the tm to covary.

sigmaDet = ((s*wm).^2); %.^size(m,2);
numValue = 1./sqrt(2*pi*sigmaDet).*exp(-.5*(s-m).^2./(s*wm).^2); % sum((s-m).^2./(s*wm).^2,2));
end

function p=pdf2d(obj,s)
s1d=s(:);
p1d=pdf(obj,s1d);
p=reshape(p1d,size(s));
end

function numValue = posteriorGaussian(s,m,wm,mu,sigma)

numValue = normpdf(m,s,wm*s).*normpdf(s,mu,sigma);
end


function numValue = posteriorGaussianMotor(s,m,wm,mu,sigma)
%% This is the simplest version of the calculation. 
%Does not allow the measurements to covary.

sigmaDet = ((wm).^2); % .^size(m,2);
sigmaDet2 = ((sigma).^2); % .^size(m,2);
numValue = 1./sqrt(2*pi*sigmaDet).*exp(-.5*(s-m).^2./(wm).^2).*...
    1./sqrt(2*pi*sigmaDet2).*exp(-.5*(s-mu).^2./(sigma).^2); % sum((s-m).^2./(s*wm).^2,2));
end

% function numValue = posteriorGaussianMotor(s,m,wm,mu,sigma)
% 
% numValue = normpdf(m,s,wm).*normpdf(s,mu,sigma);
% end