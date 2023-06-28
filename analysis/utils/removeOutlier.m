function [dataWOout,idValid,pOut,varargout]=removeOutlier(data,nSD,varargin)
% removing outliers
% [dataWOout,idValid,pOut,varargout]=removeOutlier(data,nSD,varargin)
% input: data [n x 1], nSD for criteria of SD
% output: data without outlier, idValid to indicate not outlier in the original
% data, pOut for percetage of outliers

% 2018/9/13
% change "length(data)>15", define new variable nData
% initial SD: sample SD/2
nData=2; % minimum to calculate SD; 15 before

if nSD>0 % simple SD-based outlier removeal
    idNN=(~isnan(data)); % removing NaN first
    idNO=abs(data(idNN)-mean(data(idNN)))<nSD*std(data(idNN));
    idValid=zeros(length(data),1); idValid(idNN)=idNO; idValid=logical(idValid);
    pOut=(length(data)-nnz(idValid))/length(data)*100;
    dataWOout=data(idValid);
    varargout(1)={data(~idValid)};
    
else % use fitting mixture(w*U(min(tp),max(tp))+(1-w)*N(u,s^2))
    % 3SD removel to estimate initial u,s
    nSD=3;
    idNN=(~isnan(data)); % removing NaN first
    idNO=abs(data(idNN)-mean(data(idNN)))<nSD*std(data(idNN));
    idValid=zeros(length(data),1); idValid(idNN)=idNO; idValid=logical(idValid);
    pOut=(length(data)-nnz(idValid))/length(data)*100;
    dataWOout=data(idValid);
    
    if length(varargin)==1 % uniform support
        uMin=varargin{1}(1);uMax=varargin{1}(2);
    elseif isempty(varargin)
        uMin=min(data); uMax=max(data);
    end
    
    if length(data)>nData % 15
        u0=mean(dataWOout);
        s0=std(dataWOout)/2;
        w0=pOut/100;
        p0=[w0 u0 s0];
        p0(p0==[0 0 0])=realmin; % to make sure
        if p0(1)==1, p0(1)=1-realmin; end;
        
        % 3 free parameters:
        % w as lapse rate
        % u, s for normal distribution
        % for now, use [0,max(tp)] as support of U (alternative: min and max of obsered tp, [0, 3*ts] or even fittable?)
%         likeMixUG=@(t,w,u,s) (w/(uMax-uMin)+(1-w).*normpdf(t,u,s)); % or 0
        likeMixUG=@(t,w,u,s) (w/(uMax-uMin)+(1-w).*normpdf(t,u,s));
        opts=statset('FunValCheck','on','MaxFunEvals',1000*length(p0),'MaxIter',1000*length(p0)); % ,...
        % 'Robust','on','WgtFun','logistic','Display','iter');
%         if sum(p0<=[0 0 0])|sum(p0>=[1 Inf Inf])
%             disp('');
%         end

% try
        [phat,pci]=mle(data(:)','pdf',likeMixUG,...
            'start',p0,'lowerbound',[-realmin 0 0],'upperbound',[1 Inf Inf],...
            'options',opts);
% catch
%     disp('');
% end
        u=phat(2); s=phat(3);
        %     disp(phat);
        pOut=phat(1)*100;
        
        % determine outlier based on likelihood
        L=normpdf(data,u,s);
        idValid=1/(uMax-uMin)<L; % max(data)-min(data))<L; % min(data) or 0
        dataWOout=data(idValid);
        varargout(1)={data(~idValid)};
    else
        varargout={[]};
    end
    
end

% % debug
% a=[rand(10,1); NaN];
% [d,id,p]=removeOutlier(a,3)

% data=[normrnd(500,50,90,1);unifrnd(100,1000,10,1)];
% [dataWOout,id,pOut]=removeOutlier(data,0);

function L=normpdf(t,u,s)
L=(1./sqrt(2*pi*(s.^2))).*exp(-((t-u).^2)./(2*(s.^2)));