function [mut,sdt,nt,varargout]=plotTpTs(T,tmat,varargin)

% [mut,sdt,varargout]=plotTpTs(T,tmat,varargin)
%
% varargin: idPlot) or idPlot, cmap,hFig) or idPlot, cmap,hFig,linestyle) or idPlot, cmap,hFig,linestyle,nSD)
% or idPlot, cmap,hFig,linestyle,nSD,nDataPerBin)
% varargout: nothing or h,legendLabel

% use shadedErrorBar

idPlot=varargin{1};
if nargin==5
    cmap=varargin{2};
    hFig=varargin{3};
elseif nargin==6
    cmap=varargin{2};
    hFig=varargin{3};
    linestyle=varargin{4};
elseif nargin==7
    cmap=varargin{2};
    hFig=varargin{3};
    linestyle=varargin{4};
    nSD=varargin{5};
    
elseif nargin==8
    cmap=varargin{2};
    hFig=varargin{3};
    linestyle=varargin{4};
    nSD=varargin{5};    
    nDataPerBin=varargin{6};   
end

if ~exist('nDataPerBin')
    nDataPerBin=9;
end

if ~exist('cmap')
    cmap=[0 0 0];
end
if ~exist('linestyle')
    linestyle='-';
end
if ~exist('nSD')
    nSD=10;
end
Tmat=unique(T);
if nnz(Tmat)>1 & nnz(Tmat)<10
    idDisc=true;
    dT=Tmat(end)-Tmat(end-1); % 100ms assuming uniformly sampled
elseif nnz(Tmat)==1
        dT=100;
        idDisc=true;
else
    idDisc=false;
    
    edge=binVarWidth(T,nDataPerBin);
    [nHistc,binId]=histc(T,edge);
    binIdUni=unique(binId);
    Tmat=[]; T2=zeros(size(T));
    for iEdge=1:(length(binIdUni))
        tmpId=binIdUni(iEdge)==binId;
        Tmat=[Tmat; mean(T(tmpId))];
        T2(tmpId)=mean(T(tmpId));
    end
    
    varargout{4}=edge;
    
    dT=300;
    
end
pdT=0.2; % distribution of scatter plots
nT=length(Tmat);
mut=zeros(nT,1);
sdt=zeros(nT,1);
nt=zeros(nT,1);

if idPlot
    if exist('hFig','var')
        figure(hFig);
    else
        figure;
    end
end% idPlot

for i=1:nT
    if idDisc
        id=T==Tmat(i);
    else
        id=T2==Tmat(i);
    end
    t=tmat(id);
    [tClean,idNoOut,pOut]=removeOutlier(t(t>0),nSD); % 1.5); % 10); % 3); no outlier removal assuming it's done before plotting
    if nnz(idNoOut)~=length(idNoOut) % pOut~=0,
        %         disp(['check data: ' num2str(pOut) '% data (n=' num2str(length(t)-nnz(idNoOut)) ') removed as outliers']);
    end
    [mut(i),sdt(i),nt(i)]=meanSDwoNeg(tClean); % 3SD criteria
    
    varargout{3}=Tmat;
    
    % plot T vs t
    if idPlot
        %         plot(Tmat(i)+(rand(size(t))-0.5)*dT*pdT,t,'.','color',cmap,'markerfacecolor','w','markersize',.2); drawnow; hold all
        if idDisc
            plot(Tmat(i)+(rand(size(tClean))-0.5)*dT*pdT,tClean,'.','color',cmap,'markerfacecolor','w','markersize',6);hold all; drawnow; pause(1); % 11
        else
            plot(T(T2==Tmat(i)),tmat(T2==Tmat(i)),'.','color',cmap,'markerfacecolor','w','markersize',6);hold all; drawnow; pause(1);
        end
        
        %     if i==1, pause(0.5); end;
    end % idPlot
end

if length(mut)>1 && length(sdt)>1
    stats=regstats(mut,Tmat,'linear',{'tstat'}); % beta yhat r mse rsquare tstat
%     fprintf(1, 'intercept: %d (p=%d), slope: %d (p=%d)\n',[stats.tstat.beta(1) stats.tstat.pval(1) stats.tstat.beta(2) stats.tstat.pval(2)]);
    varargout{1}=stats.tstat.beta(2);  varargout{2}=stats.tstat.beta(1);
%     legendLabel=['t_p=' num2str(tmpLegS1,2) '*t_s+' num2str(tmpLegS2,3)];
else
    varargout{1}=nan;varargout{2}=nan;% if E only one T, slope=nan;
%     legendLabel=['\mu(t_p)\pm\sigma(t_p):' num2str(mut,4) '\pm ' num2str(sdt,4)];
end

if idPlot
%     h=errorbar(Tmat,mut,sdt,linestyle,'color',cmap,'linewidth',2); drawnow;
%     h=errorbar(Tmat+0.15*(Tmat(2)-Tmat(1))*rand,mut,sdt,linestyle,'color',cmap,'linewidth',2); drawnow;
%         h=plot(Tmat,mut,linestyle,'color',cmap,'linewidth',0.5,'markerfacecolor','w','markersize',6); drawnow;

if length(Tmat)==1
    h=errorbar(Tmat,mut,sdt,linestyle,'color',cmap,'linewidth',2); drawnow;
else
    h=shadedErrorBar(Tmat,mut,sdt,{linestyle,'color',cmap,'linewidth',0.5,'markerfacecolor','w','markersize',6},1); drawnow; % 2
end
    % simple linear regression as a model-free check of prior effect (slope<1 & intercept>0)
    if length(mut)>1 && length(sdt)>1
        stats=regstats(mut,Tmat,'linear',{'tstat'}); % beta yhat r mse rsquare tstat
        fprintf(1, 'intercept: %d (p=%d), slope: %d (p=%d)\n',[stats.tstat.beta(1) stats.tstat.pval(1) stats.tstat.beta(2) stats.tstat.pval(2)]);
        tmpLegS1=stats.tstat.beta(2); tmpLegS2=stats.tstat.beta(1);
        legendLabel=['t_p=' num2str(tmpLegS1,2) '*t_s+' num2str(tmpLegS2,3)];
    else
        legendLabel=['\mu(t_p)\pm\sigma(t_p):' num2str(mut,4) '\pm ' num2str(sdt,4)];
    end
    
    axis tight; 
%     if mut(end)+sdt(end)*3>mut(1)-sdt(1)*3 
%         ylim([mut(1)-sdt(1)*3 mut(end)+sdt(end)*3]);
%     end
    plotIdentity(gca); 
%     weberF=0.2;
%     plotWeberLine(gca,weberF);
    xlabel('t_s (ms)'); ylabel('t_p (ms)');
    
%     xtickMat=get(gca,'xtick'); 
%     if ~isempty(xtickMat)
%         set(gca,'xtick',sort(xtickMat(unique(round(linspace(1,length(xtickMat),5))))));
%     end
% %     ytickMat=get(gca,'ytick'); set(gca,'ytick',sort(ytickMat(round(linspace(1,length(ytickMat),5)))));
    
    if isstruct(h)
        varargout{1}=h.mainLine;
    else
        varargout{1}=h;
    end
    varargout{2}=legendLabel;
end



%%
% function [dataWOout,id,pOut]=removeOutlier(data,nSD)
% % removing outliers
% % input: data [n x 1], nSD for criteria of SD
% % output: data without outlier, id to indicate not outlier in the original
% % data, pOut for percetage of outliers
% 
% idNN=(~isnan(data)); % removing NaN first
% idNO=abs(data(idNN)-mean(data(idNN)))<nSD*std(data(idNN));
% id=zeros(length(data),1); id(idNN)=idNO; id=logical(id);
% pOut=(length(data)-nnz(id))/length(data)*100;
% dataWOout=data(id);
% 
% % % debug
% % a=[rand(10,1); NaN];
% % [d,id,p]=removeOutlier(a,3)


%%
function [m sd n]=meanSDwoNeg(x)
n=nnz(x>0);
m=mean(x(x>0));
sd=std(x(x>0),1);

%%
function plotIdentity(hAx)
% plot identity line with black dotted
% input: hAx for gca of current plot
%
% 2013/10/19
% hansem sohn

hold on;
x=get(hAx,'xlim'); y=get(hAx,'ylim');
if strcmp(get(gca,'xscale'),'log') & strcmp(get(gca,'yscale'),'log')
    loglog([max([x(1);y(1)]);min([x(2);y(2)])],[max([x(1);y(1)]);min([x(2);y(2)])],'k:');
else
    plot([max([x(1);y(1)]);min([x(2);y(2)])],[max([x(1);y(1)]);min([x(2);y(2)])],'k:');
end
hold on;