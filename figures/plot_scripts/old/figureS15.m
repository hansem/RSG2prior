function figureS15

% simulation prior with different varaince


%% init
wm=0.08;
tm=500:1:1100;

pr1=[600 1000];
pr2=[700 900];

te1=BLS(tm,wm,pr1,'uniform');
te2=BLS(tm,wm,pr2,'uniform');

%% main

% plot
figure; ha; setFigPos(1,1);

plot(tm,te1,'-','color',0.3*ones(1,3));
plot(tm,te2,'-','color',0.7*ones(1,3));
plotIdentity(gca);
set(gca,'xtick',min(tm):100:max(tm),'ytick',(min(tm)-200):100:(max(tm)+500),'ticklength',[0.02 0.02],'tickdir','out');
xlim([600 1000]);
remLabel;


% gaussian mixture
mu=[700; 900];
mixP=[0.5; 0.5];
sigma=[25; 25];
sigma2=(reshape(sigma,[1 1 length(sigma)])).^2; 
gmobj=gmdistribution(mu,sigma2,mixP);

figure; ha; setFigPos(1,2);
teGM=BLS(tm(:),wm,mu,'gaussianmixture',mixP,sigma);
plot(tm,teGM,'-','color',0.3*ones(1,3));
plotIdentity(gca);
set(gca,'xtick',min(tm):100:max(tm),'ytick',(min(tm)-200):100:(max(tm)+500),'ticklength',[0.02 0.02],'tickdir','out');
xlim([600 1000]);
remLabel;

figure; ha; setFigPos(2,2); % GM distribution
plot(tm(:),pdf(gmobj,tm(:)),'k-');
axis tight;
set(gca,'xtick',min(tm):100:max(tm),'ytick',(min(tm)-200):100:(max(tm)+500),'ticklength',[0.02 0.02],'tickdir','out');
xlim([600 1000]);
remLabel;