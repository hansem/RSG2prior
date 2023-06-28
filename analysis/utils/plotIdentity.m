function h=plotIdentity(hAx)
% plot identity line with black dotted
% input: hAx for gca of current plot
%
% 2013/10/19
% hansem sohn

hold on;
x=get(hAx,'xlim'); y=get(hAx,'ylim');
if strcmp(get(gca,'xscale'),'log') & strcmp(get(gca,'yscale'),'log')
    h=loglog([max([x(1);y(1)]);min([x(2);y(2)])],[max([x(1);y(1)]);min([x(2);y(2)])],'k:');
else
    h=plot([max([x(1);y(1)]);min([x(2);y(2)])],[max([x(1);y(1)]);min([x(2);y(2)])],'k:');
end
hold on;