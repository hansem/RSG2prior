function setFigPos(i,j)
drawnow;
load pplot.mat;
pause(0.01);
if j<0
    set(gcf,'position',pplot.(['rect' num2str(i) '_m' num2str(-j)]));
else
    if i>2
        i=mod(i-1,2)+1;
    end
    if j>6
        j=mod(j-1,6)+1;
    end
    set(gcf,'position',pplot.(['rect' num2str(i) '_' num2str(j)]));
end
  