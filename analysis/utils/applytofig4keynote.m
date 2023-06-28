function applytofig4keynote

set(gcf,'Units','inches');
try
load('/Users/hansem/Dropbox (MIT)/timeArithmetic/pplot.mat'); % optsExpFig4keynote
catch
    load('/Users/seonminahn/Dropbox (MIT)/timeArithmetic/pplot.mat'); % optsExpFig4keynote
end
optsExpFig4keynote.LineMode='scaled';
applytofig(gcf,optsExpFig4keynote);

%%
% optsExpFig4keynote=optsExpFig;
% optsExpFig4keynote.Format='eps';
% optsExpFig4keynote.Width=4;
% optsExpFig4keynote.Height=4;
% optsExpFig4keynote.DefaultFixedLineWidth=2;
% optsExpFig4keynote.Renderer='painters';
% optsExpFig4keynote.Bounds='loose';
% optsExpFig4keynote.FontMode='fixed';
% optsExpFig4keynote.DefaultFixedFontSize=14;
% optsExpFig4keynote.FontSize=16;
% optsExpFig4keynote.LineWidth=1.5;

%% 
% /Users/hansem/Dropbox (MIT)/misc/MATLABpref_backup/jazlab2iMac_20180621/ExportSetup/4keynote.txt

% Version 1
% Format eps
% Preview none
% Width 4
% Height 4
% Units inches
% Color rgb
% Background w
% FixedFontSize 10
% ScaledFontSize 140
% FontMode scaled
% FontSizeMin 10
% FixedLineWidth 2
% ScaledLineWidth auto
% LineMode scaled
% LineWidthMin 2
% FontName Helvetica Neue
% FontWeight bold
% FontAngle auto
% FontEncoding latin1
% PSLevel 3
% Renderer painters
% Resolution 300
% LineStyleMap none
% ApplyStyle 0
% Bounds loose
% LockAxes on
% LockAxesTicks off
% ShowUI on
% SeparateText off
