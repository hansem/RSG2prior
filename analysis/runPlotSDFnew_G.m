function runPlotSDFnew_G
ccc;
% save NEV, NS3

% note
% difference from runPlotSDF: % new stat: poisson regression

% save data from around 800 trials in 12/7
% data transfer speed warning issue in 12/13? two cells in NS3 data

%%
dirName='/Users/hansem/Documents/Recording/'; cd(dirName);
fname={...
    'G_RSGprior_20170506';...
    'G_RSGprior_20170507';...
    'G_RSGprior_20170508';...
    'G_RSGprior_20170510';... % quite>2nd craniotomy
    'G_RSGprior_20170511';...
    'G_RSGprior_20170512';...
    };

behDir='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/dataMat/';
neuDir='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/neuralData/';

%% relevant cell list [electrode ID, unit ID(12345), start trial, end trial]: liberal criteria for trials
% TBD: check all units
cname={...
    % 1st craniotomy @ 5/6
    % 5/6: t 2080 (c1491); no start/end trial info
% corrected:[75 2  1 2085;2 1 1 2085]%     old:[75 1  1 2085;2 2 1 2085] % wrong note-taking of ch/unit?
    [17 4 75 2085; 21 2 1218 2085; 23 1 1 2085; 27 1 1 2085; 29 2 1205 2085; 29 4 1 2085; 29 5 1 2085; 31 1 1 2085; 31 2 55 2085; 2 1 1 2085; 4 1 1 2085; 4 2 1 2085; 6 1 1 2085; 8 4 1 2085; 8 5 1 2085; 12 1 1 2085; 12 2 1 2085; 12 4 1 2085; 16 1 1 2085; 16 2 1 2085; 16 4 1 2085;...
    33 4 1 2085; 49 1 1 2085; 49 2 1295 2085; 51 1 416 2085; 51 4 426 2085; 51 5 426 2085; 53 1 408 2085; 53 2 1275 2085; 55 4 1 2085; 57 1 1 2085; 59 4 195 2085; 59 5 195 2085; 61 1 180 2085; 61 4 1 2085; 61 5 1 2085; 34 2 1769 2085; 38 2  1 2085; 42 2 1 2085; 44 1  1 2085; 44 2  116 2085; 46 4  1 2085;...
    65 2 1 2085; 67 1  1 2085; 71 5  1 2085; 73 1  1 2085; 73 3  1 2085; 75 2  1 2085; 75 4 1 2085; 77 1 1 2085; 77 5 230 2085; 81 1 1 2085; 83 1 1 2085; 83 2 1 2085; 83 3 456 2085; 87 1 365 2085; 87 3 365 2085; 89 1 1 2085; 91 4 1 2085; 91 5 1 2085; 93 3 1 2085; 93 4 1 2085; 95 4 1 2085; 66 1 1 2085; 68 1 1 2085; 68 2 1 2085; 68 3 625 2085; 72 4 333 2085; 74 2 1 2085; 78 4 1 2085; 78 5 140 2085; 80 4 1 2085];...
    % 5/7: t 2506 (c1907)
    [65 3 1 2506;65 4 1 2506; 67 1 138 2506; 69 2  1 2506; 71 1 1 2506; 71 2 1 2506; 73 4 1 2506; 75 4 1 2506; 77 4 1 2506; 79 1 1 2506; 79 4 1 2506; 81 4 1 2506; 83 1 1 2506; 83 4 1 2506; 85 1 1 2506; 85 4 1 2506; 87 1 1 2506; 87 2 1 2506; 89 1 1 2506; 89 2 2028 2506; 89 4 1 2506;91 1 1 2506;91 2 1 2506;91 4 1 2506;93 1 1 2506;93 4 1 2506;95 1 1 2506;95 4 299 2506;66 4 1 2506;66 5 1 2506;70 2 1 2506;70 4 1 2506;70 5 1 2506;72 4 1 2506;72 1 1198 2506;74 1 1 2506;74 4 1 2506;76 1 529 2506;76 4 1 2506;78 4 1 2506;78 5 1 2506;80 1 1 2506;80 3 1 2506;80 4 1 2506;...
    43 4 1 2506;43 5 1 2506;45 1 1 2506;47 1 1 2506;49 4 1 2506;53 1 371 2506;57 2 1627 2506;59 1 1 2506;59 4 1 2506;61 1 1 2506;61 2 1 2506;63 1 1 2506;63 2 129 2506;63 3 979 2506;34 1 1 2506;34 4 1 2506;36 4 1 2506;38 1 1 2506;38 2 1 2506;40 4 1 2506;40 5 333 2506;42 1 1 2506;42 4 1 2506;42 5 294 2506;44 1 1 2506;44 2 217 2506;46 1 1 2506;46 4 1 2506;48 1 1 2506;48 2 502 2506;48 3 1 2506;48 4 627 2506;...
    1 4 1 2506;3 4 1 2506;5 1 1 2506;5 4 1 2506;7 1 1 2506;7 2 1 2506;9 4 1 2506;11 1 1 2506;11 4 1 2506;13 4 1 2506;15 1 1 2506;17 4 1 2506;19 4 1 2506;21 4 284 2506;23 1 1 2506;23 4 1 2506;25 1 1 2506;25 2 1 2506;25 4 1 2506;27 1 1 2506;27 2 1 2506;29 4 1 2506;2 4 1 2506;4 1 1 2506;4 4 1 2506;6 1 1 2506;6 2 1 2506;6 3 316 2506;8 4 1 2506;8 1 1126 2506;10 1 1 2506;10 4 1 2506;10 5 1 2506;12 3 1 2506;12 4 1 2506;12 5 409 2506;14 4 146 2506;14 5 146 2506;16 1 1 2506;16 2 1 2506;16 4 1 2506];...
    % 5/8: t2388 (1700)
    [1 4 1 2388;3 4 1 2388;5 1 1 2388;7 1 1448 2388;7 4 1 2388;7 5 1 2388;9 1 1 2388;11 1 1 2388;11 4 1 2388;13 1 1 2388;13 4 1 2388;15 4 1 2388;17 4 1922 2388;19 1 1 2388;21 4 1 2388;23 1 1 2388;27 4 1055 2388; 2 4 1 2388;4 4 1 2388;...
    33 1 1927 2388;33 4 1 2388;35 1 1 2388;41 4 286 2388;43 1 860 2388;45 4 1 2388;47 1 1 2388;47 4 1 2388;49 1 1 2388;49 2 1 2388;51 1 1 2388;51 4 1 2388;55 4 277 2388;57 4 308 2388;59 1 1 2388;...
    65 1 1 2388;67 1 1 2388;69 4 1 2388;71 1 1 2388;73 1 1 2388;73 4 1 2388;75 1 125 2388;75 4 1 2388;77 1 1 2388;77 2 313 2388;79 1 1 2388;81 1 1 2388;83 4 1 2388;85 1 397 2388;68 4 1 2388;70 4 1 2388];...
    % 5/10: t2003 (c1687)
    [7 4 1 2003;9 4 1 2003;27 4 1 2003;27 5 1 2003;31 1 1 2003;2 4 1 2003;6 1 1 2003;8 4 1 2003;8 5 1 2003;10 4 1 2003;...
    35 4 1 2003;49 4 1 2003;55 2 1 2003;...
    65 4 1 2003;75 4  1 2003];...
    % 5/11: t4338 (c3104)
    [1 4 100 4338;3 1 1763 4338;3 4 1 4338;5 4 1 4338;5 5 109 4338;7 1 3217 4338;7 4 1 4338;9 5 2443 4338;9 3 1787 4338;9 1 1 4338;9 2 1103 4338;11 4 1 4338;11 5 1 4338;13 1 1 4338;13 4 1 4338;13 5 1 4338;15 1 1 4338;17 4 1 4338;19 4 1 4338;19 1 2802 4338;21 1 1200 4338;21 2 1942 4338;21 4 1 4338;23 1 2428 4338;23 4 1 4338;25 4 1 4338;25 1 1098 4338;27 1 1 4338;27 2 1 4338;29 4 1 4338;31 4 1 4338;31 1 502 4338;2 4 1 4338;4 1 1 4338;4 4 1 4338;6 4 1 4338;6 5 1 4338;8 1 1 4338;8 2 229 4338;8 3 1 4338;8 4 1 4338;10 1 1 4338;10 2 775 4338;10 4 1 4338;12 1 1 4338;12 4 1 4338;14 4 1 4338;16 4 1 4338;...
    49 4 2940 4338;53 4 1 4338;55 4 1 4338;57 4 1 4338;59 4 1 4338;61 4 1 4338;63 4 1 4338;65 4 1 4338;34 1 1 4338;36 4 1 4338;38 4 1 4338;38 5 1 4338;40 1 1 4338;40 4 1 4338;42 1 1 4338;42 2 1 4338;42 4 1 4338;44 4 1 4338;46 4 1 4338;48 4 1 4338;...
    67 1 1277 4338;67 4 1286 4338;69 1 1 4338;71 1 1 4338;71 2 3354 4338;73 1 1 4338;73 4 1 4338;73 2 1111 4338;75 1  1 4338;75 2 1 4338;75 4 1 4338;77 1 1 4338;77 4 1 4338;77 2 1092 4338;79 1 1 4338;79 2 1 4338;79 3 1 4338;79 4 1 4338;81 4 1 4338;83 4 1324 4338;83 5 4077 4338;85 4 1 4338;87 4 1 4338;87 5 2699 4338;89 4 1 4338;89 1 3108 4338;91 5 1983 4338;91 4 1 4338;91 1 1 4338;91 2 1125 4338;93 4 1 4338;93 1 1 4338;95 4 1 4338;95 1 562 4338;66 4 1 4338;66 5 1 4338;68 1 1 4338;68 3 1 4338;68 2 2968 4338;70 3 1 4338;72 4 1 4338;74 4 1 4338;74 2 1119 4338;76 4 1 4338;78 4 1 4338;78 5 1 4338;80 1 1 4338;80 4 1 4338]; % duplicatre 72 4 1 4338>74 4 1 4338
    % 5/12: t3228 (c2349)
    [1 1 1 3228;3 1 1 3228;3 4 1 3228;5 1 1 3228;5 4 55 3228;7 4 368 3228;9 2 556 3228;9 4 1 3228;9 5 1 3228;11 5 1349 3228;11 2 1 3228;11 4 1 3228;13 1 1 3228;13 4 1 3228;15 1 124 3228;17 5 1 3228;17 2 81 3228;17 1 1 3228;17 4 1 3228;19 2 718 3228;21 1 1 3228;21 4 1 3228;23 1 1712 3228;23 4 1 3228;25 1 1 3228;27 4 1 3228;29 4 1 3228;31 4 1 3228;31 5 1 3228;2 4 1 3228;4 4 1 3228;6 5 100 3228;8 4 1 3228;10 4 1 3228;14 4 1 3228;14 5 1 3228;16 4 1 3228;16 5 331 3228;12 4 1 3228;... % duplicate 6 5 100 3228; removed
    33 1 1 3228;37 4 1 3228;37 5 1 3228;39 4 2236 3228;41 1 247 3228;45 4 1 3228;47 4 1 3228;49 4 1 3228;55 4 1748 3228;61 4 836 3228;57 4 436 3228;59 4 212 3228;63 1 1 3228;34 4 1 3228;36 4 1 3228;38 1 1 3228;40 4 1 3228;42 1 1 3228;44 4 1 3228;46 4 1 3228;48 1 1 3228;...
    65 1 1 3228;67 4 1 3228;67 5 850 3228;69 1 1 3228;71 1 1 3228;73 1 1 3228;75 1 1 3228;77 1 1544 3228;77 4 1 3228;77 5 1 3228;79 4 1 3228;79 5 242 3228;81 1 1 3228;81 2 145 3228;81 3 3115 3228;83 1 1 3228;83 4 1 3228;83 5 703 3228;85 1 937 3228;87 1 1 3228;87 4 1 3228;89 1 1 3228;91 1 1739 3228;91 4 247 3228;93 1 1 3228;95 4 3099 3228;66 4 1 3228;68 1 1 3228;70 4 877 3228;72 4 1 3228;72 5 1 3228;74 1 1 3228;76 1 1 3228;76 2 350 3228;78 1 1 3228;78 2 343 3228;80 4 318 3228;80 5 1 3228]
    };
cnamesize=cell2mat(cellfun(@size,cname,'UniformOutput',false));

% if exist([neuDir 'chkSDF.mat'])
%     load([neuDir 'chkSDF.mat']);
% else
    chkSDF=cell(size(cname));
%     % cell2mat(cellfun(@nnz,chkSDF,'UniformOutput',false)): 665
% end
% % cnameMat=cell2mat(cnameValid);

%% another manual sorting for F neuron
% 1) prior effect
%     : 1.Fixation, 2. ts, 3. both
% 2) visual response
%     : 1. target, 2. Flash (ready/set)
% 3) motor
%     : 1. Eye (directional selectivity), 2. Hand (direction), 3. Ramp (for all/hand/eye)
% 4) reward
%     : 1. magnitude

% candName={...
%     [11 1],[],[],[];... % 12/4
%     [4 1],[23 1],[],[23 1;13 1;19 1];...
%     [29 5;53 1;8 1],[11 1],[19 1;51 4],[10 4;55 1;11 1;9 5];... % 12/6
%     [],[],[],[15 4;19 4;6 4];...
%     [29 1;57 1],[57 1],[63 4;29 1;57 1],[75 1;29 1];... % 12/8
%     [67 1;67 2;71 1],[],[76 1; 81 1],[67 1;16 1;16 2];... % 12/10
%     [1 2;3 1;5 1;7 1;9 4;13 1;15 4;17 1;19 1;21 1;21 2;27 1;12 2;79 1;81 1],[11 1],[],[10 4;4 2;31 1];... % 12/11 ******
%     [],[],[],[];... % 12/12
%     [],[39 1],[48 1;36 5;34 5],[];... % 12/14
%     [27 4;19 2;3 1;1 1;8 1;6 5;4 1],[14 4],[],[];... % 12/15
%     [53 4;55 1;57 1;59 1],[],[7 4;9 1;51 1;36 1],[44 4;8 4];... % 12/18
%     [37 1;36 1;38 1;46 1],[36 1;4 4;4 1],[1 4;3 4;4 1;41 1;43 4;45 1],[47 1;8 1;6 4];... % 12/19
%     [1 1;5 1;15 1;23 1;25 1;2 4;8 1;10 4;12 1;16 1;47 4;63 2;46 1;48 2],[27 1],[1 1;13 4;15 1;2 4;41 4;43 4;45 4;47 4;49 1;51 4;55 4;63 2;34 1;40 1],[10 4;8 1];... % 12/20 ******
%     [19 2;21 1;21 5;46 4;44 4;38 4;34 5;45 1;16 2;12 4;12 2;31 1;31 2;2 1;4 1;48 1;46 2;42 4;40 4;36 4;34 4;14 4;14 2;10 2;6 2;8 1;8 2;8 3],[16 4],[21 2;21 5;23 1;31 1;44 4],[45 2;39 1;15 2];... % 12/21 ******
%     [13 2;15 1;17 1;23 1;25 1;29 4;12 2;14 4;16 2;16 3;33 4;35 1;35 2;37 1;39 1;49 4;57 1;61 1;46 2;46 3;48 1],[1 1;13 1;17 1;12 1;14 1;33 4;35 3;53 1;34 4;38 4;],[3 4;15 1;15 4;17 1;29 4;6 1;6 4;12 1;12 3;35 1;35 2;37 1;45 4;47 1;57 1;63 4;40 5],[2 4;34 4];... % 12/22 ******
%     [1 1;3 4;3 2;13 1;21 4;23 3;8 2;12 1;12 2;14 1;37 1;37 4;39 4;47 1;34 4;45 1;45 2],[15 1;23 3;37 4;39 4;41 1],[3 4;13 1;6 5;12 2;14 4;16 1;16 4;41 1;43 4;45 1;42 4],[17 4;29 1;38 4]};... % 12/23 ******

%% 
idPlot=1;
load pplot.mat;
optsExpFig.Height='13'; % 7;
optsExpFig.Width='11';
%         rmfield(optsExpFig,'Width');
%         rmfield(optsExpFig,'Height');
                optsExpFig.Format='png'; % 'png';
%         optsExpFig.Color='cmyk';
% optsExpFig.LockAxes=1;
optsExpFig.LineMode='scaled';
        
nAnal=11; % 10;
nCoeff=67;
Formula=cell(nAnal,1);
CoefficientNames=[]; % =cell(nCoeff,1);
coeff=[];coeffTmp=[];
pval=[];pvalTmp=[];

% new stat: poisson regression
%     '1)    ' prior x eyeHand 22GLM [fix fix+500]'                                  0.24361      0.06015
%     '2)    ' prior x eyeHand 22GLM [tOn-250 tOn]'                                  0.16692     0.021053
%     '3)    ' prior x eyeHand x target 222GLM [tOn tOn+250]'                        0.06015    0.0015038
%     '4)    ' prior x eyeHand x target 222GLM [ready-250 ready]'                   0.046617            0
%     '5) ' prior x eyeHand x target 222GLM   ' prior x eyeHand 22ANOVA [ready ready+480]'                              0.26617     0.082707
%     '6)   ' ts + eyeHand x target regression ' prior x eyeHand 22ANOVA [set-480 set]'                                  0.26917     0.096241
%     '7) prior x eyeHand x target 222GLM [ready+480 ready+800] for overlap ts inc. all long prior?                           0.29624      0.14436       
%     ?8) ' prior x eyeHand x target 222GLM   ' prior x eyeHand 22ANOVA [set-480 set]'
%     '9) intercept tp + eyeHand x directions regression [set set+min(tp)]'      0.8812      0.73835
%     '10) intercept tp + eyeHand x directions regression [production-min(tp) production]';
%     '11) intercept reward regression [rewardOnset onset+100]'                 0.99098      0.95489
%     ' rewarded regression [rewardOnset onset+100]'                             0.1203     0.013534
%     ' rewDur*rewarded regression [rewardOnset onset+100]'                    0.079699    0.0045113
% stat=[];

for i=1:length(fname) % length(fname):(-1):(length(fname)-2)%4:length(fname) %%%%%
    disp(['===== ' fname{i} ' =====']);
    cd([dirName fname{i}]);
    chkSDF{i}=true(size(cname{i},1),1);%     chkSDF{i}=[];
    % MAT file: RSG prior
    fid=fname{i}(end-5:end);
    if isempty(dir([fid '_1.mat'])) & ~isempty(dir([fid 'G_001.NEV'])), % RSG prior % note 'G' added
        openNEV([pwd '/' fid 'G_001.NEV'],'report','read'); % NEV
        movefile([fid 'G_001.mat'],[fid '_1.mat']); % change MAT file name from openNEV
        NS3=openNSx([pwd '/' fid 'G_001.NS3']); % NS3
        save([fid '_1.mat'],'NS3','-append');
    elseif ~isempty(dir([fid '_1.mat'])) % 001 with NIDAQ sync, 1 for photoSync
        load([fid '_1.mat']);
    end
    % extract time info
    beh=load([behDir 'G_20' fid '.mat']);
%     tInfo=extTimeInfo(NS3,beh,idNoPhoto);    save([fid '_001.mat'],'tInfo','-append');
    for j=1:size(cname{i},1)
        if chkSDF{i}(j)==1
            disp(['ch unit t(start) t(end): ' num2str(cname{i}(j,:))]);CoefficientNames=[];
            statTmp=plotSDFexport(NEV,cname{i}(j,:),tInfo,0,idPlot); % p values cell
            
%             for k=1:nAnal
%                 Formula{k}=statTmp{k}.Formula;
%                 CoefficientNames=[CoefficientNames(:); statTmp{k}.CoefficientNames(:)];
%                 coeffTmp=[coeffTmp table2array(statTmp{k}.Coefficients(:,1))']; % 1 x nCoeff
%                 pvalTmp=[pvalTmp table2array(statTmp{k}.Coefficients(:,end))'];
%             end
%             coeff=[coeff;coeffTmp];coeffTmp=[];
%             pval=[pval;pvalTmp ];pvalTmp=[];
% %             stat=[stat cell2mat(statTmp)];

        end
        if idPlot
        set(gcf,'PaperPositionMode','auto');
        saveas(gcf,[neuDir fid '_' num2str(cname{i}(j,1)) '_' num2str(cname{i}(j,2)) '.png']);
        end
%         keydown=waitforbuttonpress; %chkSDF{i}=[chkSDF{i}; keydown];
%         close all;
% exportfig(gcf,[neuDir fid '_' num2str(cname{i}(j,1)) '_' num2str(cname{i}(j,2)) '.png'],'color','rgb','Height',13,'Width',11,'format','png');

close all;
%         save([neuDir 'chkSDF.mat'],'chkSDF');
    end
%     statName={'1) prior 22ANOVA [fix fix+500]';
%     ' eyeHand 22ANOVA [fix fix+500]';
%     ' prior x eyeHand 22ANOVA [fix fix+500]';
%     '2) prior 22ANOVA [tOn-250 tOn]';
%     ' eyeHand 22ANOVA [tOn-250 tOn]';
%     ' prior x eyeHand 22ANOVA [tOn-250 tOn]';
%     '3) prior 222ANOVA [tOn tOn+250]';
%     ' eyeHand 222ANOVA [tOn tOn+250]';
%     ' target 222ANOVA [tOn tOn+250]';
%     ' prior x eyeHand 222ANOVA [tOn tOn+250]';
%     ' prior x target 222ANOVA [tOn tOn+250]';
%     ' eyeHand x target 222ANOVA [tOn tOn+250]';
%     ' prior x eyeHand x target 222ANOVA [tOn tOn+250]';
%     '4) prior 222ANOVA [ready-250 ready]';
%     ' eyeHand 222ANOVA [ready-250 ready]';
%     ' target 222ANOVA [ready-250 ready]';
%     ' prior x eyeHand 222ANOVA [ready-250 ready]';
%     ' prior x target 222ANOVA [ready-250 ready]';
%     ' eyeHand x target 222ANOVA [ready-250 ready]';
%     ' prior x eyeHand x target 222ANOVA [ready-250 ready]';
%     '5) prior 22ANOVA [ready ready+400]';
%     ' eyeHand 22ANOVA [ready ready+400]';
%     ' prior x eyeHand 22ANOVA [ready ready+400]';
%     '6) prior 22ANOVA [set-400 set]';
%     ' eyeHand 22ANOVA [set-400 set]';
%     ' prior x eyeHand 22ANOVA [set-400 set]';
%     '7) prior 22ANOVA [set-320 set] for overlap ts';
%     ' eyeHand 22ANOVA [set-320 set] for overlap ts';
%     ' prior x eyeHand 22ANOVA [set-320 set] for overlap ts';
%     '8) intercept tp + eyeHand x directions regression [set set+min(tp)]';
%     ' tp regression [set set+min(tp)]';
%     ' eyeHand regression [set set+min(tp)]';
%     ' directions regression [set set+min(tp)]';
%     ' tp x eyeHand regression [set set+min(tp)]';
%     ' tp x directions regression [set set+min(tp)]';
%     ' eyeHand x directions regression [set set+min(tp)]';
%     ' tp x eyeHand x directions regression [set set+min(tp)]';
%     '9) intercept tp + eyeHand x directions regression [production-min(tp) production]';
%     ' tp regression [production-min(tp) production]';
%     ' eyeHand regression [production-min(tp) production]';
%     ' directions regression [production-min(tp) production]';
%     ' tp x eyeHand regression [production-min(tp) production]';
%     ' tp x directions regression [production-min(tp) production]';
%     ' eyeHand x directions regression [production-min(tp) production]';
%     ' tp x eyeHand x directions regression [production-min(tp) production]';
%     '10) intercept reward regression [rewardOnset onset+100]';
%     ' rewarded regression [rewardOnset onset+100]';
%     ' rewDur*rewarded regression [rewardOnset onset+100]';};
%     save([neuDir 'stat.mat'],'stat','statName'); 
    
    
%     Generalized Linear regression model:
%     log(numSpikes) ~ 1 + idShort*idEye
%     Distribution = Poisson
% 
% Estimated Coefficients:
%                                Estimate      SE        tStat       pValue  
%                                ________    _______    _______    __________
% 
%     (Intercept)                 -1.2879      0.125    -10.303    6.8415e-25
%     idShort_Long                 0.5025    0.16428     3.0588     0.0022225
%     idEye_Hand                  0.37914     0.1583     2.3951      0.016617
%     idShort_Long:idEye_Hand    -0.35537    0.21156    -1.6798      0.092998
% 
% 
% 945 observations, 941 error degrees of freedom
% Dispersion: 1
% Chi^2-statistic vs. constant model: 14.2, p-value = 0.00262
    
    
    
    % MAT file: overlap
    if isempty(dir([fid '_002.mat'])) & ~isempty(dir([fid '_002.NEV'])), % overlap
        openNEV([pwd '/' fid '_002.NEV'],'report','read'); % NEV
        NS3=openNSx([pwd '/' fid '_002.NS3']); % NS3
        save([fid '_002.mat'],'NS3','-append');
    elseif ~isempty(dir([fid '_002.mat']))
        load([fid '_002.mat']);
    end
    
end % for i=1:length(fname)
iAnalPerCoeff=[];n=[];for i=1:length(statTmp),n=n+statTmp{i}.NumCoefficients;iAnalPerCoeff=[iAnalPerCoeff;repmat(i,statTmp{i}.NumCoefficients,1)];end
%  save([neuDir 'singleNeuronStat_smallWindow.mat'],'CoefficientNames','Formula','coeff','pval','nCoeff','nAnal','iAnalPerCoeff','statTmp'); 
 
%% checking single-neuron stats
% % index for each stat
% % iAnalPerCoeff1=5; iCoeff=24;CoefficientNames1='idShort_Long';legName1='prior(early t_s)';
% iAnalPerCoeff2=9; jCoeff=52;CoefficientNames2='t_p';legName2='production(early t_p)';
% iAnalPerCoeff1=6; iCoeff=31;CoefficientNames1='t_s';legName1='measure(late t_s)';
% % iAnalPerCoeff2=6; jCoeff=31;CoefficientNames2='t_s';legName2='measure(late t_s)';
% 
% thP=0.05;
% h=figure; load pplot.mat;kColor=1;
% legName={legName1;legName2};
% for i=1:nCoeff
%     if i==1
%         disp(['===== ']);
%         disp(Formula{iAnalPerCoeff(i)});
%     else
%         if iAnalPerCoeff(i)-iAnalPerCoeff(i-1)==1
%             disp(['===== ']);
%             disp(Formula{iAnalPerCoeff(i)});
%         end
%     end
%     pSig=mean(pval(:,i)<thP)*100;    
%     disp([CoefficientNames{i} '         :           ' num2str(pSig,3)]);
%     % hist for regression coefficients; only main effect
%     if (iAnalPerCoeff(i)==iAnalPerCoeff1 & strcmp(CoefficientNames{i},CoefficientNames1)) ||...
%             (iAnalPerCoeff(i)==iAnalPerCoeff2 & strcmp(CoefficientNames{i},CoefficientNames2))
%         subplot(2,1,kColor);
%         htmp=histStairs(coeff(:,i),40,0,h,pplot.cmap{kColor});
%         axis tight; ylabel('# neurons');legend(legName{kColor},'location','best'); legend boxoff;plotVertical(gca);
%         kColor=kColor+1;
%     end
% end
%  xlabel('regression coefficient');
% % scatter plot
% idPriorEarlyTs=(iAnalPerCoeff==iAnalPerCoeff1 & strcmp(CoefficientNames,CoefficientNames1));
% idTsLateTs=(iAnalPerCoeff==iAnalPerCoeff2 & strcmp(CoefficientNames,CoefficientNames2));
% figure; set(gcf,'position',pplot.rect2_3);
% tmpX=coeff(:,idPriorEarlyTs);
% tmpY=coeff(:,idTsLateTs);
% plot(tmpX,tmpY,'k.'); axis tight; hold on; plotVertical(gca); plotHorizon(gca);plotIdentity(gca); 
% gmfit{1}=fitgmdist(coeff(:,[iCoeff jCoeff]),1);
% % multiple cluster?
% gmfit{2}=fitgmdist(coeff(:,[iCoeff jCoeff]),2);
% % ops=statset('MaxIter',1000);
% gmfit{3}=fitgmdist(coeff(:,[iCoeff jCoeff]),3); % ,'options',ops);
% disp(num2str([gmfit{1}.AIC gmfit{2}.AIC gmfit{3}.AIC;...
%     gmfit{1}.BIC gmfit{2}.BIC gmfit{3}.BIC]));
% disp('');
% idMinBic=min([gmfit{1}.BIC gmfit{2}.BIC gmfit{3}.BIC])==[gmfit{1}.BIC gmfit{2}.BIC gmfit{3}.BIC];
% ezcontour(@(x,y)pdf(gmfit{idMinBic},[x y]),[mean(tmpX)-[1 -1]*3*std(tmpX)],[mean(tmpY)-[1 -1]*3*std(tmpY)]);
% axis tight;title([]);
% xlabel(['\beta ' legName1]); ylabel(['\beta ' legName2]);
% 
% acosd((tmpX./norm(tmpX))'*(tmpY./norm(tmpY)))
% 
% % chi2 test for independence
% contMat=[nnz(pval(:,iCoeff)>0.05&pval(:,jCoeff)>0.05) nnz(pval(:,iCoeff)<0.05&pval(:,jCoeff)>0.05);...
% nnz(pval(:,iCoeff)>0.05&pval(:,jCoeff)<0.05) nnz(pval(:,iCoeff)<0.05&pval(:,jCoeff)<0.05)];
% 
% [chi,pChi]=chi2testIndep(contMat)

%% representative cell list; check STAT

% tmp=[161218 57 1;...
% 161221 4 1;...
% 161211 7 1;...
% 161220 16 1;...
% 161222 25 1;...
% 161214 39 1;...
% 161223 23 3;...
% 161222 33 4;...
% 161208 57 1;...
% 161206 51 4;...
% 161222 15 1;...
% 161223 41 1;...
% 161222 45 4;...
% 161220 47 4;...
% 161220 15 1;...
% 161210 67 1];
% 
% iMat=[];
% for i=1:length(tmp)
%     id=find(fidMat(:,1)==tmp(i,1) &...
%         fidMat(:,2)==tmp(i,2) &...
%         fidMat(:,3)==tmp(i,3));
%     iMat=[iMat; id];
%     table(statName,stat(:,id))
%     figure; waitforbuttonpress; close;
% end
% 
% % num2str(stat(:,iMat))
% 
% table(statName,round(mean(stat<0.05,2)*100),round(mean(stat<0.001,2)*100),'variablenames',{'stat','pCellp05','pCellp001'})

%% best sessions for single-trial analysis
%% pick up common trials across cells
% listSess=[161221;...
%       161220;...
%       161222;...
%       161223;...
%       161211];
%   load pplot.mat;
%   for i=1:length(listSess)
%       idSess=fidMat(:,1)==listSess(i);
%       nC=nnz(idSess);
%       nT=max(cnameMat(fidMat(:,1)==listSess(i),end));
%       tmp=zeros(nC,nT);
%       for j=1:nC
%           startT=cnameValid{strcmp(fname,['H_RSGprior_20' num2str(listSess(i))])}(j,end-1);
%           endT=cnameValid{strcmp(fname,['H_RSGprior_20' num2str(listSess(i))])}(j,end);
%           tmp(j,startT:endT)=1;
%       end
%       figure; set(gcf,'position',pplot.(['rect1_' num2str(i)]));
%       imagesc(tmp); colormap('gray');
%       xlabel('trials'); ylabel('cell');
%       
%       edges=[1:nT]-0.5;
%       nStartT=histc(fidMat(idSess,end-1),edges);
%       nTlength=histc(fidMat(idSess,end)-fidMat(idSess,end-1),edges);
%       nEndT=histc(fidMat(idSess,end),edges);
%       figure; set(gcf,'position',pplot.(['rect2_' num2str(i)]));
%       plot(1:nT,cumsum(nStartT),'r',1:nT,cumsum(nTlength),'g',1:nT,cumsum(nEndT),'b'); 
%       axis tight; xlabel('trials'); ylabel('cumulative # cells');
% %       figure; plot([1:nT]'.*(max(cumsum(nTlength))-cumsum(nTlength)));
%   end
