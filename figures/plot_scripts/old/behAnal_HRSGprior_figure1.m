function behAnal_HRSGprior_figure1

% save NEV, NS3

% note
% save data from around 800 trials in 12/7
% data transfer speed warning issue in 12/13? two cells in NS3 data

% sumHbehRecord

% 2019/1/10: model comparison with linear model fit

%%
try
    cd('/Users/hansem/Dropbox (MIT)/figuresRSG2prior/1beh'); % SFN17/1beh'); % CCN17/beh');
    behDir='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/dataMat/';
    neuDir='/Users/hansem/Dropbox (MIT)/fileFromNHPrig/neuralData/';
    psthDir='/Users/hansem/Dropbox (MIT)/psthDataHigh/';
    dirName='/Users/hansem/Documents/Recording/'; % cd(dirName);

catch
    cd('/Users/seonminahn/Dropbox (MIT)/figuresRSG2prior/1beh'); % CCN17/beh');
    behDir='/Users/seonminahn/Dropbox (MIT)/fileFromNHPrig/dataMat/';
    neuDir='/Users/seonminahn/Dropbox (MIT)/fileFromNHPrig/neuralData/';
    psthDir='/Users/seonminahn/Dropbox (MIT)/psthDataHigh/';
    dirName='/Users/seonminahn/Documents/Recording/'; % cd(dirName);
end


fname={...
%     'H_RSGprior_20161202';... % only eye; use 003
    'H_RSGprior_20161203';... % no photosync; RSG>overlap>RSG; use 003
    'H_RSGprior_20161204';... % no photosync
    'H_RSGprior_20161205';... % no photosync; no overlap
    'H_RSGprior_20161206';...
    'H_RSGprior_20161207';...
    'H_RSGprior_20161208';... 
    % 2nd craniotomy @12/9 afterward, [5 3 20]
    'H_RSGprior_20161210';...
    'H_RSGprior_20161211';...
    'H_RSGprior_20161212';...
%     'H_RSGprior_20161213';...
    'H_RSGprior_20161214';...
    'H_RSGprior_20161215';... 
    % 3rd craniotomy @12/16
    'H_RSGprior_20161218';...
    'H_RSGprior_20161219';...
    'H_RSGprior_20161220';...
    'H_RSGprior_20161221';...
    'H_RSGprior_20161222';... % saved as 12/21
    'H_RSGprior_20161223'
    };
nS=length(fname);

idExportFig=0;

load pplot.mat;
%% relevant cell list [electrode ID, unit ID(12345), start trial, end trial]: liberal criteria for trials
% TBD: check all units
cname={...
%     [29 1 1 1658; 21 1 562 1230; 21 2 562 1230];... % %     'H_RSGprior_20161202';... % only eye; use 003
    [15 1 1 1297; 8 1 1 1297; 9 1 1 1297; 21 1 599 1297; 27 1 616 1297; 25 1 600 1297];... %     'H_RSGprior_20161203';... % no photosync; RSG>overlap>RSG; use 003
    [3 2 1 1830; 11 1 309 1628; 9 3 857 1830;29 2 804 1830;5 3 1116 1830;31 1 1 576;10 2 1490 1830;11 2 857 1830;12 3 1294 1830; 23 1 1294 1830];... %     'H_RSGprior_20161204';... % no photosync
    [4 1 208 1315; 13 4 1 2238; 17 4 572 1483;23 1 1280 2238;8 1 1 1806; 8 4 1 1806;10 1 1 1505;6 1 714 1236;4 2 902 1315;11 4 1 2238;13 1 1 2238;14 1 1 1117;16 1 1 1000;21 2 1785 2238;19 1 1785 2238;];... %     'H_RSGprior_20161205';... % no photosync; no overlap
    [6 1 1 2797;15 1 1 2029;19 1 1 2797;29 4 1 699; 29 5 1 1904;4 1 1 2797; 53 1 1 750;51 4 1 1698;4 2 439 2739;44 1 1 1309;45 4 1 2472;1 1 535 2797;29 5 579 2797;10 4 1 1725;55 1 660 1553;8 1 1 2050;11 1 830 2383; 14 2 830 1748; 19 2 661 2439; 21 1 868 1587;48 2 923 1594; 42 1 1459 2058;48 1 1637 2797;12 1 1 2394;25 1 1776 2328;2 1 1841 2575;48 3 2224 2797;9 5 1 2773;16 4 2420 2797;31 4 1 2797;57 1 1 2541];... %     'H_RSGprior_20161206';...
    [11 1 1 1993-800;15 4 1 1711;19 4 1 1711;21 4 1 1711;21 5 1 1711;23 1 1 1711;23 2 1 1711;25 4 1 1711;4 1 1 2346-800;4 4 1 1711;6 4 1 1711;14 1 1 1711;41 2 1 1711;43 1 1 1711;43 4 1 1711;43 2 1523-800 1711;45 4 926-800 1711;47 1 1 1711;47 2 1 1711;47 4 1 1711;49 4 1 1711;53 4 1 1711;55 1 1 1711;59 1 154 1711; 61 1 1 1711;63 1 1 1711;34 1 1 1711;36 4 1 1711;44 4 1 1711;48 4 1 1711];... %     'H_RSGprior_20161207';...
    [29 1 1 1711;31 1 1 1711;34 4 1 1711;34 5 1 1711;55 4 1 1711;55 5 1 1711;57 1 1 1711;57 2 1 1711;59 1 1 1711;59 2 1 1711;63 1 1 1711;63 4 1 1711;66 1 1 1711;68 1 302 1711;68 2 1425 1711;75 1 1 1711;81 1 1 1711;81 2 1 1711;81 3 1 1711;83 4 415 1711;89 1 1 1711];... %     'H_RSGprior_20161208';... 
    %     % 2nd craniotomy @12/9 afterward, [5 3 20]
    [5 4 198 1546;11 4 1 1546;11 5 1 1546;16 1 1 1546;16 2 1 1546;33 4 180 1546;33 5 180 1546;67 1 1 1420;67 2 1 1420;71 1 1 1546;73 1 1 1430;73 2 1 1430;76 1 1178 1546;79 1 1 1546;79 2 1 1546;81 1 1 1420;66 1 1 1420];... %     'H_RSGprior_20161210';...
    [1 1 1 1774;1 2 1 1774;3 1 1 1774;5 1 1 1774;5 3 1441 1774;7 1 1 1774;7 2 1 1774;9 4 1 1774;13 1 105 1774;15 4 1 1774;17 1 1 1774;19 1 1 1774;21 1 1 1774;21 2 1 1774;23 1 1 1774;27 1 1 1774;29 1 1 1774;31 1 1 1774;2 4 224 1774;4 1 1 1774; 4 2 1 1774;6 1 1 814;6 2 1 1774;8 4 1 1774;10 4 1 1774;12 2 1 1774;12 3 961 1774;14 1 1 1745; 16 1 1 854;16 4 1 1774;51 1 1 1774;53 1 1 1774;53 2 1152 1774;55 1 1 320;55 2 1 320;55 3 320 904;57 1 375 910;61 1 334 1526;61 4 1 1774;44 1 1 1774;68 1 1 1541;69 4 1 1774;69 5 1 1774;79 1 1 1774;81 1 1 1774;85 4 1 1774;85 5 1 1774;87 4 1 1774;87 5 1 1774;89 1 1 394;93 1 1 1774;95 1 1 1774;66 4 1 1774;66 5 1 1774;70 1 1 1774;72 1 1 794];... %     'H_RSGprior_20161211';...
    [7 1 1 1524;21 4 1 1524;21 5 1 1524;10 4 1 1524;33 4 1 1524;33 5 1 1524;35 4 1 1524;37 4 1 1524;37 5 1 1524;67 1 1 1367;75 4 1 1524;75 5 1 1524;79 4 1 1524;79 5 1 1524;87 4 1 1524;87 5 1 1524];... %     'H_RSGprior_20161212';...
%     [4 1 1 1287;8 4 1 1847;10 1 1 1847;43 4 1 1847;37 1 514 1847;38 1 1 730;38 2 1 1847;38 3 1 1847;40 1 1 1847;42 1 1 1847;42 2 245 1847;44 1 1 1847;48 1 1 1847;48 2 1 1847;48 3 1 1847;48 4 1 1847;85 4 1 1847;87 4 1 1847;87 5 1 1847];... % %     'H_RSGprior_20161213';...
    [48 4 1 210;46 4 221 1761;46 2 229 1761;44 1 116 1761;44 2 249 1761;44 3 249 1761;38 1 82 1761;34 1 178 1761;36 1 1 1761;47 1 338 1706;42 1 1 1761;42 4 1 1761;40 1 1 1761;37 1 1 1761;41 1 310 1761;41 2 371 1761;41 3 391 1761;48 1 1 1761;48 2 1 1761;48 3 526 1761;48 4 17 1761;46 2 229 1761;46 3 1 1761;46 4 221 1761;44 1 1 825;44 1 825 1761;44 2 1 825;44 2 825 1761;44 3 1 1761;44 4 1 1761;44 5 1 1761;42 1 1 1761;42 4 1 1761;42 5 1 1761;40 1 1 1761;40 2 549 1761;40 3 766 1761;38 4 1 1761;38 5 1 1761;38 1 1 1761;36 1 1 1222;36 2 707 1761;36 4 1 1761;36 5 1 1761;34 1 178 1761;34 4 1 1761;34 5 476 1761;34 2 996 1761;63 1 485 1761;57 1 627 1552;53 1 1 1538;49 4 1 1761;47 1 1 1495;47 4 1 1761;45 1 1 1761;45 3 1 1761;43 1 1 1761;41 1 1 1761;39 1 1 1761;37 4 1 1761;37 5 1 1761;35 1 1 1761;35 4 1 1761;33 4 581 1761;61 2 1016 1561;55 1 1422 1761;55 2 1422 1761];... %     'H_RSGprior_20161214';...
    [31 1 250 1597;31 3 1 1597;29 4 1 1597;29 5 1 1597;27 4 1 1597;27 5 262 1597;25 4 1 1193;25 1 1193 1466;25 2 1208 1597;23 4 1 1597;21 1 1 1597;21 4 1 1597;19 2 98 425;17 1 325 1597;15 1 1 1597;15 2 430 1597;13 4 436 1597;9 4 1 1597;7 1 1 1597;5 2 1 1597;5 4 1 1597;3 1 1 1597;1 2 1 475;1 1 475 1438;1 3 475 1438;16 2 1 1597;16 4 1 1597;14 2 1 1597;14 3 743 1319;14 4 1 1597;12 2 1 1597;12 3 1 1597;12 4 214 1597;10 4 1 1597;10 5 1 1597;8 1 1 1597;8 2 1 1597;6 1 185 1597;6 5 970 1597;6 4 1 970;4 1 1 1597;4 4 1 1597;4 5 1 1597;2 1 1 1597];... %     'H_RSGprior_20161215';... 
    % 3rd craniotomy @12/16
    [1 1 1 2097;5 4 1 2097;7 4 1 2097;9 1 1 2097;11 1 1 2097;13 1 1 2097;15 1 1 2097;15 2 1 314;15 4 21 2097;17 1 1268 2097;17 4 1 2097;17 5 1 2097;19 1 60 2097;21 1 1 2097;21 4 1 2097;23 4 1 2097;25 1 1 2097;25 2 1 2097;27 1 1 2097;27 2 1 2097;29 1 1 2097;29 2 1 2097;29 4 1 2097;31 4 118 2097;31 1 1 1192;4 4 1 369;4 1 369 985;12 1 1 2097;12 2 1 2097;14 1 632 2097;14 4 1 2097;14 5 1 2097;16 1 519 2097;16 2 1 2097;16 3 50 2097;2 4 1 2097;8 4 1 2097;8 1 1652 2097;10 1 1 2097;10 2 1 2097;10 3 1 2097;6 4 1670 2097;6 5 1670 2097;47 1 1 2097;47 2 1 1951;49 4 1 2097;51 1 1 2097;53 4 444 2097;55 1 1 2097;55 4 1 2097;57 1 1 2097;59 1 1 2097;61 1 1 2097;63 4 1 2097;34 1 146 2097;34 4 1 2097;36 1 164 2097;36 4 1 2097;38 1 238 2097;38 4 1 2097;42 1 249 2097;42 4 249 2097;44 1 1 2097;44 4 270 2097;44 2 672 2097;48 4 1 2097;48 5 1 2097];... %     'H_RSGprior_20161218';...
    [1 4 1 1591;3 4 108 1591;13 1 117 1591;15 4 1 1591;15 4 1 1591;17 4 148 1591;19 4 1 1591;21 1 1 1591;23 2 1386 1591;23 1 1 1591;23 4 794 1591;25 1 164 1591;27 1 1 1591;27 2 237 1591;29 4 825 1591;31 1 1 1591;2 4 1 1591;4 4 180 1591;4 5 192 643; 4 1 643 1065;4 2 643 1315;4 3 643 1591;6 4 1 1591;8 1 1 1591;10 1 360 1591;12 4 1 1591;12 5 1 1591;14 1 1 1355;14 2 1355 1591;14 4 1 1591;16 1 368 1591;16 4 1 1591;35 1 1 1591;37 1 1 1591;39 1 1 1591;41 1 1 1591;43 4 1 1591;45 1 92 1591;45 4 382 1591;47 1 1 1591;49 4 300 1197;34 4 877 1591;36 1 1 1591;36 2 1 1591;38 1 1 1591;40 1 1 1591;40 4 400 1591;44 4 1 1591;46 1 211 1591;48 4 1 1591;48 5 250 1591;];... %     'H_RSGprior_20161219';...
    [1 1 1 1312;1 2 677 1312;3 1 82 1304;3 4 1 1791;5 1 1295 1791;7 4 1 1791;9 4 1 1791;11 4 1 1791;13 4 1 1791;15 1 517 1791;15 4 1 1791;17 4 470 1791;17 5 470 1791;19 4 1 1791;23 1 1481 1791;25 1 703 1230;25 2 795 1230;25 4 1 1791;27 1 30 1791;29 4 1 1791;29 5 1215 1791;31 1 459 1204;4 1 1 1170;4 4 1 1791;4 5 1 1791;2 1 1 1791;2 2 141 1191;2 4 1 1791;6 1 1 1791;6 2 1 1791;6 3 780 1791;6 4 1327 1791;8 4 309 1103;8 5 320 1791;8 1 1 1791;10 4 1 1791;10 5 1 1791;12 1 1 1791;12 2 1 1791; 12 3 565 1791;12 4 1139 1791;12 5 1327 1791;14 1 1 1791;14 2 1 1791;14 3 530 1791;16 1 481 878;16 2 1 878;16 4 1 1791;14 5 638 1791;16 3 651 1791;10 1 1504 1791;8 2 1580 1791;8 2 1580 1791;39 4 1 1791;41 4 1 1791;43 4 1 1791;45 4 1 1791;47 4 1 1791;49 1 1 1462;49 4 1 1791;51 4 1 1565;53 4 1 1791;55 4 253 1477;57 1 352 1477; 57 4 1 1791;59 4 367 1714;61 1 1 602;61 4 1 1597;61 5 1 1597;63 1 1 1597;63 2 267 1597;34 1 1 1608;34 2 1 1608;36 4 1 1791;36 5 1 1791;40 1 219 1791;40 4 1 1791;38 1 1 1791;42 4 1 1791;42 5 1 1791;44 4 1 1791;44 5 228 1067;46 4 1 1791;46 2 932 1791;46 5 1 1791;46 1 913 1791;48 1 1 1791;48 2 1 1791;48 4 1 1791];... %     'H_RSGprior_20161220';...
    [13 5 1 2139;13 4 345 2139;15 3 379 1229;15 1 379 1229;15 2 379 1229;17 4 1 2139;19 4 1 2139;19 1 1308 2139;19 2 1308 2139;21 1 1 2139;21 2 138 1411;21 3 1411 2139;21 4 1 2139;21 5 1 2139;23 1 1 1319;23 2 1 1319;23 4 1 1319;25 4 1 2139;27 1 1 2139;29 1 1 2139;31 4 1 2139;31 1 156 2139;31 2 156 2139;2 1 1 2139;2 2 104 2139;21 4 1 2139;4 1 1 1654;4 2 1654 2139;6 1 1 2139;6 2 1 2139;8 1 1 2139;8 2 1 2139;8 3 1 2139;10 1 1 2139;10 2 1 2139;12 1 1 2139;12 2 1 2139;12 4 1 2139;14 2 1 2139;14 3 1 2139;14 4 1 2139;16 1 1 2139;16 2 412 2139;16 4 1 2139;33 4 1 2139;35 4 1 2139;37 4 1 2139;39 1 1339 2139;39 4 1 2139;39 5 1 2139;41 1 1 2139;41 2 1 2139;43 1 91 2139;43 4 1 2139;45 3 1877 2139;45 1 1 1733; 45 4 1 2139;45 2 1746 2139;47 1 1 2139;49 4 951 2139;51 1 1 2139;55 4 1 2139;57 4 1 2139;34 4 1 2139;34 5 1 2139;36 4 1 2139;38 4 1 2139;40 4 130 2139;42 4 1 2139;44 1 841 2139;44 4 1 2139;44 5 1 2139;46 1 1 2139;46 2 1600 2139;46 4 1 2139;48 2 1 2139;48 1 1050 2139;48 3 1050 2139;48 4 1 2139;];... %     'H_RSGprior_20161221';...
    [1 1 1 2158;1 2 1 2158;3 4 1 2158;13 2 139 2158;13 1 1 2158;15 1 1 2158;15 2 1839 2158;15 4 1 2158; 17 1 1 2158;21 1 1 2158;21 4 1 2158;21 5 1 2158;23 1 1 2158;23 2 1986 2158;25 1 1 2158;25 4 1 2158;27 1 617 2158;27 4 1 2158;29 4 1 2158;31 1 1 2158;31 2 1 2158;31 4 1 2158;2 1 1 2158;2 2 92  2158;2 4 1 2158;4 1 1 1706;4 4 1 2158;6 1 1 2079;6 4 1 2158;8 4 1 2158;10 1 1 2158;10 4 1 2158;12 1 1 2158;12 2 1 2158;12 3 1720 2158;14 1 1 2158;14 2 1 2158;14 4 1 2158;16 2 1 1328;16 1 1328 2158;16 3 1 2158;33 4 1 2158;35 1 1 2158;35 2 1 2158;35 3 1 2158;37 1 1 2158;39 1 1 2158;39 2 1 2158;41 1 1 2158;43 1 1 2158;45 4 196 2158;47 1 1 2158;49 4 1 2158;51 1 1 2158;51 4 1 2158;53 1 233 2158;55 4 246 2158;57 1 1 2158;61 1 1 2158;61 4 1 2158;63 4 258 2158;34 4 867 2158;36 1 1 2158;36 4 1 2158;38 4 1 2158;40 4 1 2158;40 5 1 2158;42 1 1 2158;44 1 1 2158;44 2 650 2158;46 1 209 1177;46 2 202 2158;46 4 1 2158;48 4 1 2158; 48 1 1191 2158;];... %     'H_RSGprior_20161222';... % saved as 12/21
    [1 1 1 2477;3 1 1 2477;3 4 1 2477;3 2 1588 2477;3 3 2212 2477;13 1 779 2477;13 4 1 2477;15 1 1 2477;15 4 1 2477;17 4 1763 2477;19 1 1 2477;19 4 1 2477;21 4 1 2477;21 5 630 2477;23 1 1 2477;23 4 54 2477;23 2 1558 2477;23 3 1558 2477;25 1 1 1267;25 2 1 1267;27 4 1 2477;29 4 210 2477;29 1 1 2477;29 2 1 2477;31 1 246 2477;31 4 1 2477;2 1 1 2477;2 4 1 2477;4 1 1 758;4 4 1 2477;4 5 1 2477;6 4 367 2477;6 5 1 2477;6 2 649 2477;6 3 678 2435;8 1 1 2477;8 2 967 2477;10 1 702 2477;10 2 702 1305;10 4 1 2477;12 1 1 2477;12 2 1 2477;14 1 1 2477;14 2 1 2477;14 3 166 2477;14 4 1 2477;14 5 938 2477;16 1 1 2477;16 4 823 2477;33 4 1 2477;35 4 1 2477;37 1 319 2477;37 4 1 2477;39 4 1 2477;41 1 1 2477;41 2 1 2477;43 4 1 2477;45 1 1 2477;45 2 842 2477;47 4 1 2477;51 1 1 2477;53 1 1 2477;55 1 1 2477;59 1 400 2477;59 2 1095 2477;61 4 1 2477;63 4 420 2477;34 1 1 2477;34 2 54 2477;34 3 253 2477;34 4 1 2477;36 1 1 2477;36 2 448 2477;36 4 1 2477;38 1 1 2477;38 4 1 2477;40 1 1492 2477;40 2 1492 2477;42 4 1 2477;42 5 345 2477;44 1 1 2477;44 2 1 2477;46 1 1 2477;46 4 1 2477;48 1 1 2477;48 2 1 2477;48 4 1 2477;];... %     'H_RSGprior_20161223'
    };
cnamesize=cell2mat(cellfun(@size,cname,'UniformOutput',false));
% 
% if exist([neuDir 'chkSDF.mat'])
%     load([neuDir 'chkSDF.mat']);
% else
%     chkSDF=cell(size(cname));
%     % cell2mat(cellfun(@nnz,chkSDF,'UniformOutput',false)): 665
% end

%% check across-session variability of wm,wp,offset (wFit w/o outlier removal?)
% figure; nMeas=4; 
% for iS=1:length(fname)
%     load([behDir 'H_' fname{iS}(end-7:end) '.mat']); % T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)
%     subplot(nMeas,1,1);plot(iS,wFit.w_m,'o'); ylabel('w_m');hold all;
%     subplot(nMeas,1,2);plot(iS,wFit.w_p,'o'); ylabel('w_p');hold all;
%     subplot(nMeas,1,3);plot(iS,wFit.offset1,'o'); ylabel('offset1');hold all;
%     subplot(nMeas,1,4);plot(iS,wFit.offset2,'o'); ylabel('offset2');hold all;
% end

%% define D data structure
nPr=2;    prNm={'Short','Long'};
nEH=2;ehNm={'Eye','Hand'};
nTarg=2; targNm={'Right','Left'};
nTspp=5; Tmat=[];Tmat{1}=linspace(480,800,nTspp); Tmat{2}=linspace(800,1200,nTspp); Tall=[Tmat{1} Tmat{2}];
nTp=5; 

%% outlier removal for each session
% also prior-, modality-, direction-, ts-specific

% idOut=[]; rangeOut=[0 3]; pOut=[];
% idPlot=0; % 1;
% for iS=1:nS
%     fid=fname{iS}(end-5:end);
%     if iS==1, fnmTmp=[behDir 'H_20' fid '_3.mat']; else  fnmTmp=[behDir 'H_20' fid '.mat']; end;
%     load(fnmTmp);% T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)
%     idOutS=true(size(T)); % assumig all outliers (inc. t<0)
%     for i=1:nPr
%         for j=1:nEH
%             for k=1:nTarg
%                 for m=1:length(Tmat{i})                
%                     id=idShortTrial==(2-i) &... % idShortTrial
%                         idHandEye==(2-j) &...  % idHandEye
%                         theta==(k-1)*180 &... % target location
%                         T==Tmat{i}(m) &...
%                         t>0 & t<3*Tmat{i}(m); % only valid trials for now
%                     u0=rangeOut*Tmat{i}(m); %[0.1 1.9]*overlapT; % [0.5 1.5]*overlapT; %[0 2]*overlapT; % with least outlier removal; or max(t)
%                     [tClean,idClean,pOutTmp,tOut]=removeOutlier(t(id),0,u0); %nSDremoveOutlier);
%                     idOutS(id)=~idClean;
%                     pOut=[pOut;mean(~idClean)];
%                     if idPlot
%                         figure; set(gcf,'position',pplot.(['rect' num2str(k) '_' num2str(m)]));hold all;
%                         histfit(tClean,35); if ~isempty(tOut), plot(tOut,0,'rx'); end; axis tight;
%                         title([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', ts' num2str(Tmat{i}(m)) ', p(outlier): ' num2str(mean(~idClean)*100) '%']);
%                         disp([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', ts' num2str(Tmat{i}(m)) ', p(outlier): ' num2str(mean(~idClean)*100) '%']);
%                         if mean(~idClean)>10
%                             disp('');
%                         end
%                     end                                        
%                 end % m ts
%             end % k targert
%             if idPlot, waitforbuttonpress; close all; end;
%         end % j EH
%     end % prior
%     idOut=[idOut; idOutS];
% end % for iS=1:nS
% save([behDir 'H_RSGprior_DMFC.mat'],'idOut','rangeOut','pOut','-append');
% load([behDir 'H_RSGprior_DMFC.mat'],'T','t','idShortTrial','idHandEye','theta');

load([behDir 'H_RSGprior_DMFC.mat']); % T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)

%% larger bias for long prior @ 20171117
% do linear regression between mean(tp|ts) and ts for each prior > compare slope
% two prior x two EH slopes from each session

sessUni=unique(sessId);
slope=nan(nPr,nEH,length(sessUni));
for iSess=1:length(sessUni)
    % session data
    idSess=sessId==sessUni(iSess);        
    ts=T(~idOut & idSess);
    tp=t(~idOut & idSess);
    ctidx=2-idShortTrial(~idOut & idSess); % 1 for short, 2 for long
    citdx2=2-idHandEye(~idOut & idSess); % 1 for eye, 2 for hand
    
    % getting mean(tp|ts)
    mtp=nan(nPr,nEH,nTspp);
    for i=1:nPr
        for k=1:nTspp
            for j=1:nEH
            id=ctidx==i &... % idShortTrial
                ts==Tmat{i}(k) &... % ts
                citdx2==j;            % eye hand
            mtp(i,j,k)=mean(tp(id));
            end
        end
    end
    
    % slope
    for i=1:nPr
        for j=1:nEH
            tmp=regstats(squeeze(mtp(i,j,:)),Tmat{i}(:),'linear',{'beta'}); % beta
            slope(i,j,iSess)=tmp.beta(2);
        end
    end
end

% plot
hFig=figure;H=[];
for i=1:size(slope,1)
    for j=1:size(slope,2)
        x=squeeze(slope(i,j,:));
        tmptmpcmap=pplot.cmap{(i-1)*size(slope,2)+j};
        [~,~,hTmp]=histStairs(x,15,0,hFig,tmptmpcmap);H=[H;hTmp];ha;
        plot(mean(x),max(get(gca,'ylim')),'v','color',tmptmpcmap,'markerfacecolor',tmptmpcmap);
    end
end
legend(H,'ShortEye','ShortHand','LongEye','LongHand','location','best'); legend boxoff;

% anova
matPr=repmat([1;2],[1 nEH length(sessUni)]);
matEH=repmat([1 2;1 2],[1 1 length(sessUni)]);

stat=anovan(slope(:),{matPr(:) matEH(:)},'model','interaction','display','on','varname',{'idPr','idEH'});

%% SELECT BEST SESSION FOR VISUALIZATION
idExSess=1;
iSess=161218;
if idExSess
    idSess=sessId==iSess; % or 161205
else
    idSess=true(size(sessId));
end

d.ts=T(~idOut & idSess);
d.tp=t(~idOut & idSess);
ctidx=2-idShortTrial(~idOut & idSess); % 1 for short, 2 for long
ipidx=nan(size(ctidx));

%% figure 1B: tp vs ts
plotBLS=1; plotTpMTs=1;
% data formatting: [# interval x maxRepetition]
% find max # trials 1st
for i=1:nPr
    for k=1:nTspp
        id=idShortTrial(~idOut & idSess)==(2-i) &... % idShortTrial
            d.ts==Tmat{i}(k); % ts
        if exist('nMaxTrial')
            nMaxTrial=max([nMaxTrial; nnz(id)]);
        else
            nMaxTrial=nnz(id);
        end
        ipidx(id)=k; % for each trial in d.ts specifies which time bin or ts interval it belongs to. Values: 1-5
    end
end
% filling matrix
Ts=nan(length(Tall),nMaxTrial); % [K x R]
Tp=nan(length(Tall),nMaxTrial);
for i=1:nPr
    for k=1:nTspp
        id=idShortTrial(~idOut & idSess)==(2-i) &... % idShortTrial
            d.ts==Tmat{i}(k); % ts
        Ts((i-1)*nTspp+k,1:nnz(id))=d.ts(id);
        Tp((i-1)*nTspp+k,1:nnz(id))=d.tp(id);
    end
end
d.Ts=Ts;
d.Tp=Tp;
% theory
if idExSess
    idSessTmp=find(unique(sessId)==iSess);
    wm=wFitSess(idSessTmp).w_m;
    wp=wFitSess(idSessTmp).w_p;    
    offset=[wFitSess(idSessTmp).offset2; wFitSess(idSessTmp).offset1];
else
    wm=wFit.w_m;
    wp=wFit.w_p;
    offset=[wFit.offset2; wFit.offset1];
end
for i=1:nPr
    theory{i}.tm=[min(Tmat{i}):max(Tmat{i})]'; % theory{x}.tm - Tx1 - sample corresponding to above
    theory{i}.te=BLS(theory{i}.tm,wm,[min(Tmat{i}) max(Tmat{i})],'uniform')+offset(i);theory{i}.te=theory{i}.te(:);
    theory{i}.range=1:length(theory{i}.tm);
    theory{i}.wm=wm;
end
DnPlotResponses2P(d,theory,ctidx,ipidx,plotBLS,plotTpMTs);

% paper
% set(gca,'xtick',400:400:1200,'ytick',400:400:1200,'ylim',[300 1450]); 
if plotTpMTs
    set(gca,'xtick',[480 640 800 1000 1200],'ytick',-500:100:500,'xticklabel',[],'yticklabel',[]);
    axis tight; %     ylim([350 1400]);
else
    set(gca,'xtick',[480 640 800 1000 1200],'ytick',[480 640 800 1000 1200],'xticklabel',[],'yticklabel',[]);
    ylim([350 1400]);
end
xlim([450 1250]);
xlabel([]);ylabel([]);
title([])
load pplot.mat;
savefig(gcf,'tp_ts_H.fig');
% optsExpFig.Width=5; % 2.6/2.54*2;
% optsExpFig.Height=5; % 1.9/2.54*2;
optsExpFig.Format='eps';
% exportfig(gcf,'tp_ts_H_small.eps',optsExpFig);
% poster
optsExpFig.Width=4/2.54; %6.5; % 10/2.54;
optsExpFig.Height=4/2.54; %6.5; % 7.3/2.54;
if idExportFig
exportfig(gcf,'tp_ts_H.eps',optsExpFig);
end


%% figure 1C: hist

HsHistTpOverlap(d,ctidx);
% paper
set(gca,'xtick',600:200:1000,'ytick',0:20:20,'tickdir','out'); xlabel([]);ylabel([]);
set(gca,'xticklabel',[],'yticklabel',[]);
load pplot.mat;
savefig(gcf,'histTpOverlap_H.fig');
% optsExpFig.Width=2; % 2.6/2.54*2;
% optsExpFig.Height=2; % 1.9/2.54*2;
% optsExpFig.Format='eps';exportfig(gcf,'histTpOverlap_H_small.eps',optsExpFig);
% poster
optsExpFig.Width=3; % 10/2.54;
optsExpFig.Height=2; % 7.3/2.54;
if idExportFig
exportfig(gcf,'histTpOverlap_H.eps',optsExpFig);
end

%% figure 1D: bias variance

dBLS.wm=wm; % wFit.w_m;
dBLS.wp=wp; % wFit.w_p;
dBLS.offset1=offset(2); % wFit.offset1; % for long
dBLS.offset2=offset(1); % wFit.offset2;
% with offset, discrepancy is lower; TBD: model comparison with prior-dependent wm,wp
HsPlotBiasVar(d,ctidx,dBLS);

set(gca,'xtick',0:80:160,'ytick',0:80:160); xlabel([]);ylabel([]);

load pplot.mat;
optsExpFig.Width=2.6/2.54;optsExpFig.Height=1.9/2.54;
optsExpFig.Format='eps';
if idExportFig
    exportfig(gcf,'biasVar_H.eps',optsExpFig);
end
% % poster
% optsExpFig.Width=10/2.54;optsExpFig.Height=7.3/2.54;
% exportfig(gcf,'biasVar_H_large.eps',optsExpFig);

%% stat for overlap

idAbort=t<0;
idNoResp=t>3*T;
idValid=t>0&t<3*T;

% anova with new outlier selected
overlapT=800; 
tmpId=~idOut & idSess &T==overlapT;
stat=anovan(t(tmpId),{idShortTrial(tmpId) idHandEye(tmpId) theta(tmpId)/180},'model',3,'display','on','varnames',{'idShort' 'idEye' 'idRight'});

idPlot=0; % 1;
for i=1:nPr
    for j=1:nEH
        for k=1:nTarg
            id=idShortTrial==(2-i) &... % idShortTrial
                idHandEye==(2-j) &...  % idHandEye
                theta==(k-1)*180 &... % target location
                T==overlapT &...
                ~idOut;
            
            if idPlot
                figure; set(gcf,'position',pplot.(['rect' num2str(i) '_' num2str((j-1)*nTarg+k)]));hold all;
                histfit(t(id),35); axis tight; plotHorizon(gca,mean(t(id)),[]);
                title([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, \mu&SD: ' num2str(meanSD(t(id))) ' ms']);
                disp([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, \mu&SD: ' num2str(meanSD(t(id))) ' ms']);
                xlim([400 1200]);
            end
            
        end % for k=1:nTarg
    end % for j=1:nEH
end % for i=1:nPr

%% BLS fit for each sessions
% for now, wm, wp common for all condition (prior, ts, modality, directions) & offset separately for prior
% model comparion to be done: promising (offset,wp separate for modality)

% sidUni=unique(sessId);
% bicSess=[];wFitSess=[];
% for iS=1:length(sidUni)
%     disp(['===== ' num2str(sidUni(iS)) ' =====']);
%     idS=sessId==sidUni(iS);
%     id=idS&~idOut;
%     
%     [wFitTmp,bicTmp]=estWmWpOld(t(id),T(id),idShortTrial(id),'mmse2offset',[],-1,true); % display: off
%     fieldnm=fieldnames(wFitTmp);
%     for j=1:length(fieldnm)
%         wFitSess(iS).(fieldnm{j})=wFitTmp.(fieldnm{j}); % mmse2offset: w_m w_p offset1 offset2
%     end
%     bicSess=[bicSess; bicSess];
%     disp([wFitTmp.w_m(:) wFitTmp.w_p(:) wFitTmp.offset1(:) wFitTmp.offset2(:)]);
% save([behDir 'H_RSGprior_DMFC.mat'],'wFitSess','bicSess','-append');
% 
% end
% 
% % also BLS fit for pooled data
% analHandEye2PriorRSG('H_RSGprior_DMFC.mwk'); 

%% check across-session variability of wm,wp,offset (wFit w/ outlier removal)
figure; nMeas=4; 
load([behDir 'H_RSGprior_DMFC.mat']); % T t idShortTrial idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)

for iS=1:length(wFitSess)
    subplot(nMeas,1,1);plot(iS,wFitSess(iS).w_m,'o'); ylabel('w_m');hold all;
    subplot(nMeas,1,2);plot(iS,wFitSess(iS).w_p,'o'); ylabel('w_p');hold all;
    subplot(nMeas,1,3);plot(iS,wFitSess(iS).offset1,'o'); ylabel('offset1');hold all;
    subplot(nMeas,1,4);plot(iS,wFitSess(iS).offset2,'o'); ylabel('offset2');hold all;
end

disp('check offset is needed (long,short)')
p=signrank([wFitSess.offset1])
p=signrank([wFitSess.offset2])

%% te
% load pplot.mat;
% % simulation
% nT=10000;
% wm=0.07;
% tsS=[480:80:800]; tsL=[800:100:1200]; ts={tsS;tsL};
% 
% varTe=nan(size(cell2mat(ts))); % [2 x 5]
% skewTe=nan(size(cell2mat(ts))); % [2 x 5]
% h=figure;ha;
% for iPr=1:2
%     T=repmat(ts{iPr},1,nT);T=T(:);
%     Tm=randn(size(T)) .* (wm.*T) + T; %noise
%     te=BLS(Tm,wm,[min(ts{iPr}) max(ts{iPr})],'uniform'); % BLS
%     plotTpTs(T,te,1,pplot.cmap{(iPr-1)*2+1},h,'-',1000); % last for nSD
%     for iTs=1:length(ts{iPr})
%         id=ts{iPr}(iTs)==T;
%         figure; set(gcf,'position',pplot.(['rect' num2str(iPr) '_' num2str(iTs)]));
%         histfit(te(id),35);
%         varTe(iPr,iTs)=var(te(id));
%         skewTe(iPr,iTs)=skewness(te(id));
%     end
% end
% figure(h); xlabel('t_s'); ylabel('t_e'); axis tight; plotIdentity(gca); plotWeberLine(gca,0.2);ylim([480 1250]);
% 
% figure; ha; % var(te)
% plot(tsS,varTe(1,:),'r-');plot(tsL,varTe(2,:),'b-');
% xlabel('t_s'); ylabel('var(t_e)'); axis tight;
% 
% figure; ha; % skew(te)
% plot(tsS,skewTe(1,:),'r-');plot(tsL,skewTe(2,:),'b-');
% xlabel('t_s'); ylabel('skewness(t_e)'); axis tight;


%% define D data structure
% tSort=sort(t(t>0));tB=tSort([1 round([(1/nTp):(1/nTp):1]*length(tSort))]); % 0 0.2 0.4 ... 1
% tLong=sort(t(idShortTrial==0&t>0));tLB=tLong([1 round([(1/nTp):(1/nTp):1]*length(tLong))]); % 0 0.2 0.4 ... 1
% tShort=sort(t(idShortTrial==1&t>0));tSB=tShort([1 round([(1/nTp):(1/nTp):1]*length(tShort))]); % 0 0.2 0.4 ... 1

%% session-pooled outlier removal
% % 2 prior x 2 modality x 2 direction ANOVA for overlap ts
% % after remove outliers (2 SD) for each combination
% idPlot=1;
% nSDremoveOutlier=2;
% overlapT=800; 
% tMat=[];iPr=[];iEH=[];iRL=[];
% for i=1:nPr
%     for j=1:nEH
%         for k=1:nTarg
%             id=idShortTrial==(2-i) &... % idShortTrial
%                 idHandEye==(2-j) &...  % idHandEye
%                 theta==(k-1)*180 &... % target location
%                 T==overlapT &...
%                 t>0 & t<3*overlapT; % only valid trials for now
%             u0=[0.5 1.5]*overlapT; %[0.1 1.9]*overlapT; % [0.5 1.5]*overlapT; %[0 2]*overlapT; % with least outlier removal; or max(t)
%             [tClean,idClean,pOut,tOut]=removeOutlier(t(id),0,u0); %nSDremoveOutlier);
%             tMat=[tMat;tClean(:)];
%             iPr=[iPr; repmat(prNm(i),[nnz(tClean) 1])];iEH=[iEH; repmat(ehNm(j),[nnz(tClean) 1])];iRL=[iRL; repmat(targNm(k),[nnz(tClean) 1])];
%             
%             if idPlot
%                 figure; set(gcf,'position',pplot.(['rect' num2str(i) '_' num2str((j-1)*nTarg+k)]));hold all;
%                 histfit(tClean,35); if ~isempty(tOut), plot(tOut,0,'rx'); end; axis tight;
%                 title([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, p(outlier): ' num2str(pOut) '%']);
%                 disp([prNm{i} ', ' ehNm{j} ', ' targNm{k} ', overlap, p(outlier): ' num2str(pOut) '%']);
%                 xlim([400 1200]);
%             end
%             
%         end % for k=1:nTarg
%     end % for j=1:nEH
% end % for i=1:nPr
% stat=anovan(tMat,{iPr iEH
% iRL},'model',3,'display','on','varnames',{'idShort' 'idEye' 'idRight'})
