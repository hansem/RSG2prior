function sFig_recordingSite

% plot recording site in horizontal view
% RSG2pr + RSGadapt (1pr; not 2pr-adapt)
% for G, only plot sites within chamber (AR of chamber: A27 R6; PR: A12 R65)
% 2021/1/25 hansem@mit.edu
% modified from makemap.m in dPCA
%

%% AR coordinate of each recording sessions
% (TBD): double check with pictures in H_RSGprior_recordingSites
% for now, just use estimates from labnote; 2017/4/3

coord=cell(2,2); % H/G, 2pr/1pr-adapt (red/green)

%% H
coord{1,1}={... % anterior/right to interaural stereotex zero [bank A],[B],[C]
%[23.8 7.6] %     'H_RSGprior_20161202';... % only eye; use 003
[23.5 9];...%     'H_RSGprior_20161203';... % no photosync; RSG>overlap>RSG; use 003
[24.3 7.5];... %     'H_RSGprior_20161204';... % no photosync
[24.5 9];...%     'H_RSGprior_20161205';... % no photosync; no overlap
[24.2 6.7; 24.6 7.7];...%     'H_RSGprior_20161206';...
[24.2 7.3; 24.2 7.7];... %     'H_RSGprior_20161207';...
[23.4 7.1; 23.1 7.7; 23.5 8.5];... %     'H_RSGprior_20161208';... 
%     % 2nd craniotomy @12/9 afterward, [5 3 20]
[25.3 9; 24.8 10; 26.1 10.2];...%     'H_RSGprior_20161210';...
[25.1 5.9; 25 7.3; 25.8 7];...%     'H_RSGprior_20161211';...
[25.6 8.4; 26.2 10; 26.2 8];... %     'H_RSGprior_20161212';...
% [25.4 5.5; 25.4 6; 26 5.4];...% %     'H_RSGprior_20161213';...
[26.3 6];...%     'H_RSGprior_20161214';...
[26.4 5.2];...%     'H_RSGprior_20161215';... 
%     % 3rd craniotomy @12/16
[26.7 5;26.8 5.9];...%     'H_RSGprior_20161218';...
[27.2 5.3; 27.3 6.2];...%     'H_RSGprior_20161219';...
[27.4 6.4;27.5 7.2];...%     'H_RSGprior_20161220';...
[27.3 4.3; 27.7 4.9];...%     'H_RSGprior_20161221';...
[26.8 4.6; 27 4.9];...%     'H_RSGprior_20161222';... % saved as 12/21
[27.1 6.5; 26.9 7.5]%     'H_RSGprior_20161223'
    };

% 1pr adapt (assuming inter-hole distance 0.8mm)
coord{1,2}={...
[30 -5.5-0.8;30 -5.5;30 -5.5+0.8];...%     190924*
[31 -5.5-0.8;31 -5.5;31 -5.5+0.8];... %     190925*
[32 -5.5-0.8;32 -5.5;32 -5.5+0.8];... %     190926
[29 -5.5-0.8;29 -5.5;29 -5.5+0.8];... %     190927
[34 -5.5-0.8;34 -5.5;34 -5.5+0.8];... %     200320*
[33 -5.5-0.8;33 -5.5;33 -5.5+0.8];... %     200322
[32 -5.5-0.8;32 -5.5;32 -5.5+0.8];... %     200323
[31 -5.5-0.8;31 -5.5;31 -5.5+0.8];... %     200324
[30 -7-0.8;30 -7;30 -7+0.8];... %     200326
[29 -7-0.8;29 -7;29 -7+0.8]}; %     200328*
    


%% G
xChamber=11; yChamber=16;
% for 2pr, estimate as no guide tube was used
coord{2,1}={...
[1 14; 1 15; 1 16];... %     'G_RSGprior_20170506';... % A26 R5
[1 11; 1 12; 1 13];... %     'G_RSGprior_20170507';...
[2 12; 2 13; 2 14];... %     'G_RSGprior_20170508';...
[2 14; 2 15; 2 16];... %     'G_RSGprior_20170510';... % quite>2nd craniotomy
[3 10; 3 11; 3 12];... %     'G_RSGprior_20170511';...
[2 10; 2 11; 2 12];... %     'G_RSGprior_20170512';...
%     
[1 2.5; 1 3.5; 1 4.5];... %     'G_RSGprior_20170803';...
%     
[1 4.5; 1 5.5; 1 6.5];... %     'G_RSGprior_20170817';...
[2 4.5; 2 5.5; 2 6.5];... %     'G_RSGprior_20170818';...
[3 4.5; 3 5.5; 3 6.5];... %     'G_RSGprior_20170821';...
[4 4.5; 4 5.5; 4 6.5];... %     'G_RSGprior_20170822';...
[5 4.5; 5 5.5; 5 6.5];... %     'G_RSGprior_20170823';...   
    };

coord{2,2}={...
[4 4; 4 5; 4 6];... %     180522
[3 4; 3 5; 3 6];... %     180523*
[5 4; 5 5; 5 6];... %     180524*
[2 4; 2 5; 2 6];... %     180525
[5 4; 5 5; 5 6];... %     180528
[4 4; 4 5; 4 6];... %     190425
    };

%% depth: 1st neuron


%% channel map


%% note from Hchamber docs
% 1st craniotomy
% center of chamber: 126.25, x, 216 (+R,+S,+A; 0.5mm iso MRI coordinate)
% 2.5mm drill

% 12/5: area6DR(dorsorostral, F7) or area6DC(F2), just dorsal/medial to anterior tip of PCD
% from medial wall? Unlikely as 1) gap b/t guide tube top & Vprobe tip (~1mm), 2) gap b/t guide tube bottom & dura (~1 or 2mm), 3) driving probe drag tissue, 4) depth was usually 5-6mm
% 
% 12/15: units with eye selectivity (motor) A26, R5.5 (z:212, x:121)
% 
% 12/16: 3rd craniotomy, A24.5, R6 (z:209, x:122)


%% import MRI 
% : oldest w/o artifact (grow not noticeabel from mricron overlay)
% (all in MRI coordinate: x+R, y+Superior, z+Posterior)
% (214 x 448 x 448 voxels; 0.5mm iso)
% stereotex zero
stereoZero=[110 257 264]; % [voxels]



load pplot.mat;
hFig=figure(999); set(gcf,'position',pplot.rect2_1);
plotBackground(hFig);


% H 2pr
coordMat=cell2mat(coord{1,1});
plot(coordMat(:,2),coordMat(:,1),'r.'); % recording sites
% H 1pr
coordMat=cell2mat(coord{1,2});
plot(coordMat(:,2),coordMat(:,1),'.','color',[0 0.6 0]); % recording sites

axis tight;
xlabel('right from mid-interaural'); ylabel('anterior from mid-interaural');

applytofig4keynote;

% note
% save data from around 800 trials in 12/7
% data transfer speed warning issue in 12/13? two cells in NS3 data


figure; set(gcf,'position',pplot.rect1_1); ha;
% G 2pr
coordMat=cell2mat(coord{2,1});
plot(coordMat(:,2),xChamber-coordMat(:,1),'r.'); % recording sites
% G 1pr
coordMat=cell2mat(coord{2,2});
plot(coordMat(:,2),xChamber-coordMat(:,1),'.','color',[0 0.6 0]); % recording sites

axis tight;
ylim([0 xChamber+1]);
xlim([0 yChamber+1]);
set(gca,'xtick',0:5:20,'ytick',0:5:20);
% axis equal;
% xlabel('chamber horizontal location'); 
% ylabel('chamber vertical location'); 

applytofig4keynote;



function plotBackground(hFig)
figure(hFig); hold all;
% only zx in chamber coordinate (relative to stereoZero)
% chamber anteriorLeft
% chamber anteriorRight
% chamber posteriorLeft
% chamber posteriorRight
% Chamber size (inside): 15 mm(AP) x 10 mm(ML) with 1.5 mm thickness
LA=[29 2.5]; % [z+A x+R] mm
RA=[31 11.5];
LP=[17 4.5];
RP=[19 14];

xzChamber=[LA;RA;RP;LP;LA];

C=[24 8.125];

% misc info (left-right, anterior-posterior)
AS_ua=[32 11]; % arcuate sulcus upper arm
AS_genu=[22 17.5]; % arcuate sulcus genu
AS_ua2=[19 19]; % arcuate sulcus upper arm posterior
PCD=[23 10; 15 11.5];
midline=[32 0; 15 0];

plot(midline(:,2),midline(:,1),'k-'); % midline
plot([AS_ua(2); AS_genu(2); AS_ua2(2)],[AS_ua(1); AS_genu(1); AS_ua2(1)],'k-'); % AS
plot(PCD(:,2),PCD(:,1),'k-');
plot(xzChamber(:,2),xzChamber(:,1),'k:'); % chamber
axis tight;
xlabel('right from mid-interaural'); ylabel('anterior from mid-interaural');


function normX=normMaxMin(x,X)
normX=(x-min(X))/(max(X)-min(X));

