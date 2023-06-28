function idOut=makeIdOut(tInfo,beh,fid)

% function idOut=makeIdOut(tInfo,beh,fid)
% make idOut [size(tInfo,1) x 1] for a given session's tInfo
% tInfo
% trial stage codes (2,3,4,5,6,8/9 for now)
% 1. Fix On
% 2. Fix % photodiode sync
% 3. target On % photodiode sync
% 4. Ready % photodiode sync
% 5. Set % photodiode sync
% 6. Production 
% 7. target acquired 
% 8. reward delay % photodiode sync
% 9. reward: bonusRewDur*1000+rewardDur*1000*max(0.00001,1-abs(productionInterval - interval)/interval/win_fraction ) % photodiode sync
% 10. trial end
% 0. ITI
% 8. incorrect
% 0. Bad
%
% tInfo: trial#, stageCode, t(blackRock), 
%           idShortTrial(4), idHandEye(5), theta(6), T(7), t(8), fixTimeDur(9), targetTimeDur(10), iti(11), reward(12)]
%
% fid: string e.g. 161218

%% remove outlier

% initial
idOut=false(size(tInfo,1),1); % 42602
movingBallTrials=unique(tInfo(:,1)); % ~4000
idSess=str2num(fid)==beh.sessId; % ~40000

% special treat for individual sessions
if strcmp(fid,'161207'),
    idSess(find(idSess==1,nnz(idSess)-length(movingBallTrials),'first'))=0;
end; % for H' 12/7, recorded for later 1711 trials after 800 trials w/o recording
if nnz(idSess)~=length(movingBallTrials),disp('smth wrong for trial align'); end;

iOutSess=find(beh.idOut(idSess)); % 1157

% report
disp(['outlier trials (incl. abort/wait): ' num2str(length(iOutSess)) '/' num2str(length(movingBallTrials)) ',' ...
    num2str(length(iOutSess)/length(movingBallTrials)*100) '%']);
idValid=beh.t(idSess)>0&beh.t(idSess)<3*beh.T(idSess);
disp(['outlier trials (w/o abort/wait): ' num2str(nnz(beh.idOut(idSess)&idValid)) '/' num2str(length(movingBallTrials)) ',' ...
    num2str(nnz(beh.idOut(idSess)&idValid)/length(movingBallTrials)*100) '%']);

% main
for iOut=1:length(iOutSess)
    idOut(tInfo(:,1)==movingBallTrials(iOutSess(iOut)))=true;
end