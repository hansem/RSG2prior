function [X,varargout]=estReadout(D,varargin)

% function xro=estReadout(D)
% estimating readout vector for rotation in prior support
% input: D (datahigh format; 1:5 ts) for each prior, data should end at time of Set
% output: xro [size(D(1).data,1)] normalized
%
% - use all data by averaging across ts (e.g. earlier data from longest ts)
% - also avg. across ts5-ts1 and ts4-ts2
%
% 2018/6/22
% use all mid-points, not only end points of each ts
%
% 2018/7/18
% dealing with SiL
% better estimation when data length is even
%
% 2018/9/19
% variance decrease by changing uts connecting only ts1-ts2, ts4-ts5? idUtsNew
%
% 2018/9/23
% for all ts, use ts3-ts1 for ts2, ts4-ts2 for ts3, ts5-ts3 for ts4

Dtmp=struct2mat(D,'data'); % [dim Maxtime 5conditions]
mD=squeeze(nanmean(Dtmp,3)); % [dim Maxtime]
nT=size(mD,2);
xro=[]; % for ts1-ts2

% variance decrease by changing uts connecting only ts1-ts2, ts4-ts5? idUtsNew
if ~isempty(varargin)
    idUtsNew=varargin{1};
    if idUtsNew
        % ts1-ts2
        nT1=size(D(2).data,2);
        if rem(nT1,2)==0 % even
            nTavg=nT1/2;
            for iT=1:nTavg
                xro=[xro mD(:,nT1-(iT-1))-mD(:,iT)]; % [nD x all mid points]
            end
        else % odd
            nTavg=floor(nT1/2);
            for iT=1:nTavg
                xro=[xro mD(:,nT1-(iT-1))-mD(:,iT)]; % [nD x all mid points]
            end
        end
        X=normalize(mean(xro,2),1); % [# PC x 1] normalized
        
        % 2: ts1-ts3
        xro2=[]; 
        nT22=size(D(3).data,2);
        nT21=size(D(1).data,2);
        nT2=nT22-nT21+1;
        if rem(nT2,2)==0 % even
            nTavg=nT2/2;
            for iT=1:nTavg
                xro2=[xro2 mD(:,nT22-(iT-1))-mD(:,nT21+(iT-1))]; % [nD x all mid points]
            end
        else % odd
            nTavg=floor(nT2/2);
            for iT=1:nTavg
                xro2=[xro2 mD(:,nT22-(iT-1))-mD(:,nT21+(iT-1))]; % [nD x all mid points]
            end
        end
        X2=normalize(mean(xro2,2),1); % [# PC x 1] normalized
        varargout{1}=X2;
        
        % 3: ts2-ts4
        xro3=[]; 
        nT32=size(D(4).data,2);
        nT31=size(D(2).data,2);
        nT3=nT32-nT31+1;
        if rem(nT3,2)==0 % even
            nTavg=nT3/2;
            for iT=1:nTavg
                xro3=[xro3 mD(:,nT32-(iT-1))-mD(:,nT31+(iT-1))]; % [nD x all mid points]
            end
        else % odd
            nTavg=floor(nT3/2);
            for iT=1:nTavg
                xro3=[xro3 mD(:,nT32-(iT-1))-mD(:,nT31+(iT-1))]; % [nD x all mid points]
            end
        end
        X3=normalize(mean(xro3,2),1); % [# PC x 1] normalized
        varargout{2}=X3;
        
        % 4: ts3-ts5
        xro4=[]; 
        nT42=size(D(5).data,2);
        nT41=size(D(3).data,2);
        nT4=nT42-nT41+1;
        if rem(nT4,2)==0 % even
            nTavg=nT4/2;
            for iT=1:nTavg
                xro4=[xro4 mD(:,nT42-(iT-1))-mD(:,nT41+(iT-1))]; % [nD x all mid points]
            end
        else % odd
            nTavg=floor(nT4/2);
            for iT=1:nTavg
                xro4=[xro4 mD(:,nT42-(iT-1))-mD(:,nT41+(iT-1))]; % [nD x all mid points]
            end
        end
        X4=normalize(mean(xro4,2),1); % [# PC x 1] normalized
        varargout{3}=X4;
        
        xro5=[]; % for ts4-ts5
        nT52=size(D(5).data,2);
        nT51=size(D(4).data,2);
        nT5=nT52-nT51+1;
        if rem(nT5,2)==0 % even
            nTavg=nT5/2;
            for iT=1:nTavg
                xro5=[xro5 mD(:,nT52-(iT-1))-mD(:,nT51+(iT-1))]; % [nD x all mid points]
            end
        else % odd
            nTavg=floor(nT5/2);
            for iT=1:nTavg
                xro5=[xro5 mD(:,nT52-(iT-1))-mD(:,nT51+(iT-1))]; % [nD x all mid points]
            end
        end
        X5=normalize(mean(xro5,2),1); % [# PC x 1] normalized
        varargout{4}=X5;
        
        return;
    end
end

% original readout
if rem(nT,2)==0 % even
    nTavg=nT/2;
    for iT=1:nTavg
        xro=[xro mD(:,end-(iT-1))-mD(:,iT)]; % [nD x all mid points]
    end
else % odd
    nTavg=floor(nT/2);
    for iT=1:nTavg
        xro=[xro mD(:,end-(iT-1))-mD(:,iT)]; % [nD x all mid points]
    end
end
X=normalize(mean(xro,2),1); % [# PC x 1] normalized


return;

%% old

% find prior mean
if sum(diff(nData(D,2)))==0 % SiL
    nTimeMuPr=round(size(D(1).data,2)/2);
else
    nTimeMuPr=size(D(median(1:length(D))).data,2);
end

% from prior mean, avg for each points
xro=[];
nTavg=min([nT-nTimeMuPr nTimeMuPr-1]);
for iT=1:nTavg
    xro=[xro mD(:,nTimeMuPr+iT)-mD(:,nTimeMuPr-iT)]; % [nD x all mid points]
end
X=normalize(mean(xro,2),1); % [# PC x 1] normalized



return;

%% avg. across ts5-ts1 and ts4-ts2

% get t(set)
nTime5=size(D(5).data,2); % ts(5)-ts(1)
nTime1=size(D(1).data,2);
nTime4=size(D(4).data,2); %  ts(4)-ts(2)
nTime2=size(D(2).data,2);

xro=mD(:,nTime5)-mD(:,nTime1);
xro2=mD(:,nTime4)-mD(:,nTime2);

% if not using all data
% xro=D(5).data(:,nTime5)-D(1).data(:,nTime1);
% xro2=D(4).data(:,nTime4)-D(2).data(:,nTime2);
X=normalize((xro+xro2)/2,1); % [# PC x 1] normalized

% X=normalize(mD(:,nTime5)-mD(:,nTime1),1);