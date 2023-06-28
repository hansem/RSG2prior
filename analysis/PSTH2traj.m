function [traj,varargout]=PSTH2traj(D,proj_matrix,keep_neurons,binSize,smthWidth,optimD,varargin)
% function traj=PSTH2traj(D,proj_matrix,keep_neurons,binSize,smthWidth,optimD,varargin)
%
% reduce dimensionality (for now using PCA only)
% also assuming this is for bootstrapping and proj_matrix,keep_neurons are given
% D and traj are in DataHigh format
% varargin for use_sqrt, mean_thresh
% if optimD is empty, use optimal_SVHT_coef to decide optimalD

% 2018/5/2
% if proj_matrix is provided, offset (mean) for centering should be also provided
% format of proj_matrix: [neuron x pc]

% within-function parameters
mean_thresh0=1;
thPVAF=75; % TBD: optimal_SVHT_coef

%% main
handles.binWidth=binSize;
if isempty(varargin)
    handles.use_sqrt=0;
    mean_thresh=mean_thresh0; % remove neurons with FR<1 sp/s
else    
    handles.use_sqrt=varargin{1};
    if length(varargin)>1
        mean_thresh=varargin{2};
    else
        mean_thresh=mean_thresh0; % remove neurons with FR<1 sp/s
    end
end
handles.kern=smthWidth;
if isempty(keep_neurons) %%%%% do PCA for own data    
    m = mean([D.data],2) * 1000;
    keep_neurons = m>=mean_thresh; %-Inf; %m >= mean_thresh; % -Inf; % mean_thresh; %%%%%
    handles.mean_thresh=mean_thresh; 
end
handles.keep_neurons=keep_neurons;

% Remove low firing rate neurons
for itrial = 1:length(D)
    D(itrial).data = D(itrial).data(handles.keep_neurons,:);
end
disp(['# neurons: ' num2str(size(D(1).data,1))]);

% Bin the spikes and use square root transform (use gpfa's pack)
if (handles.binWidth ~= 1 || handles.use_sqrt)
    [d(1:length(D)).spikes] = deal(D.data);
    for itrial = 1:length(D)
        d(itrial).trialId = itrial;
    end
    s = getSeq(d, handles.binWidth, 'useSqrt', handles.use_sqrt);
    [D.data] = deal(s.y);
end

% Smooth data if necessary (automatically zero if GPFA selected)
for i = 1:length(D)
    D(i).data = smoother(D(i).data, handles.kern, handles.binWidth);
end

% do PCA
alldata = [D.data]; % [neuron x time]
% figure; imagesc(alldata); % debug
[u sc lat] = pca(alldata'); % princomp(alldata'); % alldata'=sc*u' [n obs/time x p var/neuron]=[n x k]*[k x p]
% ha;plot(cumsum(lat(1:15))./sum(lat)); % debug scree

% decide optimal D (heuristic: find kink in scree plot)
if isempty(optimD)
    if size(alldata',1)<size(alldata',2), disp(['# obs ' num2str(size(alldata',1)) '< # feat ' num2str(size(alldata',2)) ': risk of overfitting noise']); end;
    
    [nRow,nCol]=size(alldata');
    sigma_known=0;
    
    if nCol<3
        traj=[];
        varargout{1}=[];
        varargout{2}=[];
        varargout{3}=[];
        return;
    end
    
    optimD=find(100*cumsum(lat)./sum(lat)>thPVAF,1,'first');
    disp(['optimD: ' num2str(optimD)]);
    if optimD<=3
        optimD=3;disp(['optimD: ' num2str(optimD)]);
    end
    
    % optimal_SVHT_coef
%     if nRow/nCol>1
%         optimD=length(lat); 
%     else
%         thSV=optimal_SVHT_coef(nRow/nCol,sigma_known)*median(sqrt(lat)); % assuming sqrt(eigenvalue)=singular value
%         optimD=find(sqrt(lat)>thSV,1,'last');
%     end
end

% C = u(:,1:optimD); % [p x k]
disp(['total variance from its own PCA: ' num2str(sum(lat))]); %     lat = (lat) ./ sum(lat);  % lat = cumsum(lat(1:dims)) ./ sum(lat(1:dims));  % eigenvalues; % include all dims
if isempty(proj_matrix) %%%%% do PCA for own data
    proj_matrix=u(:,1:optimD);% C;
    varargout{4}=mean(alldata',1)'; % [neuron x 1]
    
    disp(['% variance explained for its own PCA: ' num2str(100*sum(lat(1:optimD))./sum(lat)) ]);
else
    % projection to original PC space
    if ndims(proj_matrix)>=3 % proj_matrix=squeeze(proj_matrix(:,:,1))     proj_matrix(:,:,2) for centering offset
        meanCenter=squeeze(proj_matrix(:,1,2)); % [neuron x 1]
        proj_matrix=squeeze(proj_matrix(:,:,1)); % [
        cdata=bsxfun(@minus,alldata',meanCenter(:)'); % centering as PCA % [time x neuron]
    else % old: but not correct
        cdata=bsxfun(@minus,alldata',mean(alldata',1)); % centering as PCA % [time x neuron]
    end
    disp(['total data variance: ' num2str(sum(var(cdata,0,1)))]);
    sc =cdata/proj_matrix'; % projection
    disp(['variance explained for this projection: ' num2str(sum(var(sc,0,1)')) '(' num2str(100*sum(var(sc,0,1)')/sum(var(cdata,0,1))) '%)']);
    lat=var(sc,0,1)'; % eigen value of COV=variance of scores
%     lat2=cumsum(var(sc,0,1)') ./ sum(var(sc,0,1)');

    varargout{4}=sum(var(cdata,0,1));
end
% For each condition, store the reduced version of each data vector
index = 0;
for i=1:length(D)
    D(i).data = sc(index + (1:size(D(i).data,2)),1:optimD)';
    index = index + size(D(i).data,2);
end

% need to reformat epochStarts to match the binWidth
if (isfield(D, 'epochStarts'))
    for itrial = 1:length(D)
        if (~isempty(D(itrial).epochStarts))
            D(itrial).epochStarts = ceil(D(itrial).epochStarts / handles.binWidth);
        end
    end
end

traj = D;
    varargout{1}=proj_matrix;
    varargout{2}=keep_neurons;
    varargout{3}=lat;
