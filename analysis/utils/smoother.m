function yOut = smoother(yIn, kernSD, stepSize, varargin)
%
% yOut = smoother(yIn, kernSD, stepSize)
%
% Gaussian kernel smoothing of data across time.
%
% INPUTS:
%
% yIn      - input data (yDim x T)
% kernSD   - standard deviation of Gaussian kernel, in msec
% stepSize - time between 2 consecutive datapoints in yIn, in msec
%
% OUTPUT:
%
% yOut     - smoothed version of yIn (yDim x T)
%
% OPTIONAL ARGUMENT:
%
% causal   - logical indicating whether temporal smoothing should
%            include only past data (true) or all data (false)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

% Aug 21, 2011: Added option for causal smoothing
% Nov 25, 2019: deal with nan

  causal = false;
  assignopts(who,varargin);

  if (kernSD == 0) || (size(yIn, 2)==1)
    yOut = yIn;
    return
  end

  % Filter half length
  % Go 3 standard deviations out
  fltHL = ceil(3 * kernSD / stepSize);

  % Length of flt is 2*fltHL + 1
  flt = normpdf(-fltHL*stepSize : stepSize : fltHL*stepSize, 0, kernSD);

  if causal
    flt(1:fltHL) = 0;
  end

  [yDim, T] = size(yIn);
  yOut      = nan(yDim, T);

  % Normalize by sum of filter taps actually used
%   nm = conv(flt, ones(1, T));
  
  for i = 1:yDim
      tmp=yIn(i,:);
      idNan=isnan(tmp);
      
      tmp2=ones(size(tmp));
      tmp2(idNan)=0;  % Replace NaNs with zero
      nm = conv(flt, tmp2);
      
      tmp(idNan)=0;  % Replace NaNs with zero
      
      ys = conv(flt, tmp) ./ nm;
    % Cut off edges so that result of convolution is same length 
    % as original data
    yOut(i,:) = ys(fltHL+1:end-fltHL);
    
    yOut(i,idNan)=NaN; % replace output values with NaNs corresponding to input.
    
%   figure; 
%   subplot(4,1,1); plot(tmp); ylabel('raw');
%   subplot(4,1,2), plot(conv(flt, tmp)); ylabel('conv(flt, tmp))');
%   subplot(4,1,3), plot(nm); ylabel('conv(flt, ones))');
%   subplot(4,1,4), plot(yOut(i,:)); ylabel('output with truncation');

  end
  