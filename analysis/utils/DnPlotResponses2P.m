function []  = DnPlotResponses2P(d, theory, ctidx, ipidx,varargin)

% -----------------------------------------------------
% July 2017
% Devika Narain Aerts
% Massachusetts Institute of Technology, Cambridge, USA
% -----------------------------------------------------


% -----------------------------------------------------
%   Plot function for 2 Prior Project
% -----------------------------------------------------


% modification history
% line#115, 134, 161: changed index for Colors2P to use flipped colormap % HS 8/30/17
% line#179: histogram size normalized to inter-ts distance % HS 8/30/17

% adding option of ploting tp-ts, rather than ts 'PlotDiff' % HS 5/10/2018

%% Usage:
% Input structure: d
% d.ts - Nx1 where N: no of trials - sample
% d.tp - Nx1 where N: no of trials - responses
% d.Ts - KxR where K - no of intervals or bins (10) and R - no of reps per interval
%        or bin
% d.Tp - KxR where K - no of intervals or bins (10) and R - no of reps per interval
%        or bin

% theory{x}.te - Tx1 - continuous time BLS estimate
% theory{x}.tm - Tx1 - sample corresponding to above
% theory{x}.range - T (specifies range over which bls was calculated)
% theory{x}.wm - Wm at which BLS was calculated
% here x refers to the prior number

% ctidx - for each trial in d.ts specifies whether it belongs to prior 1 or
%         2. Enteries 1 or 2
% ipidx - for each trial in d.ts specifies which time bin or ts interval it
%         belongs to. Values: 1-5

%%

%       -----------------  Select Flags ----------------

        % 1 for Dot and 0 for Histogram        2 for shadedErrorb    
        Dot_or_Hist     =  2; % 2; % 0; % 2; % 1;
        
        % Overlay BLS? plot tp-ts?
        if ~ isempty(varargin)
            PlotDiff         =  varargin{1}; % 1;
            idStim=0;
            PlotBLS=1;
            switch length(varargin)
                case 2
                    PlotDiff=varargin{2};
                    idStim=0;
                case 3
                    PlotDiff=varargin{2};
                    idStim=varargin{3};
                case 4
                    PlotDiff=varargin{2};
                    idStim=varargin{3};
                    hFig=varargin{4};
                        
            end
%             if length(varargin)>1
%                 PlotDiff=varargin{2};
%             else
%                 PlotDiff=1;
%             end
        else
            PlotBLS=1;
            PlotDiff=0; % 1;
            idStim=0;
        end
        
        % 1 Show individual dots or 0 for only averages
        if PlotDiff==0
            PlotInd         =  1;
        else
            PlotInd         =  0;
        end
        
        markerSize=3; % 4; % 3; % avg
        markerSize2=1/10; % individual trials
        lw=0.5; % 2;
        
        % Checks
        if(~isempty(d.ts))
            Trialno     = length(d.ts);
            ts          = d.ts;
        else
            DnDisp('Input structure lacking key field: ts')
            error('Please fix structure d')
        end
        
        if(~isempty(d.tp))
            tp          = d.tp;
        else
            DnDisp('Input structure lacking key field: tp')
            error('Please fix structure d')
        end
        
        if(~isempty(d.Ts))
            [K, R]      = size(d.Ts);
            Ts          = d.Ts;
        else
            DnDisp('Input structure lacking key field: Ts')
            error('Ts must exist and be intervals x Reps')
        end
        
        if(~isempty(d.Tp))
            Tp          = d.Tp;
        else
            DnDisp('Input structure lacking key field: Tp')
            error('Tp must exist and be intervals x Reps')
        end
 
%         if(isempty(theory))
%       
%             PlotBLS     = 0
%             DnDisp('BLS structure missing')
%        
%         end
        
        
        % Add labels
        if Dot_or_Hist~=1
        f1             = figure; 
        hold on;
        end
        if ~exist('hFig')
          f2             = figure; 
          hold on;
        else
            f2=hFig; hold on; 
        end
        figure(f2);hold on; 

%         hTitle         = title('Behavior');
%         hXLabel        = xlabel('Time');
%         hYLabel        = ylabel('OP');
       
        
        % Adjust font
%         set(gca, 'FontName', 'Helvetica')
%         set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
% 
%         set([hXLabel, hYLabel], 'FontSize', 10)
%         set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')

        % Adjust axes properties
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
        'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
        'LineWidth', 1)                                    


        % Color indices
        if idStim==1
            Clo            = Colors2P(5);
        else
        Clo            = Colors2P(1); % HS 8/30/17; previously 1; 5 for hand
        end
        zo             = ctidx;
        io             = ipidx;
        
%%
  if Dot_or_Hist ==1
  % If Dot Plot
  plot([400 1250], [400 1250],'k:')
      if PlotInd
      % If plot individual points        
        for i          = 1:Trialno
            clr        = Clo(zo(i),io(i),:);
            plot(ts(i),tp(i),'.','Color',clr,'MarkerSize',markerSize2)
               
        end   
      end
      
    % Plotting Averages 
    if idStim==1
         clr     = Colors2P(6)
    else
                clr     = Colors2P(2); % HS 8/30/17; previously 2
    end
                TsA     = nanmean(Ts,2);
                TpA     = nanmean(Tp,2);
                for k = 1:K, TpSD(k)= nanstd(Tp(k,:))./sqrt(R); end

                  errorbar(TsA,TpA,TpSD,'k','LineStyle','none')
                  for k   = 1:K     
                  plot(TsA(k),TpA(k),'o','MarkerFaceColor','w','MarkerSize',markerSize,'linewidth',2,'color',clr(k,:))   
                  end   
                  
      if PlotBLS 
          for p = 1:2
              if p==1
          plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range),'-r')
              else
                  plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range),'-b')
              end
          end
          
      end
       
  elseif Dot_or_Hist==0
  %% If Hist Plot
  
   % Plotting Averages 
   figure(f2), hold on
    plot([400 1250], [400 1250],'k:')
    if idStim==1
         clr     = Colors2P(6);
    else
                clr     = Colors2P(2); % HS 8/30/17; previously 2
    end
                TsA     = nanmean(Ts,2);
                TpA     = nanmean(Tp,2);
                for k = 1:K, TpSD(k)= nanstd(Tp(k,:))./sqrt(R); end


  
        if PlotInd
            
                nbins   = 11;
            
            for k       = 1:K
               
                figure(f1), hold on; 
                h       = histogram(Tp(k,:),nbins);
                Xval    = h.BinEdges + h.BinWidth*0.5;
                Xval    = Xval(1:end-1);
                if k~=K & k~=round(K/2) % HS 8/30/2017; dealing with max ts of each prior
                    Yval    = TsA(k)+h.Values./(max(h.Values))*0.8*(TsA(k+1)-TsA(k)); % h.Values*10 + TsA(k); % normalized so that maximum is set to 0.8*betweenBinDistance
                else
                    Yval    = TsA(k)+h.Values./(max(h.Values))*0.8*(TsA(k)-TsA(k-1));  %  h.Values*10 + TsA(k); % normalized so that maximum is set to 0.8*betweenBinDistance
                end
                % setting baseline separately for each ts's bar is impossible; resolving this with 'hist' (drawing patch) is
                % only workable in old MATLAB > for now use stairs
                figure(f2), hold on
%                 h2      = stairs(Yval,Xval,'color',clr(k,:)); %  
                hBarH=barh(Xval,Yval,'FaceColor',clr(k,:),'EdgeColor','none','BaseValue',TsA(k));
                hBarH.BaseLine.Visible='off';

            end
        end
        
        
              errorbar(TsA,TpA,TpSD,'k','LineStyle','none')
              for k     = 1:K     
              plot(TsA(k),TpA(k),'ok','MarkerFaceColor',clr(k,:),'MarkerSize',markerSize)   
              end   

              
        
     if PlotBLS 
         
          for p = 1:2
              if p==1
          plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range),'-r')
              else
                  plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range),'-b')
              end
          end
          
     end

      %% actually used
  elseif Dot_or_Hist==2 % used shadedErrorBar
      if ~PlotDiff
  plot([400 1250], [400 1250],'k:')
  plot([400 1250], [800 800],'k:')
      else
          plot([400 1250], [0 0],'k:')
      end
%   plot([800 800], [400 800],'k:')
  nTspp=round(K/2);    
  tOverlap=0; % 800; % plot dots for not-overlap ts
  nbins   = 20;barWidth=0.4;
  nPr=nnz(unique(ctidx));
  if PlotInd
      % If plot individual points        
      pRand=0.1; tsList=unique(ts); dTs=tsList(2)-tsList(1); % spread a little bit around ts
      
        for i          = 1:Trialno
            clr        = Clo(zo(i),io(i),:);
            if ts(i)~=tOverlap
                if PlotDiff
            plot(ts(i)+pRand*2*(rand(1,1)-0.5)*dTs,tp(i)-ts(i),'.','Color',clr,'MarkerSize',markerSize2,'linewidth',lw/10);
                else
                    plot(ts(i)+pRand*2*(rand(1,1)-0.5)*dTs,tp(i),'.','Color',clr,'MarkerSize',markerSize2,'linewidth',lw/10);
                end
            end
               
        end
        
        % histogram for overlap
        if tOverlap~=0
            if idStim==1
                clr     = Colors2P(6)
            else
                clr     = Colors2P(2); % HS 8/30/17; previously 2
            end
            tpOveralp=tp(ts==tOverlap);
            xOverlap=linspace(min(tpOveralp),max(tpOveralp),nbins);
            % find max
            hMax=[];
            for k       = nTspp:(nTspp+1)
                
                figure(f1), hold on;
                h       = histogram(Tp(k,:),xOverlap);
                hMax    = [hMax; max(h.Values)];
                
            end
            for k       = nTspp:(nTspp+1)
                
                figure(f1), hold on;
                h       = histogram(Tp(k,:),xOverlap);
                Xval    = h.BinEdges + h.BinWidth*0.5+h.BinWidth*(k-(nTspp-1)-(nPr+1)/2)/2;
                Xval    = Xval(1:end-1);
                Yval    = tOverlap+h.Values./(max(hMax))*(1-pRand)*dTs; % h.Values*10 + TsA(k); % normalized so that maximum is set to 0.8*betweenBinDistance
                
                delete(h);
                
                figure(f2), hold on
                %                 h2      = stairs(Yval,Xval,'color',clr(k,:)); %
                hBarH=barh(Xval,Yval,'FaceColor',clr(k,:),'EdgeColor','none','BaseValue',tOverlap,'BarWidth',barWidth);
                hBarH.BaseLine.Visible='off';
                
            end
        end % if tOverlap~=0
        
        
        
  end % plotInd
  
  if PlotBLS
      if idStim==1
          clr     = Colors2P(6);
      else
          clr     = Colors2P(2); % HS 8/30/17; previously 2
      end
      for p = 1:2
          if p==1
              if PlotDiff
                  plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range)-theory{p}.tm(theory{p}.range),'-','linewidth',lw,'color',clr(1+nTspp*(p-1),:))
              else
                  plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range),'-','linewidth',lw,'color',clr(1+nTspp*(p-1),:))
              end
          else
              if PlotDiff
                  plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range)-theory{p}.tm(theory{p}.range),'-','linewidth',lw,'color',clr(1+nTspp*(p-1),:))
              else
                  plot(theory{p}.tm(theory{p}.range),theory{p}.te(theory{p}.range),'-','linewidth',lw,'color',clr(1+nTspp*(p-1),:))
              end
          end
      end
      
  end
  
  
  % Plotting Averages
  if idStim==1
      clr     = Colors2P(6);
  else
      clr     = Colors2P(2); % HS 8/30/17; previously 2
  end
  TsA     = nanmean(Ts,2);
  TpA     = nanmean(Tp,2);
  for k = 1:K, TpSD(k)= nanstd(Tp(k,:)); % ./sqrt(R);
  end
  
  % for each prior
  
  %                 hTmp=shadedErrorBar(TsA(1:nTspp),TpA(1:nTspp),TpSD(1:nTspp),{'--','color','r','linewidth',2},1); drawnow; % last input for transparent
  %                 delete(hTmp.mainLine); % delete mean line
  %                 hTmp=shadedErrorBar(TsA((1+nTspp):K),TpA((1+nTspp):K),TpSD((1+nTspp):K),{'--','color','b','linewidth',2},1); drawnow; % last input for transparent
  %                 delete(hTmp.mainLine);
  
  %                   errorbar(TsA,TpA,TpSD,'k','LineStyle','none')
  for k   = 1:K
      if PlotDiff
          plot(TsA(k),TpA(k)-TsA(k),'o','MarkerFaceColor','w','MarkerSize',markerSize,'linewidth',lw,'color',clr(k,:))
      else
          plot(TsA(k),TpA(k),'o','MarkerFaceColor','w','MarkerSize',markerSize,'linewidth',lw,'color',clr(k,:))
      end
  end
      
                  
                  
                  
                  
      
      
      
  end % End of Dot_or_Hist

% figure(f1); 
close(f1);

