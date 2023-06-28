 

function Clo = Colors2P(flag)


switch flag
    
    
    case 1
    
Clo(2,1,1:3) = rgb('RoyalBlue');
  Clo(2,2,1:3) = rgb('CornFlowerBlue');
   Clo(2,3,1:3) = rgb('DodgerBlue');
    Clo(2,4,1:3) = rgb('DeepSkyBlue');
     Clo(2,5,1:3) = rgb('SkyBlue');
     
     
     
      Clo(1,1,1:3) = rgb('FireBrick'); % short
       Clo(1,2,1:3) = rgb('OrangeRed');
        Clo(1,3,1:3) = rgb('Tomato');
         Clo(1,4,1:3) = rgb('DarkOrange');
          Clo(1,5,1:3) = rgb('Orange');
          
          
          
    case 2 
        
  Clo(6,1:3) = rgb('RoyalBlue');
  Clo(7,1:3) = rgb('CornFlowerBlue');
   Clo(8,1:3) = rgb('DodgerBlue');
    Clo(9,1:3) = rgb('DeepSkyBlue');
     Clo(10,1:3) = rgb('SkyBlue');
     
      Clo(1,1:3) = rgb('FireBrick');
       Clo(2,1:3) = rgb('OrangeRed');
        Clo(3,1:3) = rgb('Tomato');
         Clo(4,1:3) = rgb('DarkOrange');
          Clo(5,1:3) = rgb('Orange');
          
          %% flipped version of 1&2 within each prior to maximize difference at overlap
              case 3
    
Clo(1,5,1:3) = rgb('RoyalBlue');
  Clo(1,4,1:3) = rgb('CornFlowerBlue');
   Clo(1,3,1:3) = rgb('DodgerBlue');
    Clo(1,2,1:3) = rgb('DeepSkyBlue');
     Clo(1,1,1:3) = rgb('Turquoise');
     
     
     
      Clo(2,5,1:3) = rgb('YellowGreen');
       Clo(2,4,1:3) = rgb('Gold');
        Clo(2,3,1:3) = rgb('Orange');
         Clo(2,2,1:3) = rgb('DarkOrange');
          Clo(2,1,1:3) = rgb('Coral');
          
          
          
    case 4
        
  Clo(5,1:3) = rgb('RoyalBlue');
  Clo(4,1:3) = rgb('CornFlowerBlue');
   Clo(3,1:3) = rgb('DodgerBlue');
    Clo(2,1:3) = rgb('DeepSkyBlue');
     Clo(1,1:3) = rgb('Turquoise');
     
      Clo(10,1:3) = rgb('YellowGreen');
       Clo(9,1:3) = rgb('Gold');
        Clo(8,1:3) = rgb('Orange');
         Clo(7,1:3) = rgb('DarkOrange');
          Clo(6,1:3) = rgb('Coral');
          
                
                        case 5 % backup for hand % for stim
    % short
    Clo(1,5,1:3) = rgb('magenta'); % overlap
  Clo(1,4,1:3) = rgb('violet');
   Clo(1,3,1:3) = rgb('Darkviolet');
    Clo(1,2,1:3) = rgb('Darkmagenta');
     Clo(1,1,1:3) = rgb('Indigo');
     
     
     
      Clo(2,5,1:3) = rgb('cyan');
       Clo(2,4,1:3) = rgb('aquamarine');
        Clo(2,3,1:3) = rgb('mediumaquamarine');
         Clo(2,2,1:3) = rgb('mediumseagreen');
          Clo(2,1,1:3) = rgb('seagreen'); % overlap
          
% Clo(1,5,1:3) = rgb('Olive'); % overlap
%   Clo(1,4,1:3) = rgb('DarkOliveGreen');
%    Clo(1,3,1:3) = rgb('ForestGreen');
%     Clo(1,2,1:3) = rgb('Green');
%      Clo(1,1,1:3) = rgb('DarkGreen');
%      
%      
%      
%       Clo(2,5,1:3) = rgb('Purple');
%        Clo(2,4,1:3) = rgb('Indigo');
%         Clo(2,3,1:3) = rgb('DarkMagenta');
%          Clo(2,2,1:3) = rgb('DarkOrchid');
%           Clo(2,1,1:3) = rgb('Magenta'); % overlap
          
          
          
    case 6
        
          Clo(5,1:3) = rgb('magenta');
  Clo(4,1:3) = rgb('violet');
   Clo(3,1:3) = rgb('Darkviolet');
    Clo(2,1:3) = rgb('Darkmagenta');
     Clo(1,1:3) = rgb('Indigo');
     
      Clo(10,1:3) = rgb('cyan');
       Clo(9,1:3) = rgb('aquamarine');
        Clo(8,1:3) = rgb('mediumaquamarine');
         Clo(7,1:3) = rgb('mediumseagreen');
          Clo(6,1:3) = rgb('seagreen');
          
%   Clo(5,1:3) = rgb('Olive');
%   Clo(4,1:3) = rgb('DarkOliveGreen');
%    Clo(3,1:3) = rgb('ForestGreen');
%     Clo(2,1:3) = rgb('Green');
%      Clo(1,1:3) = rgb('DarkGreen');
%      
%       Clo(10,1:3) = rgb('Purple');
%        Clo(9,1:3) = rgb('Indigo');
%         Clo(8,1:3) = rgb('DarkMagenta');
%          Clo(7,1:3) = rgb('DarkOrchid');
%           Clo(6,1:3) = rgb('Magenta');
          
end

return;

%% for tmpCmap in pplot.mat
% short
tmpCmap{1,2} = flipud([rgb('magenta');... % 800
    rgb('violet');...
    rgb('Darkviolet');...
    rgb('Darkmagenta');...
    rgb('Indigo')]);

tmpCmap{2,2} = flipud([rgb('cyan');... % 1200
rgb('aquamarine');...
rgb('mediumaquamarine');...
rgb('mediumseagreen');...
rgb('seagreen')]);