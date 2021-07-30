%% Prepare Space
clear; clc;


%% Define Unit
s = sqrt(3);
Xs = [-s/2 0 s/2;...
      -1/2 1 -1/2];
Xu = [-s/2 -s/2; sqrt(s^2-s^2/4) -sqrt(s^2-s^2/4)];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
l = sqrt(sum((Xs(:,conn(:,1)) - Xu(:,conn(:,2)-3)).^2));


%% Draw
fig = figure(1); clf;
fSize = [8.5 11];

set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','inches');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','inches');

d_marg = .25;
c = [0 cumsum(l+2*d_marg)]';
cN = (c/max(c)*10.5 + .25)/11;

for i = 1:length(cN)
    annotation('line',[0 1], [1 1]*cN(i));
    if(i>1)
        annotation('line',[0 1], [1 1]*cN(i)-(d_marg/max(c)*8)/8.5,'linestyle','--');
    end
    if(i<length(cN))
        annotation('line',[0 1], [1 1]*cN(i)+(d_marg/max(c)*8)/8.5,'linestyle','--');
    end
end

for i = 1:7
    annotation('line',[1 1]*i/8, [0 1]);
end


%% Save
% fName = 'figureq5';
% set(gcf, 'Renderer', 'painters');
% fig.PaperPositionMode = 'manual';
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 fSize];
% fig.PaperSize = fSize;
% saveas(fig, ['Figures/' fName], 'pdf');