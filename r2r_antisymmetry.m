% Supplementary Figure: edge lengths
%% Prepare Space
clear; clc;
params_fig;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [8,3];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  3.00  3.00];...
        [ 5.00  0.00  3.00  3.00]];
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');


%% Quadrifolium
% Parameters
sc = 0.6;
s = sqrt(3);
Xs = [-s/2  0   s/2;...
      -1/2  1  -1/2];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
Le = 1.5;
th4 = 90; th5 = 200;
thdiff = 360-th4-th5;
th4a = th4 + thdiff;
th5a = th5 + thdiff;
Xu  = Le*[sind([th4 th5]); cosd([th4 th5])];
Xu2 = Le*[sind([th4a th5a]); cosd([th4a th5a])];

% Plot network 1
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
visualize_network(Xs*sc,Xu*sc,conn,'ucolor',[1 1 1]*.2);
line_coordinates(Xs(:,1:2)*sc,'lSh',.2,'nw',.04,'lw',.5,'style','-');
line_coordinates(Xs(:,2:3)*sc,'lSh',.2,'nw',.04,'lw',.5,'style','-');
text(0.35,subp(pInd,4)+labY,'original unit',NVTitle{:});
text(-.45,.55,'$l_1$',NVTexth{:});
text(.45,.55,'$l_2$',NVTexth{:});
axis([-sRat(pInd) sRat(pInd) -1 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;

% Plot network 2
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
visualize_network(Xs*sc,Xu2*sc,conn,'ucolor',[1 1 1]*.2);
line_coordinates(Xs(:,1:2)*sc,'lSh',.2,'nw',.04,'lw',.5,'style','-');
line_coordinates(Xs(:,2:3)*sc,'lSh',.2,'nw',.04,'lw',.5,'style','-');
text(0.3,subp(pInd,4)+labY,'mirrored unit',NVTitle{:});
text(-.45,.55,'$l_1$',NVTexth{:});
text(.45,.55,'$l_2$',NVTexth{:});
axis([-sRat(pInd) sRat(pInd) -1 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Save
fName = 'r2r_antisymmetry';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');