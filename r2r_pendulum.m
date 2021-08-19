% Supplementary Figure: edge lengths
%% Prepare Space
clear; clc;
params_fig;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [3 5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  3.00  6.00]];
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
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
CB = [100 100 255]/255;
CO = [200 100 000]/255;

x0 = [0.05, 0.8];
xf = [0.35, 0.3];
dx = (xf - x0) / norm(xf - x0);
dxp = [dx(2) -dx(1)];
th = atan2d(xf(2)-x0(2), xf(1)-x0(1));
dxa = [0, -0.25];
xa = [xf; xf+dxa];
xb = [xf; xf + (dxa*dx') * dx];
xc = [xf; xf + (dxa*dxp') * dxp];
xbo = [xf; xf - (dxa*dx') * dx];
xboo = [x0; x0 + (dxa*dx') * dx];
plot([x0(1) xf(1)], [x0(2) xf(2)],'k-','linewidth',1,'clipping',0);
plot(x0(1),x0(2),'ko','linewidth',2,'markersize',2);
plot(xf(1),xf(2),'ko','linewidth',6,'markersize',6);
arrow(xa(:,1), xa(:,2), sRat(pInd),'color',CB,'linewidth',1.5,'HeadWidth',5,'HeadLength',5);
arrow(xb(:,1), xb(:,2), sRat(pInd),'color',CB);
arrow(xc(:,1), xc(:,2), sRat(pInd),'color',CB);
arrow(xbo(:,1), xbo(:,2), sRat(pInd),'color',CO);
arrow(xboo(:,1), xboo(:,2), sRat(pInd),'color',CO);


text(x0(1)+.03,x0(2)+.02,'pin',NVTextl{:});
text(xf(1)-.06,xf(2)+.00,'mass $m$',NVTextr{:});
text(xa(2,1)-.02,xa(2,2)+.05,'$mg$',NVTextr{:},'color',CB);
text((xf(1)+x0(1))/2+.03, (xf(2)+x0(2))/2, '$F_{tension}$', NVTextl{:},'color',CO);


axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Save
fName = 'r2r_pendulum';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');