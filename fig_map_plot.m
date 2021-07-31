% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
fig = figure(3); clf;
params_fig;
fig_map_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [8.6 7.0];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .4 .0];
% Subplot position in cm [x,y,w,h]
subp = [ 0.00 -0.20  8.00  8.00];
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');


%% single module conformational motion
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .065;          % Scale drawing
aX0 = 1.28;         % Map axes start
aXF = 2.04;         % Map axes end
% Colors
CPI = interp1(DLin,CP1,D(1,:));
CPJ = interp1(DLin,CP1,D(2,:));

% Draw map ticks
line(sqrt([2 2]), [0 .04]+aX0, 'color', interp1(DLin,CP1,sqrt(2)),'linewidth', 1);
line([2 2], [0 .04]+aX0, 'color', interp1(DLin,CP1,2),'linewidth', 1);
line([0 .04]+aX0, sqrt([2 2]), 'color', interp1(DLin,CP1,sqrt(2)),'linewidth', 1);
line([0 .04]+aX0, [2 2], 'color', interp1(DLin,CP1,2),'linewidth', 1);

% Draw distance colormap
scatter(D(1,:),ones(1,size(D,2))*aX0+.012,40,interp1(DLin,CP1,D(1,:)),'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.012,D(2,:),40,interp1(DLin,CP1,D(2,:)),'filled','marker','s');
 
% Draw map axes
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF+.15], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
    
% Draw map and y=x lines
plot([aX0 2.0],[aX0 2.0], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);

% Draw spline
sCo11 = [0.00 1.00;...
         0.20 0.90] .* [.15;-.2] + [1.8;2.0];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5);

% Draw example networks
plInd = ceil([1 27 52 size(XC,3)]);
plSh = [ .18  .17 .18 .2;...
        -.06 -.06 -.06 -.045];
lSh = -.007;
for i = 1:length(plInd)
    pI = plInd(i);
    XP = XC(:,:,pI)*sc + [D(1,pI);D(2,pI)]+plSh(:,i);
    line_coordinates(XP(:,1:2), 'lSh',lSh, 'lw',1, 'color',CPI(pI,:));
    line_coordinates(XP(:,3:4), 'lSh',-lSh, 'lw',1, 'color',CPJ(pI,:));
    visualize_network(XP,[],conn);
end
scatter(D(1,plInd(2:3)),D(2,plInd(2:3)),40,zeros(2,3),...
       'filled','marker','s','linewidth',.5);
scatter(D(1,plInd([1 4])),D(2,plInd([1 4])),80,...
       [interp1(DLin,CP1,D(1,plInd(1)));interp1(DLin,CP1,D(1,plInd(4)))],...
       'filled','marker','s','linewidth',.5);
scatter(D(1,plInd([1 4])),D(2,plInd([1 4])),40,ones(2,3),'marker','s','linewidth',.1);

% Draw large template network
XP = XC(:,:,plInd(3))*sc*2 + [D(1,pI);D(2,pI)]+[.32;.51];
visualize_network(XP,[],conn,'lalpha',.3,'scolor',o*gr^.3,'ms',ms);
for i = 1:4
    text(XP(1,i), XP(2,i)-.002, num2str(i), NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
end
line_coordinates(XP(:,1:2), 'lSh',-.035, 'lw',.5,'style','-','nw',.01);
line_coordinates(XP(:,3:4), 'lSh',0.035, 'lw',.5,'style','-','nw',.01);

% Text
% Axis legend
text(sqrt(2),aX0-.04,'$\sqrt{2}$',NVTexth{:},'color',interp1(DLin,CP1,sqrt(2)));
text(2,aX0-.04,'$2$',NVTexth{:},'color',interp1(DLin,CP1,2));
text(aX0-.037,sqrt(2),'$\sqrt{2}$',NVTexth{:},'color',interp1(DLin,CP1,sqrt(2)));
text(aX0-.037,2,'$2$',NVTexth{:},'color',interp1(DLin,CP1,2));
% axis labels
text(1.01,0.04,'$l_k$',NVTextL{:});
text(0,.87,'$l_{k+1}$',NVTextL{:});
% template labels
text(XP(1,3),XP(2,3)+.15,'unit $k$',NVTextl{:},'color',o*gr);
text(.41,.79,'$l_{k+1} = f(l_k)$',NVTextL{:});
text(.145,0.71,'$l_k$',NVTextL{:});

% Axes
axis([1.24 2.2 1.24 2.2]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% Size and Save Figure
fName = 'fig_map2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');