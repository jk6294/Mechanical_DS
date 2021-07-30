%% Prepare Space
clear; clc;
fig = figure(6); clf;
params_fig;
fig_design_motivate_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [8.6 10.2];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .0 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  7.80  8.60  2.40];...
        [ 0.00  5.10  8.60  2.40];...
        [ 0.00  2.30  8.60  2.40];...
        [ 0.00 -0.50  8.60  2.40]];
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


%% Unit
fig = figure(6); clf;
% Parameters
sc = .3;
lSh1 = .08;
nw = .015;
lw = .5;

% Choose subfigure
pInd = 1;
sh = [sRat(pInd)/2-.08;.93];
subplot('position',subpN(pInd,:)); cla; hold on;
% Lengths
line_coordinates(Xs(:,1:2)*sc+sh,'lSh',lSh1,'nw',nw,'style','-','lw',lw);
line_coordinates(Xs(:,2:3)*sc+sh,'lSh',lSh1,'nw',nw,'style','-','lw',lw);
line_coordinates(Xs(:,[1 3])*sc+sh,'lSh',-lSh1,'nw',nw,'style','-','lw',lw);
% Draw network
visualize_network(Xs*sc+sh,[],[1 1]);
% Text
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'unit $k$',NVTitleH{:});
text(Xs(1,1)*sc+sh(1)-.09,Xs(2,1)*sc+sh(2)-.06,'$1$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,2)*sc+sh(1),Xs(2,2)*sc+sh(2)-.14,'$2$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,3)*sc+sh(1)+.09,Xs(2,3)*sc+sh(2)-.06,'$3$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,1)*sc+sh(1),Xs(2,1)*sc+sh(2)+.3,'$l_k$',NVTextr{:});
text(Xs(1,3)*sc+sh(1),Xs(2,3)*sc+sh(2)+.3,'$l_{k+1}$',NVTextl{:});
text(Xs(1,2)*sc+sh(1),Xs(2,1)*sc+sh(2)-.15,'$c_k$',NVTexth{:});
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% a: unit, start, and end
% Parameters
sh1 = [.45;1.12];
sh2 = [3.6;1.12];
lSh2 = .015;
tsh = 0;

% Choose subfigure
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;

% Draw start conformation
% Lengths
line_coordinates(Xs(:,1:2)*sc+sh1,'lSh',lSh2,'lw',1,'color',C1a);
line_coordinates(Xs(:,2:3)*sc+sh1,'lSh',lSh2,'lw',1,'color',C1a);
% Draw network
visualize_network(Xs*sc+sh1,[],[1 1]);
% Text
text(labX,subp(pInd,4)+labY,'\textbf{a}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-tsh,'\underline{~property 1~}',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-.5-tsh,'start and end at',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-.9-tsh,'fixed points',NVTitleH{:});
text(Xs(1,1)*sc+sh1(1)-.09,Xs(2,1)*sc+sh1(2)-.06,'$1$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,2)*sc+sh1(1),Xs(2,2)*sc+sh1(2)-.14,'$2$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,3)*sc+sh1(1)+.09,Xs(2,3)*sc+sh1(2)-.06,'$3$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,1)*sc+sh1(1),Xs(2,1)*sc+sh1(2)+.25,'$l^\bullet$',NVTexth{:},'color',C1a);
text(Xs(1,3)*sc+sh1(1)+.04,Xs(2,3)*sc+sh1(2)+.25,'$l^\bullet$',NVTexth{:},'color',C1a);

% Draw end conformation
% Lengths
line_coordinates(Xf(:,1:2)*sc+sh2,'lSh',lSh2,'lw',1,'color',C1c);
line_coordinates(Xf(:,2:3)*sc+sh2,'lSh',lSh2,'lw',1,'color',C1c);
% Draw network
visualize_network(Xf*sc+sh2,[],[1 1]);
% Text
text(Xf(1,1)*sc+sh2(1)-.09,Xf(2,1)*sc+sh2(2)-.06,'$1$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,2)*sc+sh2(1),Xf(2,2)*sc+sh2(2)-.14,'$2$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,3)*sc+sh2(1)+.09,Xf(2,3)*sc+sh2(2)-.06,'$3$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,1)*sc+sh2(1)+.1,Xf(2,1)*sc+sh2(2)+.4,'$l^\circ$',NVTexth{:},'color',C1c);
text(Xf(1,3)*sc+sh2(1)-.05,Xf(2,3)*sc+sh2(2)+.4,'$l^\circ$',NVTexth{:},'color',C1c);
text(sh1(1),1.3,'start',NVTexth{:});
text(sh2(1),1.3,'end',NVTexth{:});

% Distance colorbar
scatter(linspace(0,1,nT)+1.45, .4*ones(1,nT),15,CP1,'s','filled');
text(subp(pInd,3)/2-fMarg(1),.5,'$l$',NVTitleH{:});


% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: unit, start, end
% Colors
CP2s = interp1(DLin2,CP2,sqrt(sum(diff(Xs(:,[1,3]),1,2).^2)));
CP2e = interp1(DLin2,CP2,sqrt(sum(diff(Xf(:,[1,3]),1,2).^2)));

% Chooes subfigure
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;

% Draw start conformation
% Lengths
line_coordinates(Xs(:,[1 3])*sc+sh1,'lSh',-lSh2,'lw',1,'color',CP2s);
% Draw network
visualize_network(Xs*sc+sh1,[],[1 1]);
% Text
text(labX,subp(pInd,4)+labY,'\textbf{b}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-tsh,'\underline{~property 2~}',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-.5-tsh,'design unit shape',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-.9-tsh,'through $c$',NVTitleH{:});
text(Xs(1,1)*sc+sh1(1)-.09,Xs(2,1)*sc+sh1(2)-.06,'$1$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,2)*sc+sh1(1),Xs(2,2)*sc+sh1(2)-.14,'$2$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,3)*sc+sh1(1)+.09,Xs(2,3)*sc+sh1(2)-.06,'$3$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,2)*sc+sh1(1),Xs(2,3)*sc+sh1(2)-.13,'$c^\bullet$',NVTexth{:},'color',CP2s);

% Draw end conformation
% Lengths
line_coordinates(Xf(:,[1 3])*sc+sh2,'lSh',-lSh2,'lw',1,'color',CP2e);
% Draw network
visualize_network(Xf*sc+sh2,[],[1 1]);
% Text
text(Xf(1,1)*sc+sh2(1)-.09,Xf(2,1)*sc+sh2(2)-.06,'$1$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,2)*sc+sh2(1),Xf(2,2)*sc+sh2(2)-.14,'$2$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,3)*sc+sh2(1)+.09,Xf(2,3)*sc+sh2(2)-.06,'$3$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,2)*sc+sh2(1),Xf(2,3)*sc+sh2(2)-.13,'$c^\circ$',NVTexth{:},'color',CP2e);

% Distance colorbar
scatter(linspace(0,1,nT)+1.45, .4*ones(1,nT),15,CP2,'s','filled');
text(subp(pInd,3)/2-fMarg(1),.55,'$c$',NVTitleH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: Conformational degrees of freedom
% Choose subfigure
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;

% Draw start conformation
visualize_network(Xs*sc+sh1,Xu*sc+sh1,conn,'ucolor',CP2e);
% Text
text(labX,subp(pInd,4)+labY,'\textbf{c}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-tsh,'\underline{~property 3~}',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-.5-tsh,'move through one',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-.9-tsh,'conformational',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY-1.3-tsh,'motion',NVTitleH{:});
text(Xs(1,1)*sc+sh1(1)-.09,Xs(2,1)*sc+sh1(2)-.06,'$1$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,2)*sc+sh1(1)+.10,Xs(2,2)*sc+sh1(2),'$2$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xs(1,3)*sc+sh1(1)+.09,Xs(2,3)*sc+sh1(2)-.06,'$3$',NVTexth{:},'color',o*gr,'fontsize',8);
text(XC(1,4,1)*sc+sh1(1)-.09,XC(2,4,1)*sc+sh1(2),'$4$',NVTexth{:},'color',o*gr,'fontsize',8);
text(XC(1,5,1)*sc+sh1(1)-.09,XC(2,5,1)*sc+sh1(2),'$5$',NVTexth{:},'color',o*gr,'fontsize',8);

% Draw end conformation
visualize_network(XC(:,1:3,end)*sc+sh2,XC(:,4:5,end)*sc+sh2,conn,'ucolor',CP2e);
% Text
text(Xf(1,1)*sc+sh2(1)-.09,Xf(2,1)*sc+sh2(2)-.06,'$1$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,2)*sc+sh2(1)+.10,Xf(2,2)*sc+sh2(2),'$2$',NVTexth{:},'color',o*gr,'fontsize',8);
text(Xf(1,3)*sc+sh2(1)+.09,Xf(2,3)*sc+sh2(2)-.06,'$3$',NVTexth{:},'color',o*gr,'fontsize',8);
text(XC(1,4,end)*sc+sh2(1)-.09,XC(2,4,end)*sc+sh2(2),'$4$',NVTexth{:},'color',o*gr,'fontsize',8);
text(XC(1,5,end)*sc+sh2(1)-.05,XC(2,5,end)*sc+sh2(2)+.09,'$5$',NVTexth{:},'color',o*gr,'fontsize',8);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'fig_design_motivate1';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');