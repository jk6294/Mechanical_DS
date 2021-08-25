%% Prepare Space
% Conformational Motions and Degrees of Freedom
clear; clc;
fig = figure(1); clf;
params_fig;
fig_dof_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [8.6 6.50];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 3.25 3.50 3.25];...
        [ 0.00 0.00 3.50 3.25];....
        [ 4.00 3.25 4.60 3.25];...
        [ 4.00 0.00 4.60 3.25];...
        [ 3.10 0.00 1.30 6.50]];
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


%% a: building blocks
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sh = [.05;.04];            % Shift drawing

% Nodes
visualize_network([.3;.7]+sh,[],[1 1]);
visualize_network([.55;.7]+sh,[],[1 1],'nalpha',al,'lalpha',al);
visualize_network([.3;.45]+sh,[],[1 1],'nalpha',al,'lalpha',al);
arrow([.35 .5]+sh(1), [.7 .7]+sh(2), sRat(pInd), 'color', o*gr^2);
arrow([.3 .3]+sh(1), [.66 .50]+sh(2), sRat(pInd), 'color', o*gr^2);

% Edge
arrow([.64 .8]+sh(1), [.12 .12]+sh(2), sRat(pInd), 'color', o*gr^2);
arrow([.26 .10]+sh(1), [.12 .12]+sh(2), sRat(pInd), 'color', o*gr^2);
visualize_network([.3 .6; .12 .12]+sh,[],[1 2]);
plot(.7+sh(1),.12+sh(2),'rx','linewidth',1,'markersize',5);
plot(.2+sh(1),.12+sh(2),'rx','linewidth',1,'markersize',5);

% Text
texta = '\textbf{a}\hspace{4.5mm}building blocks';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(.45,.85,'node: $+2$ motions',NVTextH{:});
text(.45,.30,'edge: $-1$ motion',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: combining blocks
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .2;                 % Scale drawing
sh = [.32;.62];          % Shift drawing

% Draw network
visualize_network(Xf*sc+sh,[],conn,'ms',ms);

% Text
textb = '\textbf{b}\hspace{4mm}combine blocks';
text(labX,subp(pInd,4)+labY-.1,textb,NVTitle{:});
% Coordinates
text(Xf(1,1)*sc+sh(1)-.06,Xf(2,1)*sc+sh(2),'$(0,~1)$','fontsize',FS2,...
     'horizontalalignment','right','color',o*gr^3);
text(Xf(1,2)*sc+sh(1)-.06,Xf(2,2)*sc+sh(2),'$(0,-1)$','fontsize',FS2,...
     'horizontalalignment','right','color',o*gr^3);
text(Xf(1,3)*sc+sh(1)-.06,Xf(2,3)*sc+sh(2),'$(0,~0)$','fontsize',FS2,...
     'horizontalalignment','right','color',o*gr^3);
text(Xf(1,4)*sc+sh(1)+.06,Xf(2,4)*sc+sh(2),'$(0,~2)$','fontsize',FS2,...
     'horizontalalignment','left','color',o*gr^3);
% Node labels
for i = 1:4
    text(Xf(1,i)*sc+sh(1),Xf(2,i)*sc+sh(2)-.005,num2str(i),NVTexth{:},'fontsize',FS2);
end
% DOF equation
text(-.05,.28,'$4$ nodes:',NVTextL{:});
text(-.05,.16,'$4$ edges:',NVTextL{:});
text(.95,.28,'$+8$ motions',NVTextR{:});
text(.95,.16,'$-4$ motions',NVTextR{:});
text(.95,.0,'$4$ motions',NVTextR{:});
line([0 .65]+.4,[0 0]+.08,'clipping',0,'color','k','linewidth',.5);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: rigid body motions
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .15;               % Scale drawing
sh1 = [0.02;0.3];       % Shift drawing
sh2 = [0.6;0.3];
sh3 = [1.16;0.3];
Rz = rotz(20); Rz = Rz(1:2,1:2);

% Draw network: x shift
arrow([.05 .3], [.6  .6], sRat(pInd), 'color', o*gr.^2);
visualize_network(Xf*sc+sh1+[.1;0],[],conn,'scolor',o*.98,'lalpha',al);
visualize_network(Xf*sc+sh1,[],conn);

% Draw network: y shift
arrow([0 0]+.76, [.7  .5], sRat(pInd), 'color', o*gr.^2);
visualize_network(Xf*sc+sh2+[0;-.07],[],conn,'scolor',o*.98,'lalpha',al);
visualize_network(Xf*sc+sh2,[],conn);

% Draw network: rotation
th = linspace(0,3*pi/2,100);
plot(.1*cos(th)+1.32,.1*sin(th)+.6,'-','color',o*gr.^2,'linewidth',.5);
arrow([1.34 1.35], [.5 .5], sRat(pInd), 'color', o*gr.^2);
visualize_network(Rz*Xf*sc+sh3,[],conn,'scolor',o*.98,'lalpha',al);
visualize_network(Xf*sc+sh3,[],conn);

% Text
textc = '\textbf{c}\hspace{7mm}$3$ rigid body motions';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});
text(0.00,.83,'x shift',NVTextL{:});
text(0.51,.83,'y shift',NVTextH{:});
text(1.04,.83,'rotation',NVTextR{:});

% Axes limits
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% d: conformational motion
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .15;                   % Scale drawing
sh = [0.21 0.68 1.15;...    % Shift drawing
      0.30 0.30 0.30];

% Draw networks
nPl = [1 38 size(XC,3)];
for i = 1:length(nPl)
    line_coordinates(XC(:,1:2,nPl(i))*sc+sh(:,i),'lSh',-0.015,...
                     'color',interp1(DLin,CP1,D(1,nPl(i))));
    line_coordinates(XC(:,3:4,nPl(i))*sc+sh(:,i),'lSh', 0.015,...
                     'color',interp1(DLin,CP1,D(2,nPl(i))));
    visualize_network(XC(:,:,nPl(i))*sc+sh(:,i),[],conn);
end
             
% Distance
scatter(linspace(.15,sRat(pInd)-.15,nT), .64*ones(1,nT),15,CP1,'s','filled');

% Text
textd = '\textbf{d}\hspace{3.5mm}$1$ conformational motion';
text(labX,subp(pInd,4)+labY-.1,textd,NVTitle{:});
text(.5,.73,'length',NVTextH{:});
text(.1,.73,'$\sqrt{2}$',NVTextH{:},'color',CP1(1,:));
text(.9,.73,'$2$',NVTextH{:},'color',CP1(end,:));
text(-.01,.3,'$l_1$',NVTextL{:});
text(.18,.43,'$l_2$',NVTextL{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Arrow
pInd = 5;
subplot('position',subpN(pInd,:)); cla; hold on;

% Spline
sCo11 = [0.00 0.35 0.60 1.00;...
         0.00 0.04 0.80 1.00] .* [.125;.98] + [-.01;0];
sCo12 = [0.00 0.40 0.7 1.00;...
         0.00 0.08 0.8 1.00] .* [.125;.43] + [-.01;0];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'color',o*gr,'linewidth',.5);
plot_spline(sCo12,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'color',o*gr,'linewidth',.5);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'fig_dof2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');