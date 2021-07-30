% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
fig = figure(2); clf;
params_fig;
fig_motivation_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [8.6 13.00];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 9.75 8.60 3.25];...
        [ 0.00 6.50 8.60 3.25];...
        [ 0.00 3.25 8.60 3.25];...
        [ 0.00 0.00 8.60 3.25]];
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


%% a: Nodes and edges
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .05;               % Scale nodes
sh = [-.03;.08];        % Shift nodes
scE = .08;              % Scale edges
shE = [1.35;.1];       % Shift edges

% Nodes
Nrows = 13;
Ncols = 26;
[Xcoor,Ycoor] = meshgrid(1:Ncols,1:Nrows);
XN = [Xcoor(:) Ycoor(:)]';
visualize_network(XN*sc+sh,[],[1 1],'msize',2);

% Edges
Erows = 32;
Ecols = 21;
[XcoorE,YcoorE] = meshgrid((1:Ecols)*.9,(1:Erows)*.25);
XE = [XcoorE(:) YcoorE(:)]';
line([-1;1]*.025+XE(1,:)*scE+shE(1), [1;1].*XE(2,:)*scE+shE(2),...
     'color','k','linewidth',.5,'clipping',0)

% Text
texta = '\textbf{a}\hspace{2mm}how do we combine many nodes and edges...';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(.22,.83,'$N$ nodes',NVTextH{:});
text(0.75,.83,'$2N-4$ edges',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: Clover conformational motion
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;


% Parameters
sc = .028;              % Scale drawing
sh = [2.84;0.15];       % Shift drawing
Rz = rotz(45); Rz = Rz(1:2,1:2);

% start conformation
visualize_network(XCcc(:,1:size(Xscc,2),1)*sc+sh,...
                  XCcc(:,(1:size(Xucc,2))+size(Xscc,2),1)*sc+sh,...
                  conncc,'msize',2,'ucolor',CSSc,'lalpha',.1,...
                  'ucolor',CSSc.^(.05));


% end conformation
sh = [2.79; 0.18];
visualize_network(Rz*XCcc(:,1:size(Xscc,2),end)*sc+sh,...
                  Rz*XCcc(:,(1:size(Xucc,2))+size(Xscc,2),end)*sc+sh,...
                  conncc,'msize',2,'ucolor',CSSc);

% Arrow
sCo1 = [0.00 0.25 0.60 1.00;...
        0.00 0.5 0.90 1.00] .* [1.6;.35] + [.4;.3];
plot_spline(sCo1,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'color',o*gr,'linewidth',.5);

% Text
textb = '...to design precise shape changes?';
text(labX,subp(pInd,4)+labY-.1,'\textbf{b}',NVTitle{:});
text(fSize(1)-fMarg(1),subp(pInd,4)+labY-.1,textb,NVTitleR{:});
text(.12,.3,'start',NVTextR{:},'color',o*gr);
text(.68,.75,'end',NVTextR{:},'color',o*gr);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: Clover conformational motion
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;


% Parameters
sc = .028;              % Scale drawing
sh = [2.84;0.15];       % Shift drawing

% Draw networks
pI = [30 80 140 235];
pIL = [200 200 200 size(conncc,1)];
for i = 1:length(pI)
    la = .1 + .1*i;
    visualize_network(XCcc(:,unique(conncc(1:pIL(i),1)),pI(i))*sc+sh,...
                      XCcc(:,unique(conncc(1:pIL(i),2)),pI(i))*sc+sh,...
                      conncc(1:pIL(i),:) + [0 max(conncc(1:pIL(i)))-size(Xscc,2)],...
                      'lalpha',la,'msize',2,'ucolor',CSSc(1:length(unique(conncc(1:pIL(i),2))),:).^(.05*i));
end

% Arrow
sCo1 = [0.00 0.25 0.70  1.00;...
        0.00 1.00 1.00 -0.40] .* [1.55;.35] + [.0;.53];
plot_spline(sCo1,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'color',o*gr,'linewidth',.5);

% Text
textc = '...to design folding sequences?';
text(labX,subp(pInd,4)+labY-.1,'\textbf{c}',NVTitle{:});
text(fSize(1)-fMarg(1),subp(pInd,4)+labY-.1,textc,NVTitleR{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% d: Nonlinear dynamics
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;


% Parameters
ssh = 0.65;
s = sqrt(3);
sc = .04; shn = 1.4;
sh1 = [.68; .51];
sh2 = sh1 + [1;-ssh]/2*shn;
sh3 = sh1 + [1;0]*shn;
sh4 = sh1 + [3;-ssh]/2*shn;
sh0a = sh2 + [-.87;0];
sh0b = sh2 + [-.87;0]+[-.8;ssh]*.05;
sP = sRat(pInd);

% Draw networks
% Networks
visualize_network(XP0(:,1:size(Xsa4,2))*sc+sh1,...
                  XP0(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh1,conna4,...
                  'msize',2,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(XP1(:,1:size(Xsa4,2))*sc+sh2,...
                  XP1(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh2,conna4,...
                  'msize',2,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(XP2(:,1:size(Xsa4,2))*sc+sh3,...
                  XP2(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh3,conna4,...
                  'msize',2,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(XP3(:,1:size(Xsa4,2))*sc+sh4,...
                  XP3(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh4,conna4,...
                  'msize',2,'ucolor',[1 1 1]*.5,'ucolor',C3a);

% Arrows
s2 = sqrt(2);
arrow([0 .4]+1,[0 0]+.7,sRat(pInd),'linewidth',0.5,'color',o*gr);
arrow([0 .4]+1.5,[0 0]+.05,sRat(pInd),'linewidth',0.5,'color',o*gr);
arrow([0 1/10]+0.75,[0 -s/10]+.35,sRat(pInd),'linewidth',0.5,'color',o*gr);
arrow([0 1/10]+2.15,[0 -s/10]+.35,sRat(pInd),'linewidth',0.5,'color',o*gr);

% Text
textd = '...to discover exotic and nonlinear behavior?';
text(labX,subp(pInd,4)+labY-.1,'\textbf{d}',NVTitle{:});
text(fSize(1)-fMarg(1),subp(pInd,4)+labY-.1,textd,NVTitleR{:});
text(0.00,0.47,'left',NVTextL{:},'color',o*gr);
text(0.06,0.85,'top',NVTextL{:},'color',o*gr);
text(0.42,0.8,'open left',NVTextH{:},'color',o*gr);
text(0.59,0.15,'open left',NVTextH{:},'color',o*gr);
text(0.26,0.16,'open top',NVTextH{:},'color',o*gr,'rotation',-60);
text(0.79,0.32,'open top',NVTextH{:},'color',o*gr,'rotation',-60);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'fig_motivationa';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');