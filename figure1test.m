% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 6.50];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 3.25 3.50 3.25];...
        [ 0.00 0.00 3.50 3.25];....
        [ 4.50 3.25 5.00 3.25];...
        [ 4.50 0.00 5.00 3.25];...
        [10.00 3.25 9.00 3.25];...
        [10.00 0.00 9.00 3.25];...
        [ 3.10 0.00 1.80 6.50]];
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


%% Figure Parameters
FS = 10;            % Fontsize
gr = 0.9;           % Gray
al = 0.2;           % Alpha
    
% Default Name-Value Pairs
NVTitle  = {'Units','centimeters','fontsize',FS};
NVTextH  = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','right'};
NVTextR  = {'Units','Normalized','fontsize',FS};

% Color
nT = 1000;
o = [1 1 1];
% Gradient: Distance. Interpolate between 3 colors
DLin = linspace(sqrt(2)-.01,2+.01,nT);
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));


%% a: building blocks
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sh = [.1;0];            % Shift drawing

% Nodes
visualize_network([.3;.7]+sh,[],[1 1]);
visualize_network([.55;.7]+sh,[],[1 1],'nalpha',al,'lalpha',al);
visualize_network([.3;.45]+sh,[],[1 1],'nalpha',al,'lalpha',al);
arrow([.35 .5]+sh(1), [.7 .7]+sh(2), sRat(pInd), 'color', o*gr^2);
arrow([.3 .3]+sh(1), [.66 .50]+sh(2), sRat(pInd), 'color', o*gr^2);

% Edge
arrow([.64 .8]+sh(1), [.1 .1]+sh(2), sRat(pInd), 'color', o*gr^2);
arrow([.26 .10]+sh(1), [.1 .1]+sh(2), sRat(pInd), 'color', o*gr^2);
visualize_network([.3 .6; .1 .1]+sh,[],[1 2]);
plot(.7+sh(1),.1+sh(2),'rx','linewidth',1,'markersize',5);
plot(.2+sh(1),.1+sh(2),'rx','linewidth',1,'markersize',5);

% Text
texta = '\textbf{a}\hspace{4.5mm}building blocks';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(.5,.83,'node: $+2$ motions',NVTextH{:});
text(.5,.25,'edge: $-1$ motion',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: combining blocks
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .2;                % Scale drawing
sh = [.35;.62];          % Shift drawing

% Build network
L = 1;
sq = sqrt(2)*L/2;
Xs = [-sq -sq 0 sq*2;...
       sq -sq 0 0];
Xf = [0  0  0  2;...
      1 -1  0  0];
conn = [1 3; 1 4; 2 3; 2 4];

% Draw network
visualize_network(Xf*sc+sh,[],conn);

% Text
textb = '\textbf{b}\hspace{4mm}combine blocks';
text(labX,subp(pInd,4)+labY-.1,textb,NVTitle{:});
% Coordinates
text(Xf(1,1)*sc+sh(1)-.04,Xf(2,1)*sc+sh(2),'$(0,~1)$','fontsize',FS*3/4,...
     'horizontalalignment','right','color',o*gr^3);
text(Xf(1,2)*sc+sh(1)-.04,Xf(2,2)*sc+sh(2),'$(0,-1)$','fontsize',FS*3/4,...
     'horizontalalignment','right','color',o*gr^3);
text(Xf(1,3)*sc+sh(1)-.04,Xf(2,3)*sc+sh(2),'$(0,~0)$','fontsize',FS*3/4,...
     'horizontalalignment','right','color',o*gr^3);
text(Xf(1,4)*sc+sh(1)+.04,Xf(2,4)*sc+sh(2),'$(0,~2)$','fontsize',FS*3/4,...
     'horizontalalignment','left','color',o*gr^3);
% DOF equation
text(0,.28,'$4$ nodes:',NVTextR{:});
text(0,.16,'$4$ edges:',NVTextR{:});
text(1.0,.28,'$+8$ motions',NVTextRA{:});
text(1.0,.16,'$-4$ motions',NVTextRA{:});
text(1.0,.0,'$4$ motions',NVTextRA{:});
line([0 .65]+.45,[0 0]+.08,'clipping',0,'color','k','linewidth',.5);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: rigid body motions
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .15;               % Scale drawing
sh1 = [0.02;0.3];       % Shift drawing
sh2 = [0.66;0.3];
sh3 = [1.26;0.3];

% Build network
L = 1;
sq = sqrt(2)*L/2;
Rz = rotz(20); Rz = Rz(1:2,1:2);
conn = [1 3; 1 4; 2 3; 2 4];

% Draw network: x shift
arrow([.05 .3], [.6  .6], sRat(pInd), 'color', o*gr.^2);
visualize_network(Xf*sc+sh1+[.1;0],[],conn,'scolor',o*.98,'lalpha',al);
visualize_network(Xf*sc+sh1,[],conn);

% Draw network: y shift
arrow([0 0]+.8, [.7  .5], sRat(pInd), 'color', o*gr.^2);
visualize_network(Xf*sc+sh2+[0;-.07],[],conn,'scolor',o*.98,'lalpha',al);
visualize_network(Xf*sc+sh2,[],conn);

% Draw network: rotation
th = linspace(0,3*pi/2,100);
plot(.1*cos(th)+1.4,.1*sin(th)+.6,'-','color',o*gr.^2,'linewidth',.5);
arrow([1.4 1.43], [.5 .5], sRat(pInd), 'color', o*gr.^2);
visualize_network(Rz*Xf*sc+sh3,[],conn,'scolor',o*.98,'lalpha',al);
visualize_network(Xf*sc+sh3,[],conn);

% Text
textc = '\textbf{c}\hspace{7mm}$3$ rigid body motions';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});
text(0.0,.83,'x shift',NVTextR{:});
text(0.5,.83,'y shift',NVTextH{:});
text(1.0,.83,'rotation',NVTextRA{:});

% Axes limits
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% d: conformational motion
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .15;                   % Scale drawing
sh = [0.25 0.75 1.25;...    % Shift drawing
      0.30 0.30 0.30];

% Build network
L = 1;
sq = sqrt(2)*L/2;
Xs = [-sq -sq 0 sq*2;...
       sq -sq 0 0];
conn = [1 3; 1 4; 2 3; 2 4];

% Simulate network
[XC,fC] = sim_motion10(Xs,[],conn,.005,152,[0 0 0 0; 1 -1 0 0],0);
D = sqrt(squeeze(sum(diff(XC(:,:,:),1,2).^2))); D = D([1 3],:);
disp(['mean simulation error: ' num2str(mean(fC))]);

% Draw networks
nPl = [1 75 size(XC,3)];
for i = 1:length(nPl)
    line_coordinates(XC(:,1:2,nPl(i))*sc+sh(:,i),'lSh',-0.015,...
                     'color',interp1(DLin,CP1,D(1,nPl(i))));
    line_coordinates(XC(:,3:4,nPl(i))*sc+sh(:,i),'lSh', 0.015,...
                     'color',interp1(DLin,CP1,D(2,nPl(i))));
    visualize_network(XC(:,:,nPl(i))*sc+sh(:,i),[],conn);
end
             
% Distance
scatter(linspace(.15,sRat(pInd)-.15,nT), .65*ones(1,nT),5,CP1,'s','filled');

% Text
textd = '\textbf{d}\hspace{3.5mm}$1$ conformational motion';
text(labX,subp(pInd,4)+labY-.1,textd,NVTitle{:});
text(.5,.73,'length',NVTextH{:});
text(.1,.73,'$\sqrt{2}$',NVTextH{:});
text(.9,.73,'$2$',NVTextH{:});
text(0.02,.3,'$l_1$',NVTextR{:});
text(.18,.43,'$l_2$',NVTextR{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Arrow
pInd = 7;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Spline
sCo11 = [0.00 0.35 0.60 1.00;...
         0.00 0.04 0.92 1.00] .* [.18;.995] + [.01;0];
sCo12 = [0.00 0.40 0.7 1.00;...
         0.00 0.08 0.93 1.00] .* [.18;.45] + [.01;0];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd));
plot_spline(sCo12,'head',1,'headpos',1,'ratio',sRat(pInd));

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% e: Nodes and edges
pInd = 5;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .05;               % Scale nodes
sh = [-.08;.06];        % Shift nodes
scE = .08;              % Scale edges
shE = [1.35;.08];       % Shift edges

% Nodes
Nrows = 13;
Ncols = 26;
[Xcoor,Ycoor] = meshgrid(1:Ncols,1:Nrows);
XN = [Xcoor(:) Ycoor(:)]';
visualize_network(XN*sc+sh,[],[1 1],'msize',2);

% Edges
Erows = 32;
Ecols = 21;
[XcoorE,YcoorE] = meshgrid(1:Ecols,(1:Erows)*.25);
XE = [XcoorE(:) YcoorE(:)]';
line([-1;1]*.03+XE(1,:)*scE+shE(1), [1;1].*XE(2,:)*scE+shE(2),...
     'color','k','linewidth',.5,'clipping',0)

% Text
texte = '\textbf{e}\hspace{2mm}how do we combine many nodes and edges...';
text(labX,subp(pInd,4)+labY,texte,NVTitle{:});
text(.21,.83,'$N$ nodes',NVTextH{:});
text(.75,.83,'$2N-4$ edges',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% f: Clover conformational motion
pInd = 6;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Load data
load clover_data;

% Parameters
sc = .03;               % Scale drawing
sh = [2.98;0.08];       % Shift drawing
Ns = 114;               % Number of non-adde                                                                                                                                                                                                                                                                                                                                                                                                                                                        d nodes in clover example

% Correct rotation
XCcc = XCcc - XCcc(:,Ns,:);
XCcc(1,:,:) = -XCcc(1,:,:);
for i = 1:size(XCcc,3)
    Rz = rotz(atan2d(diff(XCcc(2,[-2 0]+Ns,i)),diff(XCcc(1,[-2 0]+Ns,i))));
    XCcc(:,:,i) = Rz(1:2,1:2)'*XCcc(:,:,i);
end

% Draw networks
pI = [1 260 1500];
pIL = floor(ones(1,length(pI))*size(conncc,1));
for i = 1:length(pI)
    la = .1 + .1*i;
    visualize_network(XCcc(:,unique(conncc(1:pIL(i),1)),pI(i))*sc+sh,...
                      XCcc(:,unique(conncc(1:pIL(i),2)),pI(i))*sc+sh,...
                      conncc(1:pIL(i),:) + [0 max(conncc(1:pIL(i)))-size(Xscc,2)],...
                      'lalpha',la,'msize',2,'ucolor',CSSc(1:length(unique(conncc(1:pIL(i),2))),:).^(.05*i));
end
% final conformation
sh = [2.93; 0.11];
Rz = rotz(45); Rz = Rz(1:2,1:2);
visualize_network(Rz*XCcc(:,1:size(Xscc,2),end)*sc+sh,...
                  Rz*XCcc(:,(1:size(Xucc,2))+size(Xscc,2),end)*sc+sh,...
                  conncc,'msize',2,'ucolor',CSSc);

% Text
textf = '\textbf{f}\hspace{20mm}...to design precise conformational changes?';
text(labX,subp(pInd,4)+labY-.1,textf,NVTitle{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;