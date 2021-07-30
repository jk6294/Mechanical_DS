% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 9.00];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .0 .4 .0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 5.25 3.25 3.75];...
        [ 3.75 5.25 3.25 3.75];....
        [ 8.00 5.25 5.00 3.75];...
        [13.75 5.25 5.00 3.75];...
        [ 0.00 0.10 19.0 4.75]];
% Fontsize
FS = 10;
gr = 0.8;
sc = .15;
    
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.17;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');

% Default Name-Value Pairs
NVTitle = {'Units','centimeters','fontsize',FS};
NVTextH = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center'};
NVTextR = {'Units','Normalized','fontsize',FS};



%% a: building blocks
% Plot figure
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
texta = '\textbf{a}\hspace{.5cm}building blocks';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});

% Nodes
text(0,.83,'node:',NVTextR{:});
text(.32,.83,'$+2$ motions',NVTextR{:});
visualize_network([.3;.7],[],[1 1]);
visualize_network([.5;.7],[],[1 1],'nalpha',.25);
visualize_network([.3;.5],[],[1 1],'nalpha',.25);
arrow([.35 .46], [.7 .7], sRat(pInd), 'color', [1 1 1]*gr);
arrow([.3 .3], [.66 .54], sRat(pInd), 'color', [1 1 1]*gr);

% Edge
text(0,.3,'edge:',NVTextR{:});
text(.32,.3,'$-1$ motion',NVTextR{:});
arrow([.5 .65], [.15 .15], sRat(pInd), 'color', [1 1 1]*gr);
arrow([.3 .15], [.15 .15], sRat(pInd), 'color', [1 1 1]*gr);
visualize_network([.3 .5; .15 .15],[],[1 2]);
plot(.58,.15,'rx','linewidth',1,'markersize',5);
plot(.22,.15,'rx','linewidth',1,'markersize',5);

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: combining blocks
% Plot figure
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title and text
textb = '\textbf{b}\hspace{.5cm}combine blocks';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
text(-.04,.83,'4 nodes~~~~~~4 edges',NVTextR{:});
text(-.13,.3,'$+8$ motions~$-4$ motions',NVTextR{:});
text(.15,.15,'$=4$ motions',NVTextR{:});

% Build network
L = 1;
sh = sqrt(2)*L/2;
Xs = [-sh -sh 0 sh*2;...
       sh -sh 0 0];
conn = [1 3; 1 4; 2 3; 2 4];
visualize_network(Xs*sc + [.4;.6],[],conn);


% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: rigid body motions
% Plot figure
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textb = '\textbf{c}\hspace{1cm}3 rigid body motions';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Build network
L = 1;
sh = sqrt(2)*L/2;
Xs = [-sh -sh 0 sh*2;...
       sh -sh 0 0];
Rz = rotz(20); Rz = Rz(1:2,1:2);
conn = [1 3; 1 4; 2 3; 2 4];

% x shift
text(.05,.83,'x shift',NVTextR{:});
arrow([.1 .3], [.6  .6], sRat(pInd), 'color', [1 1 1]*gr);
visualize_network(Xs*sc + [.21;.3],[],conn,'nalpha',.25,'lalpha',.25);
visualize_network(Xs*sc + [.11;.3],[],conn);

% y shift
text(.4,.83,'y shift',NVTextR{:});
arrow([.73 .73], [.7  .5], sRat(pInd), 'color', [1 1 1]*gr);
visualize_network(Xs*sc + [.67;.23],[],conn,'nalpha',.25,'lalpha',.25);
visualize_network(Xs*sc + [.67;.3],[],conn);

% rotation
th = linspace(0,3*pi/2,100);
plot(.1*cos(th)+1.2,.1*sin(th)+.6,'-','color',[1 1 1]*gr,'linewidth',.5);
arrow([1.2 1.24], [.5 .5], sRat(pInd), 'color', [1 1 1]*gr);
text(.73,.83,'rotation',NVTextR{:});
visualize_network(Rz*Xs*sc + [1.15;.3],[],conn,'nalpha',.25,'lalpha',.25);
visualize_network(Xs*sc + [1.15;.3],[],conn);

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% d: conformational motion
% Plot figure
pInd = 4;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textb = '\textbf{d}\hspace{.7cm}1 conformational motion';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Build network
L = 1;
sh = sqrt(2)*L/2;
Xs = [-sh -sh 0 sh*2;...
       sh -sh 0 0];
Rz = rotz(20); Rz = Rz(1:2,1:2);
conn = [1 3; 1 4; 2 3; 2 4];

% Simulate network
[XC,fC] = sim_motion(Xs,[],conn,.005,152,[0 0 0 0; 1 -1 0 0],0);
disp(['mean simulation error: ' num2str(mean(fC))]);

% Draw figures
C = winter(100);
scatter(linspace(.2,sRat(pInd)-.1,100), .65*ones(1,100),5,C,'s','filled');
arrow([.2 .15], [.65 .65], sRat(pInd), 'color',C(1,:),'headwidth',10,'headlength',8);
arrow([-.1 -.05]+sRat(pInd), [.65 .65], sRat(pInd), 'color',C(end,:),'headwidth',10,'headlength',8);
visualize_network(XC(:,:,1)*sc + [0.25;.3],[],conn,'lcolor',C(1,:));
visualize_network(XC(:,:,75)*sc + [0.7;.3],[],conn,'lcolor',C(50,:));
visualize_network(XC(:,:,152)*sc + [1.14;.3],[],conn,'lcolor',C(100,:));

% Lengths
line_coordinates(XC(:,1:2,1)*sc + [0.25;.3], -.06, .02, .5);
line_coordinates(XC(:,3:4,1)*sc + [0.25;.3], .12, .02, .5);
text(-.01,.3,'$d_1$',NVTextR{:});
text(.23,.5,'$d_2$',NVTextR{:});

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;



%% Size and Save Figure
fName = 'figure1c';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');