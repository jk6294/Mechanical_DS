% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 8.00];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .4 .4 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 5.25 3.75 2.75];...
        [ 4.50 5.25 3.25 2.75];....
        [ 8.25 5.25 5.25 2.75];...
        [13.75 5.25 5.25 2.75];...
        [ 0.00 0.10 19.0 4.75]];
% Fontsize
FS = 10;
lSh = .2;
nW = .04;
lw = .5;
gr = 0.8;
    
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



%% a,b: Nodes and edges
delete(findall(gca,'type','annotation'))
% Plot figure
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
hold on;

% a: Legend
% Node
visualize_network([1;2.35],[],[1 1]);
text(.01,.75,'node',NVTextR{:});
% Edge
line([0 .5]+3,[1 1]*2.35,'color','k','linewidth',1);
text(.66,.75,'edge',NVTextR{:});

% b: Examples
% arrow([0 .3]+1.9, [0 -.3]+0.8, sRat(pInd), 'color', [1 1 1]*gr);
% Node
arrow([0 0]+1, [0 .35]+.95, sRat(pInd), 'color', [1 1 1]*gr);
arrow([0 .35]+1.15, [0 0]+.8, sRat(pInd), 'color', [1 1 1]*gr);
visualize_network([1;.8],[],[1 1]);
text(0.01,.05,'$+2$',NVTextR{:});
arrow([0 .25]+1.2,[0 -.25]+0.46,sRat(pInd), 'color', [1 1 1]*gr);
% Edge
line([0 .5]+3,[1 1]*.8,'color','k','linewidth',1);
arrow([0 .45]+3.55,[0 0]+.8,sRat(pInd), 'color', [1 1 1]*gr);
arrow([0 -.45]+2.95,[0 0]+.8,sRat(pInd), 'color', [1 1 1]*gr);
plot(3.7,.8,'rx','linewidth',1,'markersize',5);
plot(2.8,.8,'rx','linewidth',1,'markersize',5);
text(.66,.05,'$-1$',NVTextR{:});
arrow([0 .25]+3.42,[0 -.25]+0.46,sRat(pInd), 'color', [1 1 1]*gr);



% Axis limits
ax = [0 sRat(pInd) 0 1]*2.3 + [[1 1]*.55 [1 1]*.2];
axis(ax);
drawnow;
set(gca,'visible',0,'xtick',[],'ytick',[]);

% Text
texta = '\textbf{a}\hspace{.39cm}building blocks';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(.275,.45,'motion',NVTextR{:});



%% b: 
% Module parameters
L = 1;
sh = sqrt(2)*L/4;
Xs = [-sh -sh 0 sh*2;...
       sh -sh 0 0];
conn = [1 3; 1 4; 2 3; 2 4];

% Plot figure
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
hold on;
visualize_network(Xs+[0;1.05],[],conn);
% Label
visualize_network([-1;0.22],[],[1,1]);
visualize_network([-.75;0.22],[],[1,1]);
visualize_network([-.5;0.22],[],[1,1]);
visualize_network([-.25;0.22],[],[1,1]);
line([0 .3]-1,[.3 0]-.4,'color','k','linewidth',1);
line([0 .3]-.8,[.3 0]-.4,'color','k','linewidth',1);
line([0 .3]-.6,[.3 0]-.4,'color','k','linewidth',1);
line([0 .3]-.4,[.3 0]-.4,'color','k','linewidth',1);
% Arrow
arrow([0 .3]+1.1,[0 -.3]+.4,sRat(pInd), 'color', [1 1 1]*gr);
arrow([0 .3]+1.1,[0 -.3]-.04,sRat(pInd), 'color', [1 1 1]*gr);
arrow([0 .3]+1.1,[0 -.3]-.49,sRat(pInd), 'color', [1 1 1]*gr);
plot([0 1]+.46,[0 0]-.45,'k-','linewidth',.5);
% plot(1.25,-.19,'rx','linewidth',1,'markersize',5);

% Axis limits
ax = [0 sRat(pInd) 0 1]*2.3 + [[1 1]*-1 [1 1]*-.8];
axis(ax);
drawnow;

% Text
textb = '~\textbf{b}~~network unit';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
text(.35,0.45,'$\Rightarrow~~8$',NVTextR{:});
text(.35,0.25,'$\Rightarrow~~4$',NVTextR{:});
text(.5,0.25,'$-$',NVTextR{:});
text( 0.63,0.05,'4',NVTextR{:});
set(gca,'visible',0)


%% d: Rigid body motions
% Module parameters
L = 1;
sh = sqrt(2)*L/4;
Xs = [-sh -sh 0 sh*2;...
       sh -sh 0 0];
conn = [1 3; 1 4; 2 3; 2 4];
rot45 = rotz(15); rot45 = rot45(1:2,1:2);
PCf = winter(2); PCf = PCf(1,:);

% Plot figure
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
hold on;
% x-translation
visualize_network(Xs+[.25;0],[],conn,'lcolor',PCf,'lalpha',.2,'nalpha',.2);
visualize_network(Xs+[0;0],[],conn,'lcolor',PCf);
% y-translation
visualize_network(Xs+[1.75;-.2],[],conn,'lcolor',PCf,'lalpha',.2,'nalpha',.2);
visualize_network(Xs+[1.75;0],[],conn,'lcolor',PCf);
% rotation
visualize_network(rot45*Xs+[3.5;0],[],conn,'lcolor',PCf,'lalpha',.2,'nalpha',.2);
visualize_network(Xs+[3.5;0],[],conn,'lcolor',PCf);
% Arrows
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[0 .75 .5 0]);
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[1.91 1.0 0 -.5]);
plot(0.2*cosd(1:270)+3.7, 0.2*sind(1:270)+.75, 'o', 'linewidth',.5/3,'markersize',.5/3,'color',[1 1 1]*gr);
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[3.7 .55 .1 0]);

% Axis limits
ax = [0 sRat(pInd) 0 1]*2.3 + [[1 1]*-.8 [1 1]*-.7];
axis(ax);

% Text
textz = '4 nodes + 4 edges $\Rightarrow$ 4 motions';
text(2.6,subp(pInd,4)+labY,textz,NVTitle{:});
text(0, .85,'\textbf{c}~~~~~~~~~3 rigid body',NVTextR{:});
drawnow;


%% d2: Conformational motion
% Define module
Xs2 = [0    0   0   L;...
       L/2 -L/2 0   0];

% Simulate
nM = 380;
CP = winter(floor(nM));
[XC, fC] = sim_motion(Xs,[],conn,.001,nM,Xs,0);
disp(max(fC));
pI = [1 140 nM];

% Plot figure
pInd = 4;
subplot('position',subpN(pInd,:)); cla;
% Modules
for i = 1:1:length(pI)
    visualize_network(XC(:,:,pI(i))+[(i-1)*1.75-.45;0],[],conn,'lcolor',CP(pI(i),:));
end
line_coordinates(XC(:,1:2,pI(1))+[-.45;0],-lSh,nW,lw);
line_coordinates(XC(:,3:4,pI(1))+[-.45;0],-lSh,nW,lw);
% Colorbar
ax = [0 sRat(pInd) 0 1]*2.3 + [[1 1]*-1.4 [1 1]*-.6];
set(annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,...
               'linewidth', 1, 'position',[2 .9 1.7 0]),...
               'color', CP(end,:),'parent', gca);
set(annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,...
               'linewidth', 1, 'position',[2 .9 -2.9 0]),...
               'color', CP(1,:), 'parent', gca);
scatter(linspace(ax(1)+0.7,ax(2)-.3,size(CP,1)), ones(1,size(CP,1))*.9,1,CP);
% Axis limits
axis(ax);
drawnow;
set(gca,'visible',0)

% Text
text(0, .85,'\textbf{d}\hspace{0.85cm}1 conformational',NVTextR{:});
text(0,.25,'$d_1$',NVTextR{:});
text(0.25,.06,'$d_2$',NVTextH{:});



%% e Modules determine map
s = sqrt(3);

% Node Colors for different curvatures
cTr1 = [126 200 255]/255;
cTr2 = [115 140 200]/255;
cTr3 = [015 082 186]/255;

l1 = [0;-s]; l2 = [0;-1.5*s];
a1 = 19.4; a2 = 38; a3 = 60;
R1 = rotz(-a1/2); R2 = rotz(-a2/2); R3 = rotz(-a3/2);
R1 = R1(1:2,1:2); R2 = R2(1:2,1:2); R3 = R3(1:2,1:2);

% 1
Xs10 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs1T = [R1*l2 [0;0] R1\l2] + [0;1];
Xu1 = [-1.5  1.5;...
        1.5  1.5];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu1,~] = construct_network(Xs10,Xs1T,Xu1,conn1,0,1);
Xu1 = Xu1(1:2,:);

% 2
Xs20 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs2T = [R2*l2 [0;0] R2\l2] + [0;1];
Xu2 = [-0.7  0.7;...
        1.5  1.5];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,~] = construct_network(Xs20,Xs2T,Xu2,conn2,0,1);
Xu2 = Xu2(1:2,:);

% 3
Xs30 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs3T = [R3*l2 [0;0] R3\l2] + [0;1];
Xu3 = [-1.0  1.0;...
        1.8  1.8];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu3,fV] = construct_network(Xs30,Xs3T,Xu3,conn3,0,1);
Xu3 = Xu3(1:2,:);


% Chain
XuC = cat(3,Xu1,Xu2,Xu3);
Xsp = [-s/2 0; -.5 1];
XuEL = [2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 3 2, 3 2, 3 2, 3 2 ,3 2, 3 2, 3 2,...
        2 2,...
        3 2, 3 2, 3 2, 3 2, 3 2,...
        2 2];
C_UNC = [cTr1; cTr2; cTr3];
[XsE,XuE,connE,C_UNE] = network_chain_x(Xsp,XuC,XuEL,C_UNC);
rot90 = [0 1; -1 0];
XsE = rot90*XsE;
XuE = rot90*XuE;

% Simulate
[XCE,fCE] = sim_motion(XsE,XuE,connE,.5,280,-[XsE,XuE],0);
XCE = XCE - XCE(:,end,:);


%% Plot e
pInd = 5;
subplot('position',subpN(pInd,:)); cla;
rot1 = rotz(-45); rot1 = rot1(1:2,1:2);
rot2 = rotz(90); rot2 = rot2(1:2,1:2);
rot3 = rotz(-46); rot3 = rot3(1:2,1:2);

% Intermediate
xshe = [0.1;.0]; xP = 50;
rote = rotz(70); rote = rote(1:2,1:2);
visualize_network(rote*XCE(:,1:size(XsE,2),xP)*.008+xshe,...
                  rote*XCE(:,size(XsE,2)+1:end,xP)*.008+xshe,connE,...
                  'scolor',[255 100 100]/255,'ucolor',[100 100 255]/255,'scale',.7,...
                  'lalpha',.3,'nalpha',.3);
xshe = [0.2;.005]; xP = 100;
rote = rotz(50); rote = rote(1:2,1:2);
visualize_network(rote*XCE(:,1:size(XsE,2),xP)*.008+xshe,...
                  rote*XCE(:,size(XsE,2)+1:end,xP)*.008+xshe,connE,...
                  'scolor',[255 100 100]/255,'ucolor',[100 100 255]/255,'scale',.7,...
                  'lalpha',.3,'nalpha',.3);
xshe = [0.28;.005]; xP = 170;
rote = rotz(0); rote = rote(1:2,1:2);
visualize_network(rote*XCE(:,1:size(XsE,2),xP)*.008+xshe,...
                  rote*XCE(:,size(XsE,2)+1:end,xP)*.008+xshe,connE,...
                  'scolor',[255 100 100]/255,'ucolor',[100 100 255]/255,'scale',.7,...
                  'lalpha',.3,'nalpha',.3);
% First and Last
xshe = [0;0]; xP = 1; 
rote = rotz(90); rote = rote(1:2,1:2);
visualize_network(rote*XCE(:,1:size(XsE,2),xP)*.008+xshe,...
                  rote*XCE(:,size(XsE,2)+1:end,xP)*.008+xshe,connE,...
                  'scolor',[255 100 100]/255,'ucolor',[100 100 255]/255,'scale',.7);
xshe = [0.38;.015]; xP = 250;
rote = rotz(-46); rote = rote(1:2,1:2);
visualize_network(rote*XCE(:,1:size(XsE,2),xP)*.008+xshe,...
                  rote*XCE(:,size(XsE,2)+1:end,xP)*.008+xshe,connE,...
                  'scolor',[255 100 100]/255,'ucolor',[100 100 255]/255,'scale',.7);

set(gca,'visible',0);
axis([0 sRat(pInd) 0 1]*.15 + [[1 1]*-.29 [1 1]*-0.01]);

% Text
texte = '\textbf{e}';
text(labX,subp(pInd,4)+labY,texte,NVTitle{:});
text(.525,1.055,'placement of \textbf{nodes} and \textbf{edges} $\longrightarrow$ trajectory of conformational \textbf{motion}',...
             NVTextH{:});
text(.48,-.08,'trajectory of conformational \textbf{motion} $\longrightarrow$ placement of \textbf{nodes} and \textbf{edges}',...
             NVTextH{:},'color',[255 100 100]/255);
text(.504,-.03,'?',NVTextH{:},'color',[255 100 100]/255);




%% Size and Save Figure
fName = 'figure1b';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');