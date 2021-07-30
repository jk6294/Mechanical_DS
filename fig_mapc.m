% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(2); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 7.25];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .4 0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 4.75 4.75 2.50];...
        [ 4.75 4.75 4.75 2.50];...
        [ 9.50 4.75 4.75 2.50];...
        [14.25 4.75 4.75 2.50];....
        [ 0.00 0.00 4.75 4.75];...
        [ 4.75 0.00 4.75 4.75];...
        [ 9.50 0.00 4.75 4.75];...
        [14.25 0.00 4.74 4.75]];
% Fontsize
FS = 10;
% Distance visualization parameters
lSh = .1;
nW = .03;
lw = .5;
gr = 0.8;
% Scaling parameters
sc = .2;
scc = .075;
    
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


%% Define module
L = 1;
Rz = rotz(45); Rz = Rz(1:2,1:2);
% Single module
Xs = [0  0  0  2*L;...
      L -L  0  0];
conn = [1 3; 1 4; 2 3; 2 4];
% Double module
Xsc = [0  0  0  2*L  L  L;...
       L -L  0  0    0 -2*L];
connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];
% Tesselate module
[Xsa, conna] = tesselate_network(Rz*Xsc, connc, [sqrt(2);0], [4;1]);


% Simulate network
Xs0 = [-1 -1 0 0; 0 0 0 0];
Xsc0 = [0 0 0 0 0  0;...
        0 0 0 0 -1 0];
Xsa0 = 0*Xsa; Xsa0(1,end-1) = 1;
[XC,fC] = sim_motion(Xs,[],conn,.001,756,Xs0,0);
[XCc,fCc] = sim_motion(Xsc,[],connc,.001,1122,Xsc0,0);
[XCa,fCa] = sim_motion(Xsa,[],conna,.001,4379,Xsa0,0);
disp(['mean simulation error: '...
     [num2str(mean(fC)) '  ' num2str(mean(fCc))] '  ' num2str(mean(fCc))]);
d1 = squeeze(sqrt(sum(diff(XC(:,[1,2],:),1,2).^2)));
d2 = squeeze(sqrt(sum(diff(XC(:,[3,4],:),1,2).^2)));
Dc = squeeze(sqrt(sum(diff(XCc,1,2).^2)));
Dc = Dc(1:2:end,:);
d1c = squeeze(sqrt(sum(diff(XCc(:,[1,2],:),1,2).^2)));
d2c = squeeze(sqrt(sum(diff(XCc(:,[3,4],:),1,2).^2)));
d3c = squeeze(sqrt(sum(diff(XCc(:,[5,6],:),1,2).^2)));
Da = squeeze(sqrt(sum(diff(XCa,1,2).^2)));
Da = Da(1:2:end,:);


%% a: single module with defined distances
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
texta = '\textbf{a}\hspace{.15cm}one module maps $d_1$ to $d_2$';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});

% Draw module
xp1 = .8;
visualize_network(Rz*Xs*sc + [xp1;.35],[],conn);
line_coordinates(Rz*Xs(:,1:2)*sc + [xp1;.35], -lSh, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(Rz*Xs(:,3:4)*sc + [xp1;.35], -lSh, nW, lw, 'color', [1 1 1]*.5);

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;

% Text
text(.28,0.2,'$d_1$',NVTextR{:});
text(.52,0.4,'$d_2$',NVTextR{:});
text(.42,0.45,'$f$',NVTextR{:});


%% a: single module conformational motion
% Plot
pInd = 5;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
texta2 = '\hspace{.4cm}changing $d_1$ changes $d_2$';
text(labX,subp(pInd,4)+labY,texta2,NVTitle{:});

% Map
CP = winter(size(XC,3));
hold on;
plot([1.3 2.1],[1.3 2.1], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(d1,d2,1,CP,'o','linewidth',.1);

% Axes
aX0 = 1.27; aXF = 2.1;
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line(sqrt([2 2]), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([2 2], [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, sqrt([2 2]), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [2 2], 'color', 'k', 'linewidth', .5);
text(.07,-.05,'$\sqrt{2}$',NVTextR{:});
text(.7,-.05,'$2$',NVTextR{:});
text(.83,0.02,'$d_1$',NVTextR{:});
text(-.01,.86,'$d_2$',NVTextR{:});
text(.31,.4,'$d_2 = d_1$',NVTextR{:},'rotation',45,'color',[1 1 1]*gr);
text(.15,.7,'$d_2 = f(d_1)$',NVTextR{:});
arrow([.01 .17]+1.8, [.18 0.08]+1.8, sRat(pInd));
axis([1.25 2.3 1.25 2.3]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);

% Examples
plInd = [1 230 480 756];
plSh = [.07 .09 .12 .14;...
        0 -.02 -.04 -.05];
for i = 1:length(plInd)
    pI = plInd(i);
    plot(d1(pI),d2(pI),'s','color',CP(pI,:),'linewidth',3,'markersize',3);
    visualize_network(XC(:,:,pI)*scc + [d1(pI);d2(pI)]+plSh(:,i),[],conn,'lcolor',CP(pI,:));
end


%% b: two modules
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textb = '\textbf{b}\hspace{.15cm}two modules map $d_1$ to $d_3$';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Draw module
xp1 = .45; xp2 = 1.4;
visualize_network(Rz*Xs*sc + [xp1;.35],[],conn);
visualize_network(Rz'*Xs*sc + [xp2;.35+1/sqrt(2)*sc],[],conn);
line_coordinates(Rz*Xs(:,1:2)*sc + [xp1;.35], -.1, .04, lw, 'color', [1 1 1]*.5);
line_coordinates(Rz*Xs(:,3:4)*sc + [xp1;.35], -.1, .04, lw, 'color', [1 1 1]*.5);
line_coordinates(Rz'*Xs(:,1:2)*sc + [xp2;.35+1/sqrt(2)*sc], -lSh, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(Rz'*Xs(:,3:4)*sc + [xp2;.35+1/sqrt(2)*sc], lSh, nW, lw, 'color', [1 1 1]*.5);
% arrow([0 .2]+xp1-.45, [.25 .25], sRat(pInd));
arrow([0 .3]+xp1+.4, [.5 .5], sRat(pInd));
% arrow([0 .2]+xp2+.35, [.25 .25], sRat(pInd));

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;

% Text
text(.11,0.2,'$d_1$',NVTextR{:});
text(.45,0.6,'$d_2$',NVTextR{:});
text(.82,0.46,'$d_3$',NVTextR{:});
text(.25,0.45,'$f$',NVTextR{:});
text(.69,0.4,'$f$',NVTextR{:});



%% b: two individual modules conformational motion
% Plot
pInd = 6;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textb2 = '\hspace{.4cm}changing $d_1$ changes $d_3$';
text(labX,subp(pInd,4)+labY,textb2,NVTitle{:});

% Map
CP = winter(size(XC,3));
hold on;
plot([1.3 2.1],[1.3 2.1], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(d1,d2,1,CP,'o','linewidth',.1);

% Axes
aX0 = 1.27; aXF = 2.1;
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line(sqrt([2 2]), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([2 2], [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, sqrt([2 2]), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [2 2], 'color', 'k', 'linewidth', .5);
text(.07,-.05,'$\sqrt{2}$',NVTextR{:});
text(.7,-.05,'$2$',NVTextR{:});
text(.83,0.02,'$d_1$',NVTextR{:});
text(-.01,.86,'$d_2$',NVTextR{:});
axis([1.25 2.3 1.25 2.3]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);

% Examples
plInds = [50 500];
plSh = [1.4 1.75;...
        1.9   1.4];
p2Sh = [1.6 1.95;...
        1.95   1.45];
for i = 1:length(plInds)
    pI = plInds(i);
    [~,pI2] = min(abs(d1 - d2(pI)));
    cobweb([d1(pI); d2(pI); d2(pI2)],'msize',3,'ncolor',CP([pI pI2],:),'mdecay',0,'linestyle',':');
%     cobweb([d2(pI); d2(pI2)],'msize',3,'ncolor',CP(pI2,:),'mdecay',0);
    visualize_network(Rz*XC(:,:,pI)*scc + plSh(:,i),[],conn,'lcolor',CP(pI,:));
    visualize_network(Rz'*XC(:,:,pI2)*scc +p2Sh(:,i),[],conn,'lcolor',CP(pI2,:));
end


%% c: combine modules
% Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textc = '\textbf{c}\hspace{.3cm}join nodes to combine';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});

% Draw module
xp1 = .45; xp2 = 1.4;
Xsp1 = Rz*Xs*sc + [xp1;.35];
Xsp2 = Rz'*Xs*sc + [xp2;.35+1/sqrt(2)*sc];
line_coordinates(Xsp1(:,1:2), -lSh, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(Xsp2(:,3:4), lSh, nW, lw, 'color', [1 1 1]*.5);
line([Xsp1(1,3) Xsp2(1,2)],[Xsp1(2,3) Xsp2(2,2)],'linestyle',':','color',[1 1 1]*gr,'linewidth',1);
line([Xsp1(1,4) Xsp2(1,1)],[Xsp1(2,4) Xsp2(2,1)],'linestyle',':','color',[1 1 1]*gr,'linewidth',1);
visualize_network(Xsp2 + [xp1-xp2+sc/sqrt(2);0],[],conn,'nalpha',.2,'lalpha',.2);
visualize_network(Xsp1,[],conn);
visualize_network(Xsp2,[],conn);

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;

% Text
text(.11,.2,'$d_1$',NVTextR{:});
text(.81,.44,'$d_3$',NVTextR{:});



%% c: two individual modules conformational motion
% Plot
pInd = 7;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textc2 = '\hspace{.4cm}changing $d_1$ changes $d_3$';
text(labX,subp(pInd,4)+labY,textc2,NVTitle{:});

% Map
hold on;
plot([1.3 2.1],[1.3 2.1], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(d1,d2,1,CP,'o','linewidth',.1);

% Axes
aX0 = 1.27; aXF = 2.1;
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line(sqrt([2 2]), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([2 2], [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, sqrt([2 2]), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [2 2], 'color', 'k', 'linewidth', .5);
text(.07,-.05,'$\sqrt{2}$',NVTextR{:});
text(.7,-.05,'$2$',NVTextR{:});
text(.83,0.02,'$d_1$',NVTextR{:});
text(-.01,.86,'$d_2$',NVTextR{:});
axis([1.25 2.3 1.25 2.3]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);

% Examples
plSh = [1.5 1.75;...
        1.9   1.4];
for i = 1:length(plInds)
    [~,plIndc] = min(abs(Dc(1,:)-d1(plInds(i))));
    [~,cInds] = min(abs(Dc(1:end-1,plIndc)' - d1));
    cobweb(Dc(:,plIndc),'msize',3,'ncolor',CP(cInds,:),'mdecay',0,'arrowind',[1 2]);
    visualize_network(Rz*XCc(:,:,plIndc)*scc + plSh(:,i),[],connc,'lcolor',CP(cInds,:));
end


%% d: combine many modules
% Plot
pInd = 4;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textd = '\textbf{d}\hspace{.3cm}combining $k$ modules';
text(labX,subp(pInd,4)+labY,textd,NVTitle{:});

% Draw module
xp1 = .3; xp2 = 1.4;
Xsp1 = Rz*Xs*sc + [xp1;.35];
Xsp2 = Rz'*Xs*sc + [xp2;.35+1/sqrt(2)*sc];
line_coordinates(Xsa(:,1:2)*sc + [xp1;.35], -lSh, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(Xsa(:,end-1:end)*sc + [xp1;.35], lSh, nW, lw, 'color', [1 1 1]*.5);
visualize_network(Xsa*sc + [xp1;.35], [], conna);

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;

% Text
text(.04,.2,'$d_1$',NVTextR{:});
text(.77,.44,'$d_{k+1}$',NVTextR{:});



%% d: two individual modules conformational motion
% Plot
pInd = 8;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textd2 = '\hspace{.3cm}changing $d_1$ changes $d_{k+1}$';
text(labX,subp(pInd,4)+labY,textd2,NVTitle{:});

% Map
hold on;
plot([1.3 2.1],[1.3 2.1], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(d1,d2,1,CP,'o','linewidth',.1);

% Axes
aX0 = 1.27; aXF = 2.1;
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line(sqrt([2 2]), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([2 2], [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, sqrt([2 2]), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [2 2], 'color', 'k', 'linewidth', .5);
text(.07,-.05,'$\sqrt{2}$',NVTextR{:});
text(.7,-.05,'$2$',NVTextR{:});
text(.83,0.02,'$d_1$',NVTextR{:});
text(-.01,.86,'$d_2$',NVTextR{:});
axis([1.25 2.3 1.25 2.3]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);

% Examples
plSh = [1.39  1.7;...
        2.01  1.38];
NC = ones(2,3) .* [.8 .7]';
plIa = [2150];
for i = 1:length(plIa)
    [~,cInds] = min(abs(Da(:,plIa(i))' - d1));
    cInds = cInds(1:end-1);
    cobweb(Da(:,plIa(i)),'msize',3,'ncolor',CP(cInds,:),'mdecay',0,'arrowind',[3 6]);
    visualize_network(XCa(:,:,plIa(i))*scc + plSh(:,i),[],conna,'lcolor',CP(cInds,:));
end


%% Size and Save Figure
fName = 'figure2c';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');