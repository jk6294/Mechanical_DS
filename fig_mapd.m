% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(2); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 7.5];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .0 .4 .0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 2.75  4.75 4.75];...
        [ 0.00 0.00  4.75 2.75];....
        [ 5.00 0.00 14.00 7.50]];
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
NVTitle = {'Units','centimeters','fontsize',FS,'Fontname','Helvetica'};
NVTextH = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center','Fontname','Helvetica'};
NVTextR = {'Units','Normalized','fontsize',FS,'Fontname','Helvetica'};


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

% Shift XC by placing the fourth node at the origin
for i = 1:size(XC,3)
    XC(:,:,i) = XC(:,:,i) - XC(:,4,i);
end
% Shift XCa by placing the first node at the origin
for i = 1:size(XCa,3)
    XCa(:,:,i) = XCa(:,:,i) - XCa(:,1,i);
end

disp(['mean simulation error: '...
     [num2str(mean(fC)) '  ' num2str(mean(fCc))] '  ' num2str(mean(fCc))]);
D = squeeze(sqrt(sum(diff(XC,1,2).^2)));
D = D([1,3],:);
Dc = squeeze(sqrt(sum(diff(XCc,1,2).^2)));
Dc = Dc(1:2:end,:);
Da = squeeze(sqrt(sum(diff(XCa,1,2).^2)));
Da = Da(1:2:end,:);


%% a: single module conformational motion
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
texta = '\textbf{a}\hspace{.15cm}each unit $k$ maps $d_k$ to $d_{k+1}$';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});

% Map
CP = winter(size(XC,3));
hold on;
plot([1.3 2.1],[1.3 2.1], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(D(1,:),D(2,:),1,CP,'o','linewidth',.1);

% Axes
aX0 = 1.27; aXF = 2.1;
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line(sqrt([2 2]), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([2 2], [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, sqrt([2 2]), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [2 2], 'color', 'k', 'linewidth', .5);
text(-.095,.15,'$\sqrt{2}$',NVTextR{:});
text(-.05,.7,'$2$',NVTextR{:});
text(.83,0.02,'$d_k$',NVTextR{:});
text(-.05,.86,'$d_{k+1}$',NVTextR{:});
text(.31,.38,'$d_{k+1} = d_k$',NVTextR{:},'rotation',45,'color',[1 1 1]*gr);
text(.47,.8,'$d_{k+1} = f(d_k)$',NVTextH{:});
arrow([1.8 1.97], [2.03 1.88], sRat(pInd));
axis([1.24 2.3 1.24 2.3]);

% Examples
plInd = [1 230 480 756];
plSh = [.22 .23 .26 .3;...
        0 -.02 -.04 -.05];
for i = 1:length(plInd)
    pI = plInd(i);
    plot(D(1,pI),D(2,pI),'s','color',CP(pI,:),'linewidth',3,'markersize',3);
    visualize_network(XC(:,:,pI)*scc + [D(1,pI);D(2,pI)]+plSh(:,i),[],conn,'lcolor',CP(pI,:)*0);
end
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);


%% b: single module with defined distances
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textb = '\textbf{b}\hspace{0.9cm}stability of one unit';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Draw module
xp1 = [.63; .38]; pSh1 = 300;
xp2 = [1.75; .38]; pSh2 = size(XC,3)-300;
% Module 1
visualize_network(XC(:,:,pSh1)*sc + xp1,[],conn,'lcolor',CP(pSh1,:),'nalpha',.5,'lalpha',.5);
visualize_network(XC(:,:,1)*sc + xp1,[],conn,'lcolor',CP(1,:));
line_coordinates(XC(:,1:2,pSh1)*sc + xp1, -lSh-.014, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(XC(:,1:2,1)*sc + xp1, -lSh-.02, nW, lw, 'color', [1 1 1]*0);
line_coordinates(XC(:,3:4,pSh1)*sc + xp1, lSh+.18, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(XC(:,3:4,1)*sc + xp1, lSh+.18, nW, lw, 'color', [1 1 1]*0);
% Module 2
visualize_network(XC(:,:,pSh2)*sc + xp2,[],conn,'lcolor',CP(pSh2,:),'nalpha',.5,'lalpha',.5);
visualize_network(XC(:,:,end)*sc + xp2,[],conn,'lcolor',CP(end,:));
line_coordinates(XC(:,1:2,pSh2)*sc + xp2, -lSh-.013, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(XC(:,1:2,end)*sc + xp2, -lSh, nW, lw, 'color', [1 1 1]*0);
line_coordinates(XC(:,3:4,pSh2)*sc + xp2, lSh+.1, nW, lw, 'color', [1 1 1]*.5);
line_coordinates(XC(:,3:4,end)*sc + xp2, lSh+.1, nW, lw, 'color', [1 1 1]*0);

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;

% Text
text(-.09,0.38,'$\Delta d_1$',NVTextR{:});
text(.17,0.75,'$\Delta d_2$',NVTextR{:});
text(.515,0.38,'$\Delta d_1$',NVTextR{:});
text(.79,0.68,'$\Delta d_2$',NVTextR{:});
text(.22,0.06,'unstable',NVTextH{:});
text(.20,-0.09,'$|\Delta d_1| < |\Delta d_2|$',NVTextH{:});
text(.83,0.06,'stable',NVTextH{:});
text(.83,-0.09,'$|\Delta d_1| > |\Delta d_2|$',NVTextH{:});


%% c-e: module combination
% Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textc = '\textbf{c}\hspace{.2cm}fixed point: if $d_1 = 2$, then $d_{k+1} = f(2) = 2$';
text(0,.975,textc,NVTextR{:});
textc = '\textbf{d}\hspace{.2cm}unstable: if $d_1 < 2$, then $d_{k+1}$ moves away from 2';
text(0,.6,textc,NVTextR{:});
textc = '\textbf{e}\hspace{.2cm}fixed point: if $d_1 = \sqrt{2}$, then $d_{k+1} = f(\sqrt{2}) = \sqrt{2}$';
text(0,.23,textc,NVTextR{:});

% Draw module
sc2 = 1/3;
pSh = [1 3000 size(XCa,3)];
pSh2 = [.25;0];
xp = [[.1;.845] [.1;.47] [.1;.10]];
sPr = {'2','2','2','2';'1.999','1.97','1.84','1.64';'$\sqrt{2}~$','$\sqrt{2}~$','$\sqrt{2}~$','$\sqrt{2}~$'};
for i = 1:3
    [~,cInds] = min(abs(Da(:,pSh(i)) - D(1,:)),[],2);
    cInds = cInds(1:end-1);

    for j = 1:3
        XCp = XCa(:,(1:6)+2*(j-1),pSh(i))*sc*sc2+xp(:,i) + pSh2*(j-1);
        line([0;pSh2(1)]+XCp(1,3:4),[0;0]+XCp(2,3:4),'linestyle',':','color',[1 1 1]*gr,'linewidth',1);
        if(j == 1)
            visualize_network(XCp(:,3:6),[],conn,'lcolor',CP(cInds(j+1),:),'nalpha',.1,'lalpha',.1);
            text(XCp(1,1)/sRat(pInd),XCp(2,1)-.09,'$d_1$',NVTextH{:});
            text(XCp(1,1)/sRat(pInd),XCp(2,1)-.005+.05*(-1)^(j-1),sPr{i,j},NVTextH{:});
            text((XCp(1,3)+pSh2(1))/sRat(pInd)-.05,XCp(2,3)+.005-.05,sPr{i,j+1},NVTextH{:});
            line_coordinates(XCp(:,1:2),(-1)^j * lSh*sc2, nW*sc2, lw*2, 'color', [1 1 1]*0);
            line_coordinates(XCp(:,3:4)+[.15;0], 0, nW*sc2, lw, 'color', [1 1 1]*0);
        else
            text((XCp(1,3)+pSh2(1))/sRat(pInd)-.05,XCp(2,3)+.05*(-1)^(j),sPr{i,j+1},NVTextH{:});
            line_coordinates(XCp(:,3:4)+[.15;0], 0, nW*sc2, lw, 'color', [1 1 1]*0);
        end
        visualize_network(XCp(:,1:4),[],conn,'lcolor',CP(cInds(j),:));
    end
end


% Lattice
text(.825,.975,'fixed point: periodic lattice',NVTextH{:});
text(.825,.6,'unstable: $|\Delta d_1| \ll |\Delta d_k|$',NVTextH{:});
texth = '\textbf{h}\hspace{.2cm}fixed point: if $d_1 = \sqrt{2}$, then $d_{k+1} = f(\sqrt{2}) = \sqrt{2}$';
text(.825,.23,'fixed point: periodic lattice',NVTextH{:});

for i = 1:3
    [~,cInds] = min(abs(Da(:,pSh(i)) - D(1,:)),[],2);
    cInds = cInds(1:end-1);
    XCp = XCa(:,:,pSh(i))*sc*sc2+xp(:,i) + [1.2;-.00];
    visualize_network(XCp,[],conna,'lcolor',CP(cInds,:));
end

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'figure2d2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');