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
NVTextH = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','right'};
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


%% Color
DLin = linspace(D(1,1)+.01,D(1,end)-.01,size(D,2));
RLin = linspace(0,1,size(D,2));
GLin = .1 * ones(1,size(D,2));
cR = .299;
cG = .587;
cB = .114;
br = 1;
BLin = sqrt(br^2 - cG*GLin.^2 - cR*RLin.^2);


CLim1 = [055  111  183]/255;
CLim2 = [183  055  096]/255;
CP = interp1([DLin(1) DLin(end)], [CLim1; CLim2], DLin);
CP = [RLin; GLin; BLin]';
% CP = flipud(winter(size(XC,3)));
CPI = interp1(DLin,CP,D(1,:));
CPJ = interp1(DLin,CP,D(2,:));


%% a: single module conformational motion
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
texta = '\textbf{a}\hspace{.2cm}each unit $k$ maps $d_k$ to $d_{k+1}$';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});

% Map
hold on;
plot([1.3 2.1],[1.3 2.1], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);

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
lSh = -.01;
for i = 1:length(plInd)
    pI = plInd(i);
    plot(D(1,pI),D(2,pI),'s','color',0*CP(pI,:),'linewidth',3,'markersize',3);
    XP = XC(:,:,pI)*scc + [D(1,pI);D(2,pI)]+plSh(:,i);
    line_coordinates(XP(:,1:2), 'lSh',lSh, 'lw',1, 'color',CPI(pI,:));
    line_coordinates(XP(:,3:4), 'lSh',lSh, 'lw',1, 'color',CPJ(pI,:));
    visualize_network(XP,[],conn);
end
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);


%% b: single module with defined distances
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textb = '\textbf{b}\hspace{.2cm}stability of one unit';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Draw module
xp1 = [.63; .38]; pSh1 = 300;
xp2 = [1.75; .38]; pSh2 = size(XC,3)-300;
lSh = .1; nW = .02; lw = 1;
% Module 1
visualize_network(XC(:,:,pSh1)*sc + xp1,[],conn,'nalpha',.5,'lalpha',.5);
visualize_network(XC(:,:,1)*sc + xp1,[],conn);
line_coordinates(XC(:,1:2,pSh1)*sc + xp1,...
    'lSh',-lSh, 'nW',nW, 'lw',lw, 'color', CPI(pSh1,:), 'style',':');
line_coordinates(XC(:,1:2,1)*sc + xp1,...
    'lSh',-lSh-.08, 'nW',nW, 'lw',lw, 'color', CPI(1,:), 'style',':');
line_coordinates(XC(:,3:4,pSh1)*sc + xp1,...
    'lSh',lSh+.12, 'nW',nW, 'lw',lw, 'color', CPJ(pSh1,:), 'style',':');
line_coordinates(XC(:,3:4,1)*sc + xp1,...
    'lSh',lSh+.18, 'nW',nW, 'lw',lw, 'color', CPJ(1,:), 'style',':');
% Module 2
visualize_network(XC(:,:,pSh2)*sc + xp2,[],conn,'nalpha',.5,'lalpha',.5);
visualize_network(XC(:,:,end)*sc + xp2,[],conn);
line_coordinates(XC(:,1:2,pSh2)*sc + xp2,...
    'lSh',-lSh-.0, 'nw',nW, 'lw',lw, 'color', CPI(pSh2,:), 'style',':');
line_coordinates(XC(:,1:2,end)*sc + xp2,...
    'lSh',-lSh-.05, 'nw',nW, 'lw',lw, 'color', CPI(end,:), 'style',':');
line_coordinates(XC(:,3:4,pSh2)*sc + xp2,...
    'lSh',lSh+.1, 'nw',nW, 'lw',lw, 'color', CPJ(pSh2,:), 'style',':');
line_coordinates(XC(:,3:4,end)*sc + xp2,...
    'lSh',lSh+.15, 'nw',nW, 'lw',lw, 'color', CPJ(end,:), 'style',':');

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;

% Text
text(-.04,0.67,'$\Delta d_1$',NVTextR{:});
text(.17,0.75,'$\Delta d_2$',NVTextR{:});
text(.55,0.65,'$\Delta d_1$',NVTextR{:});
text(.79,0.75,'$\Delta d_2$',NVTextR{:});
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
pSh = [1 3600 size(XCa,3)];
pSh2 = [.25;0];
scX = [linspace(0,pSh2(1),15); zeros(1,15)];
tySh = [.00 .00 .0];
tyWd = [-.03 -.03 -.03];
xp = [[.1;.845] [.1;.47] [.1;.08]];
lSh = .008;
sPr = {'$d_1$', '$d_2$', '$d_3$'};
for i = 1:3
    [~,cInds] = min(abs(Da(:,pSh(i)) - D(1,:)),[],2);
    cInds = cInds(1:end-1);

    for j = 1:3
        XCp = XCa(:,(1:6)+2*(j-1),pSh(i))*sc*sc2+xp(:,i) + pSh2*(j-1);
        if(j<3)
            scatter(scX(1,:)+XCp(1,3),scX(2,:)+XCp(2,3),.2,[1 1 1]*.9);
            scatter(scX(1,:)+XCp(1,4),scX(2,:)+XCp(2,4),.2,[1 1 1]*.9);
            if(i == 1)
                    visualize_network(XCp(:,3:6),[],conn,'nalpha',.1,'lalpha',.1);
                    visualize_network(XCp(:,3:6)+[.04;0],[],conn,'nalpha',.2,'lalpha',.2);
            end
        else
            scatter(scX(1,:)+XCp(1,3),scX(2,:)+XCp(2,3),.2,[1 1 1].*linspace(.9,1,size(scX,2))');
            scatter(scX(1,:)+XCp(1,4),scX(2,:)+XCp(2,4),.2,[1 1 1].*linspace(.9,1,size(scX,2))');
        end
        CPProx = interp1(DLin,CPI,Da(j,pSh(i)));
        line_coordinates(XCp(:,1:2), 'lSh',(-1)^j*lSh, 'lw',1, 'color', CPProx);
        text((XCp(1,2))/sRat(pInd)-.05,XCp(2,2)+tyWd(i)*(-1)^(j)+tySh(i),sPr{j},NVTextH{:},'color', CPProx);
        visualize_network(XCp(:,1:4),[],conn);
    end
end
text(.11,.71,'glue together',NVTextR{:},'color',[1 1 1]*.6);

% Lattice
text(1,.975,'periodic',NVTextRA{:});
text(1,.925,'lattice',NVTextRA{:});
% text(1,.6,'$|\Delta d_1|$',NVTextRA{:});
% text(.985,.55,'$\ll$',NVTextRA{:});
% text(1,.5,'$|\Delta d_9|$',NVTextRA{:});
text(1,.23,'periodic',NVTextRA{:});
text(1,.18,'lattice',NVTextRA{:});

xSh2 = [1.07 1.07 1.07; .09 .09 .09];
shL = [-.04 -.02 0];
for i = 1:3
    [~,cInds] = min(abs(Da(:,pSh(i)) - D(1,:)),[],2);
    cInds = cInds(1:end-1);
    XCp = XCa(:,:,pSh(i))*sc*sc2+xp(:,i) + xSh2(:,i);
    for j = 1:size(Da,1)
        line_coordinates(XCp(:,[-1:0]+2*j), 'lSh',(-1)^j*lSh, 'lw',1, 'color', interp1(DLin,CPI,Da(j,pSh(i))));
        if(j < 3 || j == 9)
            text(XCp(1,-1+2*j)+.03,xp(2,i)-.11+shL(i),num2str(j),'fontsize',FS,'horizontalalignment','center');
        elseif(j==3)
            text(XCp(1,-1+2*j)+.03,xp(2,i)-.11+shL(i),'...','fontsize',FS,'horizontalalignment','center');
        end
    end
    scatter(XCp(1,1:2:end)+.03,(xp(2,i)-.19+shL(i))+Da(:,pSh(i))'/11,3,interp1(DLin,CPI,Da(:,pSh(i))),...
            'marker','s','linewidth',2);
    visualize_network(XCp,[],conna);
    line([0 .6]+1.16, [0 0]+xp(2,i)-.08+shL(i),'color',[1 1 1.1]*.7,'linewidth',.5);
    line([0 0]+1.16, [0 .14]+xp(2,i)-.08+shL(i),'color',[1 1 1.1]*.7,'linewidth',.5);
    line([0 .01]+1.16, [0 0]+xp(2,i)-.19+shL(i)+2/11,'color','k');
    line([0 .01]+1.16, [0 0]+xp(2,i)-.19+shL(i)+sqrt(2)/11,'color','k');
    line([1;1].*(XCp(1,1:2:end)+.03), [0;.01].*ones(1,size(XCp(1,1:2:end),2))+xp(2,i)-.08+shL(i),'color','k');
    text(1.15,xp(2,i)-.19+shL(i)+2/11,'$2$','fontsize',FS,'horizontalalignment','right');
    text(1.15,xp(2,i)-.19+shL(i)+sqrt(2)/11,'$\sqrt{2}$','fontsize',FS,'horizontalalignment','right');
    text(1.15,xp(2,i)-.19+shL(i)+2.7/11,'$d_k$','fontsize',FS,'horizontalalignment','right','color',[1 1 1]*.7);
    text(1.78,xp(2,i)-.08+shL(i),'$k$','fontsize',FS,'color',[1 1 1]*.7);
%     visualize_network(XCp,[],conna,'lcolor',CP(cInds,:));
end

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'figure2e';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');