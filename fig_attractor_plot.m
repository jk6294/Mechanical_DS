% Attractors
%% Prepare Space
clear; clc;
fig = figure(4); clf;
params_fig;
fig_attractor_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [17.8 7.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .0 .2 .0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  7.50  7.50];...
        [ 7.30  5.00  2.50  2.50];...
        [ 7.30  2.50  2.50  2.50];...
        [ 7.30  0.00  2.50  2.50];...
        [ 9.70  5.00  9.00  2.50];...
        [ 9.70  2.50  9.00  2.50];...
        [ 9.70  0.00  9.00  2.50]];
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


%% a: single module conformational motion & big cobweb plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
aX0 = 1.28;         % Map axes start
aXF = 2.05;         % Map axes end
% Colors
CPI = interp1(DLin,CP1,D(1,:));
CPJ = interp1(DLin,CP1,D(2,:));

% Draw map ticks
line([0 .04]+aX0, sqrt([2 2]), 'color', interp1(DLin,CP1,sqrt(2)),'linewidth', 1);
line([0 .04]+aX0, [2 2], 'color', interp1(DLin,CP1,2),'linewidth', 1);

% Draw distance colormap
scatter(D(1,:),ones(1,size(D,2))*aX0+.012,40,interp1(DLin,CP1,D(1,:)),'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.012,D(2,:),40,interp1(DLin,CP1,D(2,:)),'filled','marker','s');
 
% Draw map axes
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
    
% Draw map and y=x lines
plot([aX0 2.05],[aX0 2.05], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);

% Cobweb
mInter = 1.95;
[~,mInd] = min(abs(Da(1,:) - mInter));
CPa = interp1(DLin,CP1,Da(:,1));
CPb = interp1(DLin,CP1,Da(:,mInd));
CPc = interp1(DLin,CP1,Da(:,end));
cobweb(Da(:,1),aX0,'mdecay',0,'msize',2,'lcolor',CPa,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',1);
cobweb(Da(:,mInd),aX0,'mdecay',0,'msize',2,'lcolor',CPb,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',1);
cobweb(Da(:,end),aX0,'mdecay',0,'msize',2,'lcolor',CPc,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',1);


% Text
% Title
texta = '\textbf{a}\hspace{0.4cm}each network drawn on a cobweb plot';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
% Cobweb text
text(Da(2,mInd)+.02,Da(3,mInd)-.02,'$(l_2,l_3)$');
text(Da(3,mInd)+.02,Da(4,mInd)-.02,'$(l_3,l_4)$');
text(Da(4,mInd)+.02,Da(5,mInd)-.02,'$(l_4,l_5)$');
text(Da(end,1)-.23,Da(end,1),'$l^\bullet = f(l^\bullet)$','color',CPa(1,:));
text(Da(end,end)+.015,Da(end,end)-.07,'$l^\circ = f(l^\circ)$','color',CPc(1,:));
% Axis legend
text(mInter-.03,aX0-.04,'$1.95$',NVTexth{:},'color',CPa(1,:));
text(2+.01,aX0-.04,'$2$',NVTexth{:},'color',CPb(1,:));
text(sqrt(2),aX0-.04,'$\sqrt{2}$',NVTexth{:},'color',CPc(1,:));
text(aX0-.037,sqrt(2),'$\sqrt{2}$',NVTexth{:},'color',interp1(DLin,CP1,sqrt(2)));
text(aX0-.037,2,'$2$',NVTexth{:},'color',interp1(DLin,CP1,2));
% axis labels
text(.87,.04,'$l_k$',NVTextL{:});
text(0,.89,'$l_{k+1}$',NVTextL{:});

% Axes
axis([1.24 2.2 1.24 2.2]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% Cobweb plot schematics
% Schematic 1
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;

% Plot parameters
aX0 = 1.28;         % Map axes start
aXF = 2.08;         % Map axes end

% Draw distance colormap
scatter(D(1,:),ones(1,size(D,2))*aX0+.012,10,interp1(DLin,CP1,D(1,:)),'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.012,D(2,:),10,interp1(DLin,CP1,D(2,:)),'filled','marker','s');
% Draw map axes
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
% Draw map and y=x lines
plot([aX0 2.05],[aX0 2.05], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);
% Cobweb
cobweb(Da(:,1),aX0,'mdecay',0,'msize',2,'lcolor',CPa,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',.5,...
       'headwidth',4,'headlength',3);
% Title
textb = '\textbf{b}\hspace{2.22cm}fixed point: every unit repeats because $l^\bullet = f(l^\bullet)$';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
% Axes
axis([1.24 2.2 1.24 2.2] - [0 0 1 1]*.05);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


% Schematic 2
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
% Draw distance colormap
scatter(D(1,:),ones(1,size(D,2))*aX0+.012,10,interp1(DLin,CP1,D(1,:)),'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.012,D(2,:),10,interp1(DLin,CP1,D(2,:)),'filled','marker','s');
% Draw map axes
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
% Draw map and y=x lines
plot([aX0 2.05],[aX0 2.05], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);
% Cobweb
cobweb(Da(:,mInd),aX0,'mdecay',0,'msize',2,'lcolor',CPb,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',.5,...
       'headwidth',4,'headlength',3);
% Text
textc = '\textbf{c}\hspace{2.22cm}units move from unstable fixed points to stable ones';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});
% Axes
axis([1.24 2.2 1.24 2.2] - [0 0 1 1]*.05);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


% Schematic 3
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
% Draw distance colormap
scatter(D(1,:),ones(1,size(D,2))*aX0+.012,10,interp1(DLin,CP1,D(1,:)),'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.012,D(2,:),10,interp1(DLin,CP1,D(2,:)),'filled','marker','s');
% Draw map axes
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
% Draw map and y=x lines
plot([aX0 2.05],[aX0 2.05], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);
% Cobweb
cobweb(Da(:,end),aX0,'mdecay',0,'msize',2,'lcolor',CPc,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',.5,...
       'headwidth',4,'headlength',3);
% Title
textd = '\textbf{d}\hspace{2.22cm}fixed point: every unit repeats because $l^\circ = f(l^\circ)$';
text(labX,subp(pInd,4)+labY,textd,NVTitle{:});
% Axes
axis([1.24 2.2 1.24 2.2] - [0 0 1 1]*.05);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% Actual networks
% Parameters
sc = .255;
sh1 = [.0;.5];
lSh = .025;

% Scale and shift networks for drawing
XCapa = XCa(:,:,1) * sc + sh1;
XCapb = XCa(:,:,mInd) * sc + sh1;
XCapc = XCa(:,:,end) * sc + sh1;


% Network 1
pInd = 5;
subplot('position',subpN(pInd,:)); cla; hold on;
for i = 1:size(Da,1)
    line_coordinates(XCapa(:,(-1:0)+2*i,1), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPa(i,:));
end
visualize_network(XCapa,[],conna);
% Text
text(XCapa(1,1),XCapa(2,1)-.25,'$l_1$',NVTextl{:},'color',CPa(1,:));
text(XCapa(1,end),XCapa(2,end)+.25,'$l_{15}=f^{14}(l_1)$',NVTextl{:},'color',CPa(end,:));
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


% Network 2
pInd = 6;
subplot('position',subpN(pInd,:)); cla; hold on;

for i = 1:size(Da,1)
    line_coordinates(XCapb(:,(-1:0)+2*i,1), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPb(i,:));
end
visualize_network(XCapb,[],conna);
% Text
text(XCapb(1,1),XCapb(2,1)-.25,'$l_1$',NVTextl{:},'color',CPb(1,:));
text(XCapb(1,end),XCapb(2,end)+.25,'$l_{15}$',NVTextl{:},'color',CPb(end,:));
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


% Network 3
pInd = 7;
subplot('position',subpN(pInd,:)); cla; hold on;

for i = 1:size(Da,1)
    line_coordinates(XCapc(:,(-1:0)+2*i,1), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPc(i,:));
end
visualize_network(XCapc,[],conna);
% Text
text(XCapc(1,1),XCapc(2,1)-.25,'$l_1$',NVTextl{:},'color',CPc(1,:));
text(XCapc(1,end)-.04,XCapc(2,end)+.25,'$l_{15}$',NVTextl{:},'color',CPc(end,:));
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'fig_attractor2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');