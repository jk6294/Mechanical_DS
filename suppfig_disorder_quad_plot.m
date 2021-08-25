% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
fig = figure(9); clf;
params_fig;
suppfig_disorder_quad_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [17.8 8.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  8.50  8.50];...
        [ 9.30  0.00  8.50  8.50]];
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


%% Plot errors
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
% Plot
plot(BLenDel*100, DMatMin, '.');
ticInds = [0 .2 .4 .6 .8 1];
line([1;1]*ticInds, [0;.01]*ones(1,length(ticInds)),'color','k');
line([0;.01]*ones(1,length(ticInds)), [1;1]*ticInds,'color','k');
line([0 0 1], [1 0 0], 'color', 'k');
for i = 2:length(ticInds)
    text(ticInds(i),-.03,num2str(ticInds(i)),NVTexth{:});
    text(-.01,ticInds(i),num2str(ticInds(i)),NVTextr{:});
end
text(-.01,-.03,'0',NVTextr{:});
% Title
texta = '\textbf{a}\hspace{0.4cm}node end position \textit{versus} relative bond error';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(.5,0,'average \% change in bond length',NVTextH{:});
text(0,.5,'average error of final node positions',NVTextH{:},'rotation',90);
% Axes
axis([-.13 1.1 -.1 1.1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Draw worst case network
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;

sc = .03;              % Scale drawing
sh = [.5;0.49];       % Shift drawing
[~,plInd] = max(DMatMin);
[~,mI] = min(DMat(plInd,:));
Xcce = XCcc(:,:,mI,plInd);
Xcce = Xcce - mean(Xcce(:,1:Ns),2);
[u,s,v] = svd(xPc * Xcce(:,1:Ns)');
Rp = u*v';
Xcce = Rp*Xcce;

% Get min and max values
yMin = min(Xcce(2,:));
yMax = max(Xcce(2,:));
xMin = min(Xcce(1,:));

plot((xPc(1,1:end-1)+xPc(1,2:end))/2*sc+sh(1),...
     (xPc(2,1:end-1)+xPc(2,2:end))/2*sc+sh(2),'-','linewidth',2);
% visualize_network(Xscca*sc+sh,Xucca*sc+sh,[1 1]);
visualize_network(Xcce(:,1:Ns)*sc+sh,Xcce(:,Ns+1:end)*sc+sh,conncc,'ucolor',CSSc);
% visualize_network(xPc*sc + sh,[],[1 1],'scolor',[0 0 1]);
line_coordinates([[1 1]*-15; -15 15]*sc+sh,...
                 'style','-','lw',.5,'nw',.01)
text(-15*sc+sh(1),sh(2),'30',NVTextr{:});

% Title
textb = '\textbf{b}\hspace{0.4cm}drawing of the network with largest error';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'suppfig_disorder_quad';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');