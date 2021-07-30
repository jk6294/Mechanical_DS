% Supplementary Figure: edge lengths
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 4.75];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  4.75  4.75];...
        [ 4.75  0.00  4.75  4.75];...
        [ 9.50  0.00  4.75  4.75];...
        [14.25  0.00  4.75  4.75]];
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');


%% Figure Parameters
FS = 10;            % Fontsize
gr = 0.9;           % Gray
al = 0.2;           % Alpha
    
% Default Name-Value Pairs
NVTitle  = {'Units','centimeters','fontsize',FS};
NVTextH  = {'Units','Normalized','fontsize',FS,...
           'HorizontalAlignment','center'};
NVTextHU  = {'fontsize',FS,'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,...
            'HorizontalAlignment','right'};
NVTextRAU = {'fontsize',FS,'HorizontalAlignment','right'};
NVTextR  = {'Units','Normalized','fontsize',FS};
NVTextRU  = {'fontsize',FS};

% Color
Css = [031 172 204]/255;


%% Global parameters
% Construct network
s = sqrt(3);
Xs = [-s/2 0 s/2;...                    % Initial position
      -1/2 1 -1/2];
rM = 2;
Xf = rM*Xs;                             % Final position
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];  % Network connectivity

% Parameters
R = [-2 2; -2 2];               % Range of conic visualization
mSs = 5;                        % Node size: start nodes
mSf = 3;                        % Node size: end nodes
ax = [-1 1 -1 1]*2.3 + [.3 .3 .3 .3];
o = [1 1 1];


%% a: Solution Space
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Draw solution space
visualize_conic_finite(Xs,Xf,R,'ucolori',Css,'overlay',.99);
% Draw start and end positions
visualize_network(Xs,[],[1 1],'msize',mSs);
visualize_network(Xf,[],[1 1],'msize',mSf,'scolor',o);

% Draw Legend
visualize_network([-2.05;2.2],[],[1 1],'msize',mSs);
visualize_network([-2.05;1.8],[],[1 1],'msize',mSf,'scolor',o);
line([-1 1]*.2 + 2.5, [1 1]*2.2, 'linewidth',1,'color',Css);
visualize_network([2.45;1.8],[],[1 1],'msize',mSs,'scolor',Css);

% Text Legend
text(-1.85,2.2,'start',NVTextRU{:});
text(-1.85,1.8,'end',NVTextRU{:});
text(2.2,1.8,'added',NVTextRAU{:});
text(2.2,2.2,'solution',NVTextRAU{:});

% Text Nodes
sc = 1.3;
sc2 = 1.12;
text(Xs(1,1)*sc, Xs(2,1)*sc,'1',NVTextHU{:},'color',o*.5);
text(Xs(1,2)*sc, Xs(2,2)*sc,'2',NVTextHU{:},'color',o*.5);
text(Xs(1,3)*sc, Xs(2,3)*sc,'3',NVTextHU{:},'color',o*.5);
text(Xf(1,1)*sc2, Xf(2,1)*sc2,'1',NVTextHU{:},'color',o*.5);
text(Xf(1,2)*sc2, Xf(2,2)*sc2,'2',NVTextHU{:},'color',o*.5);
text(Xf(1,3)*sc2, Xf(2,3)*sc2,'3',NVTextHU{:},'color',o*.5);

% Draw
axis(ax);
set(gca,'visible',0);
drawnow;

% Text Title
texta = '\textbf{a} \hspace{3mm}find solution space';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});


%% b: Added nodes 1
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Place extra nodes
th1 = 50; th2 = -50;
Rth1 = rotz(th1); Rth1 = Rth1(1:2,1:2);
Rth2 = rotz(th2); Rth2 = Rth2(1:2,1:2);
Xu1 = [Rth1*[0;1] Rth2*[0;1]]*rM;
D = sqrt(sum((Xs(:,conn(:,1)) - Xu1(:,conn(:,2)-3)).^2));
MD = Xu1(:,conn(:,2)-3) - Xs(:,conn(:,1));
MDA = atan2d(MD(2,:),MD(1,:));

% Draw solution space
visualize_conic_finite(Xs,Xf,R,'ucolori',Css,'overlay',.2);
% Draw network
visualize_network(Xs,Xu1,conn,'msize',mSs,'ucolor',Css);

% Draw
axis(ax);
set(gca,'visible',0);
drawnow;

% Text Angles
r = [.6 .5 .85 .85 .5 .6];
xsh = [-.4 0 0 0 0 .4];
ysh = [0 .22 .35 .35 .22 0];
for i = 1:size(conn,1)
    Rprox = rotz(MDA(i)); Rprox = Rprox(1:2,1:2);
    rprox = Rprox*[-1;0]*r(i)*D(i);
    text(rprox(1)+Xu1(1,conn(i,2)-3)+xsh(i),...
         rprox(2)+Xu1(2,conn(i,2)-3)+ysh(i),...
         sprintf('%0.3g',D(i)),NVTextHU{:});
end

% Text Nodes
sc = 1.3;
text(Xs(1,1)*sc, Xs(2,1)*sc,'1',NVTextHU{:},'color',o*.5);
text(Xs(1,2)*sc, Xs(2,2)*sc,'2',NVTextHU{:},'color',o*.5);
text(Xs(1,3)*sc, Xs(2,3)*sc,'3',NVTextHU{:},'color',o*.5);
text(Xu1(1,1)*1.15, Xu1(2,1)*1.15,'4',NVTextHU{:},'color',Css);
text(Xu1(1,2)*1.15, Xu1(2,2)*1.15,'5',NVTextHU{:},'color',Css);

% Text Title
textb = '\textbf{b} \hspace{4mm}node placement 1';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});


%% c: Added nodes 2
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Place extra nodes
th1 = 90 - acosd(s/4); th2 = 90 + acosd(s/4);
Rth1 = rotz(th1); Rth1 = Rth1(1:2,1:2);
Rth2 = rotz(th2); Rth2 = Rth2(1:2,1:2);
Xu1 = [Rth1*[0;1] Rth2*[0;1]]*rM;
D = sqrt(sum((Xs(:,conn(:,1)) - Xu1(:,conn(:,2)-3)).^2));
M = (Xs(:,conn(:,1)) + Xu1(:,conn(:,2)-3))/2;
MD = Xu1(:,conn(:,2)-3) - Xs(:,conn(:,1));
MDA = atan2d(MD(2,:),MD(1,:));

% Draw solution space
visualize_conic_finite(Xs,Xf,R,'ucolori',Css,'overlay',.2);
% Draw network
visualize_network(Xs,Xu1,conn,'msize',mSs,'ucolor',Css);

% Draw
axis(ax);
set(gca,'visible',0);
drawnow;

% Text Angles
r = [.6 .4 .7 .4 .5 .4];
xsh = [-.3 .3 .45 -.3 .4 .6];
ysh = [0 .2 0 0 0 0];
for i = 1:size(conn,1)
    Rprox = rotz(MDA(i)); Rprox = Rprox(1:2,1:2);
    rprox = Rprox*[-1;0]*r(i)*D(i);
    text(rprox(1)+Xu1(1,conn(i,2)-3)+xsh(i),...
         rprox(2)+Xu1(2,conn(i,2)-3)+ysh(i),...
         sprintf('%0.3g',D(i)),NVTextHU{:});
end

% Text Nodes
sc = 1.3;
text(Xs(1,1)*sc, Xs(2,1)*sc,'1',NVTextHU{:},'color',o*.5);
text(Xs(1,2)*sc, Xs(2,2)*sc,'2',NVTextHU{:},'color',o*.5);
text(Xs(1,3)*sc, Xs(2,3)*sc,'3',NVTextHU{:},'color',o*.5);
text(Xu1(1,1)*1.15, Xu1(2,1)*1.15,'4',NVTextHU{:},'color',Css);
text(Xu1(1,2)*1.15, Xu1(2,2)*1.15,'5',NVTextHU{:},'color',Css);

% Text Title
textb = '\textbf{c} \hspace{4mm}node placement 2';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});


%% c: Added nodes 3
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Place extra nodes
th1 = 50; th2 = 210;
Rth1 = rotz(th1); Rth1 = Rth1(1:2,1:2);
Rth2 = rotz(th2); Rth2 = Rth2(1:2,1:2);
Xu1 = [Rth1*[0;1] Rth2*[0;1]]*rM;
D = sqrt(sum((Xs(:,conn(:,1)) - Xu1(:,conn(:,2)-3)).^2));
M = (Xs(:,conn(:,1)) + Xu1(:,conn(:,2)-3))/2;
MD = Xu1(:,conn(:,2)-3) - Xs(:,conn(:,1));
MDA = atan2d(MD(2,:),MD(1,:));

% Draw solution space
visualize_conic_finite(Xs,Xf,R,'ucolori',Css,'overlay',.2);
% Draw network
visualize_network(Xs,Xu1,conn,'msize',mSs,'ucolor',Css);

% Draw
axis(ax);
set(gca,'visible',0);
drawnow;

% Text Angles
r = [.6 .3 .4 .5 .75 .4];
xsh = [-.4 .3 0 -.35 .4 .45];
ysh = [0 .15 -.43 -.2 0 0];
for i = 1:size(conn,1)
    Rprox = rotz(MDA(i)); Rprox = Rprox(1:2,1:2);
    rprox = Rprox*[-1;0]*r(i)*D(i);
    text(rprox(1)+Xu1(1,conn(i,2)-3)+xsh(i),...
         rprox(2)+Xu1(2,conn(i,2)-3)+ysh(i),...
         sprintf('%0.3g',D(i)),NVTextHU{:});
end

% Text Nodes
sc = 1.3;
text(Xs(1,1)*sc, Xs(2,1)*sc,'1',NVTextHU{:},'color',o*.5);
text(Xs(1,2)*sc, Xs(2,2)*sc,'2',NVTextHU{:},'color',o*.5);
text(Xs(1,3)*sc, Xs(2,3)*sc,'3',NVTextHU{:},'color',o*.5);
text(Xu1(1,1)*1.15, Xu1(2,1)*1.15,'4',NVTextHU{:},'color',Css);
text(Xu1(1,2)*1.15, Xu1(2,2)*1.15,'5',NVTextHU{:},'color',Css);

% Text Title
textb = '\textbf{d} \hspace{4mm}node placement 3';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});


%% Save
fName = 'suppfig_edges';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');