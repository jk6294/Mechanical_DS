% Supplementary Figure: edge lengths
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 4.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  7.25  4.50];...
        [ 7.25  0.00  7.25  4.50];...
        [14.50  0.00  4.50  4.50]];
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
NVTextHUC  = {'fontsize',FS,'HorizontalAlignment','center',...
              'VerticalAlignment','middle'};
NVTextRA = {'Units','Normalized','fontsize',FS,...
            'HorizontalAlignment','right'};
NVTextRAU = {'fontsize',FS,'HorizontalAlignment','right'};
NVTextR  = {'Units','Normalized','fontsize',FS};
NVTextRU  = {'fontsize',FS};

% Color
Css = [031 172 204]/255;
Cll = [204 120 161]/255;


%% Global parameters
% Construct network
s = sqrt(3);
Xs = [-s/2 0 s/2;...                    % Initial position
      -1/2 1 -1/2];
rM = 2;
th1 = 90 - acosd(s/4); th2 = 90 + acosd(s/4);
Rth1 = rotz(th1); Rth1 = Rth1(1:2,1:2);
Rth2 = rotz(th2); Rth2 = Rth2(1:2,1:2);
Xu = [Rth1*[0;1] Rth2*[0;1]]*rM;
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];  % Network connectivity

% Parameters
R = [-2 2; -2 2];               % Range of conic visualization
mSs = 5;                        % Node size: start nodes
mSf = 3;                        % Node size: end nodes
ax = [[-1 1]*sRat(1) -1 1]*3 + [.3 .3 .7 .7];
ax2 = [-1 1 -1 1]*3 + [.3 .3 .7 .7];
o = [1 1 1];


%% a: Two units
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Draw Network
sh1 = [-3.2;0];
sh2 = [1.5;0];
visualize_network(Xs+sh1,Xu+sh1,conn,'msize',mSs,'ucolor',Css);
visualize_network(Xs+sh2,Xu+sh2,conn,'msize',mSs,'ucolor',Css);

% Lengths
line_coordinates(Xs(:,1:2)+sh1, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs(:,2:3)+sh1, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs(:,1:2)+sh2, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs(:,2:3)+sh2, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);

% Text Units
text(sh1(1),3.2,'unit $k$',NVTextHU{:});
text(sh2(1),3.2,'unit $k+1$',NVTextHU{:});

% Text Nodes
sc = 1.4;
sc2 = 1.27;
text(Xs(1,1)*sc+sh1(1), Xs(2,1)*sc,'$1$',NVTextHUC{:},'color',o*.5);
text(Xs(1,2)*sc+sh1(1), Xs(2,2)*sc,'$2$',NVTextHUC{:},'color',o*.5);
text(Xs(1,3)*sc+sh1(1), Xs(2,3)*sc,'$3$',NVTextHUC{:},'color',o*.5);
text(Xu(1,1)+sh1(1), Xu(2,1)*sc2,'$4$',NVTextHUC{:},'color',Css);
text(Xu(1,2)+sh1(1), Xu(2,2)*sc2,'$5$',NVTextHUC{:},'color',Css);
text(Xs(1,1)*sc+sh2(1), Xs(2,1)*sc,'$1$',NVTextHUC{:},'color',o*.5);
text(Xs(1,2)*sc+sh2(1), Xs(2,2)*sc,'$2$',NVTextHUC{:},'color',o*.5);
text(Xs(1,3)*sc+sh2(1), Xs(2,3)*sc,'$3$',NVTextHUC{:},'color',o*.5);
text(Xu(1,1)+sh2(1), Xu(2,1)*sc2,'$4$',NVTextHUC{:},'color',Css);
text(Xu(1,2)+sh2(1), Xu(2,2)*sc2,'$5$',NVTextHUC{:},'color',Css);

% Text Lengths
text(mean(Xs(1,1:2))+sh1(1)-.9,mean(Xs(2,1:2)),'$l_k$',NVTextHUC{:},'color',Cll);
text(mean(Xs(1,2:3))+sh1(1)+1.1,mean(Xs(2,2:3)),'$l_{k+1}$',NVTextHUC{:},'color',Cll);
text(mean(Xs(1,1:2))+sh2(1)-1.1,mean(Xs(2,1:2)),'$l_{k+1}$',NVTextHUC{:},'color',Cll);
text(mean(Xs(1,2:3))+sh2(1)+1.1,mean(Xs(2,2:3)),'$l_{k+2}$',NVTextHUC{:},'color',Cll);

% Draw
axis(ax);
set(gca,'visible',0);
drawnow;

% Text Title
texta = '\textbf{a} \hspace{17mm}two units';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});


%% b: Flip second unit
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Draw Network
sh1 = [-3.2;0];
sh2 = [1.5;0];
Xs2 = [Xs(1,:);0.5-Xs(2,:)];
Xu2 = [Xu(1,:);0.5-Xu(2,:)];
visualize_network(Xs+sh1,Xu+sh1,conn,'msize',mSs,'ucolor',Css);
visualize_network(Xs2+sh2,Xu2+sh2,conn,'msize',mSs,'ucolor',Css);

% Lengths
line_coordinates(Xs(:,1:2)+sh1, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs(:,2:3)+sh1, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs2(:,1:2)+sh2, 'lSh',-.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs2(:,2:3)+sh2, 'lSh',-.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);

% Lengths Connecting
line_coordinates([Xs(:,2)+sh1 Xs2(:,1)+sh2], 'lSh',0, 'lw',.5, 'color',o*.5);
line_coordinates([Xs(:,3)+sh1 Xs2(:,2)+sh2], 'lSh',0, 'lw',.5, 'color',o*.5);

% Text Units
text(sh1(1),3.2,'unit $k$',NVTextHU{:});
text(sh2(1),3.2,'unit $k+1$',NVTextHU{:});

% Text Nodes
sc = 1.4;
sc2 = 1.25;
text(Xs(1,1)*sc+sh1(1), Xs(2,1)*sc,'$1$',NVTextHUC{:},'color',o*.5);
text(Xs(1,2)*sc+sh1(1), Xs(2,2)*sc,'$2$',NVTextHUC{:},'color',o*.5);
text(Xs(1,3)*sc+sh1(1), Xs(2,3)*sc,'$3$',NVTextHUC{:},'color',o*.5);
text(Xu(1,1)+sh1(1), Xu(2,1)*sc2,'$4$',NVTextHUC{:},'color',Css);
text(Xu(1,2)+sh1(1), Xu(2,2)*sc2,'$5$',NVTextHUC{:},'color',Css);
text(Xs2(1,1)*sc+sh2(1), Xs2(2,1)*sc-.2,'$1$',NVTextHUC{:},'color',o*.5);
text(Xs2(1,2)*sc+sh2(1), Xs2(2,2)*sc-.2,'$2$',NVTextHUC{:},'color',o*.5);
text(Xs2(1,3)*sc+sh2(1), Xs2(2,3)*sc-.2,'$3$',NVTextHUC{:},'color',o*.5);
text(Xu2(1,1)+sh2(1), Xu2(2,1)*sc2-.2,'$4$',NVTextHUC{:},'color',Css);
text(Xu2(1,2)+sh2(1), Xu2(2,2)*sc2-.2,'$5$',NVTextHUC{:},'color',Css);

% Text Join
text((Xs(1,2)+Xs2(1,1)+sh1(1)+sh2(1))/2,...
     (Xs(2,2)+Xs2(2,1)+sh1(2)+sh2(2))/2+.4,...
     'join 2,1', NVTextHU{:},'color',o*.5);
text((Xs(1,3)+Xs2(1,2)+sh1(1)+sh2(1))/2,...
     (Xs(2,3)+Xs2(2,2)+sh1(2)+sh2(2))/2-.35,...
     'join 3,2', NVTextHU{:},'color',o*.5);

% Text Lengths
text(mean(Xs(1,1:2))+sh1(1)-.9,mean(Xs(2,1:2)),'$l_k$',NVTextHUC{:},'color',Cll);
text(mean(Xs(1,2:3))+sh1(1)+1.1,mean(Xs(2,2:3)),'$l_{k+1}$',NVTextHUC{:},'color',Cll);
text(mean(Xs2(1,1:2))+sh2(1)-1.1,mean(Xs2(2,1:2)),'$l_{k+1}$',NVTextHUC{:},'color',Cll);
text(mean(Xs2(1,2:3))+sh2(1)+1.1,mean(Xs2(2,2:3)),'$l_{k+2}$',NVTextHUC{:},'color',Cll);

% Draw
axis(ax);
set(gca,'visible',0);
drawnow;

% Text Title
texta = '\textbf{b} \hspace{13mm}flip unit $k+1$';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});


%% c: Flip second unit
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Draw Network
sh1 = [-s/4;0];
sh2 = [s/4;0];
Xs2 = [Xs(1,:);0.5-Xs(2,:)];
Xu2 = [Xu(1,:);0.5-Xu(2,:)];
visualize_network(Xs+sh1,Xu+sh1,conn,'msize',mSs,'ucolor',Css);
visualize_network(Xs2+sh2,Xu2+sh2,conn,'msize',mSs,'ucolor',Css);

% Lengths
line_coordinates(Xs(:,1:2)+sh1, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs(:,2:3)+sh1, 'lSh',.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);
line_coordinates(Xs2(:,2:3)+sh2, 'lSh',-.4,'lw',.5,'style','-', 'nw',.1,'color',Cll);

% Lengths Connecting
line_coordinates([Xs(:,2)+sh1 Xs2(:,1)+sh2], 'lw',.5, 'color',o*.5);
line_coordinates([Xs(:,3)+sh1 Xs2(:,2)+sh2], 'lw',.5, 'color',o*.5);

% Text Lengths
text(mean(Xs(1,1:2))+sh1(1)-.9,mean(Xs(2,1:2)),'$l_k$',NVTextHUC{:},'color',Cll);
text(mean(Xs(1,2:3))+sh1(1)+.6,mean(Xs(2,2:3))+.7,'$l_{k+1}$',NVTextHUC{:},'color',Cll);
text(mean(Xs2(1,2:3))+sh2(1)+1.1,mean(Xs2(2,2:3)),'$l_{k+2}$',NVTextHUC{:},'color',Cll);

% Draw
axis(ax2);
set(gca,'visible',0);
drawnow;

% Text Title
texta = '\textbf{c} \hspace{5mm}combined network';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});


%% Save
fName = 'suppfig_combine';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');