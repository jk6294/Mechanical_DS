% Figure 6: Limit cycle
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(6); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19  9.50];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .0 .4 .0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  4.75  4.75  4.75];...
        [ 0.00  0.00  4.75  4.75];...
        [ 5.25  7.25  4.75  2.25];...
        [10.25  7.25  2.50  2.25];...
        [ 5.25  0.05  7.25  7.25];...
        [13.25  2.50  5.75  7.00];...
        [13.25  0.00  5.75  2.50]];
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
gr = 0.8;           % Gray
al = 0.2;           % Alpha

% Default Name-Value Pairs
NVTitle = {'Units','centimeters','fontsize',FS};
NVTextH = {'Units','Normalized','fontsize',FS,...
           'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,...
           'HorizontalAlignment','right'};
NVTextR = {'Units','Normalized','fontsize',FS};

% Color
nT = 1000;
o = [1 1 1];
% Gradient: Distance d1,d2. Interpolate between 3 colors
DLin = linspace(.7,1.6,nT);
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));
% Gradient: Slope. Interpolate between 2 colors
SLin = linspace(0,1,nT);
C3a = [197 066 085]/255;
C3b = [031 172 204]/255;
CP3 = interp1([0 1],[C3a;C3b],linspace(0,1,nT));


%% Some symbolic
% Positions
L = 1;
l = .5;
h = .5;
Xs = [-L  l  L;...
       0  h  0];
Xf = [-L -l  L;...
       0  h  0];
syms('x','real');
syms('y','real');

% Solve for initial position slope
x4 = [1;-l/h]*x;
x5 = [0;y];
% Motions
dx1 = diff(Xs(:,1:2),1,2); dx1 = dx1/norm(dx1);
dx2 = [0;0];
v3  = diff(Xs(:,2:3),1,2); v3 = v3/norm(v3);
% Solve for motion of nodes 3,4,5
dx4 = (Xs(:,1:2)-x4)'\(sum((Xs(:,1:2) - x4).*[dx1 dx2])');
dx5 = (Xs(:,1:2)-x5)'\(sum((Xs(:,1:2) - x5).*[dx1 dx2])');
dx3 = (Xs(:,3)-[x4 x5])'\(sum((Xs(:,3)-[x4 x5]).*[dx4 dx5])');
% Slope
sl0 = simplify(v3'*dx3);

% Solve for final position slope
x4 = [1;l/h]*x;
x5 = [0;y];
% Motions
dx1 = diff(Xf(:,1:2),1,2); dx1 = dx1/norm(dx1);
dx2 = [0;0];
v3  = diff(Xf(:,2:3),1,2); v3 = v3/norm(v3);
% Solve for motion of nodes 3,4,5
dx4 = (Xf(:,1:2)-x4)'\(sum((Xf(:,1:2) - x4).*[dx1 dx2])');
dx5 = (Xf(:,1:2)-x5)'\(sum((Xf(:,1:2) - x5).*[dx1 dx2])');
dx3 = (Xf(:,3)-[x4 x5])'\(sum((Xf(:,3)-[x4 x5]).*[dx4 dx5])');
% Slope
slf = simplify(v3'*dx3);

% Product of slopes
sl = sl0*slf;
slfun = matlabFunction(sl);

% Sample product of slopes to create phase diagram
xl = linspace(-.75,.75,1000);
yl = linspace(-1.5,0,1000);
[X,Y] = meshgrid(xl,yl);
slM = arrayfun(slfun,X',Y');


%% Network Contsruction
% Design network unit
Xu1 = [0; -1];
Xu2 = [1;-l/h]*.25;
Xu = [Xu1 Xu2];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xub,fV] = construct_network(Xs,Xf,Xu,conn,0,1);
% Initial and final distances
ds = sqrt(sum(diff(Xs,1,2).^2));
df = sqrt(sum(diff(Xf,1,2).^2));

% Combined unit
Xsc = [-L -l  L  2*L-l;...
        0  h  0  h];
Xuc = [Xub(3:4,:) [Xub(1,:)+L-l;-Xub(2,:)+h]];
connc = [1 5; 2 5; 3 5; 1 6; 2 6; 3 6;...
         2 7; 3 7; 4 7; 2 8; 3 8; 4 8];
     
% Combined network
[Xsa,Xua,conna] = tesselate_network_old(Xsc,Xuc, connc, [2*L;0], [7;1]);

% Distances
nSim = 990;
% Simulate for distances
[XCa,fCa] = sim_motion8(Xs,Xu,conn,.001,nSim,-[Xs Xu],0);
% Realign
for i = 1:nSim
    Rz = rotz(atan2d(XCa(2,3,i)-XCa(2,1,i),XCa(1,3,i)-XCa(1,1,i)));
    XCa(:,:,i) = Rz(1:2,1:2)'*XCa(:,:,i);
end
D = squeeze(sqrt(sum(diff(XCa(:,1:3,:),1,2).^2)));

% Simulate 2-module
nSim = 300;
[XCca,fCca] = sim_motion8(Xsc,Xuc,connc,.001,nSim,[Xsc Xuc],0);
Dc = squeeze(sqrt(sum(diff(XCca(:,(1:4),:),1,2).^2)));
Dc = Dc([1 3],:);

% Simulate large network
nSim = 300;
[XCaa,fCaa] = sim_motion8(Xsa,Xua,conna,.01,nSim,[Xsa Xua],0);
Da = squeeze(sqrt(sum(diff(XCaa(:,1:size(Xsa,2),:),1,2).^2)));



%% a: Design
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .175;                  % Scale network
sh = [.45;.2];              % Shift solution space
sh1 = [.15;.53];            % Shift start positions
sh2 = [.72;.53];            % Shift end positions
shl = [-.06;.84];           % Shift legend
mS = 5;                    % Size of start node markers
mF = 3;                    % Size of end node markers
lSh = .008;                 % Distance bar offset
R = [-1 1; -1 1]*1;         % Range of drawing solution space

% Draw initial network and lengths
XsP = Xs*sc;
XfP = Xf*sc;
XuP = Xu*sc;
line_coordinates(XsP(:,1:2)+sh1,'lSh',lSh,'color',interp1(DLin,CP1,ds(1)));
line_coordinates(XsP(:,2:3)+sh1,'lSh',lSh,'color',interp1(DLin,CP1,ds(2)));
visualize_network(XsP+sh1,[],[1 1],'msize',mS);
% Draw final network and lengths
line_coordinates(XfP(:,1:2)+sh2,'lSh',lSh,'color',interp1(DLin,CP1,df(1)));
line_coordinates(XfP(:,2:3)+sh2,'lSh',lSh,'color',interp1(DLin,CP1,df(2)));
visualize_network(XfP+sh2,[],[1 1],'scolor',o,'msize',mF);
% Draw solution space
visualize_conic_finite(XsP+sh,XfP+sh,R*sc+sh,'overlay',.99,'ucolori',o*gr);
visualize_network(XsP+sh,[],[1 1],'msize',mS);
visualize_network(XfP+sh,[],[1 1],'scolor',o,'msize',mF);
visualize_network(XuP+sh,[],[1 1],'scolor',interp1(SLin,CP3,0),'msize',mS);
% Draw legend
visualize_network([0;0]+shl,[],[1 1],'msize',mS);
visualize_network([0.31;0]+shl,[],[1 1],'scolor',o,'msize',mF);
line([0 .1/sc]*sc+shl(1)+.56, [0 0]*sc+shl(2),'linewidth',1,'color',o*gr);

% Spline
sCo11 = [0.00 0.4 0.50 1.00;...
         0.00 0.17 0.83 1.00] .* [.12;.107] + [.3;.05];
sCo12 = [0.00 0.30 0.60 1.00;...
         0.00 0.20 0.80 1.00] .* [.12;-.021] + [.3;.05];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),'color',o*gr^.5);
plot_spline(sCo12,'head',1,'headpos',1,'ratio',sRat(pInd),'color',o*gr^.5);
hold off;

% Text
texta = '\textbf{a} \hspace{6mm}find solution space';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(sh1(1),sh1(2)+.2,'start position',NVTextH{:});
text(sh2(1),sh2(2)+.2,'end position',NVTextH{:});
text(sh(1),sh(2)+.21,'solution space',NVTextH{:});
text(-.01,shl(2),'start',NVTextR{:},'VerticalAlignment','middle');
text(.29,shl(2),'end',NVTextR{:},'VerticalAlignment','middle');
text(.9,shl(2),'solution',NVTextRA{:},'VerticalAlignment','middle');
text(0,sh1(2)+.1,'$l_1^\circ$',NVTextR{:},'color',interp1(DLin,CP1,ds(1)));
text(.3,sh1(2)+.1,'$l_2^\circ$',NVTextR{:},'color',interp1(DLin,CP1,ds(2)));
text(.51,sh2(2)+.1,'$l_1^*$',NVTextR{:},'color',interp1(DLin,CP1,df(1)));
text(.75,sh2(2)+.1,'$l_2^*$',NVTextR{:},'color',interp1(DLin,CP1,df(2)));
text(sh(1),sh(2)-.23,'$y_1$',NVTextH{:});
text(sh(1)+.2,sh(2)-.23,'$x_2$',NVTextH{:});
text(-.08,sh(2)-.15,'add nodes',NVTextR{:},'color',o*gr);
text(-.08,sh(2)-.23,'at $y_1,x_2$',NVTextR{:},'color',o*gr);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% b: Map
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .25;                       % Scale network
sh1 = [.25;.2];                 % Shift unit 1
sh2 = [.15;.3];                 % Shift unit 2
sh3 = [.25;.25];                % Shift unit 3
ax = [.5 2.3 .5 2.3];           % Figure axis
axL = ax + [.05 -.3 .05 -.3];   % Map coordinate lines axis
lSh = .02;                      % Distance bar offset

% Draw tick marks
line([1;1].*df(1), [0;.04].*[1 1]+axL(3),'linewidth',.5,'color',...
     interp1(DLin,CP1,df(1)));
line([1;1].*df(2), [0;.04].*[1 1]+axL(3),'linewidth',.5,'color',...
     interp1(DLin,CP1,df(2)));
line([0;.04].*[1 1]+axL(3), [1;1].*df(1),'linewidth',.5,'color',...
     interp1(DLin,CP1,df(1)));
line([0;.04].*[1 1]+axL(3), [1;1].*df(2),'linewidth',.5,'color',...
     interp1(DLin,CP1,df(2)));

% Draw map
plot(axL(1:2),axL(3:4),'color',o*.8,'linewidth',.5);
plot(D(1,:),D(2,:),'color','k','linewidth',1);

% Draw axis
line([axL(1) axL(1) axL(2)],[axL(4) axL(3) axL(3)],...
     'linewidth',.5,'color','k');
arrow(axL(1:2),[1 1]*axL(3),sRat(pInd),'color','k');
arrow([1 1]*axL(1),axL(3:4),sRat(pInd),'color','k');

% Draw example units at fixed point and limit cycle points
[~,fpN] = min(abs(D(1,:)-D(2,:)));
dfp = D(1,fpN);
XCa1 = XCa(:,:,1)*sc+ds'+sh1;
plot(ds,df,'ks','linewidth',3,'markersize',3);
line_coordinates(XCa1(:,1:2),'lSh',lSh,'color',interp1(DLin,CP1,ds(1)));
line_coordinates(XCa1(:,2:3),'lSh',lSh,'color',interp1(DLin,CP1,ds(2)));
visualize_network(XCa1(:,1:3),XCa1(:,4:5),conn,...
                  'ucolor',interp1(SLin,CP3,tanh(0)));
XCa2 = XCa(:,:,end)*sc+df'+sh2;
line_coordinates(XCa2(:,1:2),'lSh',lSh,'color',interp1(DLin,CP1,df(1)));
line_coordinates(XCa2(:,2:3),'lSh',lSh,'color',interp1(DLin,CP1,df(2)));
visualize_network(XCa2(:,1:3),XCa2(:,4:5),conn,...
                  'ucolor',interp1(SLin,CP3,tanh(0)));
XCafp = XCa(:,:,fpN)*sc+dfp'+sh3;
plot(dfp,dfp,'ks','linewidth',3,'markersize',3);
line_coordinates(XCafp(:,1:2),'lSh',lSh,'color',interp1(DLin,CP1,dfp));
line_coordinates(XCafp(:,2:3),'lSh',lSh,'color',interp1(DLin,CP1,dfp));
visualize_network(XCafp(:,1:3),XCafp(:,4:5),conn,...
                  'ucolor',interp1(SLin,CP3,tanh(0)));

% Text
textb = '\textbf{b} \hspace{2mm}conformational motion';
text(labX,subp(pInd,4)+labY-.2,textb,NVTitle{:});
text(.86,.035,'$l_1$',NVTextR{:});
text(-.05,.8,'$l_2$',NVTextR{:});
text(.68,-.04,'$\sqrt{10}/2$',NVTextRA{:},'color',interp1(DLin,CP1,ds(1)));
text(.2,-.04,'$\sqrt{2}/2$',NVTextRA{:},'color',interp1(DLin,CP1,ds(2)));
text(.39,.1,'start',NVTextR{:});
text(.06,.52,'end',NVTextR{:});

axis(ax);
set(gca,'visible',0);
drawnow;


%% c: Periodicity
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .2;                            % Scale network
shL = [.45;0];                      % Shift network: iterative
sh = [.1;.4]-shL;                   % Shift network: constant offset
dotX = linspace(0,shL(1),10);       % Position of dots connecting nodes
lSh = .03;                          % Distance bar offset

% Draw 4 example units
for i = 1:4
    if(i<4)
        scatter(Xsa(1,1+i)*sc+sh(1)+shL(1)*i+dotX,...
                Xsa(2,1+i)*sc+sh(2)+shL(2)*i+zeros(1,10),.2,o*.9);
        scatter(Xsa(1,2+i)*sc+sh(1)+shL(1)*i+dotX,...
                Xsa(2,2+i)*sc+sh(2)+shL(2)*i+zeros(1,10),.2,o*.9);
    else
        scatter(Xsa(1,1+i)*sc+sh(1)+shL(1)*i+dotX,...
                Xsa(2,1+i)*sc+sh(2)+shL(2)*i+zeros(1,10),.2,...
                o.*linspace(.9,1,length(dotX))','clipping',0);
        scatter(Xsa(1,2+i)*sc+sh(1)+shL(1)*i+dotX,...
                Xsa(2,2+i)*sc+sh(2)+shL(2)*i+zeros(1,10),.2,...
                o.*linspace(.9,1,length(dotX))','clipping',0);
    end
    line_coordinates(Xsa(:,(0:1)+i)*sc+sh+shL*i,'lSh',lSh*(-1)^(i-1),...
                     'color',interp1(DLin,CP1,ds(mod(i,2)+1)));
    visualize_network(Xsa(:,(0:2)+i)*sc+sh+shL*i,...
                      Xua(:,(-1:0)+2*i)*sc+sh+shL*i,conn,...
                     'ucolor',interp1(SLin,CP3,0));
end

% Draw boxes around each 2-cycle unit combination
bD = [-.3 .95 .1 .79];
line(bD([1 1 2 2 1]), bD([3 4 4 3 3]),'color',o*.8,'linewidth',.5,...
     'clipping',0);
line(bD([1 1 2 2 1])+1.32, bD([3 4 4 3 3]),'color',o*.8,'linewidth',.5,...
     'clipping',0);

% Text
textc = '\textbf{c} \hspace{3mm}lengths repeat every 2 units';
text(labX-.5,subp(pInd,4)+labY,textc,NVTitle{:});
text(-.1,.65,'$l_1$',NVTextR{:},'color',interp1(DLin,CP1,df(1)));
text(.2,.25,'$l_2$',NVTextR{:},'color',interp1(DLin,CP1,ds(1)));
text(.46,.65,'$l_3$',NVTextR{:},'color',interp1(DLin,CP1,df(1)));
text(.76,.25,'$l_4$',NVTextR{:},'color',interp1(DLin,CP1,ds(1)));

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% d: Stability
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .2;                                % Scale network
sh = [.55;.45];                         % Shift network
Rz = rotz(3); Rz = Rz(1:2,1:2);         % Rotate network

% Draw network start and end
visualize_network(Rz*XCca(:,1:4,end)*sc+sh-[0;.01],...
                  XCca(:,5:8,end)*sc+sh,connc,'nalpha',.3,'lalpha',.3,...
                  'ucolor',interp1(SLin,CP3,0));
visualize_network(XCca(:,1:4,1)*sc+sh,XCca(:,5:8,1)*sc+sh,connc,...
                  'ucolor',interp1(SLin,CP3,0));
line_coordinates(Rz*XCca(:,1:2,end)*sc+sh,'lSh', .1,'nw',.03,...
                 'lalpha',.3,'style','-','lw',.5);
line_coordinates(Rz*XCca(:,3:4,end)*sc+sh,'lSh',-.1,'nw',.03,...
                 'lalpha',.3,'style','-','lw',.5);
line_coordinates(XCca(:,1:2,1)*sc+sh-[0;.01],'lSh',.15,'nw',.03,...
                 'style','-','lw',.5);
line_coordinates(XCca(:,3:4,1)*sc+sh-[0;.01],'lSh',-.15,'nw',.03,...
                 'style','-','lw',.5);

% Text
textd = '\textbf{d} \hspace{2.8mm}stability';
text(labX+.2,subp(pInd,4)+labY,textd,NVTitle{:});
text(.5,.05,'$|\delta l_3 / \delta l_1|$',NVTextH{:});
text(.1,.67,'$\delta l_1$',NVTextH{:});
text(.9,.28,'$\delta l_3$',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% e: Phase
pInd = 5;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .075;                  % Scale network

% Draw phase diagram
colormap(CP3);
imagesc(yl,xl,tanh(abs(slM)),'alphadata',.6);

% Draw metastable contour lines
C1 = contour(yl,xl,abs(slM),[1 1],'color',interp1(SLin,CP3,tanh(1)),...
             'linewidth',2);

% Draw superstable contour lines
plot(yl,1/4*ones(1,length(yl)),':','linewidth',1,'color',CP3(1,:)); 
plot(yl,yl./(1-yl),':','linewidth',1,'color',CP3(1,:)); 
plot(yl,-yl./(1+yl),':','linewidth',1,'color',CP3(1,:));

% Draw superstable examples
nP0 = [[  1 150 300]+333;...
       [700 750 800]+ 29;...
       [  1 150 300]+333];
nSh0 = cat(3,[-0.12 -0.12 -0.12;-0.14 -0.14 -0.14],...
             [ 0.13  0.13  0.13;-0.05 -0.05 -0.05],...
             [-0.12 -0.12 -0.12; 0.10  0.07  0.03]);
sL0 = ['k*';'k+';'kx'];
for j = 1:3
    nP = nP0(j,:);
    for i = 1:length(nP)
        if(j==1);       yV= 1/4;                      nal1 = .3; nal2 =  1;
        elseif(j==2);   yV=-yl(nP(i))./(1+yl(nP(i))); nal1 =  1; nal2 = .3;
        else;           yV= yl(nP(i))./(1-yl(nP(i))); nal1 = .3; nal2 =  1;
        end
        sh = [yl(nP(i));yV] + nSh0(:,i,j);
        Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
        Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
        Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
        plot(yl(nP(i)),yV,sL0(j,:),'linewidth',.7,'markersize',5);
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',3,'nalpha',nal1,'lalpha',nal1); 
        visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',3,'nalpha',nal2,'lalpha',nal2); 
    end
end

% Draw metastable examples
nP1 = [[  1 150 300]+333;...
       [  1 150 300]+333;...
       [600 700 800]+ 40];
nSh1 = cat(3,[-0.11 -0.11 -0.11;-0.15 -0.15 -0.15],...
             [-0.12 -0.12 -0.12; 0.12  0.11  0.10],...
             [ 0.10  0.10  0.10;-0.08 -0.08 -0.08]);
sL1 = ['ko';'ks';'kh'];
for j = 1:3
    nP = nP1(j,:);
    for i = 1:length(nP)
        if(j==1);       yV =  0;
        elseif(j==2);   yV = -yl(nP(i))./sqrt(7*yl(nP(i))^2 + 1);
        else;           yV =  yl(nP(i))./(sqrt(2-2*yl(nP(i))^2));
        end
        sh = [yl(nP(i));yV] + nSh1(:,i,j);
        Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
        Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
        Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
        plot(yl(nP(i)),yV,sL1(j,:),'linewidth',.7,'markersize',5);
        visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,'ucolor',...
            interp1(SLin,CP3,tanh(1)),'msize',3,'nalpha',1.0,'lalpha',1.0); 
    end
end

% Draw boarder
nSc = 0.0085;
mL = [min(xl) max(xl) min(yl) max(yl)];
line([mL(3) mL(3) mL(4) mL(4) mL(3) mL(3)]+[-1 -1 1 1 -1 -1]*nSc,...
     [mL(1) mL(2) mL(2) mL(1) mL(1) mL(2)]+[-1 1 1 -1 -1  1]*nSc,...
     'color', 'k', 'linewidth',2);

% Colorbar
hb = colorbar('location', 'east','position',[0.27 .05 .021 .71]);
delete(findall(gcf,'type','annotation'));
a = annotation('textbox',hb.Position,'FitBoxToText','off','FaceAlpha',0.4,...
    'EdgeColor',o*0,'BackgroundColor',o);
hb.Ticks = tanh(1); hb.TickLabels = {};

% Text
texte = '\textbf{e} \hspace{6mm}phase diagram of stability';
text(labX-.5,subp(pInd,4)+labY+.4,texte,NVTitle{:});
text(.02,.31,'unstable',NVTextR{:});
text(.02,.69,'stable',NVTextR{:});
text(.5,-.03,'$y_1$',NVTextH{:});
text(0,-.035,'$-3/2$',NVTextH{:});
text(1,-.035,'$0$',NVTextH{:});
text(1.04,.5,'$x_2$',NVTextH{:});
text(1.04,0,'$\frac{-3}{4}$',NVTextH{:});
text(1.04,1,'$\frac{3}{4}$',NVTextH{:});
text(-.1,0,'$0$',NVTextH{:},'color',interp1(SLin,CP3,0));
text(-.1,.75,'$1$',NVTextH{:},'color',interp1(SLin,CP3,tanh(1)));
text(-.105,.985,'$\infty$',NVTextH{:},'color',interp1(SLin,CP3,1));
annotation('line',[0 .0195]+.2707, [1 1]*.591, 'Units','Normalize',...
           'color',interp1(SLin,CP3,tanh(1)),'linewidth',2);
annotation('line',[0 .0195]+.2707, [1 1]*.0555, 'Units','Normalize',...
           'color',interp1(SLin,CP3,tanh(0)),'linewidth',2);
annotation(gcf,'textarrow',[1 1]*.32, [1 1]*.25,'String','$|\delta d_3 / \delta d_1|$',...
    'HeadStyle', 'none', 'LineStyle', 'none',NVTextH{:}, 'TextRotation',90,...
    'interpreter','latex','color','k');
annotation(gcf,'textarrow',[1 1]*.32, [1 1]*.5,'String','~~~stable',...
    'HeadStyle', 'none', 'LineStyle', 'none',NVTextH{:}, 'TextRotation',90,...
    'interpreter','latex','color','k');
annotation(gcf,'textarrow',[1 1]*.32, [1 1]*.675,'String','unstable',...
    'HeadStyle', 'none', 'LineStyle', 'none',NVTextH{:}, 'TextRotation',90,...
    'interpreter','latex','color','k');

% Axes
axis([min(yl) max(yl) min(xl) max(xl)] + [-1 1 -1 1]*nSc);
set(gca,'xtick',[],'ytick',[],'ydir','normal','visible',0);
drawnow;


%% f,g: Examples
pInd = 6;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .09;                   % Scale network
lSh = .0134;                % Distance bar offset

% Draw metastable examples
tY = [.84 .51 .19];
nY = [.66 .33 .06];
sLT1 = {'symmetric','asymmetric 1','asymmetric 2'};
for j = 1:3
    nP = nP1(j,:);
    for i = 1:length(nP)
        if(j==1);       yV =  0;
        elseif(j==2);   yV = -yl(nP(i))./sqrt(7*yl(nP(i))^2 + 1);
        else;           yV =  yl(nP(i))./(sqrt(2-2*yl(nP(i))^2));
        end
        sh = [.15;nY(j)];
        Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
        Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
        Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
        if(i < length(nP))
            visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,[1 1],...
                'ucolor',interp1(SLin,CP3,tanh(1)),...
                'msize',4,'nalpha',.1*i,'lalpha',.1*i); 
        else
            visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,...
                'ucolor',interp1(SLin,CP3,tanh(1)),...
                'msize',4,'nalpha',i/length(nP),'lalpha',i/length(nP)); 
        end
    end
    plot(0,tY(j),sL1(j,:),'linewidth',.7,'markersize',5);
    text(.05,tY(j),sLT1{j},NVTextR{:});
end

% Draw superstable examples
text(.67,tY(1),'1 co-linear',NVTextR{:});
plot(.5,tY(1),'k*','linewidth',.7,'markersize',5);
nP = fliplr(nP0(1,:));
for i = 1:length(nP)
    yV =  1/4;
    sh = [.62;nY(1)];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    if(i < length(nP))
        visualize_network(Xf*sc+sh,Xupf*sc+sh,[1 1],'ucolor',CP3(1,:),...
            'msize',4,'nalpha',.1*i,'lalpha',.1*i); 
    else
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',4,'nalpha',.1,'lalpha',.1); 
        visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',4,'nalpha',1,'lalpha',1); 
        XfM = mean(Xf(:,2:3),2); 
        XfD = diff(Xf(:,2:3),1,2);
        line_coordinates(([-XfD XfD]*.8+XfM)*sc+sh,'lSh',lSh,...
                         'style','-','lw',.5,'color',CP3(1,:));
    end
end
text(.67,tY(2),'2 co-linear',NVTextR{:});
plot(.5,tY(2),'kx','linewidth',.7,'markersize',5);
nP = nP0(3,:);
for i = 1:length(nP)
    yV =  yl(nP(i))./(1-yl(nP(i)));
    sh = [.62;nY(2)+.04];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    XfM = mean([Xf(:,1) Xupf(:,1)],2); 
    XfD = diff([Xf(:,1) Xupf(:,1)],1,2);
    if(i < length(nP))
        visualize_network(Xf*sc+sh,Xupf*sc+sh,[1 1],'ucolor',CP3(1,:),...
            'msize',4,'nalpha',.1*i,'lalpha',.1*i); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',-lSh,...
                         'style','-','lw',.5,'color',CP3(1,:),'lalpha',.3);
    else
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',4,'nalpha',.1,'lalpha',.1); 
        visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',4,'nalpha',1,'lalpha',1); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',-lSh,...
                         'style','-','lw',.5,'color',CP3(1,:));
    end
end
text(.67,tY(3),'2 co-linear',NVTextR{:});
plot(.5,tY(3),'k+','linewidth',.7,'markersize',5);
nP = nP0(2,:);
for i = 1:length(nP)
    yV =  -yl(nP(i))./(1+yl(nP(i)));
    sh = [.62;nY(3)-.04];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    XfM = mean([Xsp(:,1) Xup0(:,2)],2); 
    XfD = diff([Xsp(:,1) Xup0(:,2)],1,2);
    if(i < length(nP))
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,[1 1],'ucolor',CP3(1,:),...
            'msize',4,'nalpha',0.1*i,'lalpha',0.1*i); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',lSh,...
                         'style','-','lw',.5,'color',CP3(1,:),'lalpha',.3);
    else
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',4,'nalpha',1,'lalpha',1); 
        visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),...
            'msize',4,'nalpha',.1,'lalpha',.1); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',lSh,...
                         'style','-','lw',.5,'color',CP3(1,:));
    end
end

% Text
textf = '\textbf{f} \hspace{3mm}metastable';
text(labX,subp(pInd,4)+labY,textf,NVTitle{:});
textg = '\textbf{g} \hspace{3mm}superstable';
text(labX+3.1,subp(pInd,4)+labY,textg,NVTitle{:});
text(.2,.915,'$|\delta l_3/\delta l_1| = 1$',NVTextH{:});
text(.78,.915,'$|\delta l_3/\delta l_1| = 0$',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'xtick',[],'ytick',[],'visible',0);
drawnow;


%% h: Stable limit cycle
pInd = 7;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .15;                       % Scale network
sh = [.36;.6];                  % Shift network
psc = .34;                      % Scale plot
pSh = -.2;                      % Shift plot

% Draw network
plI = 200;
XCaaP = XCaa(:,:,plI)*sc+sh;
XDa = movmean(XCaaP(1,1:size(Xsa,2)),[0,1],'Endpoints','Discard');
for i = 1:size(Xsa,2)-1
    line_coordinates(XCaaP(:,(0:1)+i),'lsh',.015*(-1)^i,...
                     'color',interp1(DLin,CP1,Da(i,plI)));
end
visualize_network(XCaaP(:,1:size(Xsa,2)),...
                  XCaaP(:,(1:size(Xua,2))+size(Xsa,2)),...
                  conna,'msize',3,'ucolor',interp1(SLin,CP3,tanh(0)));
              
% Draw distance plot
lCo = [ 0.15  0.15  2.41;...
        0.45 -0.04 -0.04];
scatter(XDa,Da(:,plI)'*psc+pSh,3,interp1(DLin,CP1,Da(:,plI)),...
        'marker','s','linewidth',2);
line([1;1].*XDa, [0;.04]+ones(1,size(Xsa,2)-1)*lCo(2,2),...
     'linewidth',.5,'color','k','clipping',0);
line([0;.04]+lCo(1,1),[1;1].*ds(1)*psc+pSh,...
     'linewidth',.5,'color',interp1(DLin,CP1,ds(1)),'clipping',0);
line([0;.04]+lCo(1,1),[1;1].*ds(2)*psc+pSh,...
     'linewidth',.5,'color',interp1(DLin,CP1,ds(2)),'clipping',0);
line(lCo(1,:),lCo(2,:),'linewidth',0.5,'color','k','clipping',0);

% Text
texth = '\textbf{h} \hspace{4mm}stable units approach 2-cycle';
text(labX,subp(pInd,4)+labY+.15,texth,NVTitle{:});
text(.05,ds(1)*psc+pSh,'$\frac{\sqrt{10}}{2}$',NVTextRA{:},...
     'color',interp1(DLin,CP1,ds(1)));
text(.05,ds(2)*psc+pSh,'$\frac{\sqrt{2}}{2}$',NVTextRA{:},...
     'color',interp1(DLin,CP1,ds(2)));
text(.96,lCo(2,2),'$k$',NVTextR{:});
text(XDa(1)/sRat(pInd),-.13,'1',NVTextH{:});
text(XDa(2)/sRat(pInd),-.13,'2',NVTextH{:});
text(XDa(4)/sRat(pInd),-.13,'$\cdots$',NVTextH{:});
text(XDa(end)/sRat(pInd),-.13,num2str(length(XDa)),NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'xtick',[],'ytick',[],'visible',0);
drawnow;


%% Save
fName = 'figure6c';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');