% Figure 6: Dynamical Systems
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Some symbolic
clc;
L = 1;
l = .5;
h = .5;
x1 = [-L;0];
x2 = [l;h];
x3 = [L;0];


%% Some symbolic
clc;
syms x1 x2 x3 x4 x5 y1 y2 y3 y4 y5;
syms dx1 dx2 dx3 dx4 dx5 dy1 dy2 dy3 dy4 dy5;
assume(x1,'real'); assume(dx1,'real');
assume(x2,'real'); assume(dx2,'real');
assume(x3,'real'); assume(dx3,'real');
assume(x4,'real'); assume(dx4,'real');
assume(x5,'real'); assume(dx5,'real');
assume(y1,'real'); assume(dy1,'real');
assume(y2,'real'); assume(dy2,'real');
assume(y3,'real'); assume(dy3,'real');
assume(y4,'real'); assume(dy4,'real');
assume(y5,'real'); assume(dy5,'real');

% Equations for node 1
E1a = (x1-x4)*(dx1-dx4) + (y1-y4)*(dy1-dy4) == 0;
E1b = (x1-x5)*(dx1-dx5) + (y1-y5)*(dy1-dy5) == 0;
% Equations for node 2
E2a = (x2-x4)*(dx2-dx4) + (y2-y4)*(dy2-dy4) == 0;
E2b = (x2-x5)*(dx2-dx5) + (y2-y5)*(dy2-dy5) == 0;
% Equations for node 3
E3a = (x3-x4)*(dx3-dx4) + (y3-y4)*(dy3-dy4) == 0;
E3b = (x3-x5)*(dx3-dx5) + (y3-y5)*(dy3-dy5) == 0;

% Fix motions
dx4 = 0; dy4 = 0;           % Node 4 stationary
dy5 = 0;                    % Node 5 only moves in x direction

% Fix node positions
L = 1; 
l = .5; 
h = .5;


% Initial
x1 = -L; 
x2 =  l; 
x3 =  L; 
y1 = 0;
y2 = h;
y3 = 0;
x5 = 0;
y4 = -l/h*x4;
% Substitute
E1a0 = subs(E1a);
E1b0 = subs(E1b);
E2a0 = subs(E2a);
E2b0 = subs(E2b);
E3a0 = subs(E3a);
E3b0 = subs(E3b);
% Solve
s1 = solve([E1a0;E1b0], [dx1 dy1]);
s2 = solve([E2a0;E2b0], [dx2 dy2]);
s3 = solve([E3a0;E3b0], [dx3 dy3]);
% delta d1
v12 = [x2-x1;y2-y1]; v12 = v12/norm(v12);
v23 = [x2-x3;y2-y3]; v23 = v23/norm(v23);
deld1 = (s1.dx1 - s2.dx2)*v12(1) + (s1.dy1 - s2.dy2)*v12(2);
deld2 = (s3.dx3 - s2.dx2)*v23(1) + (s3.dy3 - s2.dy2)*v23(2);
sl0 = simplify(deld2/deld1);


% Final
x1 = -L; 
x2 = -l; 
x3 =  L; 
y1 = 0;
y2 = h;
y3 = 0;
x5 = 0;
y4 = l/h*x4;
% Substitute
E1af = subs(E1a);
E1bf = subs(E1b);
E2af = subs(E2a);
E2bf = subs(E2b);
E3af = subs(E3a);
E3bf = subs(E3b);
% Solve
s1 = solve([E1af;E1bf], [dx1 dy1]);
s2 = solve([E2af;E2bf], [dx2 dy2]);
s3 = solve([E3af;E3bf], [dx3 dy3]);
% delta d1
v12 = [x2-x1;y2-y1]; v12 = v12/norm(v12);
v23 = [x2-x3;y2-y3]; v23 = v23/norm(v23);
deld1 = (s1.dx1 - s2.dx2)*v12(1) + (s1.dy1 - s2.dy2)*v12(2);
deld2 = (s3.dx3 - s2.dx2)*v23(1) + (s3.dy3 - s2.dy2)*v23(2);
slf = simplify(deld2/deld1);

sl = simplify(sl0*slf);
slfun = matlabFunction(sl);
sl0fun = matlabFunction(sl0);
slffun = matlabFunction(slf);

pretty(simplify((sl)));

% a = solve(subs(sl,y4,.4)==0)

% Sweep
xl = linspace(-.75,.75,1000);
yl = linspace(-1.5,0,1000);
slM = zeros(length(xl),length(yl));
sl0M = slM;
slfM = slM;
for i = 1:length(xl)
    for j = 1:length(yl)
        slM(i,j) = slfun(xl(i),yl(j));
        sl0M(i,j) = sl0fun(xl(i),yl(j));
        slfM(i,j) = slffun(xl(i),yl(j));
    end
end

% figure(7); clf;
% subplot(1,3,1)
% imagesc(xl,yl,tanh((slM))); colorbar; 
% hold on; 
% % plot(xl,1/4*ones(1,length(xl)),'r-','linewidth',2); 
% % plot(xl,xl./(1-xl),'r-','linewidth',2); 
% % plot(xl,-xl./(1+xl),'r-','linewidth',2); 
% % plot(xl,-1/4*ones(1,length(xl)),'g-','linewidth',2); 
% % plot(xl,xl./(1+xl),'g-','linewidth',2); 
% % plot(xl,-xl./(1-xl),'g-','linewidth',2); 
% hold off;
% subplot(1,3,2)
% imagesc(xl,yl,tanh((sl0M))); colorbar; 
% subplot(1,3,3);
% imagesc(xl,yl,tanh((slfM))); colorbar; 


%% Figure dimensions
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
    
% Fontsize
FS = 10;
% Distance visualization parameters
nW = .08;
lw = .5;
gr = 0.8;
% Scaling parameters
scc = .25;
    
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

% Default Name-Value Pairs
NVTitle = {'Units','centimeters','fontsize',FS};
NVTextH = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','right'};
NVTextR = {'Units','Normalized','fontsize',FS};

% Coefficients
nT = 1000;

% Colormap: Distance
DLin = linspace(.7,1.6,nT);
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));

% Colormap: Slope
SLin = linspace(0,1,nT);
C3a = [197 066 085]/255;
C3b = [031 172 204]/255;
CP3 = interp1([0 1],[C3a;C3b],linspace(0,1,nT));
o = [1 1 1];


%% Network Contsruction
% Design network unit
L = 1.0;
l = .5;
h = .5;
Xs = [-L  l L;...
       0  h 0];
Xf = [-L -l L;...
       0  h 0];
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
% Da = Da([1 3],:);



%% a: Design
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;

% Parameters
sc = .175;
sh = [.45;.2];
sh1 = [.15;.53];
sh2 = [.72;.53];
shl = [-.06;.84];
mSs = 5;
mSf = 3;
lSh = .008;
R = [-1 1; -1 1]*1;

% Plot
hold on;
% Initial network
line_coordinates(Xs(:,1:2)*sc+sh1,'lSh',lSh,'color',interp1(DLin,CP1,ds(1)));
line_coordinates(Xs(:,2:3)*sc+sh1,'lSh',lSh,'color',interp1(DLin,CP1,ds(2)));
visualize_network(Xs*sc+sh1,[],[1 1],'msize',mSs);
% Final network
line_coordinates(Xf(:,1:2)*sc+sh2,'lSh',lSh,'color',interp1(DLin,CP1,df(1)));
line_coordinates(Xf(:,2:3)*sc+sh2,'lSh',lSh,'color',interp1(DLin,CP1,df(2)));
visualize_network(Xf*sc+sh2,[],[1 1],'scolor',o,'msize',mSf);
% Solution space
visualize_conic_finite(Xs*sc+sh,Xf*sc+sh,R*sc+sh,'overlay',.99,'ucolori',o*gr);
visualize_network(Xs*sc+sh,[],[1 1],'msize',mSs);
visualize_network(Xf*sc+sh,[],[1 1],'scolor',o,'msize',mSf);
visualize_network([Xu1 Xu2]*sc+sh,[],[1 1],'scolor',interp1(SLin,CP3,0),'msize',mSs);
% Legend
visualize_network([0;0]+shl,[],[1 1],'msize',mSs);
visualize_network([0.31;0]+shl,[],[1 1],'scolor',o,'msize',mSf);
line([0 .1/sc]*sc+shl(1)+.56, [0 0]*sc+shl(2),'linewidth',1,'color',o*gr);

% Spline
sCo11 = [0.00 0.4 0.50 1.00;...
         0.00 0.17 0.83 1.00] .* [.12;.107] + [.3;.05];
sCo12 = [0.00 0.30 0.60 1.00;...
         0.00 0.20 0.80 1.00] .* [.12;-.021] + [.3;.05];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),'color',o*gr^.4);
plot_spline(sCo12,'head',1,'headpos',1,'ratio',sRat(pInd),'color',o*gr^.4);
hold off;

% Text
text(labX,subp(pInd,4)+labY,'\textbf{a} \hspace{1mm}limit cycle',NVTitle{:});
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


%% b: Map
pInd = 2;
subplot('position',subpN(pInd,:)); cla;

% Axes
sc = .25;
sh1 = [.25;.2];
sh2 = [.15;.3];
sh3 = [.25;.25];
ax = [.5 2.3 .5 2.3];
axL = ax + [.05 -.3 .05 -.3];
hold on;
arrow(axL(1:2),[1 1]*axL(3),sRat(pInd),'color','k');
arrow([1 1]*axL(1),axL(3:4),sRat(pInd),'color','k');
line([1;1].*df(1), [0;.04].*[1 1]+axL(3),'linewidth',.5,'color',interp1(DLin,CP1,df(1)));
line([1;1].*df(2), [0;.04].*[1 1]+axL(3),'linewidth',.5,'color',interp1(DLin,CP1,df(2)));
line([0;.04].*[1 1]+axL(1), [1;1].*df,'linewidth',.5,'color','k');
text(.86,.035,'$l_1$',NVTextR{:});
text(-.07,.5,'$l_2$',NVTextR{:});
text(.68,-.04,'$\sqrt{10}/2$',NVTextRA{:},'color',interp1(DLin,CP1,ds(1)));
text(.2,-.04,'$\sqrt{2}/2$',NVTextRA{:},'color',interp1(DLin,CP1,ds(2)));

% Plot
lSh = .02;
% Map
plot(axL(1:2),axL(3:4),'color',o*.8,'linewidth',.5);
plot(D(1,:),D(2,:),'color','k','linewidth',1);
% Points
[~,fpN] = min(abs(D(1,:)-D(2,:)));
dfp = D(1,fpN);
plot(ds,df,'ks','linewidth',3,'markersize',3);
line_coordinates(XCa(:,1:2,1)*sc+ds'+sh1,'lSh',lSh,'color',interp1(DLin,CP1,ds(1)));
line_coordinates(XCa(:,2:3,1)*sc+ds'+sh1,'lSh',lSh,'color',interp1(DLin,CP1,ds(2)));
visualize_network(XCa(:,1:3,1)*sc+ds'+sh1,XCa(:,4:5,1)*sc+ds'+sh1,conn,...
                  'ucolor',interp1(SLin,CP3,tanh(0)));
line_coordinates(XCa(:,1:2,end)*sc+df'+sh2,'lSh',lSh,'color',interp1(DLin,CP1,df(1)));
line_coordinates(XCa(:,2:3,end)*sc+df'+sh2,'lSh',lSh,'color',interp1(DLin,CP1,df(2)));
visualize_network(XCa(:,1:3,end)*sc+df'+sh2,XCa(:,4:5,end)*sc+df'+sh2,conn,...
                  'ucolor',interp1(SLin,CP3,tanh(0)));
plot(dfp,dfp,'ks','linewidth',3,'markersize',3);
line_coordinates(XCa(:,1:2,fpN)*sc+dfp+sh3,'lSh',lSh,'color',interp1(DLin,CP1,dfp));
line_coordinates(XCa(:,2:3,fpN)*sc+dfp+sh3,'lSh',lSh,'color',interp1(DLin,CP1,dfp));
visualize_network(XCa(:,1:3,fpN)*sc+dfp+sh3,XCa(:,4:5,fpN)*sc+dfp+sh3,conn,...
                  'ucolor',interp1(SLin,CP3,tanh(0)));
hold off;

% Text
text(labX,subp(pInd,4)+labY-.2,'\textbf{b} \hspace{1mm}conformational motion',NVTitle{:});
text(.39,.1,'start',NVTextR{:});
text(.06,.52,'end',NVTextR{:});

axis(ax);
set(gca,'visible',0);


%% c: Periodicity
pInd = 3;
subplot('position',subpN(pInd,:)); cla;

lSh = .03;
sc = .2;
shL = [.45;0];
sh = [.1;.4]-shL;
dotX = linspace(0,shL(1),10);
hold on;
for i = 1:4
    if(i<4)
        scatter(Xsa(1,1+i)*sc+sh(1)+shL(1)*i+dotX,...
                Xsa(2,1+i)*sc+sh(2)+shL(2)*i+zeros(1,10),.2,o*.9,'clipping',0);
        scatter(Xsa(1,2+i)*sc+sh(1)+shL(1)*i+dotX,...
                Xsa(2,2+i)*sc+sh(2)+shL(2)*i+zeros(1,10),.2,o*.9,'clipping',0);
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
    visualize_network(Xsa(:,(0:2)+i)*sc+sh+shL*i,Xua(:,(-1:0)+2*i)*sc+sh+shL*i,conn,...
                     'ucolor',interp1(SLin,CP3,0));
end

% Box
bD = [-.3 .95 .1 .79];
line(bD([1 1 2 2 1]), bD([3 4 4 3 3]), 'color', o*.8, 'linewidth',.5,'clipping',0);
line(bD([1 1 2 2 1])+1.32, bD([3 4 4 3 3]), 'color', o*.8, 'linewidth',.5,'clipping',0);
hold off;

% Text
text(labX-.5,subp(pInd,4)+labY,'\textbf{c} \hspace{1mm}repeats every 2 units',NVTitle{:});
text(-.1,.65,'$l_1$',NVTextR{:},'color',interp1(DLin,CP1,df(1)));
text(.2,.25,'$l_2$',NVTextR{:},'color',interp1(DLin,CP1,ds(1)));
text(.46,.65,'$l_3$',NVTextR{:},'color',interp1(DLin,CP1,df(1)));
text(.76,.25,'$l_4$',NVTextR{:},'color',interp1(DLin,CP1,ds(1)));

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% d: Stability
pInd = 4;
subplot('position',subpN(pInd,:)); cla;
sc = .2;
sh = [.55;.45];
Rz = rotz(3); Rz = Rz(1:2,1:2);
visualize_network(XCca(:,1:4,1)*sc+sh,XCca(:,5:8,1)*sc+sh,connc,...
                  'nalpha',.3,'lalpha',.3,'ucolor',interp1(SLin,CP3,0));
visualize_network(Rz*XCca(:,1:4,end)*sc+sh - [0;.01],XCca(:,5:8,end)*sc+sh,connc,...
                  'ucolor',interp1(SLin,CP3,0));
line_coordinates(XCca(:,1:2,1)*sc+sh,'lSh',.1,'nw',.03,'lalpha',.3,'style','-','lw',.5);
line_coordinates(XCca(:,3:4,1)*sc+sh,'lSh',-.1,'nw',.03,'lalpha',.3,'style','-','lw',.5);
line_coordinates(Rz*XCca(:,1:2,end)*sc+sh - [0;.01],'lSh',.15,'nw',.03,'style','-','lw',.5);
line_coordinates(Rz*XCca(:,3:4,end)*sc+sh - [0;.01],'lSh',-.15,'nw',.03,'style','-','lw',.5);


% Text
text(labX+.2,subp(pInd,4)+labY,'\textbf{d} \hspace{2.8mm}stability',NVTitle{:});
text(.5,.05,'$|\delta l_3 / \delta l_1|$',NVTextH{:});
text(.1,.67,'$\delta l_1$',NVTextH{:});
text(.9,.28,'$\delta l_3$',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% e: Phase
pInd = 5;
subplot('position',subpN(pInd,:)); cla;
sc = .075;

% Contour

colormap(CP3);
imagesc(yl,xl,tanh(abs(slM)),'alphadata',.6);
hold on; 
C1 = contour(yl,xl,abs(slM),[1 1],'color',interp1(SLin,CP3,tanh(1)),'linewidth',2);
plot(yl,1/4*ones(1,length(yl)),':','linewidth',1,'color',CP3(1,:)); 
plot(yl,yl./(1-yl),':','linewidth',1,'color',CP3(1,:)); 
plot(yl,-yl./(1+yl),':','linewidth',1,'color',CP3(1,:)); 

% Networks: Line
nP2 = [1 150 300]+333;
nSh = [-0.12 -0.12 -0.12;...
       -0.14 -0.14 -0.14];
nP = nP2;
for i = 1:length(nP)
    yV = 1/4;
    sh = [yl(nP(i));yV] + nSh(:,i);
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    plot(yl(nP(i)),yV,'k*','linewidth',.7,'markersize',5);
    visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),'msize',3,...
        'nalpha',0.3,'lalpha',0.3); 
    visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),'msize',3,...
        'nalpha',1.0,'lalpha',1.0); 
end
% Networks: Curve
nP1 = [700 750 800]+29;
nSh = [ 0.13  0.13  0.13;...
       -0.05 -0.05 -0.05];
nP = nP1;
for i = 1:length(nP)
    yV = -yl(nP(i))./(1+yl(nP(i)));
    sh = [yl(nP(i));yV] + nSh(:,i);
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    plot(yl(nP(i)),yV,'k+','linewidth',.7,'markersize',5);
    visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),'msize',3,...
        'nalpha',0.3,'lalpha',0.3); 
    visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),'msize',3,...
        'nalpha',1.0,'lalpha',1.0); 
end
% Networks: Curve
nP3 = [1 150 300]+333;
nSh = [-0.12 -0.12 -0.12;...
        0.10  0.07  0.03];
nP = nP3;
for i = 1:length(nP)
    yV = yl(nP(i))./(1-yl(nP(i)));
    sh = [yl(nP(i));yV] + nSh(:,i);
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    plot(yl(nP(i)),yV,'kx','linewidth',.7,'markersize',5);
    visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),'msize',3,...
        'nalpha',0.3,'lalpha',0.3); 
    visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),'msize',3,...
        'nalpha',1.0,'lalpha',1.0); 
end


% Networks 1: Line
nP5 = [1 150 300]+333;
nSh = [-0.11 -0.11 -0.11 -0.12 -0.12 -0.12;...
       -0.15 -0.15 -0.15 -0.20 -0.20 -0.20];
nP = nP5;
for i = 1:length(nP)
    yV = 0;
    sh = [yl(nP(i));yV] + nSh(:,i);
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    plot(yl(nP(i)),yV,'ko','linewidth',.7,'markersize',5);
    visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,'ucolor',interp1(SLin,CP3,tanh(1)),...
        'msize',3,'nalpha',1.0,'lalpha',1.0); 
end
% Networks 1: Curve
nP4 = [600 700 800]+40;
nSh = [ 0.10  0.10  0.10;...
       -0.08 -0.08 -0.08];
nP = nP4;
for i = 1:length(nP)
    yV = yl(nP(i))./(sqrt(2-2*yl(nP(i))^2));
    sh = [yl(nP(i));yV] + nSh(:,i);
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    plot(yl(nP(i)),yV,'ks','linewidth',.7,'markersize',5);
    visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,'ucolor',interp1(SLin,CP3,tanh(1)),...
        'msize',3,'nalpha',1.0,'lalpha',1.0); 
end
% Networks 1: Curve
nP6 = [1 150 300]+333;
nSh = [-0.12 -0.12 -0.12;...
        0.12  0.11  0.1];
nP = nP6;
for i = 1:length(nP)
    yV = -yl(nP(i))./sqrt(7*yl(nP(i))^2 + 1);
    sh = [yl(nP(i));yV] + nSh(:,i);
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    plot(yl(nP(i)),yV,'kh','linewidth',.7,'markersize',5);
    visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,'ucolor',interp1(SLin,CP3,tanh(1)),...
        'msize',3,'nalpha',1.0,'lalpha',1.0); 
end

% Boarder
nSc = 0.0085;
line([min(yl) min(yl) max(yl) max(yl) min(yl) min(yl)]+[-1 -1 1 1 -1 -1]*nSc,...
     [min(xl) max(xl) max(xl) min(xl) min(xl) max(xl)]+[-1 1 1 -1 -1  1]*nSc,...
     'color', 'k', 'linewidth',2);
hold off;

% Colorbar
hb = colorbar('location', 'east','position',[0.27 .05 .021 .71]);
delete(findall(gcf,'type','annotation'));
a = annotation('textbox',hb.Position,'FitBoxToText','off','FaceAlpha',0.4,...
    'EdgeColor',o*0,'BackgroundColor',o);
hb.Ticks = tanh(1); hb.TickLabels = {};
annotation('line',[.7165 .7355]-.4455, [1 1]*.591, 'Units','Normalize',...
           'color',interp1(SLin,CP3,tanh(1)),'linewidth',2);
annotation('line',[.7165 .7355]-.4455, [1 1]*.054, 'Units','Normalize',...
           'color',interp1(SLin,CP3,tanh(0)),'linewidth',2);

% Text
text(labX-.5,subp(pInd,4)+labY+.4,'\textbf{e} \hspace{1mm}phase',NVTitle{:});
text(.02,.31,'unstable',NVTextR{:});
text(.02,.69,'stable',NVTextR{:});
text(.5,-.03,'$y_1$',NVTextH{:});
text(0,-.035,'$-3/2$',NVTextH{:});
text(1,-.035,'$0$',NVTextH{:});
text(1.04,.5,'$x_2$',NVTextH{:});
text(1.04,0,'$\frac{-3}{4}$',NVTextH{:});
text(1.04,1,'$\frac{3}{4}$',NVTextH{:});
text(-.1,0,'$0$',NVTextH{:});
text(-.1,.75,'$1$',NVTextH{:});
text(-.105,.98,'$\infty$',NVTextH{:});
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


%% f,g: Examples
pInd = 6;
subplot('position',subpN(pInd,:)); cla;
sc = .09;
lSh = .0134;

tY = [.84 .51 .19];
nY = [.66 .36 .04];

hold on;
text(.05,tY(1),'symmetric',NVTextR{:});
plot(0,tY(1),'ko','linewidth',.7,'markersize',5);
nP = nP5;
for i = 1:length(nP)
    yV = 0;
    sh = [.15;nY(1)];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    if(i < length(nP))
        visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,[1 1],'ucolor',interp1(SLin,CP3,tanh(1)),...
            'msize',4,'nalpha',.1*i,'lalpha',.1*i); 
    else
        visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,'ucolor',interp1(SLin,CP3,tanh(1)),...
            'msize',4,'nalpha',i/length(nP),'lalpha',i/length(nP)); 
    end
end
text(.05,tY(2),'asymmetric 1',NVTextR{:});
plot(0,tY(2),'kh','linewidth',.7,'markersize',5);
nP = nP6;
for i = 1:length(nP)
    yV = -yl(nP(i))./sqrt(7*yl(nP(i))^2 + 1);
    sh = [.15;nY(2)-.03];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    if(i < length(nP))
        visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,[1 1],'ucolor',interp1(SLin,CP3,tanh(1)),...
            'msize',4,'nalpha',.1*i,'lalpha',.1*i); 
    else
        visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,'ucolor',interp1(SLin,CP3,tanh(1)),...
            'msize',4,'nalpha',i/length(nP),'lalpha',i/length(nP)); 
    end
end
text(.05,tY(3),'asymmetric 2',NVTextR{:});
plot(0,tY(3),'ks','linewidth',.7,'markersize',5);
nP = fliplr(nP4);
for i = 1:length(nP)
    yV = yl(nP(i))./(sqrt(2-2*yl(nP(i))^2));
    sh = [.15;nY(3)+.02];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    if(i < length(nP))
        visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,[1 1],'ucolor',interp1(SLin,CP3,tanh(1)),...
            'msize',4,'nalpha',.1*i,'lalpha',.1*i); 
    else
        visualize_network(Xsc*sc+sh,[Xupf Xup0]*sc+sh,connc,'ucolor',interp1(SLin,CP3,tanh(1)),...
            'msize',4,'nalpha',i/length(nP),'lalpha',i/length(nP)); 
    end
end

% Superstable
text(.67,tY(1),'1 co-linear',NVTextR{:});
plot(.5,tY(1),'k*','linewidth',.7,'markersize',5);
nP = fliplr(nP2);
for i = 1:length(nP)
    yV =  1/4;
    sh = [.62;nY(1)];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    if(i < length(nP))
        visualize_network(Xf*sc+sh,Xupf*sc+sh,[1 1],'ucolor',CP3(1,:),'msize',4,...
            'nalpha',.1*i,'lalpha',.1*i); 
    else
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),'msize',4,...
            'nalpha',.1,'lalpha',.1); 
        visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),'msize',4,...
            'nalpha',1,'lalpha',1); 
        XfM = mean(Xf(:,2:3),2); XfD = diff(Xf(:,2:3),1,2);
        line_coordinates(([-XfD XfD]*.8+XfM)*sc+sh,'lSh',lSh,...
                         'style','-','lw',.5,'color',CP3(1,:));
    end
end
text(.67,tY(2),'2 co-linear',NVTextR{:});
plot(.5,tY(2),'kx','linewidth',.7,'markersize',5);
nP = nP3;
for i = 1:length(nP)
    yV =  yl(nP(i))./(1-yl(nP(i)));
    sh = [.62;nY(2)];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    XfM = mean([Xf(:,1) Xupf(:,1)],2); XfD = diff([Xf(:,1) Xupf(:,1)],1,2);
    if(i < length(nP))
        visualize_network(Xf*sc+sh,Xupf*sc+sh,[1 1],'ucolor',CP3(1,:),'msize',4,...
            'nalpha',.1*i,'lalpha',.1*i); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',-lSh,...
                         'style','-','lw',.5,'color',CP3(1,:),'lalpha',.3);
    else
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),'msize',4,...
            'nalpha',.1,'lalpha',.1); 
        visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),'msize',4,...
            'nalpha',1,'lalpha',1); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',-lSh,...
                         'style','-','lw',.5,'color',CP3(1,:));
    end
end
text(.67,tY(3),'2 co-linear',NVTextR{:});
plot(.5,tY(3),'k+','linewidth',.7,'markersize',5);
nP = nP1;
for i = 1:length(nP)
    yV =  -yl(nP(i))./(1+yl(nP(i)));
    sh = [.62;nY(3)-.02];
    Xup0 = [[0;-yl(nP(i))] [1;l/h]*yV]+[L-l;h];
    Xupf = [[0;yl(nP(i))] [1; l/h]*yV];
    Xsp = [Xs(1,:); -Xs(2,:)+h] + [L-l;0];
    XfM = mean([Xsp(:,1) Xup0(:,2)],2); XfD = diff([Xsp(:,1) Xup0(:,2)],1,2);
    if(i < length(nP))
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,[1 1],'ucolor',CP3(1,:),'msize',4,...
            'nalpha',0.1*i,'lalpha',0.1*i); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',lSh,...
                         'style','-','lw',.5,'color',CP3(1,:),'lalpha',.3);
    else
        visualize_network(Xsp*sc+sh,Xup0*sc+sh,conn,'ucolor',CP3(1,:),'msize',4,...
            'nalpha',1,'lalpha',1); 
        visualize_network(Xf*sc+sh,Xupf*sc+sh,conn,'ucolor',CP3(1,:),'msize',4,...
            'nalpha',.1,'lalpha',.1); 
        line_coordinates(([-XfD XfD]*.9+XfM)*sc+sh,'lSh',lSh,...
                         'style','-','lw',.5,'color',CP3(1,:));
    end
end

% Text
text(labX,subp(pInd,4)+labY,'\textbf{f} \hspace{3mm}metastable',NVTitle{:});
text(.2,.915,'$|\delta l_3/\delta l_1| = 1$',NVTextH{:});
text(labX+3.1,subp(pInd,4)+labY,'\textbf{g} \hspace{3mm}superstable',NVTitle{:});
text(.78,.915,'$|\delta l_3/\delta l_1| = 0$',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'xtick',[],'ytick',[],'visible',0);


%% h: Stable limit cycle
pInd = 7;
subplot('position',subpN(pInd,:)); cla;
sc = .15;
sh = [.36;.6];
psc = .34; pSh = -.2;

% Network
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
% Distance plot
lCo = [ 0.15  0.15  2.41;...
        0.45 -0.04 -0.04];
line(lCo(1,:),lCo(2,:),'linewidth',0.5,'color',o*.7,'clipping',0);
scatter(XDa,Da(:,plI)'*psc+pSh,3,interp1(DLin,CP1,Da(:,plI)),...
        'marker','s','linewidth',2);
line([1;1].*XDa, [0;.04]+ones(1,size(Xsa,2)-1)*lCo(2,2),...
     'linewidth',.5,'color',o*.7,'clipping',0);
line([0;.04]+[1 1]*lCo(1,1),[1;1].*ds*psc+pSh,...
     'linewidth',.5,'color',o*.7,'clipping',0);

% Text
text(labX,subp(pInd,4)+labY+.15,'\textbf{h} \hspace{1mm}stable limit cycle',NVTitle{:});
text(.05,ds(1)*psc+pSh,'$\frac{\sqrt{10}}{2}$',NVTextRA{:},'color',interp1(DLin,CP1,ds(1)));
text(.05,ds(2)*psc+pSh,'$\frac{\sqrt{2}}{2}$',NVTextRA{:},'color',interp1(DLin,CP1,ds(2)));
text(.96,lCo(2,2),'$k$',NVTextR{:});
text(XDa(1)/sRat(pInd),-.13,'1',NVTextH{:});
text(XDa(2)/sRat(pInd),-.13,'2',NVTextH{:});
text(XDa(4)/sRat(pInd),-.13,'$\cdots$',NVTextH{:});
text(XDa(end)/sRat(pInd),-.13,num2str(length(XDa)),NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'xtick',[],'ytick',[],'visible',0);


%% Save
fName = 'figure6b';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');


