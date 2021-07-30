% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(2); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 9.75];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .4 .4 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 5.00 4.75 4.75];...
        [ 4.75 5.00 4.75 4.75];...
        [ 9.50 5.00 4.75 4.75];...
        [14.25 5.00 4.74 4.75];....
        [ 0.00 0.00 4.75 4.75];...
        [ 4.75 0.00 4.75 4.75];...
        [ 9.50 0.00 4.75 4.75];...
        [14.25 0.00 4.74 4.75]];
% Fontsize
FS = 10;
lSh = .2;
nW = .04;
lw = .5;
gr = 0.7;

% Colors
CP = parula(4); CP = CP(1:3,:);
    
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.19;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');


%% a: Module: 2 fixed points
% Define module
L = 1;
Xs = [0    0   0   L;...
      L/2 -L/2 0   0];
conn = [1 3; 1 4; 2 3; 2 4];

% motion
nM = 440;
[XC, fC] = sim_motion(Xs,[],conn,.001,nM,-Xs,0);
disp(max(fC));

pI = [1 100 195 290 380];
CP = parula(floor(nM));
sC = 8/15;

pInd = 1;
subplot('position',subpN(pInd,:)); cla;
visualize_network(XC(:,:,pI(3))+[-1.5;-.2],[],conn,'lcolor',CP(pI(3),:));
line_coordinates(XC(:,1:2,pI(3))+[-1.5;-.2],-lSh,nW,lw);
line_coordinates(XC(:,3:4,pI(3))+[-1.5;-.2],-lSh,nW,lw);

% Motion
for i = 1:1:length(pI)
    visualize_network(XC(:,:,pI(i))*sC+[1.4;1.4-(i-1)*0.8],[],conn,'lcolor',CP(pI(i),:),'scale',.7);
    if(i == 1 || i == 5)
        line_coordinates(XC(:,1:2,pI(i))*sC+[1.4;1.4-(i-1)*0.8],-lSh*.7,nW*.7,lw);
        line_coordinates(XC(:,3:4,pI(i))*sC+[1.4;1.4-(i-1)*0.8],-lSh*.7,nW*.7,lw);
    end
end
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-2 [1 1]*-2]);

% Text
texta = ['\textbf{a}~~~~~~~~~fixed point: 1 conformation $i$ where $d_1 = d_2$'];
text(labX,subp(pInd,4)+labY,texta,'Units','centimeters','fontsize',FS);
text(-.03,.45,'$d_1$','Units','normalized','fontsize',FS);
text(.27,.34,'$d_2$','Units','normalized','fontsize',FS);
text(-.1,.85,'$(d_1,d_2)_{i=1}=(D_a^*,D_a^*)$','Units','normalized','fontsize',FS);
text(-.1,.05,'$(d_1,d_2)_{i=2}=(D_b^*,D_b^*)$','Units','normalized','fontsize',FS);


% a: Map: 2 fixed points
nM = 440;
% Calculate distances
[XC, fC] = sim_motion(Xs,[],conn,.001,nM,-Xs,0);
disp(max(fC));
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,3],:);
XC = XC * .08;

% Plot parameters
pF = [1 379];
pFSh = [ 0.14  0.085;...
        -0.04 -0.07];

% Begin figure
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
hold on;
% Map
plot([0 4],[0 4], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(d(1,:)',d(2,:)',1,CP(1:size(d,2),:),'o','linewidth',.1);
% Examples
for i = 1:length(pF)
    XP = (XC(:,:,pF(i))-[XC(1,4,pF(i));0])+d(:,pF(i))+pFSh(:,i);
    plot(d(1,pF(i)),d(2,pF(i)),'*','markersize',8,'linewidth',1.5,'color',CP(pF(i),:));
    visualize_network(XP,[],conn,'scale',0.7,'lcolor',CP(pF(i),:));
    line_coordinates(XP(:,[1,2]),-lSh*.7*.6/4,nW*.7*.6/4,lw);
    line_coordinates(XP(:,[3,4]),-lSh*.7*.6/4,nW*.7*.6/4,lw);
end

set(gca,'visible',1,'xtick',[],'ytick',[],'box',0);
a = [0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55];
axis(a);

% Text
% textb = '\textbf{b}~~fixed points: $D^* = f(D^*)$';
% text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(subp(pInd,3)/2,-labY,'$d_1$','Units','centimeters','fontsize',FS);
text(labX,subp(pInd,4)/2,'$d_2$','Units','centimeters','fontsize',FS,'rotation',0);
text(.35,.44,'$d_2=d_1$','Units','normalized','fontsize',FS,'color',[1 1 1]*gr,'rotation',45);
text(.54,.34,'$d_2=f(d_1)$','Units','normalized','fontsize',FS);
line([1 1], [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1],'linewidth',.7,'color','k');
line([1 1]/sqrt(2), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1]/sqrt(2), 'linewidth',.7,'color','k');
drawnow;


%% b: Module attachment and crystal structure for fixed point 1
pInd = 3;
sC = 8/15;
L = 1; sh = sqrt(2)*L/4*sC;
Xs1 = [0    0   0   L;...
       L/2 -L/2 0   0]*sC;
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);
Xs2 = rot45'*Xs1 + [0; -sh];
Xs1 = rot45*Xs1;
xSh = 15*sh; ySh = -5*sh;
conn = [1 3; 1 4; 2 3; 2 4];

C_L = parula(6);

% Plot
subplot('position',subpN(pInd,:)); cla;
hold on;
% Row 1
visualize_network(Xs2+[sh;0],[],conn,'scale',.7,'nalph',0.15,'lalph',.15,'lcolor',CP(pI(1),:));
plot([Xs1(1,3) Xs2(1,1)+xSh], [1 1]*Xs1(2,3), ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xs1(1,4) Xs2(1,2)+xSh], [1 1]*Xs1(2,4), ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xs1,[],conn,'scale',.7,'lcolor',CP(pI(1),:));
visualize_network(Xs2+[xSh;0],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
% line_coordinates(Xs1(:,3:4),lSh*7.5/4*.7,nW*7.5/4*.7,lw);
% line_coordinates(Xs2(:,1:2)+[xSh;0],-lSh*.7,nW*.7,lw);
% Row 2
visualize_network(Xs1+[2*sh;ySh],[],conn,'scale',.7,'nalph',0.15,'lalph',.15,'lcolor',CP(pI(1),:));
plot([Xs2(1,3)+sh Xs1(1,2)+xSh], [1 1]*Xs2(2,3)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xs2(1,4)+sh Xs1(1,1)+xSh], [1 1]*Xs2(2,4)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xs1+[0;ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
visualize_network(Xs2+[sh;ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
visualize_network(Xs1+[xSh;ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
% line_coordinates(Xs1(:,1:2)+[xSh;ySh],-lSh*.7,nW*.7,lw);
% Row 3
visualize_network(Xs2+[3*sh;2*ySh],[],conn,'scale',.7,'nalph',0.15,'lalph',.15,'lcolor',CP(pI(1),:));
plot([Xs1(1,3)+2*sh Xs2(1,1)+xSh], [1 1]*Xs1(2,3)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xs1(1,4)+2*sh Xs2(1,2)+xSh], [1 1]*Xs1(2,4)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xs1+[0;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
visualize_network(Xs2+[sh;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
visualize_network(Xs1+[2*sh;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
visualize_network(Xs2+[xSh;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
% line_coordinates(Xs2(:,1:2)+[xSh;2*ySh],-lSh*.7,nW*.7,lw);
% Row 4
for i = 1:6
    visualize_network(Xs1+[(2*i-2)*sh;3.3*ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
    visualize_network(Xs2+[(2*i-1)*sh;3.3*ySh],[],conn,'scale',.7,'lcolor',CP(pI(1),:));
end
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-.5 [1 1]*-3.55]);

% Text
textb = '\textbf{b}~~combine $d_1$ nodes of new module to $d_2$ nodes of network';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS);
text(.27,.24,'$D_a^*$ crystal','Units','normalized','fontsize',FS);


% b: Module attachment and crystal structure for fixed point 2
pInd = 4; sC = 8/15;
L = 1; sh = sqrt(2)*L/4*sC;
Xf1 = [-sh -sh 0 2*sh;...
        sh -sh 0 0];
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);
Xf2 = rot45'*Xf1;
Xf1 = rot45 * Xf1;
xSh = 15*sh; ySh = -5*sh;
conn = [1 3; 1 4; 2 3; 2 4];
sh = 0.5*sC;


% Plot
subplot('position',subpN(pInd,:)); cla;
hold on;
% Row 1
visualize_network(Xf2+[sh;0],[],conn,'scale',.7,'nalph',0.15,'lalph',.3,'lcolor',CP(pI(5),:));
plot([Xf1(1,3) Xf2(1,1)+xSh], [1 1]*Xf1(2,3), ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xf1(1,4) Xf2(1,2)+xSh], [1 1]*Xf1(2,4), ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xf1,[],conn,'scale',.7,'lcolor',CP(pI(5),:));
visualize_network(Xf2+[xSh;0],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
% Row 2
visualize_network(Xf1+[2*sh;ySh],[],conn,'scale',.7,'nalph',0.15,'lalph',.3,'lcolor',CP(pI(5),:));
plot([Xf2(1,3)+sh Xf1(1,2)+xSh], [1 1]*Xf2(2,3)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xf2(1,4)+sh Xf1(1,1)+xSh], [1 1]*Xf2(2,4)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xf1+[0;ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
visualize_network(Xf2+[sh;ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
visualize_network(Xf1+[xSh;ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
% Row 3
visualize_network(Xf2+[3*sh;2*ySh],[],conn,'scale',.7,'nalph',0.15,'lalph',.3,'lcolor',CP(pI(5),:));
plot([Xf1(1,3)+2*sh Xf2(1,1)+xSh], [1 1]*Xf1(2,3)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xf1(1,4)+2*sh Xf2(1,2)+xSh], [1 1]*Xf1(2,4)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xf1+[0;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
visualize_network(Xf2+[sh;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
visualize_network(Xf1+[2*sh;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
visualize_network(Xf2+[xSh;2*ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
% Row 4
for i = 1:6
    visualize_network(Xf1+[(2*i-2)*sh;3.4*ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
    visualize_network(Xf2+[(2*i-1)*sh;3.4*ySh],[],conn,'scale',.7,'lcolor',CP(pI(5),:));
end
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-.5 [1 1]*-3.55]);
text(.3,.24,'$D_b^*$ crystal','Units','normalized','fontsize',FS);


%% c: Module: limit cycle
% Define module
s = sqrt(3);
nxSh = s/2;
nySh = 0.3;
Xs = [-s/2 -nxSh  s/2;...
      -1/2  nySh -1/2];
Xf = [-s/2  nxSh  s/2;...
      -1/2  nySh -1/2];
Xus = [ 0.0000    0.3317
       -1.2691   -0.1410];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% motion
% Compute final position of limit cycle
[Q,W,v0,err] = construct_conic(Xs,Xf,1);
P = [W([1:2],:) v0(1:2);...
     zeros(1,2) 1]^-1;
Xuf = [W(3:4,:) v0(3:4)] * P * [Xus; ones(1,2)];

% Calculate distances
nM = 1700;
[XC, fC] = sim_motion(Xs,Xus,conn,.001,nM,[Xs Xus],0);
disp(max(fC));
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,2],:);
dsFP = sum(abs(d - flipud(d))); 
dsLC = sum(abs(d - flipud([sqrt(sum(diff(Xs,1,2).^2))'])));
indFP = find(dsFP == min(dsFP));
indLC = find(dsLC == min(dsLC));
XC = XC*.6;
rot3 = rotz(-15); rot3 = rot3(1:2,1:2);
for i = 1:size(XC,3)
    XC(:,:,i) = rot3*XC(:,:,i);
end
%
pF = [1 350 indFP 1200 indLC];
CP = parula(floor(nM));

%
sC = 8/15;
pInd = 5;
subplot('position',subpN(pInd,:)); cla;
visualize_network(XC(:,:,pF(3))+[-1.0;-.2],[],conn,'lcolor',CP(pF(3),:));
line_coordinates(XC(:,1:2,pF(3))+[-1.0;-.2],lSh,nW,lw);
line_coordinates(XC(:,2:3,pF(3))+[-1.0;-.2],lSh,nW,lw);


% Motion
for i = 1:1:length(pF)
    thetfp = atan2d(diff(XC(2,[1 3],pF(i))), diff(XC(1,[1 3],pF(i))));
    rotM = rotz(-thetfp); rotM = rotM(1:2,1:2);
    visualize_network(rotM*XC(:,:,pF(i))*sC+[1.6 + (i-1)*.03;1.6-(i-1)*0.8],[],conn,'lcolor',CP(pF(i),:),'scale',.7);
    if(i == 1 || i == 3 || i == 5)
        line_coordinates(rotM*XC(:,1:2,pF(i))*sC+[1.6 + (i-1)*.03;1.6-(i-1)*0.8],lSh*.7,nW*.7,lw);
        line_coordinates(rotM*XC(:,2:3,pF(i))*sC+[1.6 + (i-1)*.03;1.6-(i-1)*0.8],lSh*.7,nW*.7,lw);
    end
end
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-2 [1 1]*-2]);

% Text
textc = ['\textbf{c}~limit cycle: 2 conformations $j,k$ where $(d_1,d_2)_j = (d_2,d_1)_k$ '];
text(labX,subp(pInd,4)+labY,textc,'Units','centimeters','fontsize',FS);
text(0,.53,'$d_1$','Units','normalized','fontsize',FS);
text(.33,.53,'$d_2$','Units','normalized','fontsize',FS);
text(.44,.47,'$(D^*,D^*)$','Units','normalized','fontsize',FS);
text(-.1,.89,'$(d_1,d_2)_{j=1}=(D_a^\circ,D_b^\circ)$','Units','normalized','fontsize',FS);
text(-.1,.07,'$(d_1,d_2)_{k=2}=(D_b^\circ,D_a^\circ)$','Units','normalized','fontsize',FS);
annotation('line',[.02 .98],[1 1]*subpN(4,4)/(subpN(4,4)+subpN(1,4))*1.0,'color',[1 1 1]*.9);
annotation('line',[1 1]*subpN(4,3)/(subpN(4,3)+subpN(1,3))*0.993,[.02 .98],'color',[1 1 1]*.9);


sC = 8/15 *.6 * (1.8/4);
% Begin figure
pInd = 6;
subplot('position',subpN(pInd,:)); cla;
hold on;
% Map
scatter(d(1,:)',d(2,:)',1,CP(1:size(d,2),:),'o','linewidth',.1);
plot([0 4],[0 4], '-', 'color', [1 1 1]*gr,'linewidth',.5);
% Module 1
XP = [Xs Xus]*sC + [0.94;2.18];
visualize_network(XP,[],conn,'scale',0.7,'lcolor',CP(pF(1),:));
line_coordinates(XP(:,[1,2]),lSh*.4,nW*2/4,lw);
line_coordinates(XP(:,[2,3]),lSh*.4,nW*2/4,lw);
% Module 2
XP = XC(:,:,indFP)*sC/.6 + [1.65;1.66];
visualize_network(XP,[],conn,'scale',0.7,'lcolor',CP(pF(3),:));
line_coordinates(XP(:,[1,2]),lSh*.4,nW*2/4,lw);
line_coordinates(XP(:,[2,3]),lSh*.4,nW*2/4,lw);
% module 3
XP = [Xf Xuf]*sC + [2.05;1.02];
visualize_network(XP,[],conn,'scale',0.7,'lcolor',CP(pF(5),:));
line_coordinates(XP(:,[1,2]),lSh*.4,nW*2/4,lw);
line_coordinates(XP(:,[2,3]),lSh*.4,nW*2/4,lw);

plot(d(1,pF(1)),d(2,pF(1)),'ko','markersize',6,'linewidth',1.5,'color',CP(pF(1),:));
plot(d(1,pF(3)),d(2,pF(3)),'k*','markersize',8,'linewidth',1.5,'color',CP(pF(3),:));
plot(d(1,pF(5)),d(2,pF(5)),'ko','markersize',6,'linewidth',1.5,'color',CP(pF(5),:));

set(gca,'visible',1,'xtick',[],'ytick',[],'box',0);
a = [0 sRat(pInd) 0 1]*1.8 + [[1 1 1 1]*.6];
axis(a);

% Text
% texte = '\textbf{e}~limit cycle: $D_i^\mathrm{o} = f(f(D_i^\mathrm{o}))$';
% text(labX,subp(pInd,4)+labY,texte,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(subp(pInd,3)/2,-labY,'$d_1$','Units','centimeters','fontsize',FS);
text(labX,subp(pInd,4)/2,'$d_2$','Units','centimeters','fontsize',FS,'rotation',0);
text(.11,.20,'$d_2=d_1$','Units','normalized','fontsize',FS,...
                         'color',[1 1 1]*gr,'rotation',45);
text(.23,.18,'$d_2=f(d_1)$','Units','normalized','fontsize',FS);
line([1 1]*d(1,pF(1)), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1]*d(1,pF(1)),'linewidth',.7,'color','k');
line([1 1]*d(2,pF(3)), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1]*d(2,pF(3)), 'linewidth',.7,'color','k');
line([1 1]*d(1,pF(5)), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1]*d(1,pF(5)), 'linewidth',.7,'color','k');
drawnow;


%% d: Fixed point
% Define module
s = sqrt(3);
nxSh = s/2;
nySh = 0.3;
Xs1 = [-s/2 -nxSh  s/2;...
       -1/2  nySh -1/2]*.6;
Xf1 = [-s/2  nxSh  s/2;...
       -1/2  nySh -1/2]*.6;
Xus = [ 0.0000    0.3317
        -1.2691   -0.1410]*.6;
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Remove node offsets by setting first node to (0,0)
Xus = Xus - Xs1(:,1);
Xf1 = Xf1 - Xs1(:,1);
Xs1 = Xs1 - Xs1(:,1);

% Compute final position of limit cycle
[Q,W,v0,err] = construct_conic(Xs1,Xf1,1);
P = [W([1:2],:) v0(1:2);...
     zeros(1,2) 1]^-1;
Xuf = [W(3:4,:) v0(3:4)] * P * [Xus; ones(1,2)];

% Flip limit cycle in second module
Xf1 = [Xf1(1,:); -Xf1(2,:) + Xs1(2,2)];
Xuf = [Xuf(1,:); -Xuf(2,:) + Xs1(2,2)];

% Simulate to find fixed point
[XCa, fCa] = sim_motion(Xs1,Xus,conn,.0002,5000,[Xs1 Xus],0);
d = squeeze(sqrt(sum(diff(XCa,1,2).^2)));
ds = sum(abs(d([1,2],:) - flipud(d([1,2],:))));
fpInd = find(ds == min(ds));

% Orient fixed point
Xsfp1 = XCa(:,1:3,fpInd);
Xufp1 = XCa(:,4:5,fpInd);
thetfp = atan2d(diff(Xsfp1(2,[1 3])), diff(Xsfp1(1,[1 3])));
rotfp = rotz(-thetfp); rotfp = rotfp(1:2,1:2);
Xsfp1 = rotfp*Xsfp1;
Xufp1 = rotfp*Xufp1 - Xsfp1(:,1);
Xsfp1 = Xsfp1 - Xsfp1(:,1);
Xsfp2 = [Xsfp1(1,:); -Xsfp1(2,:) + Xsfp1(2,2)];
Xufp2 = [Xufp1(1,:); -Xufp1(2,:) + Xsfp1(2,2)];

sh = Xsfp1(1,2);
xSh = 10*sh; ySh = -3.3*sh;

pInd = 7;
subplot('position',subpN(pInd,:)); cla;
hold on;
% Row 1
visualize_network(Xsfp2+[sh;0],Xufp2+[sh;0],conn,'scale',.7,'nalph',.15,'lalph',.15,'lcolor',CP(pF(3),:));
plot([Xsfp1(1,2) Xsfp2(1,1)+xSh], [1 1]*Xsfp1(2,2), ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xsfp1(1,3) Xsfp2(1,2)+xSh], [1 1]*Xsfp1(2,3), ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xsfp1,Xufp1,conn,'scale',.7,'lcolor',CP(pF(3),:));
visualize_network(Xsfp2+[xSh;0],Xufp2+[xSh;0],conn,'scale',.7,'lcolor',CP(pF(3),:));
% Row 2
visualize_network(Xsfp1+[2*sh;ySh],Xufp1+[2*sh;ySh],conn,'scale',.7,'nalph',.15,'lalph',.15,'lcolor',CP(pF(3),:));
plot([Xsfp2(1,2)+sh Xsfp1(1,1)+xSh], [1 1]*Xsfp2(2,2)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xsfp2(1,3)+sh Xsfp1(1,2)+xSh], [1 1]*Xsfp2(2,3)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xsfp1+[0;ySh],Xufp1+[0;ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
visualize_network(Xsfp2+[sh;ySh],Xufp2+[sh;ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
visualize_network(Xsfp1+[xSh;ySh],Xufp1+[xSh;ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
% Row 3
visualize_network(Xsfp2+[3*sh;2*ySh],Xufp2+[3*sh;2*ySh],conn,'scale',.7,'nalph',.15,'lalph',.15,'lcolor',CP(pF(3),:));
plot([Xsfp1(1,2)+2*sh Xsfp2(1,1)+xSh], [1 1]*Xsfp1(2,2)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xsfp1(1,3)+2*sh Xsfp2(1,2)+xSh], [1 1]*Xsfp1(2,3)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xsfp1+[0;2*ySh],Xufp1+[0;2*ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
visualize_network(Xsfp2+[sh;2*ySh],Xufp2+[sh;2*ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
visualize_network(Xsfp1+[2*sh;2*ySh],Xufp1+[2*sh;2*ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
visualize_network(Xsfp2+[xSh;2*ySh],Xufp2+[xSh;2*ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
% Row 4
for i = 1:5
    visualize_network(Xsfp1+[(2*i-2)*sh;3.1*ySh],...
                      Xufp1+[(2*i-2)*sh;3.1*ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
    visualize_network(Xsfp2+[(2*i-1)*sh;3.1*ySh],...
                      Xufp2+[(2*i-1)*sh;3.1*ySh],conn,'scale',.7,'lcolor',CP(pF(3),:));
end
axis([0 sRat(pInd) 0 1]*7.5 + [[1 1]-1.3 [1 1]*-6.2]);

% Text
textd = ['\textbf{d}~~~~~~~limit cycles require adding alternating modules'];
text(labX,subp(pInd,4)+labY,textd,'Units','centimeters','fontsize',FS);
text(.3,.26,'$D^*$ crystal','Units','normalized','fontsize',FS);

pInd = 8;
% Define module
s = sqrt(3);
nxSh = s/2;
nySh = 0.3;
Xs1 = [-s/2 -nxSh  s/2;...
       -1/2  nySh -1/2]*.6;
Xf1 = [-s/2  nxSh  s/2;...
       -1/2  nySh -1/2]*.6;
Xus = [ 0.0000    0.3317
        -1.2691   -0.1410]*.6;
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Remove node offsets by setting first node to (0,0)
Xus = Xus - Xs1(:,1);
Xf1 = Xf1 - Xs1(:,1);
Xs1 = Xs1 - Xs1(:,1);

% Compute final position of limit cycle
[Q,W,v0,err] = construct_conic(Xs1,Xf1,1);
P = [W([1:2],:) v0(1:2);...
     zeros(1,2) 1]^-1;
Xuf = [W(3:4,:) v0(3:4)] * P * [Xus; ones(1,2)];

% Flip limit cycle in second module
sh = Xf1(1,2)/2;
Xf1 = [Xf1(1,:)-sh; -Xf1(2,:) + Xs1(2,2)];
Xuf = [Xuf(1,:)-sh; -Xuf(2,:) + Xs1(2,2)];

xSh = 11*sh; ySh = -3.5*sh;

Xsfp1 = Xs1; Xsfp2 = Xf1;
Xufp1 = Xus; Xufp2 = Xuf;


subplot('position',subpN(pInd,:)); cla;
hold on;
% Row 1
visualize_network(Xsfp2+[sh;0],Xufp2+[sh;0],conn,'scale',.7,'nalph',.3,'lalph',.3,'lcolor',CP(pF(5),:));
plot([Xsfp1(1,2) Xsfp2(1,1)+xSh], [1 1]*Xsfp1(2,2), ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xsfp1(1,3) Xsfp2(1,2)+xSh], [1 1]*Xsfp1(2,3), ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xsfp1,Xufp1,conn,'scale',.7,'lcolor',CP(pF(1),:));
visualize_network(Xsfp2+[xSh;0],Xufp2+[xSh;0],conn,'scale',.7,'lcolor',CP(pF(5),:));
% Row 2
visualize_network(Xsfp1+[2*sh;ySh],Xufp1+[2*sh;ySh],conn,'scale',.7,'nalph',.15,'lalph',.15,'lcolor',CP(pF(1),:));
plot([Xsfp2(1,2)+sh Xsfp1(1,1)+xSh], [1 1]*Xsfp2(2,2)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xsfp2(1,3)+sh Xsfp1(1,2)+xSh], [1 1]*Xsfp2(2,3)+ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xsfp1+[0;ySh],Xufp1+[0;ySh],conn,'scale',.7,'lcolor',CP(pF(1),:));
visualize_network(Xsfp2+[sh;ySh],Xufp2+[sh;ySh],conn,'scale',.7,'lcolor',CP(pF(5),:));
visualize_network(Xsfp1+[xSh;ySh],Xufp1+[xSh;ySh],conn,'scale',.7,'lcolor',CP(pF(1),:));
% Row 3
visualize_network(Xsfp2+[3*sh;2*ySh],Xufp2+[3*sh;2*ySh],conn,'scale',.7,'nalph',.3,'lalph',.3,'lcolor',CP(pF(5),:));
plot([Xsfp1(1,2)+2*sh Xsfp2(1,1)+xSh], [1 1]*Xsfp1(2,2)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
plot([Xsfp1(1,3)+2*sh Xsfp2(1,2)+xSh], [1 1]*Xsfp1(2,3)+2*ySh, ':','color',[1 1 1]*gr, 'linewidth',.7);
visualize_network(Xsfp1+[0;2*ySh],Xufp1+[0;2*ySh],conn,'scale',.7,'lcolor',CP(pF(1),:));
visualize_network(Xsfp2+[sh;2*ySh],Xufp2+[sh;2*ySh],conn,'scale',.7,'lcolor',CP(pF(5),:));
visualize_network(Xsfp1+[2*sh;2*ySh],Xufp1+[2*sh;2*ySh],conn,'scale',.7,'lcolor',CP(pF(1),:));
visualize_network(Xsfp2+[xSh;2*ySh],Xufp2+[xSh;2*ySh],conn,'scale',.7,'lcolor',CP(pF(5),:));
% Row 4
for i = 1:6
    visualize_network(Xsfp1+[(2*i-2)*sh;3.15*ySh],...
                      Xufp1+[(2*i-2)*sh;3.15*ySh],conn,'scale',.7,'lcolor',CP(pF(1),:));
    visualize_network(Xsfp2+[(2*i-1)*sh;3.15*ySh],...
                      Xufp2+[(2*i-1)*sh;3.15*ySh],conn,'scale',.7,'lcolor',CP(pF(5),:));
end
axis([0 sRat(pInd) 0 1]*7.5 + [[1 1]-1.5 [1 1]*-6.2]);
text(.32,.27,'$D_{ab}^\circ$ crystal','Units','normalized','fontsize',FS);


%% Size and Save Figure
fName = 'figure2a';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');