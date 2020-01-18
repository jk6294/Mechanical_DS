% Figure 1: Designing a Single Module
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 9.5];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .4 .4 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 0.00 4.75 9.50];...
        [ 4.75 4.75 4.75 4.75];...
        [ 4.75 0.00 4.74 4.75];....
        [ 9.50 4.75 4.75 4.75];...
        [ 9.50 0.00 4.75 4.75];...
        [14.25 4.75 4.75 4.75];...
        [14.25 0.00 4.74 4.75]];
% Fontsize
FS = 10;
lSh = .18;
nW = .04;
lw = .5;
gr = 0.7;
    
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



%% a: Modules determine map
s = sqrt(3);

% Node Colors for different curvatures
cTr1 = [126 200 255]/255;
cTr2 = [115 140 200]/255;
cTr3 = [015 082 186]/255;

l1 = [0;-s]; l2 = [0;-1.5*s];
a1 = 19.4; a2 = 38; a3 = 60;
R1 = rotz(-a1/2); R2 = rotz(-a2/2); R3 = rotz(-a3/2);
R1 = R1(1:2,1:2); R2 = R2(1:2,1:2); R3 = R3(1:2,1:2);

% 1
Xs10 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs1T = [R1*l2 [0;0] R1\l2] + [0;1];
Xu1 = [-1.5  1.5;...
        1.5  1.5];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu1,fV] = construct_network(Xs10,Xs1T,Xu1,conn1,0,1);
Xu1 = Xu1(1:2,:);

% 2
Xs20 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs2T = [R2*l2 [0;0] R2\l2] + [0;1];
Xu2 = [-0.7  0.7;...
        1.5  1.5];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,fV] = construct_network(Xs20,Xs2T,Xu2,conn2,0,1);
Xu2 = Xu2(1:2,:);

% 3
Xs30 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs3T = [R3*l2 [0;0] R3\l2] + [0;1];
Xu3 = [-1.0  1.0;...
        1.8  1.8];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu3,fV] = construct_network(Xs30,Xs3T,Xu3,conn3,0,1);
Xu3 = Xu3(1:2,:);


% Chain
XuC = cat(3,Xu1,Xu2,Xu3);
Xsp = [-s/2 0; -.5 1];
XuEL = [2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 3 2, 3 2, 3 2, 3 2 ,3 2, 3 2, 3 2,...
        2 2,...
        3 2, 3 2, 3 2, 3 2, 3 2,...
        2 2];
C_UNC = [cTr1; cTr2; cTr3];
[XsE,XuE,connE,C_UNE] = network_chain_x(Xsp,XuC,XuEL,C_UNC);
rot90 = [0 1; -1 0];
XsE = rot90*XsE;
XuE = rot90*XuE;

% Simulate
[XCE,fCE] = sim_motion(XsE,XuE,connE,.5,280,-[XsE,XuE],0);


%% Plot a
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
rot1 = rotz(-45); rot1 = rot1(1:2,1:2);
visualize_network(XCE(:,1:size(XsE,2),1)*.008+[-.017;0],[],[1,1],'scale',.7);
visualize_network(XCE(:,1:6,1)*.015+[.18;-.014],[],...
                  [1 2; 1 4; 3 2; 5 3; 3 6; 5 4; 5 6],'scale',1,'lcolor',[1 1 1]*gr);
visualize_network(XCE(:,1:size(XsE,2),80)*.008+[.054;-.13],[],[1,1],'scale',.7);
visualize_network(rot1*XCE(:,1:size(XsE,2),250)*.008+[.28;-.29],[],[1,1],'scale',.7);
line([.0 .15], [.015 .015], 'color', [1 1 1]*gr,'linewidth',1);
line([.0 .15], [-.033 -.08], 'color', [1 1 1]*gr,'linewidth',1);
line([-.03 -.03 .0 0 -.03], [-.033 .015 .015 -.033 -.033], 'color', 'k','linewidth',1);
line([.15 .15 .215 .215 .15], [.015 -.08 -.08 .015 .015], 'color', 'k','linewidth',1);

set(gca,'visible',0);
axis([0 sRat(pInd) 0 1]*.6 + [[1 1]*-.03 [1 1]*-0.58]);

% Text
textb = '\textbf{a}~~design sequence and shape';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(.15,.963,'local edge','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.15,.923,'constraints','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.18,.77,'determine the global','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.20,.73,'sequence and shape','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.54,.5,'of changes in','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.53,.46,'conformation','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.03,.12,'How do we choose local','Units','normalized','fontsize',FS);
text(.1,.08,'constraints to design','Units','normalized','fontsize',FS);
text(.02,.04,'conformational changes?','Units','normalized','fontsize',FS);


%% b: Module + RBM
% Node Positions
pInd = 2;
L = 1;
sh = sqrt(2)*L/4;
Xs = [0    0   0   L;...
      L/2 -L/2 0   0];
conn = [1 3; 1 4; 2 3; 2 4];
t = [75:-1:15]-0;
aX = 2*cosd(t)+0.3; aY = 2*sind(t)+0.2;
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);

subplot('position',subpN(pInd,:)); cla;
visualize_network(Xs+[0;2.5],[],conn);
visualize_network(Xs+[2.5;2.5],[],conn);
visualize_network(Xs,[],conn);
visualize_network(rot45*Xs+[2.5;.25],[],conn);
plot(aX,aY,'color',[1 1 1]*.7, 'linewidth', 2);
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-.5 [1 1]*-0.8]);

% Text
texta = '\textbf{b}~~~~~~rigid body motions';
text(labX,subp(pInd,4)+labY,texta,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(.38,.89,'$x$ shift','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.06,.68,'$y$ shift','Units','normalized','fontsize',FS,'color',[1 1 1]*gr,'rotation',-90);
text(.37,.65,'rotation','Units','normalized','fontsize',FS,'color',[1 1 1]*gr,'rotation',-45);
set(annotation('arrow','HeadLength',5,'HeadWidth',15,'color',[1 1 1]*gr,...
               'linewidth', 2, 'position',[1.3 2.5 1 0]), 'parent', gca);
set(annotation('arrow','HeadLength',5,'HeadWidth',15,'color',[1 1 1]*gr,...
               'linewidth', 2, 'position',[0 1.7 0 -1]), 'parent', gca);
set(annotation('arrow','HeadLength',5,'HeadWidth',15,'color',[1 1 1]*gr,...
               'linewidth', 2, 'position',[aX(end) aY(end) diff(aX(end-1:end)) diff(aY(end-1:end))]), 'parent', gca);
drawnow;

annotation('line',[0 0]+subpN(1,1)+subpN(1,3)+.015, [.02 .98],'color',[1 1 1]*.9);
annotation('line',[0 0]+subpN(3,1)+subpN(3,3)+.015, [.02 .98],'color',[1 1 1]*.9);
annotation('line',[0 0]+subpN(5,1)+subpN(5,3)+.015, [.02 .98],'color',[1 1 1]*.9);


%% c: Map
[XC, fC] = sim_motion(Xs,[],conn,.001,480,-Xs,0);
disp(max(fC));
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,3],:);

pInd = 3;
pI = [1 130 260 379 size(XC,3)];
pISh = [ 0.12  .13  .14  0.12  0.08;...
        -0.04 -.03 -.03 -.05 -0.067];

subplot('position',subpN(pInd,:)); cla;
plot(d(1,:),d(2,:),'k-','linewidth',1);
visualize_network(Xs*.6/4+[.65;1.03],[],conn);
line_coordinates(Xs(:,1:2)*.6/4+[.65;1.03],-lSh*.6/4,nW*.6/4,lw);
line_coordinates(Xs(:,3:4)*.6/4+[.65;1.03],-lSh*.6/4,nW*.6/4,lw);
plot([0 4],[0 4], '--', 'color', [1 1 1]*gr);
for i = 1:length(pI)
    plot(d(1,pI(i)),d(2,pI(i)),'ko','markersize',2.5,'linewidth',2.5);
    visualize_network((XC(:,:,pI(1))-[XC(1,4,pI(1));0])*0.08+d(:,pI(i))+pISh(:,i),[],conn,...
                      'scale',0.7,'nalpha',1-gr,'lalpha',1-gr);
    visualize_network((XC(:,:,pI(i))-[XC(1,4,pI(i));0])*0.08+d(:,pI(i))+pISh(:,i),[],conn,...
                      'scale',0.7);
end
set(gca,'visible',1,'xtick',[1/sqrt(2) 1],'ytick',[1/sqrt(2) 1],...
        'xticklabel',[],'yticklabel',[],'box',0);
axis([0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55]);

% Text
textb = '\textbf{c}~~~~conformational motion';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(subp(pInd,3)/1.07,-labY,'$d_1$','Units','centimeters','fontsize',FS);
text(labX+.18,subp(pInd,4)/1.1,'$d_2$','Units','centimeters','fontsize',FS,'rotation',90);
text(.03,.8,'$d_1$','Units','normalized','fontsize',FS);
text(.3,.69,'$d_2$','Units','normalized','fontsize',FS);
text(.62,.94,'$d_2=d_1$','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.62,.83,'$d_2=f(d_1)$','Units','normalized','fontsize',FS);
drawnow;


%% d: Combine
pInd = 4;
Xs1 = [0    0   0   L;...
       L/2 -L/2 0   0];
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);
Xs2 = rot45'*Xs1;
Xs1 = rot45*Xs1;
Xsc = [0    0   0   L   L/2  L/2;...
       L/2 -L/2 0   0   0    L];
Xsc = rot45*Xsc;
conn = [1 3; 1 4; 2 3; 2 4];
connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];
CP = parula(6);

subplot('position',subpN(pInd,:)); cla;
hold on;
plot([Xs1(1,3) Xs2(1,1)+7.4*sh], [1 1]*Xs1(2,3), '--',...
     'color',[1 1 1]*gr, 'linewidth',1);
plot([Xs1(1,4) Xs2(1,2)+7.4*sh], [1 1]*Xs1(2,4), '--',...
     'color',[1 1 1]*gr, 'linewidth',1);
visualize_network(Xs1,[],conn,'lcolor',CP(1,:));
visualize_network(Xs2+[7.4;-1]*sh,[],conn,'lcolor',CP(2,:));
visualize_network(Xsc+[1;-2.2],[],connc,'lcolor',CP(1:2,:));
% Module 1
line_coordinates(Xs1(:,1:2),-lSh,nW,lw);
line_coordinates(Xs1(:,3:4),-lSh,nW,lw);
% Module 2
line_coordinates(Xs2(:,1:2)+[7.4;-1]*sh,-lSh,nW,lw);
line_coordinates(Xs2(:,3:4)+[7.4;-1]*sh,-lSh,nW,lw);
% Combined Module
line_coordinates(Xsc(:,1:2)+[1;-2.2],-lSh,nW,lw);
line_coordinates(Xsc(:,3:4)+[1;-2.2],-lSh,nW,lw);
line_coordinates(Xsc(:,5:6)+[1;-2.2],-lSh,nW,lw);

% Text
textb = '\textbf{d}~~~~~~~combine modules';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
% Module 1
text(.03,.9,'$d_1$','Units','normalized','fontsize',FS);
text(.15,.6,'$d_2$','Units','normalized','fontsize',FS);
% Module 2
text(.73,.56,'$d_2''$','Units','normalized','fontsize',FS);
text(.94,.75,'$d_3$','Units','normalized','fontsize',FS);
% Combined Module
text(.28,.35,'$d_1$','Units','normalized','fontsize',FS);
text(.4,.04,'$d_2$','Units','normalized','fontsize',FS);
text(.61,.18,'$d_3$','Units','normalized','fontsize',FS);
text(.30,.74,'combine','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
text(.45,.67,'nodes','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
set(annotation('arrow','HeadLength',5,'HeadWidth',15,'color',[1 1 1]*gr,...
               'linewidth', 2, 'position',[1.5 -1 0 -.5]), 'parent', gca);
% text(.5,.5,'$\downarrow$','Units','normalized','fontsize',20);
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-.5 [1 1]*-3.2]);
drawnow;


%% e: Combined Map
[XCc, fCc] = sim_motion(Xsc,[],connc,.001,552,-Xsc,0);
disp(max(fCc));
dc = squeeze(sqrt(sum(diff(XCc,1,2).^2)));
dc = dc([1,3,5],:);
CP = parula(6);

pInd = 5;
pI = [1 280 size(XCc,3)];
pISh = [ 0.06  0.00 -0.02;...
         0.03 -0.06 -0.07];

subplot('position',subpN(pInd,:)); cla;
plot(d(1,:),d(2,:),'k-','linewidth',1);
hold on;
plot([0 4],[0 4], '--', 'color', [1 1 1]*gr);
for i = 1:length(pI)
    % Plot networks
    plot(dc(1,pI(i)),dc(2,pI(i)),'o','markersize',2,'linewidth',2,...
         'color',CP(i,:));
    visualize_network(XCc(:,:,pI(i))*0.08+dc(1:2,pI(i))+pISh(:,i),[],connc,...
                      'scale',0.7,'lcolor',CP(1:2,:));
    cobweb(dc(:,pI(i)),'ncolor',CP,'lcolor',[255 100 100]/255,...
           'mdecay',1,'arrowind',[1 2]);
end
set(gca,'visible',1,'xtick',[1/sqrt(2) 1],'ytick',[1/sqrt(2) 1],...
        'xticklabel',[],'yticklabel',[],'box',0);
axis([0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55]);

% Text
textb = '\textbf{e}~~~~conformational motion';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(subp(pInd,3)/1.07,-labY,'$d_k$','Units','centimeters','fontsize',FS);
text(labX+.18,subp(pInd,4)/1.2,'$d_{k+1}$','Units','centimeters','fontsize',FS,'rotation',90);
text(.72,.54,'$(d_1,d_2)$','Units','normalized','fontsize',FS,'color',CP(1,:));
text(.45,.28,'$(d_2,d_3)$','Units','normalized','fontsize',FS,'color',CP(2,:));
text(.24,.56,'$(d_2,d_2)$','Units','normalized','fontsize',FS,'color',[1 1 1]*gr);
drawnow;


%% f: Combine many
pInd = 6;
CP = parula(6);

subplot('position',subpN(pInd,:)); cla;
hold on;
plot([Xsc(1,5) Xs1(1,2)+2.5], [1 1]*Xsc(2,5), '--',...
     'color',[1 1 1]*gr, 'linewidth',1);
plot([Xsc(1,6) Xs1(1,1)+2.5], [1 1]*Xsc(2,6), '--',...
     'color',[1 1 1]*gr, 'linewidth',1);
plot([Xs2(1,3)+2*sh Xs2(1,1)+2.5], [1 1]*Xs1(2,3)-2.2, '--',...
     'color',[1 1 1]*gr, 'linewidth',1);
plot([Xs2(1,4)+2*sh Xs2(1,2)+2.5], [1 1]*Xs1(2,4)-2.2, '--',...
     'color',[1 1 1]*gr, 'linewidth',1);
visualize_network(Xsc+[ 0.0;-0.0],[],connc,'lcolor',CP(1:2,:));
visualize_network(Xs1+[ 2.5;-0.0],[],conn,'lcolor',CP(3,:));
visualize_network(Xsc+[ 0.0;-2.2],[],connc,'lcolor',CP(1:2,:));
visualize_network(Xs1+[ sh*2;-2.2],[],conn,'lcolor',CP(3,:));
visualize_network(Xs2+[ 2.5;-2.2-sh],[],conn,'lcolor',CP(4,:));
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-.5 [1 1]*-3.2]);
line_coordinates(Xsc(:,1:2),-lSh,nW,lw);
line_coordinates(Xsc(:,5:6),-lSh,nW,lw);
line_coordinates(Xs1(:,1:2)+[ 2.5;-0.0],-lSh,nW,lw);
line_coordinates(Xs1(:,3:4)+[ 2.5;-0.0],-lSh,nW,lw);
line_coordinates(Xsc(:,1:2)+[ 0.0;-2.2],-lSh,nW,lw);
line_coordinates(Xs1(:,3:4)+[ sh*2;-2.2],lSh,nW,lw);
line_coordinates(Xs2(:,1:2)+[ 2.5;-2.2-sh],-lSh,nW,lw);
line_coordinates(Xs2(:,3:4)+[ 2.5;-2.2-sh],-lSh,nW,lw);

textb = '\textbf{f}~~~combine many modules';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(.01,.89,'$d_1$','Units','normalized','fontsize',FS);
text(.38,.78,'$d_3$','Units','normalized','fontsize',FS);
text(.57,.82,'$d_3''$','Units','normalized','fontsize',FS);
text(.75,.6,'$d_4$','Units','normalized','fontsize',FS);
text(.01,.34,'$d_1$','Units','normalized','fontsize',FS);
text(.46,.19,'$d_4$','Units','normalized','fontsize',FS);
text(.59,.14,'$d_4''$','Units','normalized','fontsize',FS);
text(.9,.2,'$d_5$','Units','normalized','fontsize',FS);
drawnow;


%% g: Cobweb of many combined
[Xscc,conncc] = tesselate_network(Xsc,connc,[sh*2;0],[8;1]);
[XCc, fCc] = sim_motion(Xscc,[],conncc,.005,500,Xscc,0);
disp(max(fCc));
dc2 = squeeze(sqrt(sum(diff(XCc,1,2).^2)));
dc2 = dc2([1:2:size(dc2,1)-1],:);


%%
pInd = 7;
pI = [size(XCc,3)];
pISh = [-0.39 -0.07  0.00;...
         0.10 -0.19 -0.07];
CP = parula(18);
subplot('position',subpN(pInd,:)); cla;
plot(d(1,:),d(2,:),'k-','linewidth',1);
hold on;
plot([0 4],[0 4], '--', 'color', [1 1 1]*gr);
for i = 1:length(pI)
    % Plot networks
    plot(dc2(1,pI(i)),dc2(2,pI(i)),'o','markersize',2,'linewidth',2,...
         'color',CP(i,:));
    visualize_network(XCc(:,:,pI(i))*0.08+dc2(1:2,pI(i))+pISh(:,i),[],conncc,...
                      'scale',0.7,'lcolor',CP);
    cobweb(dc2(:,pI(i)),'ncolor',CP,'lcolor',[255 100 100]/255,...
           'mdecay',0.13,'msize',3,'arrowind',[8 11]);
end
set(gca,'visible',1,'xtick',[1/sqrt(2) 1],'ytick',[1/sqrt(2) 1],...
        'xticklabel',[],'yticklabel',[],'box',0);
axis([0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55]);

% Text
textb = '\textbf{g}~~~~conformational motion';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(subp(pInd,3)/1.07,-labY,'$d_k$','Units','centimeters','fontsize',FS);
text(labX+.18,subp(pInd,4)/1.2,'$d_{k+1}$','Units','centimeters','fontsize',FS,'rotation',90);
drawnow;


%% Size and Save Figure
fName = 'figure1a';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');