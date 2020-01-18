% Figure 1: Designing a Single Module
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(2); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 9.5];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .4 .4 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 4.75 4.75 4.75];...
        [ 4.75 4.75 4.75 4.75];...
        [ 9.50 4.75 4.75 4.75];...
        [14.25 4.75 4.74 4.75];....
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


%% a: Fixed points
% Define module
L = 1;
sh = sqrt(2)*L/4;
Xs = [0    0   0   L;...
      L/2 -L/2 0   0];
conn = [1 3; 1 4; 2 3; 2 4];

% Calculate distances
[XC, fC] = sim_motion(Xs,[],conn,.001,480,-Xs,0);
disp(max(fC));
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,3],:);
XC = XC * .08;

% Plot parameters
pI = [1 130 260 379 size(XC,3)];
pISh = [ 0.14  .14  .14  0.14  0.10;...
        -0.04 -.03 -.03 -.05 -0.067];
lAlph = [0 .7 .7 0 .7];
CC = cool(6); CC = CC([2 5],:);

% Begin figure
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
visualize_network(Xs*.6/4+[.65;1.05],[],conn);
line_coordinates(Xs(:,1:2)*.6/4+[.65;1.05],-lSh*.6/4,nW*.6/4,lw);
line_coordinates(Xs(:,3:4)*.6/4+[.65;1.05],-lSh*.6/4,nW*.6/4,lw);
plot([0 4],[0 4], '--', 'color', [1 1 1]*gr);
for i = 1:length(pI)
    XP = (XC(:,:,pI(i))-[XC(1,4,pI(i));0])+d(:,pI(i))+pISh(:,i);
    plot(d(1,pI(i)),d(2,pI(i)),'ko','markersize',2.5,'linewidth',2.5,...
         'color',[1 1 1]*lAlph(i));
    visualize_network(XP,[],conn,...
                      'nalph',1-lAlph(i),'lalph',1-lAlph(i),'scale',0.7);
end
plot(d(1,:),d(2,:),'k-','linewidth',1);
plot(d(1,pI(1)),d(2,pI(1)),'k*','markersize',8,'linewidth',1);
plot(d(1,pI(4)),d(2,pI(4)),'k*','markersize',8,'linewidth',1);
plot([1 1]*d(1,pI(1)), [0 1]*d(1,pI(1)), ':','color',CC(1,:));
plot([0 1]*d(1,pI(1)), [1 1]*d(1,pI(1)), ':','color',CC(1,:));
plot([1 1]*d(1,pI(4)), [0 1]*d(1,pI(4)), ':','color',CC(2,:));
plot([0 1]*d(1,pI(4)), [1 1]*d(1,pI(4)), ':','color',CC(2,:));

% Distance lines along fixed points
XP = (XC(:,:,pI(1))-[XC(1,4,pI(1));0])+d(:,pI(1))+pISh(:,1);
line_coordinates(XP(:,[1,2]),-lSh*.5/4,nW*.6/4,lw,'color',CC(1,:));
line_coordinates(XP(:,[3,4]),-lSh*.4/4,nW*.6/4,lw,'color',CC(1,:));
XP = (XC(:,:,pI(4))-[XC(1,4,pI(4));0])+d(:,pI(4))+pISh(:,4);
line_coordinates(XP(:,[1,2]),-lSh*.5/4,nW*.6/4,lw,'color',CC(2,:));
line_coordinates(XP(:,[3,4]),-lSh*.7/4,nW*.6/4,lw,'color',CC(2,:));

set(gca,'visible',1,'xtick',[],'ytick',[],'box',0);
a = [0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55];
axis(a);

% Text
textb = '\textbf{a}~~fixed points: $D_i^* = f(D_i^*)$';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(subp(pInd,3)/1.07,-labY,'$d_1$','Units','centimeters','fontsize',FS);
text(labX+.18,subp(pInd,4)/1.1,'$d_2$','Units','centimeters','fontsize',FS,'rotation',90);
text(subp(pInd,3)/1.4,-labY,'$D_1^*$','Units','centimeters',...
     'fontsize',FS,'color',CC(1,:));
text(subp(pInd,3)/4.5,-labY,'$D_2^*$','Units','centimeters',...
     'fontsize',FS,'color',CC(2,:));
text(labX+.18,subp(pInd,4)/1.4,'$D_1^*$','Units','centimeters',...
     'fontsize',FS,'rotation',90,'color',CC(1,:));
text(labX+.18,subp(pInd,4)/4.5,'$D_2^*$','Units','centimeters',...
     'fontsize',FS,'rotation',90,'color',CC(2,:));
text(.03,.8,'$d_1$','Units','normalized','fontsize',FS);
text(.3,.69,'$d_2$','Units','normalized','fontsize',FS);
line([1 1], [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color',CC(1,:));
line([0 diff(a(1:2))*.02]+a(1),[1 1],'linewidth',.7,'color',CC(1,:));
line([1 1]/sqrt(2), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color',CC(2,:));
line([0 diff(a(1:2))*.02]+a(1),[1 1]/sqrt(2), 'linewidth',.7,'color',CC(2,:));
drawnow;


%% b Module attachment and crystal structure
pInd = 2;
Xs1 = [0    0   0   L;...
       L/2 -L/2 0   0];
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);
Xs2 = rot45'*Xs1;
Xs1 = rot45*Xs1;
xSh = 7.5; nRep = 5;
Xsc = [0    0   0   L   L/2  L/2;...
       L/2 -L/2 0   0   0    L  ];
Xsc = rot45*Xsc;

Xf1 = [-sh -sh 0 2*sh;...
        sh -sh 0 0];
Xf2 = rot45'*Xf1;
Xf1 = rot45 * Xf1;
Xfc = [-sh -sh 0 2*sh sh sh;...
        sh -sh 0 0 sh 3*sh];
Xfc = rot45*Xfc;
ySh = 12.5;

conn = [1 3; 1 4; 2 3; 2 4];
connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];



[Xscc,conncc] = tesselate_network(Xsc,connc,[sh*2;0],[nRep;1]);
[Xfcc,connfcc] = tesselate_network(Xfc,connc,[sh*2*sqrt(2);0],[nRep;1]);
CP = lines(8); CP = CP([3:5],:);
CP = [CP; ones(2*nRep-3,3)*gr];


subplot('position',subpN(pInd,:)); cla;
hold on;
% Initial
% Combine
plot([Xs1(1,3) Xs2(1,1)+xSh*sh], [1 1]*Xs1(2,3), '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xs1(1,4) Xs2(1,2)+xSh*sh], [1 1]*Xs1(2,4), '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xs2(1,4) Xs1(1,1)+xSh*sh]+xSh*sh, [1 1]*Xs2(2,1), '--', 'color',[1 1 1]*gr, 'linewidth',1);
plot([Xs2(1,3) Xs1(1,2)+xSh*sh]+xSh*sh, [1 1]*Xs2(2,2), '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xs1(1,3) Xs2(1,1)+xSh*sh]+2*xSh*sh, [1 1]*Xs1(2,3), '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xs1(1,4) Xs2(1,2)+xSh*sh]+2*xSh*sh, [1 1]*Xs1(2,4), '--','color',[1 1 1]*gr, 'linewidth',1);
% Networks
visualize_network(Xs1,[],conn,'lcolor',CP(1,:),'scale',4/5);
visualize_network(Xs2+[xSh;-1]*sh,[],conn,'lcolor',CP(2,:),'scale',4/5);
visualize_network(Xs1+[2*xSh;0]*sh,[],conn,'lcolor',CP(3,:),'scale',4/5);
visualize_network(Xscc+[0;-1.8],[],conncc,'lcolor',CP(1:nRep*2,:),'scale',4/5);
% Distance
line_coordinates(Xs1(:,1:2),-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
line_coordinates(Xs1(:,3:4),-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
line_coordinates(Xs2(:,1:2)+[xSh;-1]*sh,-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
line_coordinates(Xs2(:,3:4)+[xSh;-1]*sh,-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
line_coordinates(Xs1(:,1:2)+[2*xSh;0]*sh,-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
line_coordinates(Xs1(:,3:4)+[2*xSh;0]*sh,-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
line_coordinates(Xscc(:,1:2)+[0;-1.8],-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
line_coordinates(Xscc(:,end-1:end)+[0;-1.8],-lSh*1.2,nW*1.5,lw,'color',CC(1,:));
% Final
% Combine
plot([Xf1(1,3) Xf2(1,1)+xSh*sh], [1 1]*Xf1(2,3)-ySh*sh, '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xf1(1,4) Xf2(1,2)+xSh*sh], [1 1]*Xf1(2,4)-ySh*sh, '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xf2(1,3) Xf1(1,2)+xSh*sh]+xSh*sh, [1 1]*Xf2(2,3)-ySh*sh, '--', 'color',[1 1 1]*gr, 'linewidth',1);
plot([Xf2(1,4) Xf1(1,1)+xSh*sh]+xSh*sh, [1 1]*Xf2(2,4)-ySh*sh, '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xf1(1,3) Xf2(1,1)+xSh*sh]+2*xSh*sh, [1 1]*Xf1(2,3)-ySh*sh, '--','color',[1 1 1]*gr, 'linewidth',1);
plot([Xf1(1,4) Xf2(1,2)+xSh*sh]+2*xSh*sh, [1 1]*Xf1(2,4)-ySh*sh, '--','color',[1 1 1]*gr, 'linewidth',1);
% Networks
visualize_network(Xf1+[0;-ySh]*sh,[],conn,'lcolor',CP(1,:),'scale',4/5);
visualize_network(Xf2+[xSh;-ySh]*sh,[],conn,'lcolor',CP(2,:),'scale',4/5);
visualize_network(Xf1+[2*xSh;-ySh]*sh,[],conn,'lcolor',CP(3,:),'scale',4/5);
visualize_network(Xfcc+[0;-6.2],[],connfcc,'lcolor',CP(1:nRep*2,:),'scale',4/5);
% Distance
line_coordinates(Xf1(:,1:2)+[0;-ySh]*sh,-lSh*1.3,nW*1.5,lw,'color',CC(2,:));
line_coordinates(Xf1(:,3:4)+[0;-ySh]*sh,-lSh*2,nW*1.5,lw,'color',CC(2,:));
line_coordinates(Xf2(:,1:2)+[xSh;-ySh]*sh,-lSh*1.3,nW*1.5,lw,'color',CC(2,:));
line_coordinates(Xf2(:,3:4)+[xSh;-ySh]*sh,-lSh*2,nW*1.5,lw,'color',CC(2,:));
line_coordinates(Xf1(:,1:2)+[2*xSh;-ySh]*sh,-lSh*1.3,nW*1.5,lw,'color',CC(2,:));
line_coordinates(Xf1(:,3:4)+[2*xSh;-ySh]*sh,-lSh*2,nW*1.5,lw,'color',CC(2,:));
line_coordinates(Xfcc(:,1:2)+[0;-6.2],-lSh*1.3,nW*1.5,lw,'color',CC(2,:));
line_coordinates(Xfcc(:,end-1:end)+[0;-6.2],-lSh*2,nW*1.5,lw,'color',CC(2,:));


axis([0 sRat(pInd) 0 1]*7.5 + [[1 1]*-.74 [1 1]*-6.7]);

% Text
textb = '\textbf{b}~~~~~~~combine modules';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(.8,.62,'$D_1^*$','Units','normalized','fontsize',FS,'color',CC(1,:));
text(.8,.04,'$D_2^*$','Units','normalized','fontsize',FS,'color',CC(2,:));


%% c: Limit Cycle
% Define module
L = 1;
s = sqrt(3);
nxSh = s/4;
nySh = 1/2;
Xs = [-s/2 -nxSh  s/2;...
      -1/2  nySh -1/2];
Xu = [-0.3861         0
      -0.6672   -0.7292];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Calculate distances
[XCa, fCa] = sim_motion(Xs,Xu,conn,.001,900,[Xs Xu],0);
[XCb, fCb] = sim_motion(Xs,Xu,conn,.001,700,-[Xs Xu],0);
XC = cat(3,flip(XCb(:,:,2:end),3),XCa);
fC = [fliplr(fCa(2:end)) fCa];
disp(max(fC));
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,2],:);
Xs = Xs/3; Xu = Xu/3;
XC = XC / 3 * .6;

% Plot parameters
pI = [1 380 size(XCb,3) 900 1041 size(XC,3)];
pISh = [-0.30 -.02  .20  .25  0.26  0.28  0.28;...
         0.15  .3   .25  .2   0.19  0.10  0.00];
lAlph = [1 1 0 1 1 1 1]*.7;
CC = cool(6); CC = CC([2 3 5],:);

% Begin figure
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
visualize_network(Xs+[.5;.5],Xu+[.5;.5],conn);
line_coordinates(Xs(:,1:2)+[.5;.5],lSh/3,nW/4,lw);
line_coordinates(Xs(:,2:3)+[.5;.5],lSh/3,nW/4,lw);
plot([0 4],[0 4], '--', 'color', [1 1 1]*gr);
for i = 1:length(pI)
    XP = (XC(:,:,pI(i))-[XC(1,4,pI(i));0])+d(:,pI(i))+pISh(:,i);
    plot(d(1,pI(i)),d(2,pI(i)),'ko','markersize',2.5,'linewidth',2.5,...
         'color',[1 1 1]*lAlph(i));
    visualize_network(XP,[],conn,...
                      'nalph',1-lAlph(i),'lalph',1-lAlph(i),'scale',0.7);
end
plot(d(1,:),d(2,:),'k-','linewidth',1);
% plot(d(2,:),d(1,:),'k-','linewidth',1);
plot(d(1,pI(2)),d(2,pI(2)),'ko','markersize',6,'linewidth',1.0);
plot(d(1,pI(4)),d(2,pI(4)),'k*','markersize',8,'linewidth',1.0);
plot(d(1,pI(6)),d(2,pI(6)),'ko','markersize',6,'linewidth',1.0);
plot([1 1]*d(1,pI(2)), [0 1]*d(2,pI(2)), ':','color',CC(1,:));
plot([0 1]*d(2,pI(2)), [1 1]*d(1,pI(2)), ':','color',CC(1,:));
plot([1 1]*d(1,pI(4)), [0 1]*d(2,pI(4)), ':','color',CC(2,:));
plot([0 1]*d(2,pI(4)), [1 1]*d(1,pI(4)), ':','color',CC(2,:));
plot([1 1]*d(1,pI(6)), [0 1]*d(2,pI(6)), ':','color',CC(3,:));
plot([0 1]*d(2,pI(6)), [1 1]*d(1,pI(6)), ':','color',CC(3,:));

XP = (XC(:,:,pI(2))-[XC(1,4,pI(2));0])+d(:,pI(2))+pISh(:,2);
line_coordinates(XP(:,[1,2]),lSh/4,nW*.6/4,lw,'color',CC(1,:));
line_coordinates(XP(:,[3,4]),-lSh/4,nW*.6/4,lw,'color',CC(1,:));
XP = (XC(:,:,pI(4))-[XC(1,4,pI(4));0])+d(:,pI(4))+pISh(:,4);
line_coordinates(XP(:,[1,2]),lSh/4,nW*.6/4,lw,'color',CC(2,:));
line_coordinates(XP(:,[3,4]),-lSh/4,nW*.6/4,lw,'color',CC(2,:));

set(gca,'visible',1,'xtick',[],'ytick',[],'box',0);
a = [0 sRat(pInd) 0 1]*2 + [[1 1]*.07 [1 1]*.0];
axis(a);

% Text
textb = '\textbf{a}~~fixed points: $D_i^* = f(D_i^*)$';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'fontweight','bold');
text(subp(pInd,3)/1.07,-labY,'$d_1$','Units','centimeters','fontsize',FS);
text(labX+.18,subp(pInd,4)/1.1,'$d_2$','Units','centimeters','fontsize',FS,'rotation',90);
% x-axis
text(subp(pInd,3)/3.5,-labY,'$D_1^o$','Units','centimeters',...
                      'fontsize',FS,'color',CC(1,:));
text(subp(pInd,3)/1.84,-labY,'$D^*$','Units','centimeters',...
                      'fontsize',FS,'color',CC(2,:));
text(subp(pInd,3)/1.47,-labY,'$D_2^o$','Units','centimeters',...
                      'fontsize',FS,'color',CC(3,:));
% y-axis
text(labX+.18,subp(pInd,4)/3.5+.15,'$D_1^o$','Units','centimeters',...
                      'fontsize',FS,'color',CC(1,:),'rotation',90);
text(labX+.18,subp(pInd,4)/1.84+.15,'$D^*$','Units','centimeters',...
                      'fontsize',FS,'color',CC(2,:),'rotation',90);
text(labX+.18,subp(pInd,4)/1.47+.15,'$D_2^o$','Units','centimeters',...
                      'fontsize',FS,'color',CC(3,:),'rotation',90);
text(.1,.36,'$d_1$','Units','normalized','fontsize',FS);
text(.3,.69,'$d_2$','Units','normalized','fontsize',FS);
line([1 1]*d(1,pI(2)), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color',CC(1,:));
line([0 diff(a(1:2))*.02]+a(1),[1 1]*d(1,pI(2)),'linewidth',.7,'color',CC(1,:));
line([1 1], [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color',CC(2,:));
line([0 diff(a(1:2))*.02]+a(1),[1 1],'linewidth',.7,'color',CC(2,:));
line([1 1]*d(2,pI(2)), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color',CC(3,:));
line([0 diff(a(1:2))*.02]+a(1),[1 1]*d(2,pI(2)), 'linewidth',.7,'color',CC(3,:));
drawnow;


%% Size and Save Figure
fName = 'figure2a';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');