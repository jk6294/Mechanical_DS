% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Parameters and Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 8.75];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .4 .4 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 5.25 19.0 3.50];...
        [ 0.00 0.00 5.25 4.75];....
        [ 5.35 0.00 4.75 4.75];...
        [10.10 0.00 8.90 4.75]];
% Fontsize
FS = 10;
lSh = .18;
nW = .04;
lw = .5;
gr = 0.8;
    
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.16;
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
[Xu1,~] = construct_network(Xs10,Xs1T,Xu1,conn1,0,1);
Xu1 = Xu1(1:2,:);

% 2
Xs20 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs2T = [R2*l2 [0;0] R2\l2] + [0;1];
Xu2 = [-0.7  0.7;...
        1.5  1.5];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,~] = construct_network(Xs20,Xs2T,Xu2,conn2,0,1);
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
XCE = XCE - XCE(:,1,1);


%% Plot a
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
rot1 = rotz(-45); rot1 = rot1(1:2,1:2);
rot2 = rotz(90); rot2 = rot2(1:2,1:2);
rot3 = rotz(70); rot3 = rot3(1:2,1:2);
visualize_network(rot2*XCE(:,1:size(XsE,2),1)*.008+[0;0],[],[1,1],'scale',.7);
visualize_network(rot2*XCE(:,1:6,1)*.015+[0;.07],[],...
                  [1 2; 1 4; 3 2; 5 3; 3 6; 5 4; 5 6],'scale',1,'lcolor',[1 1 1]*gr);
visualize_network(rot3*XCE(:,1:size(XsE,2),50)*.008+[.35;0.07],[],[1,1],'scale',.7);
visualize_network(rot1*XCE(:,1:size(XsE,2),250)*.008+[1.0;.145],[],[1,1],'scale',.7);
line([-0.007 -0.007], [.02 0.06], 'color', [1 1 1]*gr,'linewidth',1);
line([ 0.040 0.073], [.02 0.06], 'color', [1 1 1]*gr,'linewidth',1);
line([-0.007 -0.007  0.040  0.040 -0.007], [-0.01 .02  .02  -.01 -.01], 'color', 'k','linewidth',1);
line([-0.007 -0.007  0.073  0.073 -0.007], [ 0.06 .105 .105  .06  .06], 'color', 'k','linewidth',1);

set(gca,'visible',0);
axis([0 sRat(pInd) 0 1]*.15 + [[1 1]*-.01 [1 1]*-0.042]);

% Text
textb = ['~~~~~~~~~~~~local edge constraints~~~~~~~~~~~~~~~',...
         'determine the global sequence and shape~~~~~~~~~~~~~~~~~',...
         'of changes in conformation'];
text(labX,subp(pInd,4)+labY,'\textbf{a}','Units','centimeters','fontsize',FS);
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS,'color',[100 100 255]/255);
text(.285,.1,'local constraints','Units','normalized','fontsize',FS,'color',[255 100 100]/255);
text(.58,.1,'global conformations','Units','normalized','fontsize',FS,'color',[255 100 100]/255);
text(.5,.1,'?','Units','normalized','fontsize',FS,'color',[255 100 100]/255);

ah = annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,'linewidth',1);
set(ah,'parent',gca,'position',[.45 -0.027 .1 0]);
ah = annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,'linewidth',1);
set(ah,'parent',gca,'position',[.245 0.06 .07 0]);
ah = annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,'linewidth',1);
set(ah,'parent',gca,'position',[.70 0.06 .07 0]);



%% b: Module + RBM
% Node Positions
L = 1;
Xs = [0    0   0   L;...
      L/2 -L/2 0   0];
conn = [1 3; 1 4; 2 3; 2 4];
t = 75:-1:15;
aX = 2*cosd(t)+0.3; aY = 2*sind(t)+0.2;
rot45 = rotz(15); rot45 = rot45(1:2,1:2);

% motion
nM = 380;
[XC, fC] = sim_motion(Xs,[],conn,.001,nM,-Xs,0);
disp(max(fC));

pI = [380 140 1];
CP = winter(floor(nM));
PCf = winter(2); PCf = PCf(1,:);

pInd = 2;
subplot('position',subpN(pInd,:)); cla;
% x-translation
visualize_network(Xs+[-1.75;1.3],[],conn,'lcolor',PCf,'lalpha',.2,'nalpha',.2);
visualize_network(Xs+[-2;1.3],[],conn,'lcolor',PCf);
% y-translation
visualize_network(Xs+[-.25;1.1],[],conn,'lcolor',PCf,'lalpha',.2,'nalpha',.2);
visualize_network(Xs+[-.25;1.3],[],conn,'lcolor',PCf);
% rotation
visualize_network(rot45*Xs+[1.50;1.3],[],conn,'lcolor',PCf,'lalpha',.2,'nalpha',.2);
visualize_network(Xs+[1.50;1.3],[],conn,'lcolor',PCf);
ax = [0 sRat(pInd) 0 1]*4 + [[1 1]*-2 [1 1]*-2];
axis(ax);
% Motion
for i = 1:1:length(pI)
    visualize_network(XC(:,:,pI(i))+[(i-1)*1.55-1.65;-1.15],[],conn,'lcolor',CP(pI(i),:));
end
line_coordinates(XC(:,1:2,pI(1))+[-1.65;-1.15],-lSh,nW,lw);
line_coordinates(XC(:,3:4,pI(1))+[-1.65;-1.15],-lSh,nW,lw);

% Text
textb = '\textbf{b}\hspace{0.8cm}3 rigid body motions';
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS);
textc = '\hspace{0.7cm}1 conformational motion';
text(labX,subp(pInd,4)/2.4+labY,textc,'Units','centimeters','fontsize',FS);
text(-.06,.21,'$d_1$','Units','normalized','fontsize',FS);
text(.18,.11,'$d_2$','Units','normalized','fontsize',FS);
set(annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,...
               'linewidth', 1, 'position',[-1.9 -0.4 diff(ax(1:2))-.1 0]),...
               'color', CP(1,:),'parent', gca);
set(annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,...
               'linewidth', 1, 'position',[ax(2)-.5 -0.4 -diff(ax(1:2))+.5 0]),...
               'color', CP(end,:), 'parent', gca);
scatter(linspace(ax(2)-.1,ax(1)+.1,size(CP,1)), ones(1,size(CP,1))*-.4,1,CP);
drawnow;
% 
% annotation('line',[.02 .98],[1 1]*subpN(4,4)/(subpN(4,4)+subpN(1,4))*1.0,'color',[1 1 1]*.9);


%% c: Map
L = 1;
sh = sqrt(2)*L/4;
Xs = [0    0   0   L;...
      L/2 -L/2 0   0];
conn = [1 3; 1 4; 2 3; 2 4];
nM = 380;
[XC, ~] = sim_motion(Xs,[],conn,.001,nM,-Xs,0);
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,3],:);

% Compute intermediate distance
nInt = 140;
ddiff = abs(d(1,:)-d(2,nInt));
nInt2 = find(ddiff == min(ddiff));
ddiff = abs(d(1,:)-d(2,nInt2));
nInt3 = find(ddiff == min(ddiff));

% pInd = 4;
pI = [1 nInt nInt2 380];
pISh = [ 0.130  .11  0.075  0.075  0.075;...
        -0.025 -.06 -0.067 -0.067 -0.067];
CP = winter(floor(nM));


pInd = 3;
subplot('position',subpN(pInd,:)); cla;
hold on;
plot([0 4],[0 4], '-', 'color', [1 1 1]*gr,'linewidth',.5);
scatter(d(1,:)',d(2,:)',1,CP(1:size(d,2),:),'o','linewidth',.1);
for i = 1:length(pI)
    plot(d(1,pI(i)),d(2,pI(i)),'s','markersize',3,'linewidth',3,'color',CP(pI(i),:));
    visualize_network((XC(:,:,pI(i))-[XC(1,4,pI(i));0])*0.08+d(:,pI(i))+pISh(:,i),[],conn,...
                      'scale',0.7,'lcolor',CP(pI(i),:));
end
set(gca,'visible',1,'xtick',[1/sqrt(2) 1],'ytick',[1/sqrt(2) 1],...
        'xticklabel',[],'yticklabel',[],'box',0);
axis([0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55]);

% Text
textc = '\textbf{c}~~plot $(d_1,d_2)$ along motion';
text(labX,subp(pInd,4)+labY,textc,'Units','centimeters','fontsize',FS);
text(subp(pInd,3)/2,-labY,'$d_1$','Units','centimeters','fontsize',FS);
text(labX,subp(pInd,4)/2,'$d_2$','Units','centimeters','fontsize',FS);
text(.68,.78,'$2$','Units','normalized','fontsize',FS,'rotation',45);
text(.17,.27,'$\sqrt{2}$','Units','normalized','fontsize',FS,'rotation',45);
text(.36,.46,'$d_2=d_1$','Units','normalized','fontsize',FS,...
                         'color',[1 1 1]*gr,'rotation',45);
text(.73,.55,'$d_2=\textbf{f}(d_1)$','Units','normalized','fontsize',FS);

drawnow;


%% d: 
pInd = 4;
sC = 1; fSc = 1;
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);
rot90 = rotz(90); rot90 = rot90(1:2,1:2);
xA1 = 0.15; xA2 = xA1 + 1.5;
yA1 = -0.6; yA2 = -2.4; yA3 = -2.8;
delyA = -.5;
xN1 = 0.83;
lA = 0.3;
XSh = 2.28;
% % Plot
subplot('position',subpN(pInd,:)); cla;
% Row 1
for i = 1:2
    XP = rot90^mod(i-1,2)*rot45*XC(:,:,pI(1))*sC+[xN1+(i-1)*XSh;-.7-sh*mod(i-1,2)];
    visualize_network(XP,[],conn,'scale',fSc,'lcolor',CP(pI(1),:));
    line_coordinates(XP(:,1:2),-lSh,nW,lw);
    line_coordinates(XP(:,3:4),(-1)^(i-1)*lSh*fSc,nW,lw);
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                    'color',[1 1 1]*gr,'linewidth',.5);
    set(ah,'parent',gca,'position',[xA1+XSh*(i-1) yA1+delyA*mod(i-1,2) lA 0.0]);
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                    'color',[1 1 1]*gr,'linewidth',.5);
    set(ah,'parent',gca,'position',[xA2+XSh*(i-1) yA1+delyA*mod(i,2) lA 0.0]);
end
XP = rot45*XC(:,:,pI(1))*sC+[xN1+2.38*XSh;-.7];
visualize_network(XP,[],conn,'scale',fSc,'lcolor',CP(pI(1),:));
line_coordinates(XP(:,1:2),-lSh,nW,lw);
line_coordinates(XP(:,3:4),lSh*fSc,nW,lw);
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[xA1+XSh*2.38 yA1+delyA*mod(2,2) lA 0.0]);
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[xA2+XSh*2.38 yA1+delyA*mod(3,2) lA 0.0]);


% Network 2
pInd4 = zeros(1,3); pInd4(1) = pI(2);
for i = 1:2
    XP = rot90^mod(i-1,2)*rot45*XC(:,:,pInd4(i))*sC+[xN1+(i-1)*XSh;-2.5-sh*mod(i-1,2)];
    visualize_network(XP,[],conn,'scale',fSc,'lcolor',CP(pInd4(i),:));
    line_coordinates(XP(:,1:2),-lSh,nW,lw);
    line_coordinates(XP(:,3:4),(-1)^(i-1)*lSh*fSc,nW,lw);
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                    'color',[1 1 1]*gr,'linewidth',.5);
    set(ah,'parent',gca,'position',[xA1+XSh*(i-1) yA2+delyA*mod(i-1,2) lA 0.0]);
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                    'color',[1 1 1]*gr,'linewidth',.5);
    set(ah,'parent',gca,'position',[xA2+XSh*(i-1) yA2+delyA*mod(i,2) lA 0.0]);
    
    % calculate next index
    pdiff = abs(d(1,:) - d(2,pInd4(i)));
    pInd4(i+1) = find(pdiff == min(pdiff));
end
XP = rot45*XC(:,:,pI(4))*sC+[xN1+2.42*XSh;-2.5];
visualize_network(XP,[],conn,'scale',fSc,'lcolor',CP(pI(4),:));
line_coordinates(XP(:,1:2),-lSh,nW,lw);
line_coordinates(XP(:,3:4),lSh*fSc,nW,lw);
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[xA1+XSh*2.41 yA2+delyA*mod(2,2) lA 0.0]);
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[xA2+XSh*2.39 yA2+delyA*mod(3,2) lA 0.0]);

% hold on;
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-.5 [1 1]*-3.4]);

% Text
textd = ['~~~~~\textbf{d}~~~~~~~~~each unit is a map $\textbf{f}$ from $d_k$'...
         ' to $d_{k+1}$'];
text(labX,subp(pInd,4)+labY,textd,'Units','centimeters','fontsize',FS);
text(.03,.9,['$d_1\rightarrow~~~ \textbf{f}~~~\rightarrow d_2'...
         '\rightarrow~~~ \textbf{f}~~~\rightarrow d_3~...~d_k'...
         '\rightarrow~~~ \textbf{f}~~~\rightarrow d_{k+1}$'],...
         'Units','normalized','fontsize',FS);
% Main schematic
delX = 0.28;
xT = 0.047;
text(xT+0*delX,.70,'$2$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+1*delX,.57,'$2$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+2*delX,.70,'$2$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+2.2*delX,.70,'...','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+2.38*delX,.70,'$2$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+3.38*delX,.57,'$2$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+0*delX,.24,'$1.8$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+1*delX,.12,'$1.6$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+2*delX,.24,'$1.5$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+2.2*delX,.24,'...','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+2.38*delX,.24,'$\sqrt{2}$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');
text(xT+3.38*delX,.12,'$\sqrt{2}$','Units','normalized','fontsize',FS,'HorizontalAlignment','center');



% %% e: Combine
% rot45 = rotz(-45); rot45 = rot45(1:2,1:2);
% Xsc = [0    0   0   L   L/2  L/2;...
%        L/2 -L/2 0   0   0    L];
% Xsc = rot45*Xsc;
% connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];
% [Xscc, conncc] = tesselate_network(Xsc,connc,[L/sqrt(2);0],[6;1]);
% 
% [XCc,fC] = sim_motion(Xscc,[],conncc,.005,920,Xscc,0);
% 
% %%
% pInd = 5;
% subplot('position',subpN(pInd,:)); cla;
% pSh = [6 4 2 0]; pI = [1 200 500 745];
% for pV = 1:length(pI)
%     % Closest color pairings for each module
%     disV = sqrt(sum(diff(XCc(:,:,pI(pV)),1,2).^2)); 
%     disV = disV(1:2:end);
%     disV = [disV(1:end-1); disV(2:end)];
%     C_LM = zeros(size(disV,2),3);
%     for i = 1:size(disV,2)
%         fdiff = sum(abs( d-disV(:,i)));
%         fdInd = find(fdiff == min(fdiff));
%         C_LM(i,:) = CP(fdInd,:);
%     end
%     visualize_network(XCc(:,:,pI(pV))+[-XCc(1,2,pI(pV));pSh(pV)],[],...
%                       conncc,'scale',.7,'lcolor',C_LM);
% end
% axis([0 sRat(pInd) 0 1]*7.5 + [[1 1]*-1,[1 1]*-0.8]);
% 
% % Text
% texte = '\textbf{e}~~motion of combined network';
% text(labX,subp(pInd,4)+labY,texte,'Units','centimeters','fontsize',FS);
% set(annotation('arrow','HeadLength',7,'HeadWidth',7,'color',[1 1 1]*gr,...
%                'linewidth', 1, 'position',[-.5 6.3 0 -7]), 'parent', gca);


%% Size and Save Figure
fName = 'figure1a';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');