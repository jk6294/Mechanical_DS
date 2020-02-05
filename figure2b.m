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
        [ 5.15 5.00 8.90 4.75];...
        [14.25 5.00 4.75 4.75];....
        [ 0.00 0.00 4.75 4.75];...
        [ 4.75 0.00 4.75 4.75];...
        [ 9.50 0.00 4.75 4.75];...
        [14.25 0.00 4.74 4.75]];
% Fontsize
FS = 10;
lSh = .2;
nW = .04;
lw = .5;
gr = 0.8;

% Colors
CP = parula(4); CP = CP(1:3,:);
    
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.17;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');

% Default Name-Value Pairs
NVTitle = {'Units','centimeters','fontsize',FS};
NVTextH = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center'};
NVTextR = {'Units','Normalized','fontsize',FS};


%% a: Module: 2 fixed points
% Define module
L = 1;
Xs = [0    0   0   L;...
      L/2 -L/2 0   0];
conn = [1 3; 1 4; 2 3; 2 4];

% Map: 2 fixed points
nM = 380;
% Calculate distances
[XC, fC] = sim_motion(Xs,[],conn,.001,nM,-Xs,0);
disp(max(fC));
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,3],:);
% Compute intermediate distance
nInt = 140;
ddiff = abs(d(1,:)-d(2,nInt));
nInt2 = find(ddiff == min(ddiff));

% Plot parameters
pF = [1 380];
pFSh = [ 0.130  0.075;...
        -0.025 -0.067];
pI = [nInt nInt2];
pISh = [  .11  0.075;...
         -.06 -0.077];
CP = winter(floor(nM));
sC = 8/15;

% Begin figure
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
hold on;
% Map
plot([0 4],[0 4], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(d(1,:)',d(2,:)',1,CP(1:size(d,2),:),'o','linewidth',.1);

% Fixed points
for i = 1:length(pF)
    XP = (XC(:,:,pF(i))-[XC(1,4,pF(i));0])*.08+d(:,pF(i))+pFSh(:,i);
    plot(d(1,pF(i)),d(2,pF(i)),'s','markersize',3,'linewidth',3,'color',CP(pF(i),:));
    visualize_network(XP,[],conn,'scale',0.7,'lcolor',CP(pF(i),:));
end
line_coordinates(XP(:,1:2),-lSh*.12,nW*.12,lw);
line_coordinates(XP(:,3:4),-lSh*.12,nW*.12,lw);
% Examples
for i = 1:length(pI)
    XP = (XC(:,:,pI(i))-[XC(1,4,pI(i));0])*.08+d(:,pI(i))+pISh(:,i);
    plot(d(1,pI(i)),d(2,pI(i)),'s','markersize',3,'linewidth',3,'color',CP(pI(i),:));
    visualize_network(XP,[],conn,'scale',0.7,'lcolor',CP(pI(i),:));
end
hold off;
set(gca,'visible',1,'xtick',[],'ytick',[],'box',0);
a = [0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55];
axis(a);

% Text
texta = '\textbf{a}~~motion of unit: $d_2 = \textbf{f}(d_1)$';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(subp(pInd,3)/2,-labY,'$d_1$',NVTitle{:});
text(labX,subp(pInd,4)/2,'$d_2$',NVTitle{:},'rotation',0);
text(.36,.46,'$d_2=d_1$',NVTextR{:},'color',[1 1 1]*gr,'rotation',45);
text(.60,.37,'$\textbf{f}$',NVTextR{:});
text(.69,.79,'$2$',NVTextR{:},'rotation',45);
text(.16,.26,'$\sqrt{2}$',NVTextR{:},'rotation',45);
text(.11,.06,'$d_1$',NVTextR{:});
text(.39,.06,'$d_2$',NVTextR{:});

line([1 1], [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1],'linewidth',.7,'color','k');
line([1 1]/sqrt(2), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1]/sqrt(2), 'linewidth',.7,'color','k');
drawnow;


%% b: Module is a map
pInd = 2;
sC = 1; fSc = 1; sh = sqrt(2)/4;
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);
rot90 = rotz(90); rot90 = rot90(1:2,1:2);
% Arrow parameters
xA1 = -0.18;
yA1 = -0.6; yA2 = -2.4; yA3 = -2.8;
delyA = -.5;
lA = 0.5;
% node parameters
xN1 = 0.8;
XSh = 2.4;
pI = [1 nInt nInt2 380];
pISh = [ 0.130  .11  0.075  0.075;...
        -0.025 -.06 -0.077 -0.067];
CP = winter(floor(nM));

% Plot
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
end
XP = rot45*XC(:,:,pI(1))*sC+[xN1+2.38*XSh;-.7];
visualize_network(XP,[],conn,'scale',fSc,'lcolor',CP(pI(1),:));
line_coordinates(XP(:,1:2),-lSh,nW,lw);
line_coordinates(XP(:,3:4),lSh*fSc,nW,lw);
ah = annotation('arrow','HeadLength',3,'HeadWidth',3,...
                'color',[1 1 1]*gr,'linewidth',.5);
set(ah,'parent',gca,'position',[xA1+XSh*2.38 yA1+delyA*mod(2,2) lA 0.0]);


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

% hold on;
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-.4 [1 1]*-3.4]);

% Text
textb = ['~\textbf{b}~~~~~~~~~~~each unit is a map $\textbf{f}$ from $d_k$'...
         ' to $d_{k+1}$'];
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
text(-.03,.9,['$d_1\longrightarrow~\textbf{f}(d_1)= d_2'...
         '\longrightarrow~\textbf{f}(d_2)= d_3~...~d_k'...
         '\longrightarrow\textbf{f}(d_k)=d_{k+1}$'],...
         NVTextR{:});
% Main schematic
delX = 0.297;
xT = -.016;
text(xT+0*delX,.70,'$2$',NVTextH{:});
text(xT+1*delX,.57,'$2$',NVTextH{:});
text(xT+2*delX,.70,'$2$',NVTextH{:});
text(xT+2.2*delX,.70,'...',NVTextH{:});
text(xT+2.38*delX,.70,'$2$',NVTextH{:});
text(xT+3.38*delX,.57,'$2$',NVTextH{:});
text(xT+0*delX,.24,'$1.8$',NVTextH{:});
text(xT+1*delX,.12,'$1.6$',NVTextH{:});
text(xT+2*delX,.24,'$1.5$',NVTextH{:});
text(xT+2.2*delX,.24,'...',NVTextH{:});
text(xT+2.38*delX,.24,'$\sqrt{2}$',NVTextH{:});
text(xT+3.38*delX,.12,'$\sqrt{2}$',NVTextH{:});


%% c: Combine networks
pInd = 3;

% Placement parameters
xSh = 2.4;
rot45 = rotz(-45); rot45 = rot45(1:2,1:2);

% Module positions
XP1 = rot45*XC(:,:,pI(1));
XP2 = rot45'*XC(:,:,pI(1)) + [xSh;-sh];
XP3 = rot45*XC(:,:,pI(2)); XP3 = XP3 + [0;-1.8-XP3(2,3)];
XP4 = rot45'*XC(:,:,pI(3)); XP4 = XP4 + [xSh;-1.8-XP4(2,1)];

% Plot
subplot('position',subpN(pInd,:)); cla; hold on;
% Lines 1
plot([XP1(1,3) XP2(1,1)],[XP1(2,3) XP2(2,1)],':','linewidth',1,'color',[1 1 1]*gr);
plot([XP1(1,4) XP2(1,2)],[XP1(2,4) XP2(2,2)],':','linewidth',1,'color',[1 1 1]*gr);
% Lines 2
plot([XP3(1,3) XP4(1,1)],[XP3(2,3) XP4(2,1)],':','linewidth',1,'color',[1 1 1]*gr);
plot([XP3(1,4) XP4(1,2)],[XP3(2,4) XP4(2,2)],':','linewidth',1,'color',[1 1 1]*gr);
% Combination 1
visualize_network(XP2+[XP1(1,3)-XP2(1,1);0],[],conn,'scale',fSc,...
                 'lcolor',CP(pI(1),:),'nalpha',.15,'lalpha',.15);
visualize_network(XP1,[],conn,'scale',fSc,'lcolor',CP(pI(1),:));
visualize_network(XP2,[],conn,'scale',fSc,'lcolor',CP(pI(1),:));
line_coordinates(XP1(:,1:2),-lSh,nW,lw);
line_coordinates(XP1(:,3:4)+[xSh/2;0],0,nW,lw);
line_coordinates(XP2(:,3:4),-lSh,nW,lw);
% Combination 2
visualize_network(XP4+[XP3(1,3)-XP4(1,1);0],[],conn,'scale',fSc,...
                 'lcolor',CP(pI(3),:),'nalpha',.15,'lalpha',.2);
visualize_network(XP3,[],conn,'scale',fSc,'lcolor',CP(pI(2),:));
visualize_network(XP4,[],conn,'scale',fSc,'lcolor',CP(pI(3),:));
line_coordinates(XP3(:,1:2),-lSh,nW,lw);
line_coordinates(XP3(:,3:4)+[xSh/2;0],0,nW,lw);
line_coordinates(XP4(:,3:4),-lSh,nW,lw);

% Set axis
axis([0 sRat(pInd) 0 1]*4 + [[1 1]*-0.7 [1 1]*-2.7]);

% Text
textb = ['~~~~~\textbf{c}~~combine nodes of units'];
text(labX,subp(pInd,4)+labY,textb,'Units','centimeters','fontsize',FS);
text(.08,.9,'$d_1$',NVTextH{:});
text(.08,.73,'2',NVTextH{:});
text(.08,.33,'1.8',NVTextH{:});
text(.5,.9,'$\textbf{f}(d_1)$',NVTextH{:});
text(.5,.55,'2',NVTextH{:});
text(.47,.12,'1.6',NVTextH{:});
text(0.95,.9,'$$\textbf{f}(\textbf{f}(d_1))$$',NVTextH{:});
text(0.95,.55,'2',NVTextH{:});
text(0.95,.12,'1.5',NVTextH{:});

set(gca,'visible',0);


%% d: 2 modules
% Define module
L = 1;
Xs = [0    0   0   L;...
      L/2 -L/2 0   0];
conn = [1 3; 1 4; 2 3; 2 4];

% Define combined module
Xsc = [0    0    0    L    L/2  L/2;...
       L/2 -L/2  0    0    0    L];
connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];


% Map: 2 fixed points
nM = 380;
% Calculate distances
[XC, fC] = sim_motion(Xs,[],conn,.001,nM,-Xs,0);
disp(max(fC));
d = squeeze(sqrt(sum(diff(XC,1,2).^2)));
d = d([1,3],:);
% Compute intermediate distance
nInt = 140;
ddiff = abs(d(1,:)-d(2,nInt));
nInt2 = find(ddiff == min(ddiff));

% Combined
nMc = 760;
[XCc, fC] = sim_motion(Xsc,[],connc,.001,nMc,-Xsc,0);
dc = squeeze(sqrt(sum(diff(XCc,1,2).^2)));
dc = dc([1,3,5],:);
pF = [1 760];
pFSh = [1 .7; 1.01 1.01];

ddiff1 = sum(abs([dc(1:2,:);dc(2:3,:)] - [d(:,nInt); d(:,nInt2)]));
nIntc = find(ddiff1 == min(ddiff1));
ddiff2 = sum(abs([dc(1:2,:);dc(2:3,:)] - [1 1 1 1]'*sqrt(2)));
nIntc2 = find(ddiff2 == min(ddiff2));


% Plot parameters: 1 module
pI = [1 nInt nInt2 380];
pISh = [ 0.130  .11  0.075  0.075;...
        -0.025 -.06 -0.077 -0.067];
CP = winter(floor(nM));
sC = 8/15;



% Begin figure
pInd = 4;
subplot('position',subpN(pInd,:)); cla;
hold on;
% Map
plot([0 4],[0 4], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(d(1,:)',d(2,:)',1,CP(1:size(d,2),:),'o','linewidth',.1);


% Examples
for i = 1:length(pI)
    XP = (XC(:,:,pI(i))-[XC(1,4,pI(i));0])*.08+d(:,pI(i))+pISh(:,i);
    plot(d(1,pI(i)),d(2,pI(i)),'s','markersize',3,'linewidth',3,'color',CP(pI(i),:));
    visualize_network(XP,[],conn,'scale',0.7,'lcolor',CP(pI(i),:));
end
for i = 1:length(pF)
    XP = (XCc(:,:,pF(i))-[XCc(1,4,pF(i));0])*.08+pFSh(:,i);
    visualize_network(XP,[],connc,'scale',0.7,'lcolor',CP(round(pF(i)/2),:));
end
hold off;
set(gca,'visible',1,'xtick',[],'ytick',[],'box',0);
a = [0 sRat(pInd) 0 1]*.6 + [[1 1]*.55 [1 1]*.55];
axis(a);

% Text
textd = '\textbf{d}~~motion of combined units';
text(labX,subp(pInd,4)+labY,textd,NVTitle{:});
text(subp(pInd,3)/2,-labY,'$d_1$',NVTitle{:});
text(labX,subp(pInd,4)/2,'$d_2$',NVTitle{:},'rotation',0);
text(.36,.46,'$d_2=d_1$',NVTextR{:},'color',[1 1 1]*gr,'rotation',45);

line([1 1], [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1],'linewidth',.7,'color','k');
line([1 1]/sqrt(2), [0 diff(a(3:4))*.02]+a(3),'linewidth',.7,'color','k');
line([0 diff(a(1:2))*.02]+a(1),[1 1]/sqrt(2), 'linewidth',.7,'color','k');
drawnow;


%% Size and Save Figure
fName = 'figure2b';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');