% Figure 5: Dynamical Systems
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(5); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19  10];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .0 .0 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  6.50  3.75 3.50];...
        [ 4.00  6.50  4.75 3.50];...
        [ 9.75  6.50  9.25 3.50];...
        [ 0.00  0.00 19.00 6.00]];
    
% Fontsize
FS = 10;
% Distance visualization parameters
nW = .08;
lw = .5;
gr = 0.9;
dgr = .6;
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

% Colormap: Distance
nT = 720;
o = [1 1 1];
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));
C3a = [197 066 085]/255;


%% Superstability simulations
% Define Unit
s = sqrt(3);
Xs = [-s/2 0 s/2;...
      -1/2 1 -1/2];
Xu = [-s/2 -s/2; sqrt(s^2-s^2/4) -sqrt(s^2-s^2/4)];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Simulation
[XC,fC] = sim_motion(Xs,Xu,conn,.001,1902,[Xs,Xu],0);
D = squeeze(sqrt(sum(diff(XC(:,[1 2 3],:),1,2).^2)));

% Combine
nR = 10;
[Xsa,Xua,conna] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR));

% Conformation
[XCa,fCa] = sim_motion8(Xsa,Xua,conna,.01,860,[Xsa,Xua],0);
XCa(1,:,:) = XCa(1,:,:) - XCa(1,12,:);
Da = squeeze(sqrt(sum(diff(XCa(:,(1:12),:),1,2).^2)));
DLin = linspace(min(Da(:))-.01,max(Da(:))+.01,nT);


%% Combine 2
nR2 = 10;
[Xsa2,Xua2,conna2] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR2*2));
[Xsa3,Xua3,conna3] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR2));
% Connectivity
% conna3(:,2) = conna3(:,2) + max(conna2(:,2));
% conna2(:,2) = conna2(:,2) + max(conna3(:,1));
% conna3(:,1) = conna3(:,1) + max(conna2(:,1));
% Shift
Xua2 = Xua2 - Xsa2(:,nR2+4); Xsa2 = Xsa2 - Xsa2(:,nR2+4);
Xua3 = Xua3 - Xsa3(:,nR2+2); Xsa3 = Xsa3 - Xsa3(:,nR2+2);
% Rotate
Rz = rotz(-60);
Xua3 = Rz(1:2,1:2)*Xua3; Xsa3 = Rz(1:2,1:2)*Xsa3; 
% Combine
[Xsa4,Xua4,conna4] = tesselate_network_old([Xsa2 Xsa3], [Xua2 Xua3],...
                     [conna2(:,1) conna2(:,2)+max(conna3(:,1));...
                      conna3(:,1)+max(conna2(:,1)) conna3(:,2)+max(conna2(:,2))], [0;0], [1;1]);
                  
% Simulate
X041 = zeros(2,size(Xsa4,2)+size(Xua4,2));
X042 = X041;
X041(1,1) = -1;
X042(:,8) = [-1/2;sqrt(3)/2];
[XCa4a,fCa4a] = sim_motion8(Xsa4,Xua4,conna4,.01,1221,X041,0);
[XCa4b,fCa4b] = sim_motion8(Xsa4,Xua4,conna4,.01,995,X042,0);
[XCa4c,fCa4c] = sim_motion8(XCa4a(:,1:size(Xsa4,2),end),...
                            XCa4a(:,(1:size(Xua4,2))+size(Xsa4,2),end),...
                            conna4,.025,2158,X042,0);
XCa4 = cat(3,XCa4a, XCa4c(:,:,2:end));
                        

%% a: Cobweb and superstability
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Plot Stability
sc = .095;
sh1 = [.25;.85];
sh3 = [.89;.17];

nV = 300;
XCp1 = (XC(:,:,1) - XC(:,2,1))*sc;
XCp2 = (XC(:,:,nV) - XC(:,2,nV))*sc;
Rz = rotz(-atan2d(XCp2(2,3)-XCp2(2,2),XCp2(1,3)-XCp2(1,2))-60);
XCp2 = Rz(1:2,1:2)*XCp2;
visualize_network(XCp2(:,1:3)+sh1,XCp2(:,4:5)+sh1,conn,...
                  'scolor',o*gr.^.2,'lalpha',.2,'ucolor',C3a);
visualize_network(XCp1(:,1:3)+sh1,XCp1(:,4:5)+sh1,conn,'ucolor',C3a);
line_coordinates(XCp1(:,1:2)+sh1,'style','-','lw',.7,'color',...
                 interp1(DLin,CP1,D(1,1)),'nw',.015,'lSh',.2);
line_coordinates(XCp1(:,2:3)+sh1,'style','-','lw',.7,'color',...
                 interp1(DLin,CP1,D(2,1)),'nw',.015,'lSh',.12);
line_coordinates(XCp2(:,1:2)+sh1,'style','-','lw',.7,'color',...
                 interp1(DLin,CP1,D(1,nV)).^.2,'nw',.015,'lSh',.14);
line_coordinates(XCp2(:,2:3)+sh1,'style','-','lw',.7,'color',...
                 interp1(DLin,CP1,D(2,nV)).^.2,'nw',.015,'lSh',.08);

visualize_network(XCa(:,1:(nR+2),1)*sc+sh3,XCa(:,(nR+3):end,1)*sc+sh3,conna,'ucolor',C3a);

axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);

% Text
texta = '\textbf{a}\hspace{3mm}superstable unit';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(0.03,.67,'$\delta l_1$',NVTextH{:});
text(0.41,.67,'$\delta l_2$',NVTextH{:});
text(.84,.85,'$\frac{\delta l_2}{\delta l_1} = 0$',NVTextRA{:});
text(.39,.45,'superstable network',NVTextH{:});


%% b: Conformational motion
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
scd = .25;
sh4 = [1.37;.65];
axsh = [.08 .05];
axL = [1.3 .45];


for i = 1:(nR+1)
    line_coordinates(XCa(:,i:i+1,end)*sc+sh4,'lsh',(-1)^(i<4)*(-1)^i*.015,'color',interp1(DLin,CP1,Da(i,end)));
end
visualize_network(XCa(:,1:(nR+2),end)*sc+sh4,XCa(:,(nR+3):end,end)*sc+sh4,conna,'ucolor',C3a);
% Ticks
line([1;1].*movmean(XCa(1,1:12,end),2,'endpoints','discard')*sc+sh4(1),...
     [0;1].*ones(2,nR+1)*.015+axsh(2),'color','k','linewidth',.5);
line([0;1].*ones(2,2)*.015+axsh(1),[1;1].*(s-Da(end,end))*scd+axsh(2)+.04,...
     'color',interp1(DLin,CP1,sqrt(3)),'linewidth',.5);
line([0;1].*ones(2,2)*.015+axsh(1),[1;1].*(3-Da(end,end))*scd+axsh(2)+.04,...
     'color',interp1(DLin,CP1,3),'linewidth',.5);
line([0 0 axL(1)]+axsh(1), [axL(2) 0 0]+axsh(2),'color','k','linewidth',.5);
% Points
scatter(movmean(XCa(1,1:12,end),2,'endpoints','discard')*sc+sh4(1),...
       (Da(:,end)'-Da(end,end))*scd+axsh(2)+.04, 20, interp1(DLin,CP1,Da(:,end)),'s','filled');

% Text
xSh = .05;
yShi = .1;
textb = '\textbf{b}\hspace{7mm}conformational motion';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
text(.3,.95,'open',NVTextH{:});
text(.75,.95,'closed',NVTextH{:});
text((axsh(1)+axL(1)+.005)/sRat(pInd),axsh(2),'$k$',NVTextR{:});
text(axsh(1)-.05,axsh(2)+axL(2)+.03,'$d_k$',NVTextRA{:});
text((mean(XCa(1,1:2,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.07,'$1$',NVTextH{:});
text((mean(XCa(1,2:3,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.07,'$2$',NVTextH{:});
text((mean(XCa(1,3:4,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.07,'$3$',NVTextH{:});
text((mean(XCa(1,4:5,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.07,'$\cdots$',NVTextH{:});
text((mean(XCa(1,nR+1:nR+2,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.07,num2str(nR+1),NVTextH{:});
text(axsh(1)-.05,(3-Da(end,end))*scd+axsh(2)+.03,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(axsh(1)-.05,(0)*scd+axsh(2)+.04,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(axsh(1)+.02,axsh(2)-.03+3*yShi,'decays',NVTextR{:},'color',[1 1 1]*.7);
text(axsh(1)+.02,axsh(2)-.03+2*yShi,'faster than',NVTextR{:},'color',[1 1 1]*.7);
text(axsh(1)+.02,axsh(2)-.03+yShi,'exponential',NVTextR{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+xSh)/sRat(pInd),axsh(2)+.12+2*yShi,'$d_{11} = \sqrt{3}$ to',NVTextRA{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+xSh)/sRat(pInd),axsh(2)+.12+yShi,'numerical',NVTextRA{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+xSh)/sRat(pInd),axsh(2)+.12,'error',NVTextRA{:},'color',[1 1 1]*.7);

axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% c: Combine chains
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc  = .065;
sh  = [.65; .2];
sh1 = [.65; .4];
sh2 = [2.350; .2];

% Networks
visualize_network(Xsa2*sc+sh,Xua2*sc+sh,conna2,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(Xsa3*sc+sh1,Xua3*sc+sh1,conna3,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(XCa4b(:,1:size(Xsa4,2),1)*sc+sh2,...
                  XCa4b(:,(1:size(Xua4,2))+size(Xsa4,2),1)*sc+sh2,conna4,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);

Xsp1 = Xsa2*sc+sh;
Xsp2 = Xsa3*sc+sh1;
scX1 = linspace(Xsp2(1,end-1),Xsp1(1,12),10);
scX2 = linspace(Xsp2(1,end),Xsp1(1,14),10);
scY1 = linspace(Xsp2(2,end-1),Xsp1(2,12),10);
scY2 = linspace(Xsp2(2,end),Xsp1(2,14),10);
scatter(scX1,scY1,1,o*dgr,'filled','marker','o','linewidth',1);
scatter(scX2,scY2,1,o*dgr,'filled','marker','o','linewidth',1);

arrow([1.2 1.5],[0 0]+.5,sRat(pInd),'linewidth',1.5,'headwidth',10,...
      'headlength',10,'lp',.9,'color',o*gr^2);

% Text
textc = '\textbf{c}\hspace{10mm}combine two networks by joining nodes';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});
text(.465,.62,'glue together',NVTextH{:},'color',o*gr^2);
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
    
    
%% d: Mechanical AND gate. Make references to Newton Raphson optimization, and deadbeat control
% Plot
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% scale and shift parameters
sc = .038;
sh1 = [.8; .59];
sh2 = [2.05; .59];
sh3 = [1.35; .115];
sh4 = [2.83; .127];
sP = sRat(pInd);

% Correct position and rotation
% Initial
XP0 = XCa4b(:,:,1); XP0 = XP0 - XP0(:,24);
% Top extension
XP1 = XCa4b(:,:,end); XP1 = XP1 - XP1(:,24);
Rz = rotz(atan2d(diff(XP1(2,[1 3])),diff(XP1(1,[1 3]))));
XP1 = Rz(1:2,1:2)'*XP1;
% Left extension
XP2 = XCa4a(:,:,end); XP2 = XP2 - XP2(:,24);
Rz = rotz(atan2d(diff(XP2(2,[1 3])),diff(XP2(1,[1 3]))));
XP2 = Rz(1:2,1:2)'*XP2;
% Full extension
XP3 = XCa4(:,:,end); XP3 = XP3 - XP3(:,24);
Rz = rotz(atan2d(diff(XP3(2,[1 3])),diff(XP3(1,[1 3]))));
XP3 = Rz(1:2,1:2)'*XP3;

% Lines
lSh = .062; nW = .01;
C_3 = interp1(DLin,CP1,3);
C_s3 = interp1(DLin,CP1,sqrt(3));
% Initial
line_coordinates(XP0(:,1:2)*sc+sh1,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP0(:,[8 13])*sc+sh1,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP0(:,[31 32])*sc+sh1,'lsh',-lSh/2,'nw',nW,'color',C_s3);
% Open d_in,2
line_coordinates(XP1(:,1:2)*sc+sh2,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP1(:,[8 13])*sc+sh2,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP1(:,[31 32])*sc+sh2,'lsh',-lSh/2,'nw',nW,'color',C_s3);
% Open d_in,1
line_coordinates(XP2(:,1:2)*sc+sh3,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP2(:,[8 13])*sc+sh3,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP2(:,[31 32])*sc+sh3,'lsh',-lSh/2,'nw',nW,'color',C_s3);
% Open both d_in,1 and d_in,2
line_coordinates(XP3(:,1:2)*sc+sh4,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP3(:,[8 13])*sc+sh4,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP3(:,[31 32])*sc+sh4,'lsh',-lSh/2,'nw',nW,'color',C_3);

% Networks
visualize_network(XP0(:,1:size(Xsa4,2))*sc+sh1,...
                  XP0(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh1,conna4,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(XP1(:,1:size(Xsa4,2))*sc+sh2,...
                  XP1(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh2,conna4,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(XP2(:,1:size(Xsa4,2))*sc+sh3,...
                  XP2(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh3,conna4,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(XP3(:,1:size(Xsa4,2))*sc+sh4,...
                  XP3(:,(1:size(Xua4,2))+size(Xsa4,2))*sc+sh4,conna4,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
              
% Arrows and annotations
arrow([0 .33]+.35+sh1(1),[0 0]+sh1(2)-.03,sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',10,'lp',.9,'color',o*gr^2);
arrow([0 .33]+.35+sh3(1),[0 0]+sh3(2)-.03,sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',10,'lp',.9,'color',o*gr^2);
arrow([0 .13]+.16+sh1(1),[0 -.13*s]-.125+sh1(2),sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',10,'lp',.9,'color',o*gr^2);
arrow([0 .13]+.16+sh2(1),[0 -.13*s]-.125+sh2(2),sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',10,'lp',.9,'color',o*gr^2);
line([.35 1]*.56-.04,-[.35 1]*s*.56+1,'color',o*gr,'linewidth',.7);

% Text
textd = '\textbf{d}\hspace{2mm} mechanical AND gate';
text(labX,subp(pInd,4)+labY-.1,textd,NVTitle{:});
% Network
text(.083,.64,'$d_{\mathrm{in},1}$',NVTextH{:});
text(.15,.92,'$d_{\mathrm{in},2}$',NVTextH{:});
text(.32,.66,'$d_{\mathrm{out}}$',NVTextH{:});
% Arrows and annotations
text(.273,.46,'open $d_{\mathrm{in},1}$',NVTextR{:},'color',[1 1 1]*.8,...
     'rotation',-60);
text(.65,.46,'open $d_{\mathrm{in},1}$',NVTextR{:},'color',[1 1 1]*.8,...
     'rotation',-60);
text(.345,.61,'open $d_{\mathrm{in},2}$',NVTextR{:},'color',[1 1 1]*.8);
text(.511,.135,'open $d_{\mathrm{in},2}$',NVTextR{:},'color',[1 1 1]*.8);
text(.02,.5,'inputs',NVTextH{:});
text(.02,.42,'$d_{\mathrm{in},1}$',NVTextH{:});
text(.02,.34,'$d_{\mathrm{in},2}$',NVTextH{:});
text(.02,.2,'output',NVTextH{:});
text(.02,.12,'$d_{\mathrm{out}}$',NVTextH{:});
% Networks
text(sh1(1)/sP+.025,sh1(2)+.15,'i',NVTextH{:});
text(sh2(1)/sP+.025,sh2(2)+.15,'ii $\alpha$',NVTextH{:});
text(sh3(1)/sP+.025,sh3(2)+.15,'ii $\beta$',NVTextH{:});
text(sh4(1)/sP+.025,sh4(2)+.15,'iii',NVTextH{:});
scatter([sh1(1) sh2(1) sh3(1) sh4(1)]+.025*sP,...
        [sh1(2) sh2(2) sh3(2) sh4(2)]+.15,400,'k','linewidth',.7);
% Table
xS = 2.95; yS = .85;
xSh = .18; ySh = -.09; 
txSh = -.02; tySh = .04;
tM = {'', '$\sqrt{3}~$', '$3~$';...
      '$\sqrt{3}~$', '$\sqrt{3}~$', '$\sqrt{3}~$';...
      '$3~$', '$\sqrt{3}~$', '$3~$'};
lM = {'i', 'ii $\alpha$';...
      'ii $\beta$', 'iii'};
tC = {[0 0 0], interp1(DLin,CP1,s), interp1(DLin,CP1,3);...
      interp1(DLin,CP1,s),interp1(DLin,CP1,s),interp1(DLin,CP1,s);...
      interp1(DLin,CP1,3),interp1(DLin,CP1,s),interp1(DLin,CP1,3)};
for i = 1:3
    for j = 1:3
        text((xS+(i-1)*xSh+txSh)/sP,...
             yS+(j-1)*ySh+tySh,tM{j,i},NVTextRA{:},...
             'color', tC{j,i});
    end
end
for i = 1:2
    for j = 1:2
        text((xS+(i-1)*xSh+.008)/sP,...
             yS+j*ySh+tySh+.025,lM{j,i},NVTextR{:},'fontsize',6);
    end
end
text((-xSh+xS+txSh)/sP,ySh+yS,'$d_{\mathrm{in},1}$',NVTextRA{:});
text((1.5*xSh+xS)/sP,-ySh+tySh+yS,'$d_{\mathrm{in},2}$',NVTextRA{:});
text((xS-xSh)/sP,yS-1*ySh+tySh,'truth table',NVTextH{:});
text((xS-xSh)/sP,yS-0*ySh+tySh,'$d_{\mathrm{out}}$',NVTextH{:});
% Table lines
line([-1;2]*xSh+xS,[1;1]*ySh+yS,'color','k');
line([1;1]*xSh+xS,[-1;2]*ySh+yS,'color','k');
line([-1;-1]*xSh+xS,[0;2]*ySh+yS,'color',o*gr^2);
line([0;2]*xSh+xS,[-1;-1]*ySh+yS,'color',o*gr^2);
line([-2;2;2;-2;-2;2]*xSh+xS,[2;2;-2;-2;2;2]*ySh+yS,...
     'color','k','clipping',0);
line([-2;2]*xSh+xS,[0;0]*ySh+yS,'color','k','linewidth',1,'clipping',0);
line([0;0]*xSh+xS,[-2;2]*ySh+yS,'color','k','linewidth',1,'clipping',0);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% Save
fName = 'figure5d';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
% fig.PaperUnits = 'centimeters';

fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');


