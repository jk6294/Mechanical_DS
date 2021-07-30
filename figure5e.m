% Figure 5: Dynamical Systems
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(5); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19  7.0];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .0 .1 .3];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  5.00 7.00]
        [ 5.00  0.00 14.00 7.00]];
    
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
nR2 = 12;
[Xsa2,Xua2,conna2] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR2));
[Xsa3,Xua3,conna3] = network_chain_x(Xs(:,1:2), Xu, ones(1,nR2-1));

% Center
Xua2 = Xua2 - Xsa2(:,end);
Xsa2 = Xsa2 - Xsa2(:,end);
Xua3 = Xua3 - Xsa3(:,end);
Xsa3 = Xsa3 - Xsa3(:,end);

% % Flip
Xua3(2,:) = -Xua3(2,:);
Xsa3(2,:) = -Xsa3(2,:);

% Rotate
Rz = rotz(-60);
Xua3 = Rz(1:2,1:2)*Xua3;
Xsa3 = Rz(1:2,1:2)*Xsa3;

% figure(6); clf;
% visualize_network(Xsa2,Xua2,conna2,'ucolor',[0 0 1]);
% visualize_network(Xsa3,Xua3,conna3,'ucolor',[1 0 0]);

% Combine
[Xsa4,Xua4,conna4] = tesselate_network_old([Xsa2 Xsa3], [Xua2 Xua3],...
                     [conna2(:,1) conna2(:,2)+max(conna3(:,1));...
                      conna3(:,1)+max(conna2(:,1)) conna3(:,2)+max(conna2(:,2))], [0;0], [1;1]);

% Simulate
X041 = zeros(2,size(Xsa4,2)+size(Xua4,2));
X042 = X041;
X041(1,1) = -1;
X042(:,11) = [-1/2;sqrt(3)/2];
[XCa4a,fCa4a] = sim_motion8(Xsa4,Xua4,conna4,.01,1161,X041,0);
[XCa4b,fCa4b] = sim_motion8(Xsa4,Xua4,conna4,.01,1168,X042,0);
[XCa4c,fCa4c] = sim_motion8(XCa4a(:,1:size(Xsa4,2),end),...
                            XCa4a(:,(1:size(Xua4,2))+size(Xsa4,2),end),...
                            conna4,.025,1595,X042,0);
XCa4 = cat(3,XCa4a, XCa4c(:,:,2:end));
                        

%% a: Cobweb and superstability
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Plot Stability
sc = .045;
sh1 = [.05;.91];
sh3 = [.44;.65];

% Unit
line_coordinates(Xs(:,1:2)*sc+sh1,'color',interp1(DLin,CP1,D(1,1)),...
                 'lSh',.007);
line_coordinates(Xs(:,2:3)*sc+sh1,'color',interp1(DLin,CP1,D(2,1)),...
                 'lSh',.007);
visualize_network(Xs*sc+sh1,Xu*sc+sh1,conn,'ucolor',C3a);
% Network
visualize_network(XCa(:,1:(nR+2),1)*sc+sh3,XCa(:,(nR+3):end,1)*sc+sh3,conna,'ucolor',C3a);

% Conformational motion
% Parameters
scd = .12;
sh4 = [0.63;.32];
axsh = [.02 .035];
axL = [0.62 .20];

for i = 1:(nR+1)
    line_coordinates(XCa(:,i:i+1,end)*sc+sh4,'lsh',(-1)^(i<4)*(-1)^i*.007,'color',interp1(DLin,CP1,Da(i,end)));
end
visualize_network(XCa(:,1:(nR+2),end)*sc+sh4,XCa(:,(nR+3):end,end)*sc+sh4,conna,'ucolor',C3a);
% Ticks
line([1;1].*movmean(XCa(1,1:12,end),2,'endpoints','discard')*sc+sh4(1),...
     [0;1].*ones(2,nR+1)*.015+axsh(2),'color','k','linewidth',.5);
line([0;1].*ones(2,2)*.015+axsh(1),[1;1].*(s-Da(end,end))*scd+axsh(2)+.03,...
     'color',interp1(DLin,CP1,sqrt(3)),'linewidth',.5);
line([0;1].*ones(2,2)*.015+axsh(1),[1;1].*(3-Da(end,end))*scd+axsh(2)+.03,...
     'color',interp1(DLin,CP1,3),'linewidth',.5);
line([0 0 axL(1)]+axsh(1), [axL(2) 0 0]+axsh(2),'color','k','linewidth',.5);
% Points
scatter(movmean(XCa(1,1:12,end),2,'endpoints','discard')*sc+sh4(1),...
       (Da(:,end)'-Da(end,end))*scd+axsh(2)+.03, 20, interp1(DLin,CP1,Da(:,end)),'s','filled');
% Text
xSh = .05;
yShi = .05;
text(.2,.445,'open',NVTextH{:},'color',interp1(DLin,CP1,3));
text(.7,.445,'closed',NVTextH{:},'color',interp1(DLin,CP1,s));
text((axsh(1)+axL(1)+.005)/sRat(pInd),axsh(2),'$k$',NVTextR{:});
text(axsh(1)-.1,axsh(2)+axL(2)+.03,'$l_k$',NVTextR{:});
text((mean(XCa(1,1:2,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.03,'$1$',NVTextH{:});
text((mean(XCa(1,2:3,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.03,'$2$',NVTextH{:});
text((mean(XCa(1,3:4,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.03,'$3$',NVTextH{:});
text((mean(XCa(1,4:5,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.03,'$\cdots$',NVTextH{:});
text((mean(XCa(1,nR+1:nR+2,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.03,num2str(nR+1),NVTextH{:});
text(axsh(1)-.01,(3-Da(end,end))*scd+axsh(2)+.025,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(axsh(1)-.01,(0)*scd+axsh(2)+.025,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(axsh(1)+.05,axsh(2)-.01+3*yShi,'decays',NVTextR{:},'color',[1 1 1]*.7);
text(axsh(1)+.05,axsh(2)-.01+2*yShi,'faster than',NVTextR{:},'color',[1 1 1]*.7);
text(axsh(1)+.05,axsh(2)-.01+yShi,'exponential',NVTextR{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+xSh)/sRat(pInd)-.02,axsh(2)+.07+2*yShi,'$d_{11} = \sqrt{3}$ to',NVTextRA{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+xSh)/sRat(pInd)-.02,axsh(2)+.07+yShi,'numerical',NVTextRA{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+xSh)/sRat(pInd)-.02,axsh(2)+.07,'error',NVTextRA{:},'color',[1 1 1]*.7);
   
line([0 1]*.6+.3,[0;-s]*.6+1.03,'color',o*gr,'linewidth',.7,'clipping',0)

axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);

% Text
texta = '\textbf{a}\hspace{2mm}unit';
textb = '\textbf{b}\hspace{2mm}network';
textc = '\textbf{c}\hspace{2mm}conformational motion';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(labX,subp(pInd,4)+labY-1.55,textb,NVTitle{:});
text(labX,subp(pInd,4)+labY-3.45,textc,NVTitle{:});
text(-.02,.95,'$l_1$',NVTextH{:},'color',interp1(DLin,CP1,D(2,1)));
text(0.15,.95,'$l_2$',NVTextH{:},'color',interp1(DLin,CP1,D(2,1)));

    
%% d: Mechanical AND gate. Make references to Newton Raphson optimization, and deadbeat control
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% scale and shift parameters
sc = .03; shn = .74;
sh1 = [.85; .73];
sh2 = sh1 + [1;-s]/2*shn;
sh3 = sh1 + [1;0]*shn;
sh4 = sh1 + [3;-s]/2*shn;
sh0a = sh1 + [-.68;-.05]-[-.8;s]*.05;
sh0b = sh1 + [-.68;-.05];
sP = sRat(pInd);

% Connect networks
lp1 = (sh0b-sh0a).*linspace(0,1,10)+sh0a;
scatter(lp1(1,:),lp1(2,:),1,o*gr,'filled');
scatter(lp1(1,:)+Xsa3(1,end-1)*sc,lp1(2,:)+Xsa3(2,end-1),1,o*gr,'filled');
visualize_network(Xsa2*sc+sh0a,Xua2*sc+sh0a,conna2,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(Xsa3*sc+sh0b,Xua3*sc+sh0b,conna3,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
% Arrow
plot_spline([0 .4 .6 1;0 .15 .85 1].*[.27;.19]+[.23;.65],'head',1,...
             'headpos',1,'color',o*gr);
% Text
text(-.34,.9,'combine','fontsize',FS,'color',o*dgr);
text(-.34+.03,.9-.03*s,'two','fontsize',FS,'color',o*dgr);
text(-.34+.06,.9-.06*s,'networks','fontsize',FS,'color',o*dgr);


% Correct position and rotation
% Initial
XP0 = XCa4b(:,:,1); XP0 = XP0 - XP0(:,size(Xsa4,2));
% Top extension
XP1 = XCa4b(:,:,end); XP1 = XP1 - XP1(:,size(Xsa4,2));
Rz = rotz(atan2d(diff(XP1(2,[1 3])),diff(XP1(1,[1 3]))));
XP1 = Rz(1:2,1:2)'*XP1;
% Left extension
XP2 = XCa4a(:,:,end); XP2 = XP2 - XP2(:,size(Xsa4,2));
Rz = rotz(atan2d(diff(XP2(2,[1 3])),diff(XP2(1,[1 3]))));
XP2 = Rz(1:2,1:2)'*XP2;
% Full extension
XP3 = XCa4(:,:,end); XP3 = XP3 - XP3(:,size(Xsa4,2));
Rz = rotz(atan2d(diff(XP3(2,[1 3])),diff(XP3(1,[1 3]))));
XP3 = Rz(1:2,1:2)'*XP3;

% Lines
lSh = .055; nW = .008;
C_3 = interp1(DLin,CP1,3);
C_s3 = interp1(DLin,CP1,sqrt(3));
% Initial
line_coordinates(XP0(:,1:2)*sc+sh1,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP0(:,[8 11])*sc+sh1,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP0(:,[23 25])*sc+sh1,'lsh',-lSh/2,'nw',nW,'color',C_s3);
% Open d_in,2
line_coordinates(XP1(:,1:2)*sc+sh2,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP1(:,[8 11])*sc+sh2,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP1(:,[23 25])*sc+sh2,'lsh',-lSh/2,'nw',nW,'color',C_s3);
% Open d_in,1
line_coordinates(XP2(:,1:2)*sc+sh3,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP2(:,[8 11])*sc+sh3,'lsh',lSh,'nw',nW,'color',C_s3);
line_coordinates(XP2(:,[23 25])*sc+sh3,'lsh',-lSh/2,'nw',nW,'color',C_s3);
% Open both d_in,1 and d_in,2
line_coordinates(XP3(:,1:2)*sc+sh4,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP3(:,[8 11])*sc+sh4,'lsh',lSh/2,'nw',nW,'color',C_3);
line_coordinates(XP3(:,[23 25])*sc+sh4,'lsh',-lSh/2,'nw',nW,'color',C_3);

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
aL1 = .2; aL2 = aL1/2;
aSh1 = [.03; .11];
aSh2 = [-.11; -.13];
arrow([0 aL1]+aSh1(1)+sh1(1),[0 0]+aSh1(2)+sh1(2),sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',7,'lp',.9,'color',o*gr^2);
arrow([0 aL1]+aSh1(1)+sh2(1),[0 0]+aSh1(2)+sh2(2),sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',7,'lp',.9,'color',o*gr^2);
arrow([0 aL2]+aSh2(1)+sh1(1),[0 -aL2*s]+aSh2(2)+sh1(2),sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',7,'lp',.9,'color',o*gr^2);
arrow([0 aL2]+aSh2(1)+sh3(1),[0 -aL2*s]+aSh2(2)+sh3(2),sRat(pInd),'linewidth',1.5,...
      'headwidth',10,'headlength',7,'lp',.9,'color',o*gr^2);

% Text
textd = '\textbf{d}\hspace{3mm}a mechanical AND gate';
text(labX-2.35,subp(pInd,4)+labY,textd,NVTitle{:});
% Arrows and annotations
text(.34,.58,'open $l_{\mathrm{in},b}$',NVTextR{:},'color',[1 1 1]*.8,...
     'rotation',-60);
text(.7,.58,'open $l_{\mathrm{in},b}$',NVTextR{:},'color',[1 1 1]*.8,...
     'rotation',-60);
text(.426,.892,'open $l_{\mathrm{in},a}$',NVTextR{:},'color',[1 1 1]*.8);
text(.606,.25,'open $l_{\mathrm{in},a}$',NVTextR{:},'color',[1 1 1]*.8);
sSh = [1;-s]*.033;
tXs = [.106;.42];
text(tXs(1),tXs(2)+0*sSh(2),'inputs',NVTextRA{:});
text(tXs(1),tXs(2)+1*sSh(2),'$l_{\mathrm{in},a}$',NVTextRA{:});
text(tXs(1),tXs(2)+2*sSh(2),'$l_{\mathrm{in},a}$',NVTextRA{:});
text(tXs(1),tXs(2)+4*sSh(2),'output',NVTextRA{:});
text(tXs(1),tXs(2)+5*sSh(2),'$l_{\mathrm{out}}$',NVTextRA{:});
% Networks
text(.188,.7,'$l_{\mathrm{in},a}$',NVTextR{:});
text(.262,.96,'$l_{\mathrm{in},b}$',NVTextR{:});
text(.43,.7,'$l_{\mathrm{out}}$',NVTextR{:});
xSh = -.12; ySh = .11;
text(sh1(1)/sP+xSh,sh1(2)+ySh,'i',NVTextH{:});
text(sh2(1)/sP+xSh,sh2(2)+ySh,'ii$b$',NVTextH{:});
text(sh3(1)/sP+xSh,sh3(2)+ySh,'ii$a$',NVTextH{:});
text(sh4(1)/sP+xSh,sh4(2)+ySh,'iii',NVTextH{:});
scatter([sh1(1) sh2(1) sh3(1) sh4(1)]+xSh*sP,...
        [sh1(2) sh2(2) sh3(2) sh4(2)]+ySh,400,'k','linewidth',.7);
    
% Table
xS = 1.77; yS = 0.96;
xSi = .092; ySi = -.075;
txSh = -.006; tySh = .033;
line([0;1].*[1 1 1 1 1]*3*xSi+xS, [0;0]+(0:4)*ySi+yS,'color','k');
line([0;0]+(0:3)*xSi+xS, [0;4].*[1 1 1 1]*ySi+yS,'color','k','clipping',0);
line([0;3]*xSi+xS, [0;0]+yS,'color','k','linewidth',1);
line([2;2]*xSi+xS, [0;4]*ySi+yS,'color','k','linewidth',1,'clipping',0);
line([0 3 3 0 0 3]*xSi+xS, [0 0 4 4 0 0]*ySi+yS,'color','k');
tM = {'',      '$l_{\mathrm{in},a}$','$l_{\mathrm{in},b}$','$l_{\mathrm{out}}$';...
      'i',     '$\sqrt{3}$',         '$\sqrt{3}$',         '$\sqrt{3}$';...
      'ii$a$', '$3$',                '$\sqrt{3}$',         '$\sqrt{3}$';...
      'ii$b$', '$\sqrt{3}$',         '$3$',                '$\sqrt{3}$';...
      'iii',   '$3$',                '$3$',                '$3$'};
tC = {[0 0 0], [0 0 0], [0 0 0], [0 0 0];...
      [0 0 0],interp1(DLin,CP1,s),interp1(DLin,CP1,s),interp1(DLin,CP1,s);...
      [0 0 0],interp1(DLin,CP1,3),interp1(DLin,CP1,s),interp1(DLin,CP1,s);...
      [0 0 0],interp1(DLin,CP1,s),interp1(DLin,CP1,3),interp1(DLin,CP1,s);...
      [0 0 0],interp1(DLin,CP1,3),interp1(DLin,CP1,3),interp1(DLin,CP1,3)};
for i = 1:5
    for j = 1:4
        text((xS+(j-1)*xSi+txSh)/sP,...
              yS+(i-1)*ySi+tySh,tM{i,j},NVTextRA{:},...
              'color', tC{i,j});
    end
end

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% Save
fName = 'figure5e';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');


