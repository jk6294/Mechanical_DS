% Figure 5: Dynamical Systems
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(5); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19  6.0];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .0 .0 .4];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  4.75 6.00];...
        [ 4.75  0.00 14.25 6.00]];
    
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

% Colormap: Distance
nT = 720;
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
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
hold on;

% Plot Stability
sc = .05;
sh1 = [.04;.81];
sh3 = [.72;.81];
sh4 = [.72;.45];
axsh = [.02 .06];
axL = [.7 .3];

% Unit
visualize_network(XC(:,1:3,300)*sc+sh1,XC(:,4:5,300)*sc+sh1,conn,...
                  'nalpha',.5,'lalpha',.5,'ucolor',C3a);
visualize_network(XC(:,1:3,1)*sc+sh1,XC(:,4:5,1)*sc+sh1,conn,'ucolor',C3a);

% Networks
scd = .16;
DLin = linspace(min(Da(:)),max(Da(:)),nT);
for i = 1:(nR+1)
    line_coordinates(XCa(:,i:i+1,end)*sc+sh4,'lsh',(-1)^i*.005,'color',interp1(DLin,CP1,Da(i,end)));
end
visualize_network(XCa(:,1:(nR+2),1)*sc+sh3,XCa(:,(nR+3):end,1)*sc+sh3,conna,'ucolor',C3a);
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
text((axsh(1)+axL(1)+.005)/sRat(pInd),axsh(2),'$k$',NVTextR{:});
text(axsh(1)-0,axsh(2)+axL(2)+.03,'$d_k$',NVTextRA{:});
text((mean(XCa(1,1:2,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.035,'$1$',NVTextH{:});
text((mean(XCa(1,2:3,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.035,'$2$',NVTextH{:});
text((mean(XCa(1,3:4,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.035,'$3$',NVTextH{:});
text((mean(XCa(1,5:6,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.035,'$\cdots$',NVTextH{:});
text((mean(XCa(1,nR+1:nR+2,end))*sc+sh4(1))/sRat(pInd),axsh(2)-.035,num2str(nR+1),NVTextH{:});
text(axsh(1)-.005,(3-Da(end,end))*scd+axsh(2)+.03,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(axsh(1)-.005,(0)*scd+axsh(2)+.04,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(axsh(1)+.05,axsh(2)+.18,'decays',NVTextR{:},'color',[1 1 1]*.7);
text(axsh(1)+.05,axsh(2)+.12,'faster than',NVTextR{:},'color',[1 1 1]*.7);
text(axsh(1)+.05,axsh(2)+.06,'exponential',NVTextR{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+.03)/sRat(pInd),axsh(2)+.25,'$d_{11} = \sqrt{3}$ to',NVTextRA{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+.03)/sRat(pInd),axsh(2)+.19,'numerical',NVTextRA{:},'color',[1 1 1]*.7);
text((axsh(1)+axL(1)+.03)/sRat(pInd),axsh(2)+.13,'error',NVTextRA{:},'color',[1 1 1]*.7);
text(.08,.97,'unit',NVTextH{:});
text(.6,.97,'combined',NVTextH{:});
text(.44,.62,'conformational motion',NVTextH{:});

hold off;
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);

% Text
text(labX,subp(pInd,4)+labY,'\textbf{a} \hspace{4mm}superstable network',NVTitle{:});


%% b: Mechanical AND gate. Make references to Newton Raphson optimization, and deadbeat control
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;

% scale and shift parameters
sc = .038;
sh = [.45; .45];
sh1 = [.45; .55];
sh2 = [1.15; .65];
sh3 = [1.15; .115];
sh4 = [2; .45];

% Correct position and rotation
% Top extension
XP1 = XCa4b(:,:,end); XP1 = XP1 - XP1(:,24);
Rz = rotz(atan2d(diff(XP1(2,[-2 0]+size(Xsa4,2))),diff(XP1(1,[-2 0]+size(Xsa4,2)))));
XP1 = Rz(1:2,1:2)'*XP1;
% Left extension
XP2 = XCa4a(:,:,end); XP2 = XP2 - XP2(:,24);
Rz = rotz(atan2d(diff(XP2(2,[-2 0]+size(Xsa4,2))),diff(XP2(1,[-2 0]+size(Xsa4,2)))));
XP2 = Rz(1:2,1:2)'*XP2;
% Full extension
XP3 = XCa4(:,:,end); XP3 = XP3 - XP3(:,24);
Rz = rotz(atan2d(diff(XP3(2,[-2 0]+size(Xsa4,2))),diff(XP3(1,[-2 0]+size(Xsa4,2)))));
XP3 = Rz(1:2,1:2)'*XP3;

% arrows
PSC1 = [0.00  0.20  1.00;...
        0.00  0.80  1.00] .*[.3;.12] + [.5;.7];
PSC2 = [0.00  0.20  1.00;...
        0.00 -0.80 -1.00] .*[.3;.12] + [.5;.3];
PSC3 = [0.00  0.80  1.00;...
        1.00  0.80  0.00] .*[.3;.12] + [1.25;.7];
PSC4 = [0.00  0.80  1.00;...
       -1.00 -0.80  0.00] .*[.3;.12] + [1.25;.3];
plot_spline(PSC1,'head',1,'headpos',1,'color',[1 1 1]*.8,'ratio',sRat(pInd));
plot_spline(PSC2,'head',1,'headpos',1,'color',[1 1 1]*.8,'ratio',sRat(pInd));
plot_spline(PSC3,'head',1,'headpos',1,'color',[1 1 1]*.8,'ratio',sRat(pInd));
plot_spline(PSC4,'head',1,'headpos',1,'color',[1 1 1]*.8,'ratio',sRat(pInd));

% Lines
% Initial
line_coordinates(Xsa2(:,1:2)*sc+sh,'lsh',0.005,'color',interp1(DLin,CP1,sqrt(3)));
line_coordinates(Xsa3(:,1:2)*sc+sh1,'lsh',0.005,'color',interp1(DLin,CP1,sqrt(3)));
line_coordinates(Xsa2(:,end-1:end)*sc+sh,'lsh',-0.005,'color',interp1(DLin,CP1,sqrt(3)));
% Open top
line_coordinates(XP1(:,1:2)*sc+sh2,'lsh',0.005,'color',interp1(DLin,CP1,sqrt(3)));
line_coordinates(XP1(:,[8 13])*sc+sh2,'lsh',0.005,'color',interp1(DLin,CP1,3));
line_coordinates(XP1(:,[31 32])*sc+sh2,'lsh',-0.005,'color',interp1(DLin,CP1,sqrt(3)));
% Open left
line_coordinates(XP2(:,1:2)*sc+sh3,'lsh',0.005,'color',interp1(DLin,CP1,3));
line_coordinates(XP2(:,[8 13])*sc+sh3,'lsh',0.005,'color',interp1(DLin,CP1,sqrt(3)));
line_coordinates(XP2(:,[31 32])*sc+sh3,'lsh',-0.005,'color',interp1(DLin,CP1,sqrt(3)));
% Final
line_coordinates(XP3(:,1:2)*sc+sh4,'lsh',0.005,'color',interp1(DLin,CP1,3));
line_coordinates(XP3(:,[8 13])*sc+sh4,'lsh',0.005,'color',interp1(DLin,CP1,3));
line_coordinates(XP3(:,[31 32])*sc+sh4,'lsh',-0.005,'color',interp1(DLin,CP1,3));

% Networks
visualize_network(Xsa2*sc+sh,Xua2*sc+sh,conna2,...
                  'msize',3,'ucolor',[1 1 1]*.5,'ucolor',C3a);
visualize_network(Xsa3*sc+sh1,Xua3*sc+sh1,conna3,...
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

% Text
text(0,.55,'$d_{\mathrm{in},1}$',NVTextR{:});
text(.05,.9,'$d_{\mathrm{in},2}$',NVTextR{:});
text(.27,.33,'$d_{\mathrm{out}}$',NVTextR{:});
% Text arrows
text(.22,.87,'open $d_{\mathrm{in},2}$',NVTextR{:},'color',[1 1 1]*.8);
text(.22,.25,'open $d_{\mathrm{in},1}$',NVTextR{:},'color',[1 1 1]*.8);
text(.50,.87,'open $d_{\mathrm{in},1}$',NVTextR{:},'color',[1 1 1]*.8);
text(.50,.25,'open $d_{\mathrm{in},2}$',NVTextR{:},'color',[1 1 1]*.8);
% Enumeration
text(.15,.86,'i',NVTextR{:});
text(.43,.42,'ii',NVTextR{:});
text(.40,1.03,'iii',NVTextR{:});
text(.715,0.93,'iv',NVTextR{:});
% Text table
xS = .88; yS = .95;
xSh = .055; ySh = .09;
text(xS+0*xSh,yS-0*ySh,'$d_{\mathrm{in},1}$',NVTextRA{:});
text(xS+1*xSh,yS-0*ySh,'$d_{\mathrm{in},2}$',NVTextRA{:});
text(xS+2*xSh,yS-0*ySh,'$d_{\mathrm{out}}$',NVTextRA{:});
text(xS+0*xSh,yS-1*ySh,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(xS+1*xSh,yS-1*ySh,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(xS+2*xSh,yS-1*ySh,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(xS+0*xSh,yS-2*ySh,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(xS+1*xSh,yS-2*ySh,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(xS+2*xSh,yS-2*ySh,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(xS+0*xSh,yS-3*ySh,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(xS+1*xSh,yS-3*ySh,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(xS+2*xSh,yS-3*ySh,'$\sqrt{3}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(3)));
text(xS+0*xSh,yS-4*ySh,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(xS+1*xSh,yS-4*ySh,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(xS+2*xSh,yS-4*ySh,'$3$',NVTextRA{:},'color',interp1(DLin,CP1,3));
text(xS-xSh,yS-1*ySh,'i',NVTextRA{:});
text(xS-xSh,yS-2*ySh,'ii',NVTextRA{:});
text(xS-xSh,yS-3*ySh,'iii',NVTextRA{:});
text(xS-xSh,yS-4*ySh,'iv',NVTextRA{:});
% Table lines
line(([1;1].*(-1:2)*xSh + xS)*sRat(pInd) + .014,...
     [.55;-4.45].*[1 1 1 1]*ySh + yS,'color','k','clipping',0);
line(([-1;2].*[1 1 1 1 1 1]*xSh + xS)*sRat(pInd) + .014,...
     [1;1].*(.55:-1:-4.45)*ySh + yS,'color','k','clipping',0);
line(([-1 2 2 -1 -1 2]*xSh+xS)*sRat(pInd)+.014,...
      [.55 .55 -4.45 -4.45 .55 .55]*ySh+yS,'color','k','clipping',0);

axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);

% Text
text(labX,subp(pInd,4)+labY,'\textbf{b} \hspace{1.5mm} mechanical AND gate',NVTitle{:});


%% Save
fName = 'figure5b';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');


