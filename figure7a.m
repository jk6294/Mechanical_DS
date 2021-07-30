% Figure 7: Period doubling route to chaos
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(7); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19  9.50];
% Margins in cm, [l,r,d,u]
fMarg = [.4 0 .4 0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  4.75  4.75  4.75];...
        [ 0.00  0.00  4.75  4.75];...
        [ 5.00  3.50 14.00  6.00];...
        [ 5.20  0.00 13.80  3.50]];
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


%% Figure parameters
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
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));
% Gradient: Slope. Interpolate between 2 colors
SLin = linspace(-.6,-1.605,nT);
C3a = [197 066 085]/255;
C3b = [031 172 204]/255;
CP3 = interp1([0 1],[C3a;C3b],linspace(0,1,nT));


%% a: Velocities
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .2;
sh = [.25;.85];
mS = 5;

% Construct
L = 1;
th = 45;
Rz = rotz(th); Rz = Rz(1:2,1:2);
sh2 = -Rz*[.012;0];
nVal = linspace(-1.6,-0.6,4);
d1dot = -Rz*[1;0];
Xs = [-L*cosd(th)  0  L*cosd(th);...
      -L*sind(th)  0 -L*sind(th)];

% Draw arrows
for i = 1:length(nVal)
    d2dot = [nVal(i);0];
    dX = [d1dot zeros(2,1) Rz'*d2dot];
    % Plot
    arrow((Xs(1,1)+[0;dX(1,1)/2])*sc+sh(1)+sh2(1),...
          (Xs(2,1)+[0;dX(2,1)/2])*sc+sh(2)+sh2(2),sRat(pInd));
    arrow((Xs(1,3)+[0;dX(1,3)/2])*sc+sh(1)+sh2(1)*(i-2.5),...
          (Xs(2,3)+[0;dX(2,3)/2])*sc+sh(2)+sh2(2)*(i-2.5),sRat(pInd),...
          'color',interp1(SLin,CP3,nVal(i)));
end
% Network
visualize_network(Xs*sc+sh,[],[1 1],'msize',mS);
line_coordinates(Xs(:,1:2)*sc+sh,'lSh',.05,'style','-',...
                 'color',0*o,'nw',.01,'lw',.5);
line_coordinates(Xs(:,2:3)*sc+sh,'lSh',.05,'style','-',...
                 'color',0*o,'nw',.01,'lw',.5);

% Phase
syms({'a','b'},'real');
sh = [.25;.4];
R = [-1.25 1.25; -2.00 0.50];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
Xu = zeros(2,2,length(nVal));
for i = length(nVal):-1:1
    d2dot = [nVal(i);0];
    dX = [d1dot zeros(2,1) Rz'*d2dot];
    % Plot
    visualize_conic(Xs*sc+sh,dX*sc,R*sc+sh,'ucolori',...
                    interp1(SLin,CP3,nVal(i)),'overlay',.99);
    % Solve
    [Q,W,v,err] = construct_conic(Xs,dX,0);
    P = [W(1:2,:), v(1:2);...
         zeros(1,2), 1]^-1;
    Q = P'*Q*P;
    E1 = [a b 1]*Q*[a;b;1]==0;
    if(nVal(i) <-1)
        E2 = a == -b - sind(th);
        nP = 1;
    else
        E2 = a == b + sind(th);
        nP = 2;
    end
    an = solve([E1;E2]);
    Xu(:,:,i) = [0 double(an.a(nP)); -1/sind(th) double(an.b(nP))];
end
for i = length(nVal):-1:1
    if(i == 1)
        visualize_network(Xs*sc+sh,Xu(:,:,i)*sc+sh,conn,...
            'ucolor',interp1(SLin,CP3,nVal(i)),'msize',mS);
    else
        visualize_network(Xs*sc+sh,Xu(:,:,i)*sc+sh,[1 1],...
            'ucolor',interp1(SLin,CP3,nVal(i)),'msize',mS);
    end
end
colormap(CP3);
cNr = mS/(2*72)*2.54/subp(pInd,4);
cNx = -cNr*cosd(linspace(0,360,120));
cNy = cNr*sind(linspace(0,360,120));
patch(cNx+sh(1),cNy+Xu(2,1,1)*sc+sh(2),...
      1./(1+exp(20*[linspace(1,-.7,120/2) linspace(-.7,1,120/2)])),...
      'clipping',0,'linewidth',.5);

% Legend
xSh = .75;
ySh = .885; 
ySh2 = .05;
yShi = .164;
lS = linspace(-1,1,100)*.135;

visualize_network([xSh;ySh-ySh2-.01],[],[1 1],'msize',mS);
arrow([1 -1]*.05+xSh,[1 1]*ySh-ySh2-1*yShi,sRat(pInd));
scatter(lS+xSh,ones(1,length(lS))*ySh-ySh2-2*yShi,.3,...
        CP3(floor(linspace(1,nT,length(lS))),:),'filled');
scatter(lS+xSh,ones(1,length(lS))*ySh-ySh2-3*yShi,1,...
        CP3(floor(linspace(1,nT,length(lS))),:),'filled','marker','s');
patch(cNx+xSh,cNy+ySh-ySh2-4*yShi-.01,...
      1./(1+exp(20*[linspace(1,-.7,120/2) linspace(-.7,1,120/2)])),...
      'clipping',0,'linewidth',.5);
visualize_network([-.115 .115;0 0]+[xSh;ySh-ySh2-4*yShi-.01],[],[1 1],...
                  'msize',mS,'scolor',CP3([1 end],:));
scatter(lS+xSh,ones(1,length(lS))*ySh-ySh2-5*yShi,1,...
        CP1(floor(linspace(1,nT,length(lS))),:),'filled','marker','s');

% Text
texta = '\textbf{a}\hspace{2mm}design motion';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(0,.53,'solution space',NVTextR{:});
text( .1,.32,'$1$',NVTextH{:},'color',o*gr^2);
text( .2,.39,'$2$',NVTextH{:},'color',o*gr^2);
text( .4,.32,'$3$',NVTextH{:},'color',o*gr^2);
text( .15,.645,'$dl_1$',NVTextH{:});
text( .35,.645,'$dl_2$',NVTextH{:});
text( .10,.85,'$l_1$',NVTextH{:});
text( .40,.85,'$l_2$',NVTextH{:});
text(xSh,ySh,'position',NVTextH{:});
text(xSh,ySh-1*yShi,'motion',NVTextH{:});
text(xSh,ySh-2*yShi,'$dl_2 / dl_1$',NVTextH{:});
text(xSh,ySh-3*yShi,'solution',NVTextH{:});
text(xSh,ySh-4*yShi,'added',NVTextH{:});
text(xSh,ySh-5*yShi,'distance',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% b: Simulate
% Simulate
nSim = 60;
[XC1a,fC1a] = sim_motion10(Xs,Xu(:,:,1),conn,.001,nSim, [Xs Xu(:,:,1)],0);
[XC1b,fC1b] = sim_motion10(Xs,Xu(:,:,1),conn,.001,nSim,-[Xs Xu(:,:,1)],0);
[XC2a,fC2a] = sim_motion10(Xs,Xu(:,:,2),conn,.001,nSim, [Xs Xu(:,:,2)],0);
[XC2b,fC2b] = sim_motion10(Xs,Xu(:,:,2),conn,.001,nSim,-[Xs Xu(:,:,2)],0);
[XC3a,fC3a] = sim_motion10(Xs,Xu(:,:,3),conn,.001,nSim, [Xs Xu(:,:,3)],0);
[XC3b,fC3b] = sim_motion10(Xs,Xu(:,:,3),conn,.001,nSim,-[Xs Xu(:,:,3)],0);
[XC4a,fC4a] = sim_motion10(Xs,Xu(:,:,4),conn,.001,nSim, [Xs Xu(:,:,4)],0);
[XC4b,fC4b] = sim_motion10(Xs,Xu(:,:,4),conn,.001,nSim,-[Xs Xu(:,:,4)],0);

% Combine
XC1 = cat(3,flip(XC1a(:,:,2:end),3),XC1b);
XC2 = cat(3,flip(XC2a(:,:,2:end),3),XC2b);
XC3 = cat(3,flip(XC3a(:,:,2:end),3),XC3b);
XC4 = cat(3,flip(XC4a(:,:,2:end),3),XC4b);

% Distances
D1 = sqrt(squeeze(sum(diff(XC1(:,1:3,:),1,2).^2)));
D2 = sqrt(squeeze(sum(diff(XC2(:,1:3,:),1,2).^2)));
D3 = sqrt(squeeze(sum(diff(XC3(:,1:3,:),1,2).^2)));
D4 = sqrt(squeeze(sum(diff(XC4(:,1:3,:),1,2).^2)));


%% b: Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;

% Axes
ax1 = min([min(D1(:)) min(D2(:)) min(D3(:)) min(D4(:))])-.03;
ax2 = max([max(D1(:)) max(D2(:)) max(D3(:)) max(D4(:))])+.03;
ax = [ax1 ax2 ax1 ax2] + [-1 2 -1 2]*.03;
axL = ax + [.003 -.03 .003 -.03];
% Plot
hold on;
% Axes
arrow(axL(1:2),[1 1]*axL(3),sRat(pInd),'color','k');
arrow([1 1]*axL(1),axL(3:4),sRat(pInd),'color','k');
line([1;1].*1, [0;.005].*[1 1]+axL(3),'linewidth',.5,'color','k');
line([0;.005].*[1 1]+axL(1), [1;1].*1,'linewidth',.5,'color','k');
text(.85,-.05,'$l_1$',NVTextR{:});
text(-.08,.86,'$l_2$',NVTextR{:});
text(.43,-.05,'$1$',NVTextR{:});
text(-.05,.45,'$1$',NVTextR{:});
% Map
plot(D1(1,:),D1(2,:),'color',interp1(SLin,CP3,nVal(1)),'linewidth',1);
plot(D2(1,:),D2(2,:),'color',interp1(SLin,CP3,nVal(2)),'linewidth',1);
plot(D3(1,:),D3(2,:),'color',interp1(SLin,CP3,nVal(3)),'linewidth',1);
plot(D4(1,:),D4(2,:),'color',interp1(SLin,CP3,nVal(4)),'linewidth',1);
plot([.88 1.12], [1.12 .88],'-','linewidth',.5,...
     'color',interp1(SLin,CP3,-1));
% Networks
sc = .036;
sh = [1.09 0.99 0.915 0.915;...
      1.12 1.12 1.040 0.940];
for i = 1:4
    visualize_network(Xs*sc+sh(:,i),Xu(:,:,i)*sc+sh(:,i),conn,...
                      'ucolor',interp1(SLin,CP3,nVal(i)));
end

% Text
textb = '\textbf{b}\hspace{3.5mm}from stable to unstable';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
text(.7,.27,'$l_2 = -l_1$',NVTextR{:},'rotation',-45,...
     'color',interp1(SLin,CP3,-1));
text(.21,.42,'stable', NVTextR{:},'color',CP3(1,:));
text(.37,.63,'unstable', NVTextR{:},'color',CP3(end,:));
text(.72,.53,'slope:',NVTextH{:});
text(.72,.44,'$s=dl_2 / dl_1$',NVTextH{:});

% Axes
axis(ax);
set(gca,'visible',0);


%% Phase
% Added node positions
s = .5/sind(th);
Xupd = @(vV) [0 vV; -2*s -L*sind(th)-vV];

% Period double
nS = 5000;
nW = linspace(100000,100000,nS);
% nValpd = .275547;
% nValpd = .2175;       % Chaos
nValpd = .214659;
nValpdc = .217;

load period_double_params.mat;
% Uncomment to get period doubling bifurcation values in real time
% [cV,lc,LC,dlc] = period_double(Xs,Xupd(.01),Xupd(nValpd),nS,1.2,nW);

% Find boundaries
nStr = 1; nEnd = 0;
cVp = cV(nStr:(end-nEnd));
lcp = lc(nStr:(end-nEnd));
LCp = LC(nStr:(end-nEnd));
lcDI = [0 find(diff(lcp)>0)];
CP2 = parula(length(lcDI)-0);

% Construct  and simulate networks
pLC = 400;
Xsc1 = Xs; Xuc1 = Xupd(0); connc1 = conn;
[Xsc2,Xuc2,connc2] = network_chain_x(Xs(:,1:2),Xupd(cVp(pLC)*nValpd),ones(1,2));
[Xsc3,Xuc3,connc3] = network_chain_x(Xs(:,1:2),Xupd(cVp(lcDI(4))*nValpd),ones(1,4));
[Xsc4,Xuc4,connc4] = network_chain_x(Xs(:,1:2),Xupd(cVp(lcDI(5))*nValpd),ones(1,8));
[Xsc5,Xuc5,connc5] = network_chain_x(Xs(:,1:2),Xupd(nValpdc),ones(1,22));
% Simulate
XCc1 = [Xsc1, Xuc1];
[XCc2,fC2] = sim_motion10(Xsc2,Xuc2,connc2,.005,200,[Xsc2,Xuc2],0);
[XCc3,fC3] = sim_motion10(Xsc3,Xuc3,connc3,.005,400,[Xsc3,Xuc3],0);
[XCc4,fC4] = sim_motion10(Xsc4,Xuc4,connc4,.005,1000,[Xsc4,Xuc4],0);
% Uncomment to run simulation in real-time
% [XCc5,fC5] = sim_motion10(Xsc5,Xuc5,connc5,.01,12500,[Xsc5,Xuc5],0);
load chaos_network.mat;
disp(['mean simulation error: ' [num2str(mean(fC2)) '  ',...
       num2str(mean(fC3)) ' ' num2str(mean(fC4)) ' ' num2str(mean(fC5))]]);
% Distances
D1c = sqrt(squeeze(sum(diff(XCc1(:,1:(size(Xsc1,2)),:),1,2).^2)));
D2c = sqrt(squeeze(sum(diff(XCc2(:,1:(size(Xsc2,2)),:),1,2).^2)));
D3c = sqrt(squeeze(sum(diff(XCc3(:,1:(size(Xsc3,2)),:),1,2).^2)));
D4c = sqrt(squeeze(sum(diff(XCc4(:,1:(size(Xsc4,2)),:),1,2).^2)));
D5c = sqrt(squeeze(sum(diff(XCc5(:,1:(size(Xsc5,2)),:),1,2).^2)));
% Find minima
nSt = 20;
[mV2,mI2] = min(abs(D2c(1,nSt:end)-D2c(end,nSt:end))); mI2 = mI2 + nSt-1;
[mV3,mI3] = min(abs(D3c(1,nSt:end)-D3c(end,nSt:end))); mI3 = mI3 + nSt-1;
[mV4,mI4] = min(abs(D4c(1,nSt:end)-D4c(end,nSt:end))); mI4 = mI4 + nSt-1;
disp([mV2 mV3 mV4]);
disp([mI2 mI3 mI4]);


%% c: Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;

cVP = -(linspace(.26,.475,length(cVp))).^(1);
DLin = linspace(.4,1.5,size(CP1,1));
slC = @(xP) 1./(8.*xP.^2 - 1);

lcDI = [0 find(diff(lcp)>0)];
% Compute Feigenbaum constant
fc = (slC(cV(lcDI(3:end-1))*nValpd)-slC(cV(lcDI(2:end-2))*nValpd)) ./...
     (slC(cV(lcDI(4:end))*nValpd)-slC(cV(lcDI(3:end-1))*nValpd));


hold on;
nF = [1 1.85 1.5 1.57*ones(1,length(lcDI))];
% nF = ones(1,length(lcDI));
PM0 = LCp{lcDI(2)};
mLS = ['k^';'ks';'ko'];
for i = 1:length(lcDI)-1
% for i = 1:1
    plI = (lcDI(i)+1):lcDI(i+1);
    PM = [LCp{plI}];
    if(i==1)
        cVPl = ones(size(PM,1),1).*cVP(plI);
        LCpp = [LCp{plI}];
        scatter(PM(:),cVPl(:),1,interp1(DLin,CP1,LCpp(:)),'filled','clipping',0);
        plot(LCp{pLC},[1;1]*cVP(pLC),mLS(1,:),'linewidth',1,'markersize',4);
    else
        PM0 = reshape([1;1].*PM0',[2*lcp(lcDI(i)-1),1]);
        PM0r = reshape([1;1].*LCp{lcDI(i)-1}',[2*lcp(lcDI(i)-1),1]);
        PM = (PM-PM0r)*nF(i)^(i-1) + PM0;
        cVPl = ones(size(PM,1),1).*cVP(plI);
        LCpp = [LCp{plI}];
        if(i < 5)
            for j = 1:(size(LCpp,1)/2)
                line_coordinates([PM((1:2)+2*(j-1),end)';[1 1]*cVP(plI(i-1))],...
                                 'lSh',0,'style','-','lw',.5,'nw',.005,'color',[1 1 1]*.8);
            end
        end
        scatter(PM(:),cVPl(:),1,interp1(DLin,CP1,LCpp(:)),'filled','clipping',0);
        if(i>2 && i<5)
            plot(PM0(1:2:end),ones(2^(i-1),1)*cVP(lcDI(i)),mLS(i-1,:),...
                 'linewidth',1,'clipping',0,'markersize',4);
        end
        PM0 = PM(:,end);
    end
end

% Plot networks
sc = .05; 
sh1 = [0.47;-.10];
sh2 = [0.71;-.11];
sh3 = [0.97;-.11];
sh4 = [1.26;-.11];
sh =  [.540;-.15];
Rz2 = rotz(-18);
Rz3 = rotz(-11);
Rz4 = rotz(-07);
XCc1p = XCc1(:,:)*sc+sh1;
XCc2p = Rz2(1:2,1:2)*XCc2(:,:,mI2)*sc+sh2;
XCc3p = Rz3(1:2,1:2)*XCc3(:,:,mI3)*sc+sh3;
XCc4p = Rz4(1:2,1:2)*XCc4(:,:,mI4)*sc+sh4;
XCc1pp = XCc1(:,:)*sc*.75+sh1+sh;
nL2 = [1 0 -1];
nL3 = [1 0 -1 -1 -1];
nL4 = [1 0 1 1 1 1 -1 0 -1];
% Distance measurements
for i = 1:size(Xsc1,2)-1
    line_coordinates(XCc1p(:,[0 1]+i),'lSh',0.004,'color',...
                     interp1(DLin,CP1,D1c(i)));
end
for i = 1:size(Xsc2,2)-1
    line_coordinates(XCc2p(:,[0 1]+i),'lSh',.004*nL2(i),'color',...
                     interp1(DLin,CP1,D2c(i,mI2)));
end
for i = 1:size(Xsc3,2)-1
    line_coordinates(XCc3p(:,[0 1]+i),'lSh',.004*nL3(i),'color',...
                     interp1(DLin,CP1,D3c(i,mI3)));
end
for i = 1:size(Xsc4,2)-1
    line_coordinates(XCc4p(:,[0 1]+i),'lSh',.004*nL4(i),'color',...
                     interp1(DLin,CP1,D4c(i,mI4)));
end
for i = 1:size(Xsc1,2)-1
    line_coordinates(XCc1pp(:,[0 1]+i),'lSh',0.004,'color',...
                     interp1(DLin,CP1,D1c(i)));
end
% Networks
visualize_network(XCc1p(:,1:size(Xsc1,2),end),...
                  XCc1p(:,(1:size(Xuc1,2))+size(Xsc1,2),end),connc1,...
                  'ucolor',interp1(SLin,CP3,-1));
visualize_network(XCc2p(:,1:size(Xsc2,2)),...
                  XCc2p(:,(1:size(Xuc2,2))+size(Xsc2,2)),connc2,...
                  'ucolor',interp1(SLin,CP3,slC(cVp(pLC)*nValpd)));
visualize_network(XCc3p(:,1:size(Xsc3,2)),...
                  XCc3p(:,(1:size(Xuc3,2))+size(Xsc3,2)),connc3,...
                  'ucolor',interp1(SLin,CP3,slC(cVp(lcDI(4))*nValpd)));
visualize_network(XCc4p(:,1:size(Xsc4,2)),...
                  XCc4p(:,(1:size(Xuc4,2))+size(Xsc4,2)),connc4,...
                  'ucolor',interp1(SLin,CP3,slC(cVp(lcDI(5))*nValpd)));
visualize_network(XCc1pp(:,1:size(Xsc1,2),end),...
                  XCc1pp(:,(1:size(Xuc1,2))+size(Xsc1,2),end),connc1,...
                  'ucolor',interp1(SLin,CP3,-1));
% Symbols
plot(.66,sh2(2)+.052,mLS(1,:),'linewidth',1);
plot(.95,sh3(2)+.056,mLS(2,:),'linewidth',1);
plot(1.33,sh4(2)+.056,mLS(3,:),'linewidth',1);
% Line
arrow([1;1]*.957,cVP([1 end])+.01,sRat(pInd),'linewidth',.5);
line([-1;1]*.005+ones(1,length(lcDI)-1)*.957,[0;0]+cVP(lcDI(1:end-1)+1),...
     'color','k','linewidth',.5);
for i = 1:length(lcDI)-1
    text(.922,cVP(lcDI(i)+1),['$s_' num2str(i-1) '$'],'fontsize',FS);
end

% Text
textc = '\textbf{c} \hspace{42mm}period doubling route to chaos';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});
textc2 = ['$2^n$--cycles lose stability at $s_n$'];
text(.47,.52,textc2,NVTextH{:});
text(.000,.88,'fixed point',NVTextR{:});
text(.250,.88,'2-cycle',NVTextR{:});
text(.500,.88,'4-cycle',NVTextR{:});
text(.830,.88,'8-cycle',NVTextR{:});
text(.23,.78,'$l_1$',NVTextR{:});
text(.31,.67,'$l_3$',NVTextR{:});
text(.465,.78,'$l_1$',NVTextR{:});
text(.615,.67,'$l_5$',NVTextR{:});
text(.725,.78,'$l_1$',NVTextR{:});
text(.98,.67,'$l_9$',NVTextR{:});
text(.127,.37,['zoom $\times$' num2str(nF(2))],NVTextR{:},'color',[1 1 1]*.7);
text(.332,.305,['$\times$' num2str(nF(3))],NVTextR{:},'color',[1 1 1]*.7);
text(.067,.251,['$\times$' num2str(nF(4))],NVTextR{:},'color',[1 1 1]*.7);


hold off;
axis([0 sRat(pInd) -1 0]*abs(min(cVP)) + [1 1 0 0]*.4);
set(gca,'visible',0);


%% d: Chaos
pInd = 4;
subplot('position',subpN(pInd,:)); cla;

sc = .19;
scd = .22;
sh = [.05;.66];
shd = [0;-.04];
pld = 10850;
XCc5p = XCc5(:,:,pld)*sc+sh;
D5cp = D5c(:,pld);
for i = 1:size(Xsc5,2)-1
    line_coordinates(XCc5p(:,(0:1)+i),'lSh',-.007*(-1)^i,'color',interp1(DLin,CP1,D5cp(i)));
end
visualize_network(XCc5p(:,1:size(Xsc5,2)),...
                  XCc5p(:,(1:size(Xuc5,2))+size(Xsc5,2)),connc5,...
                  'ucolor',interp1(SLin,CP3,slC(cVp(lcDI(5))*nValpd)));

% Axes
XDa = movmean(XCc5p(1,1:size(Xsc5,2)),[0,1],'Endpoints','Discard');
axV = [0 2.85 0 .4];
text(XDa(1)/sRat(pInd),-.08,'$1$',NVTextH{:});
text(XDa(2)/sRat(pInd),-.08,'$2$',NVTextH{:});
text(XDa(3)/sRat(pInd),-.08,'$3$',NVTextH{:});
text(XDa(5)/sRat(pInd),-.08,'$\cdots$',NVTextH{:});
text(XDa(end)/sRat(pInd),-.08,num2str(length(D5cp)),NVTextH{:});
text(axV(2)/sRat(pInd)+.005,-.0,'$k$',NVTextR{:});
text(axV(1)/sRat(pInd)-.003,1*scd+shd(2),'$1$',NVTextRA{:},'color',interp1(DLin,CP1,1));
text(axV(1)/sRat(pInd)-.003,.4*scd+shd(2),'$0.4$',NVTextRA{:},'color',interp1(DLin,CP1,.4));
text(axV(1)/sRat(pInd)-.003,1.5*scd+shd(2),'$1.5$',NVTextRA{:},'color',interp1(DLin,CP1,1.5));
text(axV(1)/sRat(pInd),axV(4)+.06,'$l_k$',NVTextH{:});

% Distances
scatter(XDa,D5cp*scd+shd(2),20,interp1(DLin,CP1,D5cp),'filled','marker','s');
line([1;1].*XDa, [0;.03].*ones(1,size(XDa,2)) + axV(3),'color','k','linewidth',.5);
line([0;.03]+axV(1),[1;1]*1*scd+shd(2),'color',interp1(DLin,CP1,1));
line([0;.03]+axV(1),[1;1]*min(DLin)*scd+shd(2),'color',interp1(DLin,CP1,min(DLin)));
line([0;.03]+axV(1),[1;1]*max(DLin)*scd+shd(2),'color',interp1(DLin,CP1,max(DLin)));
line(axV([1 1 2]),axV([4 3 3]),'linewidth',.5,'color','k');

% e: Feigenbaum constant
fcInf = 4.669201609;
axV2 = [0 .98 0 .7] + [3.24 3.24 0 0];
scf = .22;
shf = [3.16;-.95];
plx = (1:length(fc))*.19+shf(1);
line(axV2(1:2)',[1;1]*fcInf*scf+shf(2),'color',o*gr,...
     'linewidth',.5,'linestyle','--');
text(.79,.15,'$\delta_{\infty}$',NVTextH{:},'color',o*gr);
scatter(plx,fc*scf+shf(2),20,'k','filled','clipping',0);
% Ticks
line([1;1].*plx, [0;.03].*ones(1,length(plx)) + axV2(3),'color','k','linewidth',.5);
line([0;.03].*[1 1 1]+axV2(1),[1;1].*[5 6 7]*scf+shf(2),'color','k');
line(axV2([1 1 2]),axV2([4 3 3]),'linewidth',.5,'color','k');
% Text
text(axV2(1)/sRat(pInd)-.015,5*scf+shf(2),'$5$',NVTextH{:});
text(axV2(1)/sRat(pInd)-.015,6*scf+shf(2),'$6$',NVTextH{:});
text(axV2(1)/sRat(pInd)-.015,7*scf+shf(2),'$7$',NVTextH{:});
text(axV2(1)/sRat(pInd)-.01,axV2(4)+.06,...
     '$\delta_n = \frac{s_{n-1}-s_{n-2}}{s_n - s_{n-1}}$',NVTextR{:});
text(axV2(2)/sRat(pInd)+.005,0,'$n$',NVTextR{:});
text(plx(1)/sRat(pInd),-.08,'$2$',NVTextH{:});
text(plx(2)/sRat(pInd),-.08,'$3$',NVTextH{:});
text(plx(3)/sRat(pInd),-.08,'$4$',NVTextH{:});
text(plx(4)/sRat(pInd),-.08,'$5$',NVTextH{:});
text(plx(5)/sRat(pInd),-.08,'$6$',NVTextH{:});


% Text
textd = '\textbf{d}\hspace{3mm}chaotic modules yield non-repeating lattice';
texte = '\textbf{e}\hspace{2mm}Feigenbaum constant';
text(labX-.2,subp(pInd,4)+labY,textd,NVTitle{:});
text(labX+10.1,subp(pInd,4)+labY,texte,NVTitle{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);



%% Save
fName = 'figure7a';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');