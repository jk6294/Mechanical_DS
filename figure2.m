% Figure 2: Iterated Maps and 1D Discrete Dynamical Systems
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');

% Figure Colors
cTr1 = [126 200 255]/255;
cTr2 = [115 140 200]/255;
cTr3 = [015 082 186]/255;
cSr  = [255 100 100]/255;
% Plot Parameters
lw_d = .5;
pSc = 0.6;
pSc2 = 0.5;
lSh = .5;
nW = .1;
% Label Parameters
labX = -.1;
labY = 1.0;
labRowX = -0.35;
labColY = 1.3;

% Subplot Indices
nRow = [6 10 12];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [7 7 7 7];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-2), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(2); clf;
set(gcf, 'Renderer', 'painters', 'Position', [10 10 24.5 13.1], 'Units', 'centimeters'); 
lwL = 1;
arrHL = 5;
arrHW = 5;

annotation('line','linewidth',.5,'position',[.298 .15 0 .8],'color',[.9 .9 .9]);
annotation('line','linewidth',.5,'position',[.493 .15 0 .8],'color',[.9 .9 .9]);
annotation('line','linewidth',.5,'position',[.689 .15 0 .8],'color',[.9 .9 .9]);


%% a: Two Modules
subplot(NRow,NCol,cellM{1}); cla;
tSx = -.05;
tSy = -.15;
Xs1 = [-1 -1  1  1;...
       -1  1 -1  1];
conn1 = [1 3; 1 4; 2 3; 2 4];

% Individual modules
visualize_network(Xs1-[4;0],[],conn1,pSc);
visualize_network(Xs1-[0;0],[],conn1,pSc);
text(tSx+.04,tSy+.2, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.20,tSy+.2, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.33,tSy+.2, '$\mathrm{d_2^\prime}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.48,tSy+.2, '$\mathrm{d_3}$', 'Units', 'Normalized', 'fontsize',10);
hold on;
line_coordinates(Xs1(:,1:2)-[4;0], lSh, nW, lw_d);
line_coordinates(Xs1(:,3:4)-[4;0], -lSh, nW, lw_d);
line_coordinates(Xs1(:,1:2)-[0;0], lSh, nW, lw_d);
line_coordinates(Xs1(:,3:4)-[0;0], -lSh, nW, lw_d);
hold off;

% Combined network
[Xs1a,conn1a] = tesselate_network(Xs1,conn1,[2;0],[2;1]);
visualize_network(Xs1a+[5;0],[],conn1a,pSc);
text(tSx+.68,tSy+.2, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.82,tSy+.2, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+0.96,tSy+.2, '$\mathrm{d_3}$', 'Units', 'Normalized', 'fontsize',10);

% Labels
text(.09, 1, 'Modules', 'Units', 'Normalized', 'fontsize', 10);
text(.65, 1, 'Combine', 'Units', 'Normalized', 'fontsize', 10);

axis([0 14.0 0 5.2] + [-5.6 -5.6 -2.7 -2.7]);
text(labX,labY,'\textbf{a}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.2,labColY,'\textbf{4-Bar Linkage}','Units','Normalized','fontsize',10,'fontweight','bold');
drawnow;


%% b: Coupled Modules
subplot(NRow,NCol,cellM{2}); cla;
% Simulate Combined Network for Propagation
[XC1a,~] = sim_motion(Xs1a,[],conn1a,.001,1200,Xs1a,0);
D1 = sqrt(squeeze(sum(diff(XC1a,1,2).^2)));
D1 = D1(1:2:end,:); 

% Indices to plot
pInd1 = 1184;
pInd2 = 592;
pInd3 = 1;
% Networks
visualize_network(XC1a(:,:,pInd1)+[0.0;-4.5],[],conn1a,pSc2,cTr1,cTr1);
visualize_network(XC1a(:,:,pInd2)+[0.0;0.0],[],conn1a,pSc2,cTr2,cTr2);
visualize_network(XC1a(:,:,pInd3)+[0.0;4.5],[],conn1a,pSc2,cTr3,cTr3);
hold on;
line_coordinates(XC1a(:,1:2,pInd3)+[0;+4.5], lSh, nW, lw_d);
line_coordinates(XC1a(:,5:6,pInd3)+[0;+4.5], -lSh, nW, lw_d);
line_coordinates(XC1a(:,1:2,pInd2)+[0;+0.0], lSh, nW, lw_d);
line_coordinates(XC1a(:,5:6,pInd2)+[0;+0.0], -lSh, nW, lw_d);
line_coordinates(XC1a(:,1:2,pInd1)+[0;-4.5], lSh, nW, lw_d);
line_coordinates(XC1a(:,5:6,pInd1)+[0;-4.5], -lSh, nW, lw_d);
hold off;
axis([0 17.5 0 14] + [-7 -7 -6.8 -6.8])
% Text
tSx = -8;
text(tSx, 4.5, '$\mathrm{d_1=2}$', 'fontsize', 10);
text(-tSx-3, 4.5, '$\mathrm{d_3=2}$', 'fontsize', 10);
text(tSx, 0, '$\mathrm{d_1=2.5}$', 'fontsize', 10);
text(-tSx-3, 0, '$\mathrm{d_3=2.5}$', 'fontsize', 10);
text(tSx, -4.5, '$\mathrm{d_1=3}$', 'fontsize', 10);
text(-tSx-3, -4.5, '$\mathrm{d_3=3}$', 'fontsize', 10);
text(labX,labY,'\textbf{b}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.32, 1, 'Motion', 'Units', 'Normalized', 'fontsize', 10);
drawnow;


%% c: Cobweb Plot
subplot(NRow,NCol,cellM{3}); cla;

% Simulate Single Module for Cobweb
dX1 = [0 0 0 0; -1 1 1 -1];
[XCa,~] = sim_motion(Xs1,[],conn1,.01,170,dX1,0);   % Simulate
[XCb,~] = sim_motion(Xs1,[],conn1,.01,170,-dX1,0);  % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,3:4,:),[],2)).^2)));

plot(d1,d2,'k-','linewidth',1);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
caF = .8;
% Cobweb 1
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr1,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
plot(dP(1),1,'ko','markersize',2.5,'linewidth',2.5);
plot(dP(2),dP(3),'ko','markersize',2.5,'linewidth',2.5);
plot(dP(4),dP(5),'ko','markersize',2.5,'linewidth',2.5);
% Cobweb 2
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr2,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 3
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr3,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
plot(D1(1,pInd3),D1(2,pInd3),'k*','linewidth',1,'markersize',6.5);
plot(d1,d2,'k-','linewidth',1);
hold off;
% Formatting
set(gca,'visible',1,'XTick',[2 2.5 3],'YTick',[2 2.5 3],'XTickLabel',[],'YTickLabel',[]);
axis([1 3.4 1 3.4]);

% Text
text(.6,.5,'$\mathrm{d_{k+1}=d_k}$','Units','normalized','fontsize',10,'color',[200 200 200]/255);
text(.3,.7,'$\mathrm{d_{k+1}=f(d_k)}$','Units','normalized','fontsize',10);
text(.47,-.1,'$\mathrm{d_k}$','Units','normalized','fontsize',10);
text(-.1,.4,'$\mathrm{d_{k+1}}$','Units','normalized','rotation',90,'fontsize',10);
text(.85,.05,'$\mathrm{d_1}$','Units','normalized','fontsize',10);
text(.71,.27,'$\mathrm{(d_1,d_2)}$','Units','normalized','fontsize',10);
text(.17,.88,'$\mathrm{(d_2,d_3)}$','Units','normalized','fontsize',10);
text(.38,-.06,'2','Units','normalized','fontsize',10);
text(.58,-.06,'2.5','Units','normalized','fontsize',10);
text(.86,-.06,'3','Units','normalized','fontsize',10);
text(labX,labY*1.1,'\textbf{c}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.31,labY*1.1,'Cobweb','Units','Normalized','fontsize',10);
drawnow;


%% d: 2 FP + Super Stability
subplot(NRow,NCol,cellM{4}); cla;
tSx = -.03;
tSy = -.15;

s = sqrt(3);
Xs2 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
Xu2 = [-0.86 -0.86;...
       -1.45  1.47];
Xs2p = [Xs2(1,:); -Xs2(2,:)];
Xu2p = [Xu2(1,:); -Xu2(2,:)];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Modules
visualize_network(Xs2,Xu2,conn2,pSc);
visualize_network(Xs2p+[3.5;.5],Xu2p+[3.5;.5],conn2,pSc);
text(tSx-.05,tSy+.8, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.22,tSy+.8, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.35,tSy+.04, '$\mathrm{d_2^\prime}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.5,tSy+.04, '$\mathrm{d_3}$', 'Units', 'Normalized', 'fontsize',10);
hold on;
line_coordinates(Xs2(:,1:2)+[0;0], lSh, nW, lw_d);
line_coordinates(Xs2(:,2:3)+[0;0], lSh, nW, lw_d);
line_coordinates(Xs2p(:,1:2)+[3.5;.5], -lSh, nW, lw_d);
line_coordinates(Xs2p(:,2:3)+[3.5;.5], -lSh, nW, lw_d);
hold off;

% Combined
Xs2c = [-s/2  0  s/2  s;...
        -1/2  1 -1/2  1];
Xu2c = [Xu2 Xu2p+[s/2;.5]];
conn2c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[1;1]);
visualize_network(Xs2a+[8;0],Xu2a+[8;0],conn2a,pSc);

% Labels
text(.08, 1, 'Modules', 'Units', 'Normalized', 'fontsize', 10);
text(.66, 1, 'Combine', 'Units', 'Normalized', 'fontsize', 10);

axis([0 11.5 0 4.5]+[-1.4 -1.4 -1.5 -1.5]);
text(labX,labY,'\textbf{d}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.101,labColY,'\textbf{2 Fixed Points}','Units','Normalized','fontsize',10,'fontweight','bold');
drawnow;


%% e: Tesselate
subplot(NRow,NCol,cellM{5}); cla;
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[4;1]);
% Simulate Combined Network for Propagation
[XC2a,fC] = sim_motion(Xs2a,Xu2a,conn2a,.02,570,[Xs2a Xu2a],0);

pInd3 = 1;
pInd2 = 250;
pInd1 = 570;

% Networks
XC2a(1,:,:) = XC2a(1,:,:) - XC2a(1,1,:);
cC1 = [cSr;cSr;repmat(cTr1,[6,1]);cSr;cSr];
cC2 = [cSr;cSr;repmat(cTr2,[6,1]);cSr;cSr];
cC3 = [cSr;cSr;repmat(cTr3,[6,1]);cSr;cSr];
visualize_network(XC2a(:,1:10,pInd1)+[0;-4.5],...
                  XC2a(:,11:end,pInd1)+[0;-4.5],conn2a,pSc2,cC1,cTr1);
visualize_network(XC2a(:,1:10,pInd2)+[0;0.0],...
                  XC2a(:,11:end,pInd2)+[0;0.0],conn2a,pSc2,cC2,cTr2);
visualize_network(XC2a(:,1:10,pInd3)+[0;4.5],...
                  XC2a(:,11:end,pInd3)+[0;4.5],conn2a,pSc2,cC3,cTr3);
hold on;
line_coordinates(XC2a(:,1:2,pInd3)+[0;+4.5], lSh, nW, lw_d);
line_coordinates(XC2a(:,9:10,pInd3)+[0;+4.5], -lSh, nW, lw_d);
line_coordinates(XC2a(:,1:2,pInd2)+[0;0], lSh, nW, lw_d);
line_coordinates(XC2a(:,9:10,pInd2)+[0;0], -lSh, nW, lw_d);
line_coordinates(XC2a(:,1:2,pInd1)+[0;-4.5], lSh, nW, lw_d);
line_coordinates(XC2a(:,9:10,pInd1)+[0;-4.5], -lSh, nW, lw_d);
hold off;

axis([0 17.5 0 14] + [-2.3 -2.3 -5.8 -5.8])
% Text
tSx = -3;
text(tSx,4.9,'$\mathrm{D_1^*}$','fontsize',10);
text(-tSx+5.4,4.1,'$\mathrm{D_1^*}$','fontsize',10);
text(tSx-1.7,0.5,'$\mathrm{\rightarrow D_2^*}$','fontsize',10);
text(-tSx+8.3,-.4,0,'$\mathrm{\rightarrow D_1^*}$','fontsize',10);
text(tSx,-4.1,'$\mathrm{D_2^*}$','fontsize',10);
text(-tSx+10.5,-4.9,'$\mathrm{D_2^*}$','fontsize',10);
text(labX,labY,'\textbf{e}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.32, 1, 'Motion', 'Units', 'Normalized', 'fontsize', 10);
drawnow;


%% f: Simulate
subplot(NRow,NCol,cellM{6}); cla;% Simulate Single Module for Cobweb
[XCa,~] = sim_motion(Xs2,Xu2,conn2,.01,185,[Xs2 Xu2],0);   % Simulate
[XCb,~] = sim_motion(Xs2,Xu2,conn2,.01,15,-[Xs2 Xu2],0);   % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));

% Combined Network Distances
D1 = sqrt(squeeze(sum(diff(XC2a,1,2).^2)));
D1 = D1(1:size(Xs2a,2)-2,:);
plot(d1,d2,'k-','linewidth',1);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
% Cobweb 3
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr3,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 2
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr2,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 1
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr1,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF*.8 diff(dPb(i:i+1))*caF*.8]);
end
plot(D1(1,pInd3),D1(2,pInd3),'k*','linewidth',1,'markersize',6.5);
plot(D1(1,pInd1),D1(2,pInd1),'k*','linewidth',1,'markersize',6.5);
plot(d1,d2,'k-','linewidth',1);
hold off;

% Formatting
set(gca,'visible',1,'XTick',[1.75 2.95],'YTick',[1.75 2.95],'XTickLabel',[],'YTickLabel',[]);
axis([1.6 3.1 1.6 3.1]);

% Text
text(.47,-.1,'$\mathrm{d_k}$','Units','normalized','fontsize',10);
text(-.1,.4,'$\mathrm{d_{k+1}}$','Units','normalized','rotation',90,'fontsize',10);
text(.05,.25,'$\mathrm{D_1^*}$','Units','normalized','fontsize',10);
text(.70,.92,'$\mathrm{D_2^*}$','Units','normalized','fontsize',10);
text(labX,labY*1.1,'\textbf{f}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.32,labY*1.1,'Cobweb','Units','Normalized','fontsize',10);
drawnow;


%% g: Limit Cycle
subplot(NRow,NCol,cellM{7}); cla;
Xs3 = [-s/2  0.0  s/2;...
        -0.5  1.0 -0.5];
Xu3 = [ 0.10 -0.30;...
       -0.25 -0.90];
Xs3p = [Xs3(1,:); -Xs3(2,:)];
Xu3p = [Xu3(1,:); -Xu3(2,:)];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

xSh = .4;
ySh = xSh/sqrt(3);
tSx = -.03;
tSy = -.15;

visualize_network(Xs3,Xu3,conn3,pSc);
visualize_network(Xs3p+[3.5;.5],Xu3p+[3.5;.5],conn3,pSc);
text(tSx-.01,tSy+.8, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.22,tSy+.8, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.35,tSy+.1, '$\mathrm{d_2^\prime}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.5,tSy+.1, '$\mathrm{d_3}$', 'Units', 'Normalized', 'fontsize',10);

% Distances
hold on;
line_coordinates(Xs3(:,1:2),lSh,nW,lw_d);
line_coordinates(Xs3(:,2:3),lSh,nW,lw_d);
line_coordinates(Xs3p(:,1:2)+[3.5;.5],-lSh,nW,lw_d);
line_coordinates(Xs3p(:,2:3)+[3.5;.5],-lSh,nW,lw_d);
hold off;

% Combined
Xs3c = [-s/2  0  s/2  s;...
        -1/2  1 -1/2  1];
Xu3c = [Xu3 Xu3p+[s/2;.5]];
conn3c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[1;1]);
visualize_network(Xs3a+[8;0],Xu3a+[8;0],conn3a,pSc);

% Labels
text(.08, 1, 'Modules', 'Units', 'Normalized', 'fontsize', 10);
text(.66, 1, 'Combine', 'Units', 'Normalized', 'fontsize', 10);

text(labX,labY,'\textbf{g}','Units','Normalized','fontsize',10);
text(-.02,labColY,'\textbf{Isolated Limit Cycle}','Units','Normalized','fontsize',10,'fontweight','bold');
axis([0 11.5 0 4.5]+[-1.4 -1.4 -1.5 -1.5]);
drawnow;


%% h: Tesselate Limit Cycle
subplot(NRow,NCol,cellM{8}); cla;
Xs3c = [-s/2  0  s/2  s;...
        -1/2  1 -1/2  1];
Xu3c = [Xu3 Xu3p+[s/2;.5]];
conn3c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[8;1]);
% Simulate Combined Network for Propagation
[XC3a,fC] = sim_motion(Xs3a,Xu3a,conn3a,.02,700,-[Xs3a Xu3a],0);

% Networks
pInd3 = 1;
pInd2 = 280;
pInd1 = 603;
XC3a(1,:,:) = XC3a(1,:,:) - XC3a(1,2,:);
cC1 = [cSr;cSr;cSr;repmat(cTr1,[12,1]);cSr;cSr;cSr];
cC2 = [cSr;cSr;cSr;repmat(cTr2,[12,1]);cSr;cSr;cSr];
cC3 = [cSr;cSr;cSr;repmat(cTr3,[12,1]);cSr;cSr;cSr];
visualize_network(XC3a(:,1:18,pInd1)+[0;-3.4],...
                  XC3a(:,19:end,pInd1)+[0;-3.4],conn3a,pSc2,cC1,cTr1);
visualize_network(XC3a(:,1:18,pInd2)+[0;0.5],...
                  XC3a(:,19:end,pInd2)+[0;0.5],conn3a,pSc2,cC2,cTr2);
visualize_network(XC3a(:,1:18,pInd3)+[0;4.5],...
                  XC3a(:,19:end,pInd3)+[0;4.5],conn3a,pSc2,cC3,cTr3);
hold on;
line_coordinates(XC3a(:,1:2,pInd3)+[0;+4.5], lSh, nW, lw_d);
line_coordinates(XC3a(:,17:18,pInd3)+[0;+4.5], -lSh, nW, lw_d);
line_coordinates(XC3a(:,1:2,pInd2)+[0;.5], lSh, nW, lw_d);
line_coordinates(XC3a(:,17:18,pInd2)+[0;.5], -lSh, nW, lw_d);
line_coordinates(XC3a(:,1:2,pInd1)+[0;-3.5], lSh, nW, lw_d);
line_coordinates(XC3a(:,17:18,pInd1)+[0;-3.5], -lSh, nW, lw_d);
hold off;

% Text
tSx = -2.8;
text(tSx-.2,4.9,'$\mathrm{D^*}$','fontsize',10);
text(-tSx+11.5,4.1,'$\mathrm{D^*}$','fontsize',10);
text(tSx-1.9,0,'$\mathrm{\rightarrow D_1^o}$','fontsize',10);
text(-tSx+9.5,-2.2,'$\mathrm{\rightarrow D_1^o}$','fontsize',10);
text(tSx+1.5,-6.5,'$\mathrm{D_1^o}$','fontsize',10);
text(-tSx+9.2,-6.5,'$\mathrm{D_1^o}$','fontsize',10);
text(labX,labY,'\textbf{h}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.32, 1, 'Motion', 'Units', 'Normalized', 'fontsize', 10);

axis([0 17.5 0 14] + [-2 -2 -6.2 -6.2])
drawnow;


%% i: Simulate
subplot(NRow,NCol,cellM{9}); cla;% Simulate Single Module for Cobweb
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[7;1]);
[XCa,~] = sim_motion(Xs3,Xu3,conn3,.01,90,[Xs3 Xu3],0);   % Simulate
[XCb,~] = sim_motion(Xs3,Xu3,conn3,.01,100,-[Xs3 Xu3],0);   % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
D1 = sqrt(squeeze(sum(diff(XC3a,1,2).^2)));
D1 = D1(1:size(Xs3a,2)-2,:);

plot(d1,d2,'k-','linewidth',1);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
% Cobweb 1
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr1,'linewidth',lwL);
dLC1 = dP(1);                   % First Limit Cycle
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF*.8 diff(dPb(i:i+1))*caF*.8]);
end
% Cobweb 2
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr2,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 3
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr3,'linewidth',lwL);
dFP1 = dP(1);                   % First Fixed Point
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
plot(D1(1,pInd3),D1(2,pInd3),'k*','linewidth',1,'markersize',6.5);
plot(D1(1,pInd1),D1(2,pInd1),'ko','linewidth',1,'markersize',6.5);
plot(D1(2,pInd1),D1(3,pInd1),'ko','linewidth',1,'markersize',6.5);
plot(d1,d2,'k-','linewidth',1);
hold off;

% Formatting
set(gca,'visible',1,'XTick',[1.192 2.088],'YTick',[1.192 2.088],...
                    'XTickLabel',[],'YTicklabel',[],'fontsize',10);
axis([1.08 2.2 1.08 2.2]);
text(.47,-.1,'$\mathrm{d_k}$','Units','normalized','fontsize',10);
text(-.1,.4,'$\mathrm{d_{k+1}}$','Units','normalized','rotation',90,'fontsize',10);
text(.55,.93,'$\mathrm{D^*}$','Units','normalized','fontsize',10);
text(.025,.93,'$\mathrm{D_1^o}$','Units','normalized','fontsize',10);
text(.89,.24,'$\mathrm{D_2^o}$','Units','normalized','fontsize',10);
text(labX,labY*1.1,'\textbf{i}','Units','Normalized','fontsize',10);
text(.32,labY*1.1,'Cobweb','Units','Normalized','fontsize',10);
drawnow;


%% j: Chaos
subplot(NRow,NCol,cellM{10}); cla;

xSh = .4;
ySh = 2*xSh/sqrt(3);
tSx = .02;
tSy = 0;

Xs4 = [-s/2  0.0  s/2;...
        -0.5  0.0 -0.5]*1.5;
Xu4 = [ 0.15 -0.22;...
       -0.50 -0.90]*1.5;
Xs4p = [Xs4(1,:); -Xs4(2,:)];
Xu4p = [Xu4(1,:); -Xu4(2,:)];
conn4 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_network(Xs4,Xu4,conn4,pSc);
visualize_network(Xs4p+[3.5;.5],Xu4p+[3.5;.5],conn4,pSc);
text(tSx+.02, 0.6+tSy, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(tSx+.16, 0.6+tSy, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(.34+tSx, .1+tSy, '$\mathrm{d_2^\prime}$', 'Units', 'Normalized', 'fontsize',10);
text(.47+tSx, .1+tSy, '$\mathrm{d_3}$', 'Units', 'Normalized', 'fontsize',10);


% Distances
hold on;
line_coordinates(Xs4(:,1:2),lSh,nW,lw_d);
line_coordinates(Xs4(:,2:3),lSh,nW,lw_d);
line_coordinates(Xs4p(:,1:2)+[3.5;.5],-lSh,nW,lw_d);
line_coordinates(Xs4p(:,2:3)+[3.5;.5],-lSh,nW,lw_d);
hold off;

% Combined
Xs4c = [-s/2  0  s/2  s;...
        -1/2  0 -1/2  0]*1.5;
Xu4c = [Xu4 Xu4p+[s;-1]*1.5/2];
conn4c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs4a,Xu4a,conn4a] = tesselate_network_old(Xs4c,Xu4c,conn4c,[2*s;0],[1;1]);
visualize_network(Xs4a+[7.1;0],Xu4a+[7.1;0],conn4a,pSc);


% Labels
text(.08, 1, 'Modules', 'Units', 'Normalized', 'fontsize', 10);
text(.66, 1, 'Combine', 'Units', 'Normalized', 'fontsize', 10);

text(labX,labY,'\textbf{j}','Units','Normalized','fontsize',10);
text(.32,labColY,'\textbf{Chaos}','Units','Normalized','fontsize',10,'fontweight','bold');
axis([0 11.5 0 4.5]+[-1.8 -1.8 -1.4 -1.4]);
drawnow;


%% k: Tesselate Chaos
subplot(NRow,NCol,cellM{11}); cla;

% Simulate Combined Network for Propagation
[Xs4a,Xu4a,conn4a] = tesselate_network_old(Xs4c,Xu4c,conn4c,[1.5*s;0],[4;1]);
[XC4a,fC] = sim_motion(Xs4a,Xu4a,conn4a,.002,8500,-[Xs4a Xu4a],0);

% Networks
cC1 = [cSr;cSr;cSr;repmat(cTr1,[4,1]);cSr;cSr;cSr];
cC2 = [cSr;cSr;cSr;repmat(cTr2,[4,1]);cSr;cSr;cSr];
cC3 = [cSr;cSr;cSr;repmat(cTr3,[4,1]);cSr;cSr;cSr];
XC4a(1,:,:) = XC4a(1,:,:) - XC4a(1,2,:);
pInd3 = 1;
pInd2 = 3750;
pInd1 = 8445;
visualize_network(XC4a(:,1:10,pInd1)+[0;-3.5],...
                  XC4a(:,11:end,pInd1)+[0;-3.5],conn4a,pSc2,cC1,cTr1);
visualize_network(XC4a(:,1:10,pInd2)+[0;.5],...
                  XC4a(:,11:end,pInd2)+[0;.5],conn4a,pSc2,cC2,cTr2);
visualize_network(XC4a(:,1:10,pInd3)+[0;4.5],...
                  XC4a(:,11:end,pInd3)+[0;4.5],conn4a,pSc2,cC3,cTr3);
hold on;
line_coordinates(XC4a(:,1:2,pInd3)+[0;+4.5], lSh, nW, lw_d);
line_coordinates(XC4a(:,9:10,pInd3)+[0;+4.5], -lSh, nW, lw_d);
line_coordinates(XC4a(:,1:2,pInd2)+[0;.5], lSh, nW, lw_d);
line_coordinates(XC4a(:,9:10,pInd2)+[0;.5], -lSh, nW, lw_d);
line_coordinates(XC4a(:,1:2,pInd1)+[0;-3.5], lSh, nW, lw_d);
line_coordinates(XC4a(:,9:10,pInd1)+[0;-3.5], -lSh, nW, lw_d);
hold off;

% Text
tSx = -2.8;
text(tSx,5.2,'$\mathrm{D^*}$','fontsize',10);
text(-tSx+7.5,2.7,'$\mathrm{D^*}$','fontsize',10);
text(tSx-2.4,0,'$\mathrm{\rightarrow D_1^o}$','fontsize',10);
text(-tSx+7.5,-1.5,'$\mathrm{\rightarrow ?}$','fontsize',10);
text(tSx-0.1,-4,'$\mathrm{D_1^o}$','fontsize',10);
text(-tSx+7.5,-5.5,'$\mathrm{D_1^o}$','fontsize',10);
text(labX,labY,'\textbf{k}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.32, 1, 'Motion', 'Units', 'Normalized', 'fontsize', 10);
              
axis([0 17.5 0 14] + [-3 -3 -5.6 -5.6])
drawnow;


%% l: Simulate
subplot(NRow,NCol,cellM{12}); cla;
% Simulate Single Module for Cobweb
[XCa,~] = sim_motion(Xs4,Xu4,conn4,.01,125,[Xs4 Xu4],0);   % Simulate
[XCb,~] = sim_motion(Xs4,Xu4,conn4,.01,70,-[Xs4 Xu4],0);   % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
D1 = sqrt(squeeze(sum(diff(XC4a,1,2).^2)));
D1 = D1(1:9,:);

plot(d1,d2,'k-','linewidth',1);
hold on;
plot([.5 4],[.5 4], '--', 'color', [200 200 200]/255);
% Cobweb 1
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [.5;dP(3:end)];
line(dP(1:end-1),[.5;dP(3:end)],'color',cTr1,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF*.8 diff(dPb(i:i+1))*caF*.8]);
end
% Cobweb 2
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [.5;dP(3:end)];
line(dP(1:end-1),[.5;dP(3:end)],'color',cTr2,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 3
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:); dP = dP(1:end-1);
dPa = dP(1:end-1); dPb = [.5;dP(3:end)];
line(dP(1:end-1),[.5;dP(3:end)],'color',cTr3,'linewidth',lwL);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end

% Markers
plot(D1(1,pInd3),D1(2,pInd3),'k*','linewidth',1,'markersize',6.5);
plot(D1(1,pInd1),D1(2,pInd1),'ko','linewidth',1,'markersize',6.5);
plot(D1(2,pInd1),D1(3,pInd1),'ko','linewidth',1,'markersize',6.5);
plot(d1,d2,'k-','linewidth',1);
hold off;

% Formatting
set(gca,'visible',1,'XTick',[.75 1.95],'YTick',[.75 1.95],...
                    'XTickLabel',[],'YTickLabel',[],'fontsize',10);
axis([.6 2.1 .6 2.1]);
text(.47,-.1,'$\mathrm{d_k}$','Units','normalized','fontsize',10);
text(-.1,.4,'$\mathrm{d_{k+1}}$','Units','normalized','rotation',90,'fontsize',10);
text(.58,.90,'$\mathrm{D^*}$','Units','normalized','fontsize',10);
text(.05,.90,'$\mathrm{D_1^o}$','Units','normalized','fontsize',10);
text(.87,.185,'$\mathrm{D_2^o}$','Units','normalized','fontsize',10);
text(labX,labY*1.1,'\textbf{l}','Units','Normalized','fontsize',10);
text(.32,labY*1.1,'Cobweb','Units','Normalized','fontsize',10);
drawnow;


%% Size and Save Figure
fName = 'figure2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters'; 
fig.PaperPosition = [-2.6 -1.65 24.5 13.1];
fig.PaperSize = [19 11.1];
% saveas(fig, ['Figures/' fName], 'pdf');