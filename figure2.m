% Figure 2: Iterated Maps and 1D Discrete Dynamical Systems
%% Prepare Space
clear; clc;

% Plot
cTr1 = [126 240 240]/255;
cTr2 = [115 164 211]/255;
cTr3 = [015 082 186]/255;
pSc = 0.6;
labX = -.2;
labY = 0.94;
labColY = 1.3;

% Subplot Indices
nRow = [7 7 16];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [5 5 5 5];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(2); clf;

% Global Annotation
annotation('line','linewidth',.2,'position',[.285 .1 0 .83],'color',[.9 .9 .9]);
annotation('line','linewidth',.2,'position',[.482 .1 0 .83],'color',[.9 .9 .9]);
annotation('line','linewidth',.2,'position',[.679 .1 0 .83],'color',[.9 .9 .9]);


%% a: Two Modules
subplot(NRow,NCol,cellM{1}); cla;
tSx = .7;
tSy = .1;
Xs1 = [-2 -2  2  2;...
       -1  1 -1  1];
conn1 = [1 3; 1 4; 2 3; 2 4];
visualize_network(Xs1-[3;0],[],conn1,pSc);
visualize_network(Xs1+[3;0],[],conn1,pSc);
text(-5-tSx, tSy, 'd_1','fontsize',10);
text(-1-tSx, tSy, 'd_2','fontsize',10);
text(1-tSx, tSy, 'd_1','fontsize',10);
text(5-tSx, tSy, 'd_2','fontsize',10);
axis([-5 5 -2 2]);
text(labX,labY,'a','Units','Normalized','fontsize',10,'fontweight','bold');
text(.15,labColY,'4-bar linkage','Units','Normalized','fontsize',10);


%% b: Coupled Modules
subplot(NRow,NCol,cellM{2}); cla;
[Xs1a,conn1a] = tesselate_network(Xs1,conn1,[4;0],[2;1]);
visualize_network(Xs1a,[],conn1a,pSc);
text(-2-tSx, tSy, 'd_1','fontsize',10);
text(2-tSx, tSy, 'd_2','fontsize',10);
text(6-tSx, tSy, 'd_3','fontsize',10);
axis([-5 5 -2 2]+[2 2 0 0]);
text(labX,labY,'b','Units','Normalized','fontsize',10,'fontweight','bold');


%% c: Cobweb Plot
subplot(NRow,NCol,cellM{3}); cla;
caF = .8;
% Simulate Single Module for Cobweb
dX1 = [0 0 0 0; -1 1 1 -1];
[XCa,~] = sim_motion(Xs1,[],conn1,.01,170,dX1,0);   % Simulate
[XCb,~] = sim_motion(Xs1,[],conn1,.01,170,-dX1,0);  % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,3:4,:),[],2)).^2)));
% Simulate Combined Network for Propagation
[XC1a,~] = sim_motion(Xs1a,[],conn1a,.01,140,Xs1a,0);
D1 = sqrt(squeeze(sum(diff(XC1a,1,2).^2)));
D1 = D1(1:2:end,:);
plot(d1,d2,'k-');
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
% Cobweb 1
pInd1 = 140;
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr1);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 2
pInd2 = 70;
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr2);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 3
pInd3 = 1;
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr3);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
plot(d1,d2,'k-');
hold off;
% Networks
visualize_network(XC1a(:,:,pInd1)/8+[1.7;3.6],[],conn1a,.33,cTr1);
visualize_network(XC1a(:,:,pInd2)/8+[1.85;2.9],[],conn1a,.33,cTr2);
visualize_network(XC1a(:,:,pInd3)/8+[2.6;2.0],[],conn1a,.33,cTr3);
% Formatting
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3],'fontsize',10);
axis([1 4 1 4]);
text(.47,-.18,'d_k','Units','normalized');
text(-.18,.4,'d_{k+1}','Units','normalized','rotation',90);
text(labX,labY,'c','Units','Normalized','fontsize',10,'fontweight','bold');


%% d: 2 FP + Super Stability
subplot(NRow,NCol,cellM{4});
s = sqrt(3);
Xs2 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
Xu2 = [-0.86 -0.86;...
       -1.45  1.47];
Xs2p = [Xs2(1,:); -Xs2(2,:)];
Xu2p = [Xu2(1,:); -Xu2(2,:)];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_network(Xs2,Xu2,conn2,pSc);
visualize_network(Xs2p+[3.5;.5],Xu2p+[3.5;.5],conn2,pSc);
visualize_network(Xs2+[7;0],Xu2+[7;0],conn2,pSc);
text(-1.2-tSx, .4+tSy, 'd_1','fontsize',10);
text(1-tSx, .4+tSy, 'd_2','fontsize',10);
text(2.3-tSx,-.1+tSy, 'd_1','fontsize',10);
text(4.5-tSx,-.1+tSy, 'd_2','fontsize',10);
text(5.8-tSx, .4+tSy, 'd_1','fontsize',10);
text(8.0-tSx, .4+tSy, 'd_2','fontsize',10);
axis([-1 8 -1.5 2]);
text(labX,labY,'d','Units','Normalized','fontsize',10,'fontweight','bold');
text(-.15,labColY,'Super-Stable Fixed Point','Units','Normalized','fontsize',10);


%% e: Tesselate
subplot(NRow,NCol,cellM{5});
Xs2c = [-s/2  0  s/2  s;...
        -1/2  1 -1/2  1];
Xu2c = [Xu2 Xu2p+[s/2;.5]];
conn2c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[4;1]);
visualize_network(Xs2a,Xu2a,conn2a,pSc);
axis([-1 8 -1.5 2]);
text(-1.2-tSx, .4+tSy, 'd_1','fontsize',10);
text(7.3-tSx,-.1+tSy, 'd_9','fontsize',10);
text(labX,labY,'e','Units','Normalized','fontsize',10,'fontweight','bold');


%% f: Simulate
subplot(NRow,NCol,cellM{6});% Simulate Single Module for Cobweb
[XCa,~] = sim_motion(Xs2,Xu2,conn2,.01,185,[Xs2 Xu2],0);   % Simulate
[XCb,~] = sim_motion(Xs2,Xu2,conn2,.01,15,-[Xs2 Xu2],0);   % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
% Simulate Combined Network for Propagation
[XC2a,fC] = sim_motion(Xs2a,Xu2a,conn2a,.02,570,[Xs2a Xu2a],0);
D1 = sqrt(squeeze(sum(diff(XC2a,1,2).^2)));
D1 = D1(1:9,:);
plot(d1,d2,'k-');
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
% Cobweb 3
pInd3 = 1;
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr3);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 2
pInd2 = 250;
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr2);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 1
pInd1 = 570;
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr1);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF*.8 diff(dPb(i:i+1))*caF*.8]);
end
plot(d1,d2,'k-');
hold off;
% Networks
visualize_network(XC2a(:,1:10,pInd1)/10+[1.51;3.05],...
                  XC2a(:,11:end,pInd1)/10+[1.51;3.05],conn2a,.33,cTr1);
visualize_network(XC2a(:,1:10,pInd2)/10+[1.43;2.6],...
                  XC2a(:,11:end,pInd2)/10+[1.43;2.6],conn2a,.33,cTr2);
visualize_network(XC2a(:,1:10,pInd3)/10+[1.23;2.07],...
                  XC2a(:,11:end,pInd3)/10+[1.23;2.07],conn2a,.33,cTr3);
% Formatting
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3],'fontsize',10);
axis([1 3.3 1 3.3]);
text(.47,-.18,'d_k','Units','normalized');
text(-.18,.4,'d_{k+1}','Units','normalized','rotation',90);
text(labX,labY,'f','Units','Normalized','fontsize',10,'fontweight','bold');


%% g: Limit Cycle
subplot(NRow,NCol,cellM{7});
Xs3 = [-s/2  0.0  s/2;...
        -0.5  1.0 -0.5];
Xu3 = [ 0.10 -0.30;...
       -0.25 -0.90];
Xs3p = [Xs3(1,:); -Xs3(2,:)];
Xu3p = [Xu3(1,:); -Xu3(2,:)];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_network(Xs3,Xu3,conn3,pSc);
visualize_network(Xs3p+[3.5;.5],Xu3p+[3.5;.5],conn3,pSc);
visualize_network(Xs3+[7;0],Xu3+[7;0],conn3,pSc);
text(-0.7-tSx, .4+tSy, 'd_1','fontsize',10);
text(1.2-tSx, .4+tSy, 'd_2','fontsize',10);
text(2.8-tSx,-.1+tSy, 'd_1','fontsize',10);
text(4.7-tSx,-.1+tSy, 'd_2','fontsize',10);
text(6.3-tSx, .4+tSy, 'd_1','fontsize',10);
text(8.2-tSx, .4+tSy, 'd_2','fontsize',10);
text(labX,labY,'g','Units','Normalized','fontsize',10,'fontweight','bold');
text(-.02,labColY,'Isolated Limit Cycle','Units','Normalized','fontsize',10);
axis([-1 8 -1.5 2]);


%% h: Tesselate Limit Cycle
subplot(NRow,NCol,cellM{8});
Xs3c = [-s/2  0  s/2  s;...
        -1/2  1 -1/2  1];
Xu3c = [Xu3 Xu3p+[s/2;.5]];
conn3c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[4;1]);
visualize_network(Xs3a,Xu3a,conn3a,pSc);
axis([-1 8 -1.5 2]);
text(-0.7-tSx, .4+tSy, 'd_1','fontsize',10);
text(7.5-tSx,-.1+tSy, 'd_9','fontsize',10);
text(labX,labY,'h','Units','Normalized','fontsize',10,'fontweight','bold');


%% i: Simulate
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[7;1]);
subplot(NRow,NCol,cellM{9});% Simulate Single Module for Cobweb
[XCa,~] = sim_motion(Xs3,Xu3,conn3,.01,90,[Xs3 Xu3],0);   % Simulate
[XCb,~] = sim_motion(Xs3,Xu3,conn3,.01,100,-[Xs3 Xu3],0);   % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
% Simulate Combined Network for Propagation
[XC3a,fC] = sim_motion(Xs3a,Xu3a,conn3a,.02,570,-[Xs3a Xu3a],0);
D1 = sqrt(squeeze(sum(diff(XC3a,1,2).^2)));
D1 = D1(1:9,:);
plot(d1,d2,'k-');
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
% Cobweb 3
pInd3 = 1;
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr3);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 2
pInd2 = 300;
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr2);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 1
pInd1 = 462;
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];
line(dP(1:end-1),[1;dP(3:end)],'color',cTr1);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF*.8 diff(dPb(i:i+1))*caF*.8]);
end
plot(d1,d2,'k-');
hold off;
% Networks
visualize_network(XC3a(:,1:10,pInd1)/10+[1.6;3.03],...
                  XC3a(:,11:end,pInd1)/10+[1.6;3.03],conn3a,.33,cTr1);
visualize_network(XC3a(:,1:10,pInd2)/10+[1.6;2.63],...
                  XC3a(:,11:end,pInd2)/10+[1.6;2.63],conn3a,.33,cTr2);
visualize_network(XC3a(:,1:10,pInd3)/10+[1.6;2.23],...
                  XC3a(:,11:end,pInd3)/10+[1.6;2.23],conn3a,.33,cTr3);
% Formatting
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3],'fontsize',10);
axis([1 3.3 1 3.3]);
text(.47,-.18,'d_k','Units','normalized');
text(-.18,.4,'d_{k+1}','Units','normalized','rotation',90);
text(labX,labY,'i','Units','Normalized','fontsize',10,'fontweight','bold');


%% j: Chaos
subplot(NRow,NCol,cellM{10});
Xs4 = [-s/2  0.0  s/2;...
        -0.5  0.0 -0.5]*2;
Xu4 = [ 0.15 -0.22;...
       -0.50 -0.90]*2;
Xs4p = [Xs4(1,:); -Xs4(2,:)];
Xu4p = [Xu4(1,:); -Xu4(2,:)];
conn4 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_network(Xs4,Xu4,conn4,pSc);
visualize_network(Xs4p+[3.5;.5],Xu4p+[3.5;.5],conn4,pSc);
visualize_network(Xs4+[7;0],Xu4+[7;0],conn4,pSc);
text(-0.8-tSx, .0+tSy, 'd_1','fontsize',10);
text(1.1-tSx, .0+tSy, 'd_2','fontsize',10);
text(2.7-tSx, .5+tSy, 'd_1','fontsize',10);
text(4.6-tSx, .5+tSy, 'd_2','fontsize',10);
text(6.2-tSx, .0+tSy, 'd_1','fontsize',10);
text(8.1-tSx, .0+tSy, 'd_2','fontsize',10);
axis([-1 8 -1.5 2]*1.2 - .7*[1 1 0 0]);
text(labX,labY,'j','Units','Normalized','fontsize',10,'fontweight','bold');
text(.25,labColY,'Chaos','Units','Normalized','fontsize',10);


%% k: Tesselate Chaos
subplot(NRow,NCol,cellM{11});
Xs4c = [-s/2  0  s/2  s;...
        -1/2  0 -1/2  0]*2;
Xu4c = [Xu4 Xu4p+[s;-1]];
conn4c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs4a,Xu4a,conn4a] = tesselate_network_old(Xs4c,Xu4c,conn4c,[2*s;0],[4;1]);
visualize_network(Xs4a,Xu4a,conn4a,pSc);
axis([-1 8 -1.5 2]*1.8 + 0*[1 1 0 0]);
text(-1.5-tSx, .4+tSy, 'd_1','fontsize',10);
text(13.5-tSx,-1.3+tSy, 'd_9','fontsize',10);
text(labX,labY,'k','Units','Normalized','fontsize',10,'fontweight','bold');


%% l: Simulate
subplot(NRow,NCol,cellM{12});% Simulate Single Module for Cobweb
[XCa,~] = sim_motion(Xs4,Xu4,conn4,.01,150,[Xs4 Xu4],0);   % Simulate
[XCb,~] = sim_motion(Xs4,Xu4,conn4,.01,90,-[Xs4 Xu4],0);   % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
% Simulate Combined Network for Propagation
[XC4a,fC] = sim_motion(Xs4a,Xu4a,conn4a,.002,12000,-[Xs4a Xu4a],0);
D1 = sqrt(squeeze(sum(diff(XC4a,1,2).^2)));
D1 = D1(1:9,:);
plot(d1,d2,'k-');
hold on;
plot([.5 4],[.5 4], '--', 'color', [200 200 200]/255);
% Cobweb 3
pInd3 = 1;
dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [.5;dP(3:end)];
line(dP(1:end-1),[.5;dP(3:end)],'color',cTr3);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr3);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 2
pInd2 = 5000;
dP = [D1(:,pInd2)';D1(:,pInd2)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [.5;dP(3:end)];
line(dP(1:end-1),[.5;dP(3:end)],'color',cTr2);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr2);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
% Cobweb 1
pInd1 = 11260;
dP = [D1(:,pInd1)';D1(:,pInd1)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [.5;dP(3:end)];
line(dP(1:end-1),[.5;dP(3:end)],'color',cTr1);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',3,'HeadWidth',3,'color',cTr1);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF*.8 diff(dPb(i:i+1))*caF*.8]);
end
plot(d1,d2,'k-');
hold off;
% Networks
visualize_network(XC4a(:,1:10,pInd1)/10+[1.6;3.8],...
                  XC4a(:,11:end,pInd1)/10+[1.6;3.8],conn4a,.33,cTr1);
visualize_network(XC4a(:,1:10,pInd2)/10+[1.6;3.3],...
                  XC4a(:,11:end,pInd2)/10+[1.6;3.3],conn4a,.33,cTr2);
visualize_network(XC4a(:,1:10,pInd3)/10+[1.6;2.83],...
                  XC4a(:,11:end,pInd3)/10+[1.6;2.83],conn4a,.33,cTr3);
% Formatting
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3],'fontsize',10);
axis([.8 4 .8 4]);
text(.47,-.18,'d_k','Units','normalized');
text(-.18,.4,'d_{k+1}','Units','normalized','rotation',90);
text(labX,labY,'l','Units','Normalized','fontsize',10,'fontweight','bold');


%% Size and Save Figure
fName = 'figure2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters'; 
fig.PaperPosition = [-2.28 -0.3 24.4 8.9];
fig.PaperSize = [19 8.5];
saveas(fig, ['Figures/' fName], 'pdf');