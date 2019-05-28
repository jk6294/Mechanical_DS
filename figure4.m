% Figure 4: Large-Scale Conformational Changes
%% Prepare Space
clear; clc;

% Subplot Indices
NCol = 14;
colBal = {[1 1 1 4], [1 1 1 4], [3 4]};
rowBal = [9 9 9]; 
NRow = sum(rowBal);
ru = cumsum(rowBal)-1;
rd = [0 ru(1:end-1)+1]+1;
cellL = [0 cumsum(cellfun('length',colBal))];
cellM = cell(1,cellL(end));
for i = 1:length(colBal)
    s = colBal{i}; 
    su = cumsum(s/sum(s)*(NCol+1))-1;
    sd = [1 su(1:end-1)+2];
    for j = 1:length(s)
        cellM(j+cellL(i)) = {[[sd(j) su(j)]+rd(i)*NCol,...
                              [sd(j) su(j)]+ru(i)*NCol]};
    end
end

% Figure Axis Bounds
axM = [-1.5 1.5 -1.0 2.8];
labX = -.1;
labY = 0.90;
labColY = 1.3;

fig = figure(4); clf;

% Sizing
s = sqrt(3);
pSc = 0.6;

% Node Colors for different curvatures
cTr1 = [126 240 240]/255;
cTr2 = [115 164 211]/255;
cTr3 = [015 082 186]/255;

% Module Parameters
l1 = [0;-s]; l2 = [0;-1.5*s];
a1 = 19.4; a2 = 38; a3 = 60;
R1 = rotz(-a1/2); R2 = rotz(-a2/2); R3 = rotz(-a3/2);
R1 = R1(1:2,1:2); R2 = R2(1:2,1:2); R3 = R3(1:2,1:2);


%% a: Decrease Distance
subplot(NRow,NCol,cellM{1}); cla;
% Construct Network
Xs10 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs1T = [R1*l2 [0;0] R1\l2] + [0;1];
Xu1 = [-1.5  1.5;...
        1.5  1.5];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu1,fV] = construct_network(Xs10,Xs1T,Xu1,conn1,0,1);
Xu1 = Xu1(1:2,:);
visualize_network(Xs10,Xu1,conn1,1,[255 100 100]/255, cTr1);
axis(axM);
text(labX,labY,'a','Units','Normalized','fontsize',10,'fontweight','bold');

annotation('arrow','HeadLength',6,'HeadWidth',8,'color',[.7 .7 .7],...
           'linewidth',2,'position',[.155 .68 0 -.08]);
annotation('arrow','HeadLength',6,'HeadWidth',8,'color',[.7 .7 .7],...
           'linewidth',2,'position',[.276 .68 0 -.08]);
annotation('arrow','HeadLength',6,'HeadWidth',8,'color',[.7 .7 .7],...
           'linewidth',2,'position',[.397 .68 0 -.08]);

subplot(NRow,NCol,cellM{5}); cla;
% Simulate Network
[XC1,fC] = sim_motion(Xs10,Xu1,conn1,.01,220,-[Xs10 Xu1],0);
d1 = sqrt(squeeze(sum((XC1(:,1,:)-XC1(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC1(:,2,:)-XC1(:,3,:)).^2)));
d3 = sqrt(squeeze(sum((XC1(:,1,:)-XC1(:,3,:)).^2)));
distv = sqrt(sum(diff([Xs1T Xs1T(:,1)],1,2).^2))
distL = sum(abs(distv - [d1 d2 d3]),2);
% Find Correct Point in Simulation
pInd = find(distL==min(distL),1);
visualize_network(XC1(:,1:3,pInd),XC1(:,4:5,pInd),conn1,...
                  1,[255 100 100]/255, cTr1);
% plot(d1); hold on; plot(d2); plot(d3); hold off;
axis(axM);
drawnow;


%% b: Maintain Distance
subplot(NRow,NCol,cellM{2}); cla;
% Construct Network
Xs20 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs2T = [R2*l2 [0;0] R2\l2] + [0;1];
Xu2 = [-0.7  0.7;...
        1.5  1.5];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,fV] = construct_network(Xs20,Xs2T,Xu2,conn2,0,1);
Xu2 = Xu2(1:2,:);
% Show Network
visualize_network(Xs20,Xu2,conn2,1,[255 100 100]/255, cTr2);
axis(axM);
text(labX,labY,'b','Units','Normalized','fontsize',10,'fontweight','bold');

subplot(NRow,NCol,cellM{6}); cla;
% Simulate Network
[XC2,fC] = sim_motion(Xs20,Xu2,conn2,.01,200,[Xs20 Xu2],0);
d1 = sqrt(squeeze(sum((XC2(:,1,:)-XC2(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC2(:,2,:)-XC2(:,3,:)).^2)));
d3 = sqrt(squeeze(sum((XC2(:,1,:)-XC2(:,3,:)).^2)));
distv = sqrt(sum(diff([Xs2T Xs2T(:,1)],1,2).^2))
distL = sum(abs(distv - [d1 d2 d3]),2);
% Find Correct Point in Simulation
pInd = find(distL==min(distL),1);
visualize_network(XC2(:,1:3,pInd),XC2(:,4:5,pInd),conn1,...
                  1,[255 100 100]/255, cTr2);
% plot(d1); hold on; plot(d3); hold off;
axis(axM);
drawnow;


%% c: Increase Distance
subplot(NRow,NCol,cellM{3}); cla;
% Construct Network
Xs30 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs3T = [R3*l2 [0;0] R3\l2] + [0;1];
Xu3 = [-1.0  1.0;...
        1.8  1.8];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu3,fV] = construct_network(Xs30,Xs3T,Xu3,conn3,0,1);
Xu3 = Xu3(1:2,:);
% Show Network
visualize_network(Xs30,Xu3,conn1,1,[255 100 100]/255, cTr3);
axis(axM);
text(labX,labY,'c','Units','Normalized','fontsize',10,'fontweight','bold');

subplot(NRow,NCol,cellM{7}); cla;
% Simulate Network
[XC3,fC] = sim_motion(Xs30,Xu3,conn3,.01,150,[Xs30 Xu3],0);
d1 = sqrt(squeeze(sum((XC3(:,1,:)-XC3(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC3(:,2,:)-XC3(:,3,:)).^2)));
d3 = sqrt(squeeze(sum((XC3(:,1,:)-XC3(:,3,:)).^2)));
distv = sqrt(sum(diff([Xs3T Xs3T(:,1)],1,2).^2))
distL = sum(abs(distv - [d1 d2 d3]),2);
% Find Correct Point in Simulation
pInd = find(distL==min(distL),1);
visualize_network(XC3(:,1:3,pInd),XC3(:,4:5,pInd),conn1,...
                  1,[255 100 100]/255, cTr3);
% plot(d1); hold on; plot(d3); hold off;
axis(axM);
drawnow;


%% d: Example Line
% Collect Unspecified Node Posiitons and Colors
XuC = cat(3,Xu1,Xu2,Xu3);
C_UNC = [cTr1; cTr2; cTr3];
% Decide Placement Order of Unspecified Nodes
XuNL = [2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        1 3, 1 3, 1 3, 1 3,...
        2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        3 1 ,3 1, 3 1, 3 1,...
        2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2];
[XsN,XuN,connN,C_UNN] = network_chain_x([-s/2 0; -.5 1],XuC,XuNL,C_UNC);

% Simulate
[XCN,fC] = sim_motion(XsN,XuN,connN,.5,300,-[XsN,XuN],0);

subplot(NRow,NCol,cellM{9}); cla;
% Visualize
visualize_network(XsN,XuN,connN,.3,[255 100 100]/255,C_UNN);
da = mean(squeeze(sqrt(sum(diff(XCN(:,1:size(XsN,2),:),1,2).^2))) + l2(2));
pInd = find(abs(da) == min(abs(da)),1);
RV = rotz(6.5); RV = RV(1:2,1:2);
XSh = 54; YSh = -3;
visualize_network(RV*XCN(:,1:size(XsN,2),round(pInd/2))+[XSh;YSh],...
                  RV*XCN(:,[1:size(XuN,2)]+size(XsN,2),round(pInd/2))+[XSh;YSh],connN,...
                  .3,[255 100 100]/255,C_UNN);
XSh = 92; YSh = -3;
visualize_network(RV*XCN(:,1:size(XsN,2),pInd)+[XSh;YSh],...
                  RV*XCN(:,[1:size(XuN,2)]+size(XsN,2),pInd)+[XSh;YSh],connN,.3,...
                  [255 100 100]/255,C_UNN);
axis([-1,124,-10,13]);
annotation('arrow','HeadLength',6,'HeadWidth',8,'color',[.7 .7 .7],...
           'linewidth',2,'position',[.255 .21 .02 0]);
annotation('arrow','HeadLength',6,'HeadWidth',8,'color',[.7 .7 .7],...
           'linewidth',2,'position',[.355 .21 .02 0]);
text(-.02,labY,'d','Units','Normalized','fontsize',10,'fontweight','bold');


%% All Letters: NETWORKS
% Collect Unspecified Node Posiitons and Colors
XuC = cat(3,Xu1,Xu2,Xu3);
C_UNC = [cTr1; cTr2; cTr3];
Xsp = [-s/2 0; -.5 1];

% N
XuNL = [2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        1 3, 1 3, 1 3, 1 3,...
        2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        3 1 ,3 1, 3 1, 3 1,...
        2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2];
[XsN,XuN,connN,C_UNN] = network_chain_x(Xsp,XuC,XuNL,C_UNC);

% E
XuEL = [2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 3 2, 3 2, 3 2, 3 2 ,3 2, 3 2, 3 2,...
        2 2,...
        3 2, 3 2, 3 2, 3 2, 3 2,...
        2 2];
[XsE,XuE,connE,C_UNE] = network_chain_x(Xsp,XuC,XuEL,C_UNC);

% T
XuTL = [2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        1 3, 1 3, 1 3, 1 3, 1 3, 1 3, 1 2,...
        2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 3 1];
[XsT,XuT,connT,C_UNT] = network_chain_x(Xsp,XuC,XuTL,C_UNC);

% W
XuWL = [2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        1 3, 1 3, 1 3, 1 3,...
        2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 3 1,...
        2 2, 2 2, 2 2,...
        1 3, 1 3, 1 3, 1 3,...
        2 2, 2 2, 2 2, 2 2, 2 2, 2 2];
[XsW,XuW,connW,C_UNW] = network_chain_x(Xsp,XuC,XuWL,C_UNC);

% O
XuOL = [3 2, 3 2, 3 2, 3 2, 3 2, 3 2, 3 2, 3 2,...
        3 2, 3 2, 3 2, 3 2, 3 2, 3 2, 3 2, 3 2, 3];
[XsO,XuO,connO,C_UNO] = network_chain_x(Xsp,XuC,XuOL,C_UNC);

% R
XuRL = [2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 3 1, 3 1, 3 1, 3 1, 3 2, 3 1,...
        1 3, 1 3, 1 3, 2 3, 2 3];
[XsR,XuR,connR,C_UNR] = network_chain_x(Xsp,XuC,XuRL,C_UNC);

% K
XuKL = [2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 2 1, 2 1,...
        2 2, 2 2, 2 2, 2 2, 2 2,...
        3 1, 3 1, 3 1, 3 1, 2 2, 3 1, 3 1,...
        2 2, 2 2, 2 2, 2 2, 2 2, 2 2, 2 2];
[XsK,XuK,connK,C_UNK] = network_chain_x(Xsp,XuC,XuKL,C_UNC);


%% Simulate
% N
[XCN,fCN] = sim_motion(XsN,XuN,connN,.5,300,-[XsN,XuN],0);
da = mean(squeeze(sqrt(sum(diff(XCN(:,1:size(XsN,2),:),1,2).^2))) + l2(2));
pIndN = find(abs(da) == min(abs(da)),1);
% E
[XCE,fCE] = sim_motion(XsE,XuE,connE,1,140,-[XsE,XuE],0);
da = mean(squeeze(sqrt(sum(diff(XCE(:,1:size(XsE,2),:),1,2).^2))) + l2(2));
pIndE = find(abs(da) == min(abs(da)),1);
% T
[XCT,fCT] = sim_motion(XsT,XuT,connT,1,200,-[XsT,XuT],0);
da = mean(squeeze(sqrt(sum(diff(XCT(:,1:size(XsT,2),:),1,2).^2))) + l2(2));
pIndT = find(abs(da) == min(abs(da)),1);
% W
[XCW,fCW] = sim_motion(XsW,XuW,connW,1,150,-[XsW,XuW],0);
da = mean(squeeze(sqrt(sum(diff(XCW(:,1:size(XsW,2),:),1,2).^2))) + l2(2));
pIndW = find(abs(da) == min(abs(da)),1);
% O
[XCO,fCO] = sim_motion(XsO,XuO,connO,1,90,-[XsO,XuO],0);
da = mean(squeeze(sqrt(sum(diff(XCO(:,1:size(XsO,2),:),1,2).^2))) + l2(2));
pIndO = find(abs(da) == min(abs(da)),1);
% R
[XCR,fCR] = sim_motion(XsR,XuR,connR,1,200,-[XsR,XuR],0);
da = mean(squeeze(sqrt(sum(diff(XCR(:,1:size(XsR,2),:),1,2).^2))) + l2(2));
pIndR = find(abs(da) == min(abs(da)),1);
% K
[XCK,fCK] = sim_motion(XsK,XuK,connK,1.5,240,-[XsK,XuK],0);
da = mean(squeeze(sqrt(sum(diff(XCK(:,1:size(XsK,2),:),1,2).^2))) + l2(2));
pIndK = find(abs(da) == min(abs(da)),1);


%% e: Visualize Original and Transition
subplot(NRow,NCol,[cellM{4} cellM{8}]); cla;
nM = 60;
nI = 2;
nSh = 25;
axis([-1,174,0 nM-nI+3]);
R90 = rotz(-90); R90 = R90(1:2,1:2);

% N
visualize_network(R90*XsN+[0;nM-nI],R90*XuN+[0;nM-nI],connN,.3,[],C_UNN);
RV = rotz(-90); RV = RV(1:2,1:2);
XSh = 10; YSh = 45;
visualize_network(RV*XCN(:,1:size(XsN,2),round(pIndN/2))+[XSh;YSh],...
                  RV*XCN(:,[1:size(XuN,2)]+size(XsN,2),round(pIndN/2))+[XSh;YSh],connN,.3,...
                  [255 100 100]/255,C_UNN);
drawnow;

% E
visualize_network(R90*XsE+[1*nSh;nM-nI],R90*XuE+[1*nSh;nM-nI],connE,.3,[],C_UNE);
RV = rotz(-100); RV = RV(1:2,1:2);
XSh = 36; YSh = 34;
visualize_network(RV*XCE(:,1:size(XsE,2),round(pIndE/2))+[XSh;YSh],...
                  RV*XCE(:,[1:size(XuE,2)]+size(XsE,2),round(pIndE/2))+[XSh;YSh],connE,.3,...
                  [255 100 100]/255,C_UNE);
drawnow;

% T
visualize_network(R90*XsT+[2*nSh;nM-nI],R90*XuT+[2*nSh;nM-nI],connT,.3,[],C_UNT);
RV = rotz(-37); RV = RV(1:2,1:2);
XSh = 44; YSh = 28;
visualize_network(RV*XCT(:,1:size(XsT,2),round(0.7*pIndT))+[XSh;YSh],...
                  RV*XCT(:,[1:size(XuT,2)]+size(XsT,2),round(0.7*pIndT))+[XSh;YSh],connT,.3,...
                  [255 100 100]/255,C_UNT);
drawnow;

% W
visualize_network(R90*XsW+[3*nSh;nM-nI],R90*XuW+[3*nSh;nM-nI],connW,.3,[],C_UNW);
RV = rotz(270); RV = RV(1:2,1:2);
XSh = 85; YSh = 48;
visualize_network(RV*XCW(:,1:size(XsW,2),round(pIndW/2))+[XSh;YSh],...
                  RV*XCW(:,[1:size(XuW,2)]+size(XsW,2),round(pIndW/2))+[XSh;YSh],connW,.3,...
                  [255 100 100]/255,C_UNW);
drawnow;
              
% O
visualize_network(R90*XsO+[4*nSh;nM-nI],R90*XuO+[4*nSh;nM-nI],connO,.3,[],C_UNO);
RV = rotz(-90); RV = RV(1:2,1:2);
XSh = 103; YSh = 30;
visualize_network(RV*XCO(:,1:size(XsO,2),round(pIndO/2))+[XSh;YSh],...
                  RV*XCO(:,[1:size(XuO,2)]+size(XsO,2),round(pIndO/2))+[XSh;YSh],connO,.3,...
                  [255 100 100]/255,C_UNO);
drawnow;

% R
visualize_network(R90*XsR+[5*nSh;nM-nI],R90*XuR+[5*nSh;nM-nI],connR,.3,[],C_UNR);
RV = rotz(-21); RV = RV(1:2,1:2);
XSh = 108; YSh = 16;
visualize_network(RV*XCR(:,1:size(XsR,2),round(pIndR/2))+[XSh;YSh],...
                  RV*XCR(:,[1:size(XuR,2)]+size(XsR,2),round(pIndR/2))+[XSh;YSh],connR,.3,...
                  [255 100 100]/255,C_UNR);
drawnow;

% K
visualize_network(R90*XsK+[6*nSh-5;nM-nI],R90*XuK+[6*nSh-5;nM-nI],connK,.3,[],C_UNK);
RV = rotz(-10); RV = RV(1:2,1:2);
XSh = 133; YSh = 15;
visualize_network(RV*XCK(:,1:size(XsK,2),round(pIndK/2))+[XSh;YSh],...
                  RV*XCK(:,[1:size(XuK,2)]+size(XsK,2),round(pIndK/2))+[XSh;YSh],connK,.3,...
                  [255 100 100]/255,C_UNK);
drawnow;

text(nSh/175-.1,0.955,'N','Units','Normalized','fontsize',10);
text(2*nSh/175-.1,0.955,'E','Units','Normalized','fontsize',10);
text(3*nSh/175-.1,0.955,'T','Units','Normalized','fontsize',10);
text(4*nSh/175-.1,0.955,'W','Units','Normalized','fontsize',10);
text(5*nSh/175-.1,0.955,'O','Units','Normalized','fontsize',10);
text(6*nSh/175-.1,0.955,'R','Units','Normalized','fontsize',10);
text(7*nSh/175-.1,0.955,'K','Units','Normalized','fontsize',10);
text(-.05,0.955,'e','Units','Normalized','fontsize',10,'fontweight','bold');


%% f: Full Word: NETWORK
subplot(NRow,NCol,cellM{10}); cla;
axis([15,180,-10,13] - [20 20 0 0]);

% N
RV = rotz(6.5); RV = RV(1:2,1:2);
XSh = 0; YSh = -3;
visualize_network(RV*XCN(:,1:size(XsN,2),pIndN)+[XSh;YSh],...
                  RV*XCN(:,[1:size(XuN,2)]+size(XsN,2),pIndN)+[XSh;YSh],connN,.3,...
                  [255 100 100]/255,C_UNN);
drawnow;

% E
RV = rotz(-133); RV = RV(1:2,1:2);
XSh = 52; YSh = 12;
visualize_network(RV*XCE(:,1:size(XsE,2),pIndE)+[XSh;YSh],...
                  RV*XCE(:,[1:size(XuE,2)]+size(XsE,2),pIndE)+[XSh;YSh],connE,.3,...
                  [255 100 100]/255,C_UNE);
drawnow;

% T
RV = rotz(-21); RV = RV(1:2,1:2);
XSh = 36; YSh = 10;
visualize_network(RV*XCT(:,1:size(XsT,2),pIndT+5)+[XSh;YSh],...
                  RV*XCT(:,[1:size(XuT,2)]+size(XsT,2),pIndT+5)+[XSh;YSh],connT,.3,...
                  [255 100 100]/255,C_UNT);
drawnow;
              
% W
RV = rotz(180); RV = RV(1:2,1:2);
XSh = 102; YSh = -1.5;
visualize_network(RV*XCW(:,1:size(XsW,2),pIndW)+[XSh;YSh],...
                  RV*XCW(:,[1:size(XuW,2)]+size(XsW,2),pIndW)+[XSh;YSh],connW,.3,...
                  [255 100 100]/255,C_UNW);
drawnow;

% O
RV = rotz(0); RV = RV(1:2,1:2);
XSh = 82; YSh = -1.5;
visualize_network(RV*XCO(:,1:size(XsO,2),pIndO)+[XSh;YSh],...
                  RV*XCO(:,[1:size(XuO,2)]+size(XsO,2),pIndO)+[XSh;YSh],connO,.3,...
                  [255 100 100]/255,C_UNO);
drawnow;
              
% R
RV = rotz(222); RV = RV(1:2,1:2);
XSh = 123; YSh = 14;
visualize_network(RV*XCR(:,1:size(XsR,2),pIndR)+[XSh;YSh],...
                  RV*XCR(:,[1:size(XuR,2)]+size(XsR,2),pIndR)+[XSh;YSh],connR,.3,...
                  [255 100 100]/255,C_UNR);
drawnow;
              
% K
RV = rotz(67); RV = RV(1:2,1:2);
XSh = 116; YSh = -26;
visualize_network(RV*XCK(:,1:size(XsK,2),pIndK)+[XSh;YSh],...
                  RV*XCK(:,[1:size(XuK,2)]+size(XsK,2),pIndK)+[XSh;YSh],connK,.3,...
                  [255 100 100]/255,C_UNK);
drawnow;
text(-.05,labY,'f','Units','Normalized','fontsize',10,'fontweight','bold');


%% Animate Network Alone
% fig = figure(5); clf;
% fName = 'animation_net.gif';
% dT = 0.03;
% nSV = 1;
% Xs1a = XsN;
% Xu1a = XuN;
% conn1a = connN;
% XC = XCN;
% 
% RV = rotz(6.5); RV = RV(1:2,1:2);
% for i = 1:size(XC,3)
%     XC(:,:,i) = RV*XC(:,:,i);
% end
% 
% 
% for i = 1:5:size(XC,3)-10
%     cla;
%     visualize_network(XC(:,1:size(Xs1a,2),i),...
%                       XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),conn1a,.7);
%     axis([min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))]);
%     set(gca,'visible',0);
%     drawnow;
% 
%     % Capture the plot as an image
%     frame = getframe(fig);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if nSV == 1
%       imwrite(imind,cm,fName,'gif', 'Loopcount',inf,'DelayTime',dT);
%       nSV = 0;
%     else
%       imwrite(imind,cm,fName,'gif','WriteMode','append','DelayTime',dT);
%     end
% end


%% Size and Save Figure
fName = 'figure4';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters'; 
fig.PaperPosition = [-3.0 -0.75 24.3 6.4];
fig.PaperSize = [19 5.0];
saveas(fig, ['Figures/' fName], 'pdf');

