% Figure 4: Large-Scale Conformational Changes
%% Prepare Space
clear; clc;

% Subplot Indices
NCol = 14;
colBal = {[1 1 1 3], [1 1 1 3], [3 3]};
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
axM = [-1.5 1.5 -0.8 2];
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
% Flipped
Xs1f = [Xs10(1,:); -Xs10(2,:)+.5];
Xu1f = [Xu1(1,:); -Xu1(2,:)+.5];
% Show Network
% visualize_conic_finite(Xs10,Xs1T,[-2 2;-2 2],[100;100],0,1,1);
visualize_network(Xs10,Xu1,conn1,1,[255 100 100]/255, cTr1);
% axis(axM);
text(labX,labY,'a','Units','Normalized','fontsize',10,'fontweight','bold');

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
% Flipped
Xs2f = [Xs20(1,:); -Xs20(2,:)+.5];
Xu2f = [Xu2(1,:); -Xu2(2,:)+.5];
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
% Flipped
Xs3f = [Xs30(1,:); -Xs30(2,:)+.5];
Xu3f = [Xu3(1,:); -Xu3(2,:)+.5];
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


%% Simulate
[XCN,fC] = sim_motion(XsN,XuN,connN,.5,300,-[XsN,XuN],0);


%% d: Folded Network
subplot(NRow,NCol,cellM{9}); cla;
% Visualize
visualize_network(XsN,XuN,connN,.3,[255 100 100]/255,C_UNN);
da = mean(squeeze(sqrt(sum(diff(XCN(:,1:size(XsN,2),:),1,2).^2))) + l2(2));
pInd = find(abs(da) == min(abs(da)),1);
RV = rotz(7); RV = RV(1:2,1:2);
XSh = 42; YSh = -3;
visualize_network(RV*XCN(:,1:size(XsN,2),round(pInd/2))+[XSh;YSh],...
                  RV*XCN(:,[1:size(XuN,2)]+size(XsN,2),round(pInd/2))+[XSh;YSh],connN,...
                  .3,[255 100 100]/255,C_UN4);
XSh = 70; YSh = -3;
visualize_network(RV*XCN(:,1:size(XsN,2),pInd)+[XSh;YSh],...
                  RV*XCN(:,[1:size(XuN,2)]+size(XsN,2),pInd)+[XSh;YSh],connN,.3,...
                  [255 100 100]/255,C_UNN);
axis([-1,102,-10,10]);


%% e: All Letters: NETWORKS
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
[XsE,XuE,connE] = network_chain_x(Xsp,XuC,XuEL,C_UNC);



axis([-1,102,-10,10]);


%% Simulate
% N
[XCN,fCN] = sim_motion(XsN,XuN,connN,.5,300,-[XsN,XuN],0);
% E
[XCE,fCE] = sim_motion(XsE,XuN,connE,1,140,-[XsE,XuN],0);


%% T





%% Visualize
% Visualize
subplot(NRow,NCol,cellM{8}); cla;
visualize_network(XsE,XuN,connE,.3,[255 100 100]/255,C_UEN);
da = mean(squeeze(sqrt(sum(diff(XC(:,1:nSE+2,:),1,2).^2))) + l2(2));
pInd = find(abs(da) == min(abs(da)),1);
[pInd da(pInd)]
RV = rotz(106); RV = RV(1:2,1:2);
XSh = 58; YSh = -16;
visualize_network(RV*XC(:,1:nSE+2,round(pInd/2))+[XSh;YSh],...
                  RV*XC(:,[1:2*nSE]+nSE+2,round(pInd/2))+[XSh;YSh],connE,...
                  .3,[255 100 100]/255,C_UEN);
RV = rotz(-134); RV = RV(1:2,1:2);
XSh = 90; YSh = 12;
visualize_network(RV*XC(:,1:nSE+2,pInd)+[XSh;YSh],...
                  RV*XC(:,[1:2*nSE]+nSE+2,pInd)+[XSh;YSh],connE,.3,...
                  [255 100 100]/255,C_UEN);
axis([-1,102,-10,10]);




%% Animate Network Alone
% fig = figure(5); clf;
% fName = 'animation_net.gif';
% dT = 0.03;
% nSV = 1;
% Xs1a = Xs4;
% Xu1a = Xu4;
% conn1a = conn4;
% XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);
% 
% for i = 1:3:size(XC,3)-10
%     cla;
%     visualize_network(XC(:,1:size(Xs1a,2),i),...
%                       XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),conn1a,.5);
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
fig.PaperPosition = [-2.28 -0.0 24.4 7.2];
fig.PaperSize = [19 7.0];
saveas(fig, ['Figures/' fName], 'pdf');
