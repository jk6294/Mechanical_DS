%% Prepare Space
clear; clc;


%% Simple Example: 4 Nodes
clear; clf;
Xs = [ 0  0;...
       0  1];
Xu = [ 1.5  1.5;...
       0  1];
conn = [1 1; 1 2; 2 1; 2 2];

% Simulate
[XC, fC] = sim_motion(Xs, Xu, conn, .005, 2000, ones(size([Xs, Xu])), 0);

% Configuration 1
pInd1 = 1;
fig = figure(1); clf;
subplot(1,3,1);
visualize_network(XC(:,1:2,pInd1),XC(:,3:4,pInd1),conn);
hold on;
plot(XC(1,1:2,pInd1), XC(2,1:2,pInd1), 'r--');
plot(XC(1,3:4,pInd1), XC(2,3:4,pInd1), 'b--');
hold off;
axis([-.5 2.0 -.5 2.0]);

% Configuration 2
pInd2 = 100;
subplot(1,3,2);
visualize_network(XC(:,1:2,pInd2),XC(:,3:4,pInd2),conn);
hold on;
plot(XC(1,1:2,pInd2), XC(2,1:2,pInd2), 'r--');
plot(XC(1,3:4,pInd2), XC(2,3:4,pInd2), 'b--');
hold off;
axis([-.5 2.0 -.5 2.0]);

% Distances
LVal = sqrt(sum((Xs(:,conn(:,1)) - Xu(:,conn(:,2))).^2));
LV1 = squeeze(sqrt(sum((XC(:,1,:) - XC(:,2,:)).^2)));
LV2 = squeeze(sqrt(sum((XC(:,3,:) - XC(:,4,:)).^2)));
subplot(1,3,3);
plot(LV1,LV2,'k-');
hold on;
plot([0 3], [0 3], 'k--');
plot(LV1(pInd1), LV2(pInd1), 'ko');
plot(LV1(pInd2), LV2(pInd2), 'kx');
hold off;
axis([0 3.5 0 3.5]);
set(gca,'XColor','r', 'YColor','b');
xlabel('r_1', 'color', 'r');
ylabel('r_2', 'color', 'b');

% Size and Save
fName = '4bar';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 2.3];
fig.PaperSize = [8.2 2.3];
saveas(fig, ['Figures\' fName], 'pdf');


%% Combinations
Xs = [ 0  0;...
       0  1];
Xu = [ 1.5  1.5;...
       0  1];
conn = [1 1; 1 2; 2 1; 2 2];

% Multiple Modules
fig = figure(2); clf;
subplot(2,3,1);
visualize_network(Xs,Xu,conn);
visualize_network(Xs + [2;0],Xu + [2;0],conn);
hold on;
plot([Xs(1,1)+2, Xu(1,1)], [Xs(2,1) Xu(2,1)], 'k--');
plot([Xs(1,2)+2, Xu(1,2)], [Xs(2,2) Xu(2,2)], 'k--');
hold off;
axis([-.5 3.5 -.5 2.0]);

% Combined Modules
subplot(2,3,4);
visualize_network(Xu+[1.5;0],Xs+[1.5;0],conn);
visualize_network(Xs,Xu,conn);
axis([-.5 3.5 -.5 2.0]);
hold on;
plot(Xs(1,:), Xs(2,:), 'r--');
plot(Xu(1,:), Xu(2,:), 'm--');
plot(Xu(1,:)+1.5, Xu(2,:), 'b--');
hold off;

% Combined: 1
subplot(2,3,2);
pInd = 50;
Xs = [ 0.0 0.0 3.0 3.0;...
       0.0 1.0 0.0 1.0];
Xu = [ 1.5 1.5;...
       0.0 1.0];
conn = [1 1; 1 2; 2 1; 2 2;...
        3 1; 3 2; 4 1; 4 2];
[XC, fC] = sim_motion(Xs, Xu, conn, 0.01, 200, ones(size([Xs Xu])), 0);
visualize_network(XC(:,1:size(Xs,2),pInd), XC(:,[1:size(Xu,2)]+size(Xs,2),pInd), conn);
axis([-.5 3.5 -.5 2.0]);

% Map 1
dV = repmat(sqrt(sum((XC(:,[1 5 3],pInd) - XC(:,[2 6 4],pInd)).^2)),[2,1]);
dV = dV(:);
subplot(2,3,3);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
xlabel('r_{n}');
ylabel('r_{n+1}');
axis([0 3.5 0 3.5]);

% Combined: 2
subplot(2,3,5);
pInd = 100;
Xs = [ 0.0 0.0 3.0 3.0;...
       0.0 1.0 0.0 1.0];
Xu = [ 1.5 1.5;...
       0.0 1.0];
conn = [1 1; 1 2; 2 1; 2 2;...
        3 1; 3 2; 4 1; 4 2];
[XC2, fC] = sim_motion(Xs, Xu, conn, 0.01, 200, -ones(size([Xs Xu])), 0);
visualize_network(XC2(:,1:size(Xs,2),pInd), XC2(:,[1:size(Xu,2)]+size(Xs,2),pInd), conn);
axis([-.5 3.5 -.5 2.0]);

% Map 2
dV = repmat(sqrt(sum((XC2(:,[1 5 3],pInd) - XC2(:,[2 6 4],pInd)).^2)),[2,1]);
dV = dV(:);
subplot(2,3,6);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
xlabel('r_{n}');
ylabel('r_{n+1}');
axis([0 3.5 0 3.5]);

% Size and Save
fName = 'combined';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 5];
fig.PaperSize = [8.2 5];
saveas(fig, ['Figures\' fName], 'pdf');


%% Construct Super Stable Fixed Point
Xs = [-sqrt(3)/2 0 sqrt(3)/2;...
            -1/2 1 -1/2];
Us = [0 0 0;...
      .3 0 0];
conn = [1 1; 2 1; 3 1; 1 2; 2 2; 3 2];

% Instantaneous Motion
fig = figure(3); clf;
subplot(2,3,1);
visualize_conic(Xs,Us,[-2 2; -2 2], [100; 100], 0, 1, 1);
hold on;
plot(Xs(1,1:2), Xs(2,1:2), 'g--');
plot(Xs(1,2:3), Xs(2,2:3), 'c--');
hold off;
axis([-1.2 1.2 -1.0 1.4]);

% Motion
subplot(2,3,2);
Xu = [ 0    0.3434;...
      -1/2  0.3838];
LVal = sqrt(sum((Xs(:,conn(:,1)) - Xu(:,conn(:,2))).^2));
[XC, fC] = sim_motion3D_congrad(Xs, Xu, conn, LVal, 0.002, 500, [1 2], Xs(:,1:2) + [0 0; -.24 .24], 1);
[XC2, fC2] = sim_motion3D_congrad(Xs, Xu, conn, LVal, 0.002, 500, [1 2], Xs(:,1:2) + [0 0; .7 -.7], 1);
visualize_network(Xs,Xu,conn);
axis([-1.2 1.2 -1.0 1.4]);

% Diagram
subplot(2,3,3);
LV1 = squeeze(sqrt(sum((XC(:,1,:) - XC(:,2,:)).^2)));
LV2 = squeeze(sqrt(sum((XC(:,2,:) - XC(:,3,:)).^2)));
LV12 = squeeze(sqrt(sum((XC2(:,1,:) - XC2(:,2,:)).^2)));
LV22 = squeeze(sqrt(sum((XC2(:,2,:) - XC2(:,3,:)).^2)));
LV1 = [flipud(LV1); LV12];
LV2 = [flipud(LV2); LV22];
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(LV12(1), LV22(1), 'ko');
hold off;
xlabel('r_1');
ylabel('r_2');
set(gca,'XColor','g', 'YColor','c');

% Combine 1
subplot(2,3,4);
visualize_network(Xs,Xu,conn);
visualize_network([Xs(1,:) + 1; -Xs(2,:)+.8],[Xu(1,:) + 1; -Xu(2,:)+.8],conn);
axis([-1.2 2 -1.4 1.8]);

% Full Combined
subplot(2,3,[5,6]);
xSh = (Xs(1,3) - Xs(1,1))/2;
Xsa = [Xs [Xs(1,2)+2*xSh; Xs(2,2)]];
Xua = [Xu [Xu(1,:)+xSh; -Xu(2,:)+(Xs(2,2)+Xs(2,1))]];
conn = [1 1; 2 1; 3 1;...
        1 2; 2 2; 3 2;...
        2 3; 3 3; 4 3;...
        2 4; 3 4; 4 4];
[Xsa, Xua, conn] = tesselate_network(Xsa,Xua,conn,[xSh*2;0], [4,1]);
visualize_network(Xsa,Xua,conn);
axis([-1.2 7.5 -1.4 1.8]);

% Size and Save
fName = 'fixed_point';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 5];
fig.PaperSize = [8.2 5];
saveas(fig, ['Figures\' fName], 'pdf');


%% Simulate Super-Stable Fixed-Point
Xs = [-sqrt(3)/2 0 sqrt(3)/2;...
            -1/2 1 -1/2];
Xu = [ 0    0.3434;...
      -1/2  0.3838];
Xsa = [Xs [Xs(1,2)+2*xSh; Xs(2,2)]];
Xua = [Xu [Xu(1,:)+xSh; -Xu(2,:)+(Xs(2,2)+Xs(2,1))]];
conn = [1 1; 2 1; 3 1;...
        1 2; 2 2; 3 2;...
        2 3; 3 3; 4 3;...
        2 4; 3 4; 4 4];
[Xsa, Xua, conn] = tesselate_network(Xsa,Xua,conn,[xSh*2;0], [8,1]);

% Simulate
Xv0 = zeros(size([Xsa Xua]));
Xv0(2,1) = 1;
Xv0(2,2) = -1;
[XC, fC] = sim_motion(Xsa,Xua,conn,0.05,900,Xv0,0);

% Plot
fig = figure(4); clf;

subplot(5,10,[1:8]);
pInd = 1;
visualize_network(XC(:,1:size(Xsa,2),pInd),...
                  XC(:,[1:size(Xua,2)]+size(Xsa,2),pInd),conn);
axis([-1 14.2 -1 1.5]);
subplot(5,10,[9.5 10]);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
axis([min([LV1;LV2]) max([LV1;LV2]) min([LV1;LV2]) max([LV1;LV2])]);

subplot(5,10,[1:8]+10);
pInd = 100;
visualize_network(XC(:,1:size(Xsa,2),pInd),...
                  XC(:,[1:size(Xua,2)]+size(Xsa,2),pInd),conn);
axis([-1 14.2 -1 1.5]);
subplot(5,10,[9.5 10]+10);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
axis([min([LV1;LV2]) max([LV1;LV2]) min([LV1;LV2]) max([LV1;LV2])]);

subplot(5,10,[1:8]+20);
pInd = 200;
visualize_network(XC(:,1:size(Xsa,2),pInd),...
                  XC(:,[1:size(Xua,2)]+size(Xsa,2),pInd),conn);
axis([-1 14.2 -1 1.5]);
subplot(5,10,[9.5 10]+20);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
axis([min([LV1;LV2]) max([LV1;LV2]) min([LV1;LV2]) max([LV1;LV2])]);

subplot(5,10,[1:8]+30);
pInd = 300;
visualize_network(XC(:,1:size(Xsa,2),pInd),...
                  XC(:,[1:size(Xua,2)]+size(Xsa,2),pInd),conn);
axis([-1 14.2 -1 1.5]);
subplot(5,10,[9.5 10]+30);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
axis([min([LV1;LV2]) max([LV1;LV2]) min([LV1;LV2]) max([LV1;LV2])]);

subplot(5,10,[1:8]+40);
pInd = 392;
visualize_network(XC(:,1:size(Xsa,2),pInd),...
                  XC(:,[1:size(Xua,2)]+size(Xsa,2),pInd),conn);
axis([-1 14.2 -1 1.5]);
subplot(5,10,[9.5 10]+40);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
axis([min([LV1;LV2]) max([LV1;LV2]) min([LV1;LV2]) max([LV1;LV2])]);

% Size and Save
fName = 'fixed_point_sim';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 6.4];
fig.PaperSize = [8.2 6.4];
saveas(fig, ['Figures\' fName], 'pdf');


%% Aperiodic
fig = figure(5); clf;
Xs = [-sqrt(3)/1.5    0.0  sqrt(3);...
            -1/2  1.0       -1/2];

Xu = [ 0.551 -0.020;...
      -0.400 -1.200];
conn = [1 1; 2 1; 3 1; 1 2; 2 2; 3 2];

subplot(2,3,1);
visualize_network(Xs,Xu,conn);
axis([-2 3 -1.4 2.2]);

[XC, fC] = sim_motion(Xs,Xu,conn,0.001,1800,[Xs Xu],0);
[XC2, fC] = sim_motion(Xs,Xu,conn,0.001,1600,-[Xs Xu],0);

subplot(2,3,2);
LV1 = squeeze(sqrt(sum((XC(:,1,:) - XC(:,2,:)).^2)));
LV2 = squeeze(sqrt(sum((XC(:,2,:) - XC(:,3,:)).^2)));
LV12 = squeeze(sqrt(sum((XC2(:,1,:) - XC2(:,2,:)).^2)));
LV22 = squeeze(sqrt(sum((XC2(:,2,:) - XC2(:,3,:)).^2)));
LV1 = [flipud(LV1); LV12];
LV2 = [flipud(LV2); LV22];
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(LV12(1), LV22(1), 'ko');
hold off;
xlabel('r_n');
ylabel('r_{n+1}');
axis([0 3 0 3]);

subplot(2,3,3);
dA = abs(LV12 - LV22);
ii = find(dA == min(dA));
Xs = XC2(:,1:3,ii);
Xu = XC2(:,4:5,ii);
th1 = atan2d(Xs(2,1)-Xs(2,2), Xs(1,1)-Xs(1,2));
th2 = atan2d(Xs(2,3)-Xs(2,2), Xs(1,3)-Xs(1,2));
d1p = sqrt(sum((Xs(:,1)-Xs(:,2)).^2));
d2p = sqrt(sum((Xs(:,2)-Xs(:,3)).^2));
R = rotz((abs(180 + th1) - abs(th2))/2)^-1; R = R(1:2,1:2);
Xs = R*Xs;
Xu = R*Xu;
visualize_network(Xs,Xu,conn);
visualize_network([Xs(1,:) + 1.7; -Xs(2,:)+.8],[Xu(1,:) + 1.7; -Xu(2,:)+.8],conn);
axis([-1.5 3.5 -1.4 2.2]);

subplot(2,3,[4:6]);
xSh = (Xs(1,3) - Xs(1,1))/2;
Xsa = [Xs [Xs(1,2)+2*xSh; Xs(2,2)]];
Xua = [Xu [Xu(1,:)+xSh; -Xu(2,:)+(Xs(2,2)+Xs(2,1))]];
conn = [1 1; 2 1; 3 1;...
        1 2; 2 2; 3 2;...
        2 3; 3 3; 4 3;...
        2 4; 3 4; 4 4];
[Xsa, Xua, conn] = tesselate_network(Xsa,Xua,conn,[xSh*2;0], [5,1]);
visualize_network(Xsa,Xua,conn);
axis([-2 15 -1.5 2]);

% Size and Save
fName = 'aperiod';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 5];
fig.PaperSize = [8.2 5];
saveas(fig, ['Figures\' fName], 'pdf');


%% Simulate
xSh = (Xs(1,3) - Xs(1,1))/2;
Xsa = [Xs [Xs(1,2)+2*xSh; Xs(2,2)]];
Xua = [Xu [Xu(1,:)+xSh; -Xu(2,:)+(Xs(2,2)+Xs(2,1))]];
conn = [1 1; 2 1; 3 1;...
        1 2; 2 2; 3 2;...
        2 3; 3 3; 4 3;...
        2 4; 3 4; 4 4];
[Xsa, Xua, conn] = tesselate_network(Xsa,Xua,conn,[xSh*2;0], [10,1]);

LVal = sqrt(sum((Xsa(:,conn(:,1)) - Xua(:,conn(:,2))).^2));
[XC, fC] = sim_motion3D_congrad(Xsa, Xua, conn, LVal, 0.0002, 2000, [1 2], Xs(:,1:2) + [0 0; -.15 .15], 0);


%% Plot
fig = figure(6); clf;
subplot(2,2,1);
pInd = 1;
visualize_network(XC(:,[1:size(Xsa,2)], pInd), XC(:,[1:size(Xua,2)]+size(Xsa,2), pInd), conn);
axis([min(min(XC(1,:,:))), max(max(XC(1,:,:))), min(min(XC(2,:,:))), max(max(XC(2,:,:)))]);
subplot(2,2,2);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
xlabel('r_n');
ylabel('r_{n+1}');

subplot(2,2,3);
pInd = 200;
visualize_network(XC(:,[1:size(Xsa,2)], pInd), XC(:,[1:size(Xua,2)]+size(Xsa,2), pInd), conn);
axis([min(min(XC(1,:,:))), max(max(XC(1,:,:))), min(min(XC(2,:,:))), max(max(XC(2,:,:)))]);
subplot(2,2,4);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
xlabel('r_n');
ylabel('r_{n+1}');

% Size and Save
fName = 'aperiod_sim';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 7.5];
fig.PaperSize = [8.2 7.5];
saveas(fig, ['Figures\' fName], 'pdf');


fig = figure(7); clf;
subplot(2,2,1);
pInd = 400;
visualize_network(XC(:,[1:size(Xsa,2)], pInd), XC(:,[1:size(Xua,2)]+size(Xsa,2), pInd), conn);
axis([min(min(XC(1,:,:))), max(max(XC(1,:,:))), min(min(XC(2,:,:))), max(max(XC(2,:,:)))]);
subplot(2,2,2);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
xlabel('r_n');
ylabel('r_{n+1}');

subplot(2,2,3);
pInd = 600;
visualize_network(XC(:,[1:size(Xsa,2)], pInd), XC(:,[1:size(Xua,2)]+size(Xsa,2), pInd), conn);
axis([min(min(XC(1,:,:))), max(max(XC(1,:,:))), min(min(XC(2,:,:))), max(max(XC(2,:,:)))]);
subplot(2,2,4);
dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],pInd) - XC(:,[2:size(Xsa,2)],pInd)).^2)),[2,1]);
dV = dV(:);
plot(LV1, LV2, 'k-');
hold on;
plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
plot(dV(2:end-1), dV(3:end), 'r-');
plot([dV(1) dV(1)], [0 dV(3)], 'r-');
plot(dV(2),min([LV1;LV2]), 'ko');
hold off;
xlabel('r_n');
ylabel('r_{n+1}');

% Size and Save
fName = 'aperiod_sim2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 7.5];
fig.PaperSize = [8.2 7.5];
saveas(fig, ['Figures\' fName], 'pdf');



%%
fig = figure(8); clf;
fName = 'aperiodic_cobweb.gif';
nSV = 1;
dT = 0.03;
for i = 1:10:size(XC,3)
    dV = repmat(sqrt(sum((XC(:,[1:size(Xsa,2)-1],i) - XC(:,[2:size(Xsa,2)],i)).^2)),[2,1]);
    dV = [dV(:)];
    plot(LV1, LV2, 'k-');
    hold on;
    plot([min([LV1;LV2]) max([LV1;LV2])], [min([LV1;LV2]) max([LV1;LV2])], 'k--');
    plot(dV(2:end-1), dV(3:end), 'r-');
    plot([dV(1) dV(1)], [0 dV(3)], 'r-');
    plot(dV(2),min([LV1;LV2]), 'ko');
    hold off;
    xlabel('r_n');
    ylabel('r_{n+1}');

    % Capture the plot as an image
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if nSV == 1
      imwrite(imind,cm,fName,'gif', 'Loopcount',inf,'DelayTime',dT);
      nSV = 0;
    else
      imwrite(imind,cm,fName,'gif','WriteMode','append','DelayTime',dT);
    end
end
