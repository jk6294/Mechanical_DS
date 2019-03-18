%% Prepare Space
clear; clc;


%% Simple Example: 4 Nodes
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