%% Figure 1: Designing a Single Module
% Prepare Space
clear; clc;

% Plot Preliminaries
cArrow = [100 235 100]/255;
cDist = [200 200 200]/255;

% Subplot Indices
nRow = [4 4 7];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [5 5 5 5];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nCol)
    for j = 1:length(nRow)
        cP = figM((nR(j)+1):(nR(j+1)-1), (nC(i)+1):(nC(i+1)-1));
        cellM(j+(i-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(1); clf;

% Rigid Body Motion
subplot(NRow,NCol,cellM{1});
Xs1 = [-2 -2  2  2;...
       -1  1 -1  1];
conn1 = [1 3; 1 4; 2 3; 2 4];
visualize_network(Xs1,[],conn1);
hold on;
h = annotation('arrow');
set(h,'parent',gca,'position',[2.5 0 1 0], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
h = annotation('arrow');
set(h,'parent',gca,'position',[0 1.2 0 1], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
t = 80:-2:10;
aX = cosd(t)+2.2; aY = sind(t)+1;
plot(aX,aY,'color',cArrow, 'linewidth', 1);
h = annotation('arrow');
set(h,'parent',gca,'position',[aX(end) aY(end) diff(aX(end-1:end)) diff(aY(end-1:end))], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
hold off;
axis([-3.5 3.5 -1.7 1.7] + [0 0 .6 .6]);

% Define Distances
subplot(NRow,NCol,cellM{2});
hold on;
line(Xs1(1,1:2)-.8,Xs1(2,1:2),'color','k','LineWidth',.5);
line(Xs1(1,1:2)-[.9 .7],[Xs1(2,1) Xs1(2,1)],'color','k','LineWidth',.5);
line(Xs1(1,1:2)-[.9 .7],[Xs1(2,2) Xs1(2,2)],'color','k','LineWidth',.5);
line(Xs1(1,3:4)+.8,Xs1(2,3:4),'color','k','LineWidth',.5);
line(Xs1(1,3:4)+[.9 .7],[Xs1(2,3) Xs1(2,3)],'color','k','LineWidth',.5);
line(Xs1(1,3:4)+[.9 .7],[Xs1(2,4) Xs1(2,4)],'color','k','LineWidth',.5);
visualize_network(Xs1,[],conn1);
text(Xs1(1,1)-1.8, mean(Xs1(2,1:2)), 'd_1');
text(Xs1(1,3)+1, mean(Xs1(2,1:2)), 'd_2');
hold off;
axis([-3.5 3.5 -1.7 1.7] + [0 0 0 0]);

% Motion Plot
subplot(NRow,NCol,cellM{3});
dX1 = [0 0 0 0; -1 1 1 -1];
[XCa,~] = sim_motion(Xs1,[],conn1,.01,170,dX1,0);   % Simulate
[XCb,~] = sim_motion(Xs1,[],conn1,.01,170,-dX1,0);  % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,3:4,:),[],2)).^2)));
plot(d1,d2);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
hold off;
visualize_network(XCb(:,:,end)/8+[1.6,3.6]',[],conn1,.3);
visualize_network(XCa(:,:,end)/8+[3.6,1.6]',[],conn1,.3);
visualize_network(Xs1/8+[2.3 2.3]',[],conn1,.3);
set(gca,'visible',1,'XTick',[1 4],'YTick',[1 4]);
xlabel('d_1');
ylabel('d_2');
axis([1 4 1 4]);

subplot(NRow,NCol,cellM{4});
s = sqrt(3);
Xs2 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
dXs2 = [];



% Size and Save Figure
fName = 'Test';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 3.6];
fig.PaperSize = [8.2 3.6];
saveas(fig, ['Figures/' fName], 'pdf');


