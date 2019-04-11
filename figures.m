%% Figure 1: Designing a Single Module
% Prepare Space
clear; clc;

% Plot Preliminaries
cArrow = [100 235 100]/255;
cDist = [200 200 200]/255;

% Subplot Indices
nRow = [8 8 8];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
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

% a: Rigid Body Motion
subplot(NRow,NCol,cellM{1});
Xs1 = [-2 -2  2  2;...
       -1  1 -1  1];
conn1 = [1 3; 1 4; 2 3; 2 4];
visualize_network(Xs1,[],conn1);
hold on;
h = annotation('arrow');
set(h,'parent',gca,'position',[2.5 0 1 0], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
h = annotation('arrow');
set(h,'parent',gca,'position',[0 1.5 0 1], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
t = 90:-2:0;
aX = 1.5*cosd(t)+1.7; aY = 1.5*sind(t)+0.8;
plot(aX,aY,'color',cArrow, 'linewidth', 1);
h = annotation('arrow');
set(h,'parent',gca,'position',[aX(end) aY(end) diff(aX(end-1:end)) diff(aY(end-1:end))], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
hold off;
axis([-3.5 3.5 -3.5 3.5] + [0 0 .6 .6]);

% b: Define Distances
nW = .1;
subplot(NRow,NCol,cellM{2});
hold on;
line(Xs1(1,1:2)-.8,Xs1(2,1:2),'color','k','LineWidth',.5);
line(Xs1(1,1:2)-.8 + [nW -nW],[Xs1(2,1) Xs1(2,1)],'color','k','LineWidth',.5);
line(Xs1(1,1:2)-.8 + [nW -nW],[Xs1(2,2) Xs1(2,2)],'color','k','LineWidth',.5);
line(Xs1(1,3:4)+.8,Xs1(2,3:4),'color','k','LineWidth',.5);
line(Xs1(1,3:4)+.8 + [nW -nW],[Xs1(2,3) Xs1(2,3)],'color','k','LineWidth',.5);
line(Xs1(1,3:4)+.8 + [nW -nW],[Xs1(2,4) Xs1(2,4)],'color','k','LineWidth',.5);
visualize_network(Xs1,[],conn1);
text(Xs1(1,1)-1.8, mean(Xs1(2,1:2)), 'd_1');
text(Xs1(1,3)+1, mean(Xs1(2,1:2)), 'd_2');
hold off;
axis([-3.5 3.5 -3.5 3.5] - [0 0 .6 .6]);

% c: Motion Plot
subplot(NRow,NCol,cellM{3});
dX1 = [0 0 0 0; -1 1 1 -1];
[XCa,~] = sim_motion(Xs1,[],conn1,.01,170,dX1,0);   % Simulate
[XCb,~] = sim_motion(Xs1,[],conn1,.01,170,-dX1,0);  % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,3:4,:),[],2)).^2)));
diff1 = abs(d1-2) + abs(d2-2);       dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-1.143) + abs(d2-3.5); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-3.5) + abs(d2-1.143); dInd3 = find(diff3==min(diff3),1);
plot(d1,d2);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'bo','markersize',5);
plot(d1(dInd2),d2(dInd2),'bo','markersize',5);
plot(d1(dInd3),d2(dInd3),'bo','markersize',5);
hold off;
visualize_network(XC(:,:,dInd1)/8+[d1(dInd1)+.3 d2(dInd1)+.3]',[],conn1,.3);
visualize_network(XC(:,:,dInd2)/8+[d1(dInd2)+.5 d2(dInd2)+0]',[],conn1,.3);
visualize_network(XC(:,:,dInd3)/8+[d1(dInd3)+0 d2(dInd3)+.5]',[],conn1,.3);
set(gca,'visible',1,'XTick',[1 4],'YTick',[1 4]);
xlabel('d_1');
ylabel('d_2');
axis([1 4 1 4]);

% d: Solution Space Example
subplot(NRow,NCol,cellM{4});
s = sqrt(3);
Xs2 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
dXs2 = [-s/2 0 s/2;...
        -1/2 1 -1/2]/2.5;
visualize_conic(Xs2,dXs2,[-2 2; -2 2],[100;100],3,1,1);
axis([-1 1 -1 1]*1.5 + [0 0 .2 .2]);

% e: Constructed Network
nW = .02;
xSh = .6;
ySh = xSh/sqrt(3);
subplot(NRow,NCol,cellM{5});
Xu0 = [-s/2 0;...
        1/2 -1];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,~] = construct_network(Xs2,dXs2,Xu0,conn2,1,0);
Xu2 = Xu2(1:2,:);
hold on;
line(Xs2(1,1:2)-xSh,Xs2(2,1:2)+ySh,'color','k','LineWidth',.5);
line(Xs2(1,1)-xSh+nW*[s -s], Xs2(2,1)+ySh+nW*[-1 1],'color','k','LineWidth',.5);
line(Xs2(1,2)-xSh+nW*[s -s], Xs2(2,2)+ySh+nW*[-1 1],'color','k','LineWidth',.5);
line(Xs2(1,2:3)+xSh,Xs2(2,2:3)+ySh,'color','k','LineWidth',.5);
line(Xs2(1,2)+xSh+nW*[s -s], Xs2(2,2)+ySh+nW*[1 -1],'color','k','LineWidth',.5);
line(Xs2(1,3)+xSh+nW*[s -s], Xs2(2,3)+ySh+nW*[1 -1],'color','k','LineWidth',.5);
text(Xs2(1,1)-.7, mean(Xs2(2,1:2))+.4, 'd_1');
text(Xs2(1,3)+.25, mean(Xs2(2,2:3))+.4, 'd_2');
hold off;
axis([-1 1 -1 1]*1.5 + [0 0 .2 .2]);

% f: Simulate
subplot(NRow,NCol,cellM{6});
[XC2,~] = sim_motion(Xs2,Xu2,conn2,.005,920,[Xs2,Xu2],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC2(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC2(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-s) + abs(d2-s);       dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-1.96) + abs(d2-2.74); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-1.24) + abs(d2-2.48); dInd3 = find(diff3==min(diff3),1);
plot(d1,d2);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'bo','markersize',5);
plot(d1(dInd2),d2(dInd2),'bo','markersize',5);
plot(d1(dInd3),d2(dInd3),'bo','markersize',5);
hold off;
visualize_network(XC2(:,1:3,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',...
                  XC2(:,4:5,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',conn2,.3);
visualize_network(XC2(:,1:3,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',...
                  XC2(:,4:5,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',conn2,.3);
visualize_network(XC2(:,1:3,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)+.3]',...
                  XC2(:,4:5,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)+.3]',conn2,.3);
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3]);
xlabel('d_1');
ylabel('d_2');
axis([.9 3.1 .9 3.1]);





% Size and Save Figure
fName = 'figure1';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 5];
fig.PaperSize = [8.2 5];
saveas(fig, ['Figures/' fName], 'pdf');


