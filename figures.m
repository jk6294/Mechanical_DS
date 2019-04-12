%% Figure 1: Designing a Single Module
% Prepare Space
clear; clc;

%% Plot Preliminaries
cArrow = [76 187 23]/255;
cDist = [200 200 200]/255;

% Subplot Indices
nRow = [8 8 8];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [4 4 4 4];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
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
text(Xs1(1,1)-2.2, mean(Xs1(2,1:2)), 'd_1', 'fontsize', 10);
text(Xs1(1,3)+1.2, mean(Xs1(2,1:2)), 'd_2', 'fontsize', 10);
hold off;
axis([-3.5 3.5 -3.5 3.5] + [0 0 .6 .6]);

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
visualize_network(XC(:,:,dInd1)/8+[d1(dInd1)+.5 d2(dInd1)+.3]',[],conn1,.3);
visualize_network(XC(:,:,dInd2)/8+[d1(dInd2)+.5 d2(dInd2)+0]',[],conn1,.3);
visualize_network(XC(:,:,dInd3)/8+[d1(dInd3)+0 d2(dInd3)+.5]',[],conn1,.3);
set(gca,'visible',1,'XTick',[1 2.5 4],'YTick',[1 2.5 4],'fontsize',10,...
        'xticklabel',{'1','d_1','4'},'yticklabel',{'1','d_2','4'});
axis([1 4 1 4]);

% d: Solution Space Example
subplot(NRow,NCol,cellM{4});
s = sqrt(3);
Xs2 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
dXs2 = [-s/2 0 s/2;...
        -1/2 1 -1/2]/2.5;
visualize_conic(Xs2,dXs2,[-2 2; -2 2],[100;100],3,1,1);
axis([-1 1 -1 1]*1.5 + [0 0 0 0]);

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
text(Xs2(1,1)-.9, mean(Xs2(2,1:2))+.4, 'd_1','fontsize',10);
text(Xs2(1,3)+.3, mean(Xs2(2,2:3))+.4, 'd_2','fontsize',10);
hold off;
axis([-1 1 -1 1]*1.5 + [0 0 0 0]);

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
plot(d1(dInd1),d2(dInd1),'o','markersize',5,'color',cArrow);
plot(d1(dInd2),d2(dInd2),'bo','markersize',5);
plot(d1(dInd3),d2(dInd3),'bo','markersize',5);
plot([d1(dInd1)-.2,d1(dInd1)+.2],[d2(dInd1)-.2,d1(dInd1)+.2],'-','color',cArrow);
hold off;
visualize_network(XC2(:,1:3,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',...
                  XC2(:,4:5,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',conn2,.3);
visualize_network(XC2(:,1:3,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',...
                  XC2(:,4:5,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',conn2,.3);
visualize_network(XC2(:,1:3,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)+.3]',...
                  XC2(:,4:5,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)+.3]',conn2,.3);
set(gca,'visible',1,'XTick',[1 2 3],'YTick',[1 2 3],'fontsize',10,...
        'xticklabel',{'1','d_1','3'},'yticklabel',{'1','d_2','3'});
axis([.9 3.1 .9 3.1]);

% g: Finite Solution Space
subplot(NRow,NCol,cellM{7})
s = sqrt(3);
Xs30 = [-s/2 0 s/2;...
        -1/2 1 -1/2];
Xs3T = [-0.5  0.0  0.5;...
        -1.0  1.5 -1.0];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs30,Xs3T,[-1 1; -1 1]*1.8,[100;100],3,.7,.85);
axis([-1 1 -1 1]*1.5 + [0 0 .4 .4]);

% h: Finite Network
subplot(NRow,NCol,cellM{8})
Xu0 = [-0.38  0.00;...
        1.58 -0.75];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu3,~] = construct_network(Xs30,Xs3T,Xu0,conn3,1,1);
Xu3 = Xu3(1:2,:);
axis([-1 1 -1 1]*1.5 + [0 0 .4 .4]);

% i: Simulate
subplot(NRow,NCol,cellM{9});
[XC3,~] = sim_motion(Xs30,Xu3,conn3,.005,500,[Xs30,Xu3],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC3(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC3(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-s) + abs(d2-s);       dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-2.27) + abs(d2-2.00); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-2.55) + abs(d2-2.55); dInd3 = find(diff3==min(diff3),1);
plot(d1,d2);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'o','markersize',5,'color',cArrow);
plot(d1(dInd2),d2(dInd2),'bo','markersize',5);
plot(d1(dInd3),d2(dInd3),'o','markersize',5,'color',cArrow);
hold off;
visualize_network(XC3(:,1:3,dInd1)/8+[d1(dInd1)+.0 d2(dInd1)-.4]',...
                  XC3(:,4:5,dInd1)/8+[d1(dInd1)+.0 d2(dInd1)-.4]',conn3,.3);
visualize_network(XC3(:,1:3,dInd2)/8+[d1(dInd2)+.2 d2(dInd2)-.2]',...
                  XC3(:,4:5,dInd2)/8+[d1(dInd2)+.2 d2(dInd2)-.2]',conn3,.3);
visualize_network(XC3(:,1:3,dInd3)/8+[d1(dInd3)+.3 d2(dInd3)-.1]',...
                  XC3(:,4:5,dInd3)/8+[d1(dInd3)+.3 d2(dInd3)-.1]',conn3,.3);
set(gca,'visible',1,'XTick',[1 2 3],'YTick',[1 2 3],'fontsize',10,...
        'xticklabel',{'1','d_1','3'},'yticklabel',{'1','d_2','3'});
axis([.9 3.1 .9 3.1]);

% j: Infinitesimal and Finite
subplot(NRow,NCol,cellM{10});
Xs40 = [-s/2 0 s/2;...
        -1/2 1 -1/2];
Xs4T = Xs40*1.7;
dXs4 = [0 0 .5;...
        0 0 -1]/2;
visualize_conic_finite(Xs40,Xs4T,[-2 2;-2 2],[100;100],0,.7,0);
visualize_conic(Xs40,dXs4,[-2 2;-2 2],[100;100],0,1,0);
axis([-1 1 -1 1]*2);

% k: Construct Infinitesimal and Finite Network
subplot(NRow,NCol,cellM{11});
Xu4 = [-0.950  1.690;...
       -1.394 -0.080];
conn4 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_network(Xs40,Xu4,conn4);
axis([-1 1 -1 1]*2);

% l: Simulate Infinitesimal and Finite
subplot(NRow,NCol,cellM{12});
[XC4,~] = sim_motion(Xs40,Xu4,conn4,.005,400,[Xs40,Xu4],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC4(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC4(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-s) + abs(d2-s);       dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-2.00) + abs(d2-2.62); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-2.93) + abs(d2-2.93); dInd3 = find(diff3==min(diff3),1);
plot(d1,d2);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'o','markersize',5,'color',cArrow);
plot(d1(dInd2),d2(dInd2),'bo','markersize',5);
plot(d1(dInd3),d2(dInd3),'o','markersize',5,'color',cArrow);
plot([d1(dInd1),d1(dInd1)],[d2(dInd1)-.3,d1(dInd1)+.3],'-','color',cArrow);
hold off;
visualize_network(XC4(:,1:3,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',...
                  XC4(:,4:5,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',conn4,.3);
visualize_network(XC4(:,1:3,dInd2)/8+[d1(dInd2)-.4 d2(dInd2)+0]',...
                  XC4(:,4:5,dInd2)/8+[d1(dInd2)-.4 d2(dInd2)+0]',conn4,.3);
visualize_network(XC4(:,1:3,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)-.4]',...
                  XC4(:,4:5,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)-.4]',conn4,.3);
set(gca,'visible',1,'XTick',[1 2 3],'YTick',[1 2 3],'fontsize',10,...
        'xticklabel',{'1','d_1','3'},'yticklabel',{'1','d_2','3'});
axis([.9 3.1 .9 3.1]);


% Size and Save Figure
fName = 'figure1';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.7 -0.8 18 10.4];
fig.PaperSize = [14 8.8];
saveas(fig, ['Figures/' fName], 'pdf');



%% Figure 2
% Prepare Space
clear; clc;

% Plot
cCob = [200 200 100]/255;

% Subplot Indices
nRow = [6 6 12];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [4 4 4 4];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(2); clf;

% Two Modules
subplot(NRow,NCol,cellM{1});
tSx = .7;
tSy = .1;
Xs1 = [-2 -2  2  2;...
       -1  1 -1  1];
conn1 = [1 3; 1 4; 2 3; 2 4];
visualize_network(Xs1-[3;0],[],conn1);
visualize_network(Xs1+[3;0],[],conn1);
text(-5-tSx, tSy, 'd_1','fontsize',10);
text(-1-tSx, tSy, 'd_2','fontsize',10);
text(1-tSx, tSy, 'd_2*','fontsize',10);
text(5-tSx, tSy, 'd_3','fontsize',10);
axis([-5 5 -2 2]);

subplot(NRow,NCol,cellM{2});
[Xs1a,conn1a] = tesselate_network(Xs1,conn1,[4;0],[2;1]);
visualize_network(Xs1a,[],conn1a);
text(-2-tSx, tSy, 'd_1','fontsize',10);
text(2-tSx, tSy, 'd_2','fontsize',10);
text(6-tSx, tSy, 'd_3','fontsize',10);
axis([-5 5 -2 2]+[2 2 0 0]);

% Cobweb Plot
subplot(NRow,NCol,cellM{3});
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
plot(d1,d2);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
% Cobweb 1
pInd = 140;
dP = [D1(:,pInd)';D1(:,pInd)']; dP = dP(:);
line(dP,[0;dP(1:end-1)],'color',cCob);
visualize_network(XC1a(:,:,pInd)/8+[2.5;3.6],[],conn1a,.3);
% Cobweb 2
pInd = 20;
dP = [D1(:,pInd)';D1(:,pInd)']; dP = dP(:);
line(dP,[0;dP(1:end-1)],'color',cCob);
visualize_network(XC1a(:,:,pInd)/8+[2.2;2.4],[],conn1a,.3);
hold off;
set(gca,'visible',1,'XTick',[1 2.5 4],'YTick',[1 2.5 4],'fontsize',10,...
        'xticklabel',{'1','d_k','4'},'yticklabel',{'1','d_{k+1}','4'});
axis([1 4 1 4]);




subplot(NRow,NCol,cellM{4});

subplot(NRow,NCol,cellM{5});

subplot(NRow,NCol,cellM{6});

subplot(NRow,NCol,cellM{7});




subplot(NRow,NCol,cellM{12});




% Size and Save Figure
fName = 'figure2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-2.3 -0.6 24.7 9];
fig.PaperSize = [19 8];
saveas(fig, ['Figures/' fName], 'pdf');





