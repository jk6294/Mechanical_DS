-1% Figure 1: Designing a Single Module
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');

% Plot Preliminaries
cArrow = [255 100 100]/255;
cDist = [200 200 200]/255;
ms = 1.5;
lw = 1.5;
lw_d = .5;
labX = -.2;
labY = 0.90;
netSX = .2;
netSc = .4;
labRowX = -0.35;
labColY = 1.1;

% Subplot Indices
nRow = [12 12 12 8];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [5 5 5 5];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(1); clf;
set(gcf, 'Renderer', 'painters', 'Position', [10 10 [23.7 17.1]*14.5/19], 'Units', 'centimeters'); 

% Global Annotation
annotation('arrow','HeadLength',8,'HeadWidth',12,'color',[.7 .7 .7],...
           'linewidth',4,'position',[.078 .765 0 -.07]);
annotation('arrow','HeadLength',8,'HeadWidth',12,'color',[.7 .7 .7],...
           'linewidth',4,'position',[.078 .53 0 -.07]);
annotation('line','linewidth',.5,'position',[.285 .25 0 .68],'color',[.9 .9 .9]);
annotation('line','linewidth',.5,'position',[.48 .25 0 .68],'color',[.9 .9 .9]);
annotation('line','linewidth',.5,'position',[.675 .25 0 .68],'color',[.9 .9 .9]);


%% a: Rigid Body Motion
subplot(NRow,NCol,cellM{1}); cla;
Xs1 = [-1 -1  1  1;...
       -1  1 -1  1];
conn1 = [1 3; 1 4; 2 3; 2 4];
visualize_network(Xs1,[],conn1);
hold on;
h = annotation('arrow');
set(h,'parent',gca,'position',[2.0 0 1 0], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
h = annotation('arrow');
set(h,'parent',gca,'position',[0 1.5 0 1], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
t = 90:-2:0;
aX = 1.5*cosd(t)+1.2; aY = 1.5*sind(t)+0.8;
plot(aX,aY,'color',cArrow, 'linewidth', 1);
h = annotation('arrow');
set(h,'parent',gca,'position',[aX(end) aY(end) diff(aX(end-1:end)) diff(aY(end-1:end))], 'HeadLength',5,'HeadWidth',5,'color',cArrow,'LineWidth',1);
hold off;
axis([-3.5 3.5 -3.5 3.5] + [netSX+.2 netSX+.2 -.2 -.2]);
text(labX,labY,'\textbf{a}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.02,labColY,'4-Bar Linkage','Units','Normalized','fontsize',10);
text(labRowX,.3,'Design','rotation',90,'Units','Normalized','fontsize',10);


%% b: Define Distances
subplot(NRow,NCol,cellM{2}); cla;
visualize_network(Xs1,[],conn1);
line_coordinates(Xs1(:,1:2),.8,.1,lw_d);
line_coordinates(Xs1(:,3:4),-.8,.1,lw_d);
text(Xs1(1,1)-2.0, mean(Xs1(2,1:2)), '$\mathrm{d_1}$', 'fontsize', 10);
text(Xs1(1,3)+1.2, mean(Xs1(2,1:2)), '$\mathrm{d_2}$', 'fontsize', 10);
axis([-3.5 3.5 -3.5 3.5] + [netSX+.2 netSX+.2 -.2 -.2]);
text(labX,labY,'\textbf{b}','Units','Normalized','fontsize',10,'fontweight','bold');
text(labRowX,.23,'Construct','rotation',90,'Units','Normalized','fontsize',10);


%% c: Motion Plot
subplot(NRow,NCol,cellM{3}); cla;
dX1 = [0 0 0 0; -1 1 1 -1];
[XCa,~] = sim_motion(Xs1,[],conn1,.01,140,dX1,0);   % Simulate
[XCb,~] = sim_motion(Xs1,[],conn1,.01,140,-dX1,0);  % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,3:4,:),[],2)).^2)));
diff1 = abs(d1-2) + abs(d2-2);       dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-1.143) + abs(d2-3.5); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-3.5) + abs(d2-1.143); dInd3 = find(diff3==min(diff3),1);
diff4 = abs(d1-1.5) + abs(d2-2.66); dInd4 = find(diff4==min(diff4),1);
diff5 = abs(d1-2.66) + abs(d2-1.5); dInd5 = find(diff5==min(diff5),1);
plot(d1,d2,'k-','linewidth',.5);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd4),d2(dInd4),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd5),d2(dInd5),'ko','markersize',ms,'linewidth',lw);
hold off;
visualize_network(XC(:,:,dInd1)/8+[d1(dInd1)+.3 d2(dInd1)+.4]',[],conn1,netSc);
visualize_network(XC(:,:,dInd2)/8+[d1(dInd2)+.5 d2(dInd2)+.05]',[],conn1,netSc);
visualize_network(XC(:,:,dInd3)/8+[d1(dInd3)+.1 d2(dInd3)+.4]',[],conn1,netSc);
visualize_network(XC(:,:,dInd4)/8+[d1(dInd4)+.4 d2(dInd4)+.25]',[],conn1,netSc);
visualize_network(XC(:,:,dInd5)/8+[d1(dInd5)+.2 d2(dInd5)+.35]',[],conn1,netSc);
set(gca,'visible',1,'XTick',[1.2 3.6],'YTick',[1.2 3.6],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);
axis([.9 3.9 .9 3.9]);
text(labX,labY+.2,'\textbf{c}','Units','Normalized','fontsize',10,'fontweight','bold');
text(labRowX,.18,'Trajectory','rotation',90,'Units','Normalized','fontsize',10);
text(.35,.9,'$\mathrm{d_2=d_1}$','Units','Normalized','fontsize',10,'color',[200 200 200]/255);
text(.03,.09,'$\mathrm{d_2=f(d_1)}$','Units','Normalized','fontsize',10)


%% d: Solution Space Example
subplot(NRow,NCol,cellM{5}); cla;
s = sqrt(3);
Xs2 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
dXs2 = [-s/2 0 s/2;...
        -1/2 1 -1/2]/2.5;
visualize_conic(Xs2,dXs2,[-2 2; -2 2],[100;100],3,1,1);
axis([-1 1 -1 1]*1.5 + [netSX netSX 0 0]);
text(labX,labY,'\textbf{d}','Units','Normalized','fontsize',10,'fontweight','bold');
text(.16,labColY,'Velocity','Units','Normalized','fontsize',10);


%% e: Constructed Network
subplot(NRow,NCol,cellM{6}); cla;
Xu0 = [-s/2 0;...
        1/2 -1];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,~] = construct_network(Xs2,dXs2,Xu0,conn2,0,0);
Xu2 = Xu2(1:2,:);
visualize_conic(Xs2,dXs2,[-2 2; -2 2],[100;100],0,0,0,.15);
visualize_network(Xs2,Xu2,conn2);
line_coordinates(Xs2(:,1:2),.25,.05,lw_d);
line_coordinates(Xs2(:,2:3),.25,.05,lw_d);
text(Xs2(1,1)-.70, mean(Xs2(2,1:2))+.0, '$\mathrm{d_1}$','fontsize',10);
text(Xs2(1,3)+.25, mean(Xs2(2,2:3))+.0, '$\mathrm{d_2}$','fontsize',10);
axis([-1 1 -1 1]*1.5 + [netSX netSX 0 0]);
text(labX,labY,'\textbf{e}','Units','Normalized','fontsize',10,'fontweight','bold');


%% f: Simulate
subplot(NRow,NCol,cellM{7}); cla;
[XC2,~] = sim_motion(Xs2,Xu2,conn2,.005,920,[Xs2,Xu2],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC2(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC2(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-s) + abs(d2-s);        dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-1.99) + abs(d2-2.27); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-1.94) + abs(d2-2.81);  dInd3 = find(diff3==min(diff3),1);
diff4 = abs(d1-1.45) + abs(d2-2.81);  dInd4 = find(diff4==min(diff4),1);
plot(d1,d2,'k-','linewidth',.5);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot([d1(dInd1)-.2,d1(dInd1)+.2],[d2(dInd1)-.2,d1(dInd1)+.2],'-','color',cArrow);
plot(d1(dInd1),d2(dInd1),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd1),d2(dInd1),'ko','markersize',2*ms);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd4),d2(dInd4),'ko','markersize',ms,'linewidth',lw);
hold off;
visualize_network(XC2(:,1:3,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',...
                  XC2(:,4:5,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',conn2,netSc);
visualize_network(XC2(:,1:3,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',...
                  XC2(:,4:5,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',conn2,netSc);
visualize_network(XC2(:,1:3,dInd3)/8+[d1(dInd3)+.3 d2(dInd3)+.0]',...
                  XC2(:,4:5,dInd3)/8+[d1(dInd3)+.3 d2(dInd3)+.0]',conn2,netSc);
visualize_network(XC2(:,1:3,dInd4)/8+[d1(dInd4)-.3 d2(dInd4)+.0]',...
                  XC2(:,4:5,dInd4)/8+[d1(dInd4)-.3 d2(dInd4)+.0]',conn2,netSc);
set(gca,'visible',1,'XTick',[1.2 2.8],'YTick',[1.2 2.8],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);
axis([.9 3.1 .9 3.1]);
text(labX,labY+.2,'\textbf{f}','Units','Normalized','fontsize',10,'fontweight','bold');


%% g: Finite Solution Space
subplot(NRow,NCol,cellM{9}); cla;
s = sqrt(3);
Xs30 = [-s/2 0 s/2;...
        -1/2 1 -1/2];
Xs3T = [-0.8  0.0  0.8;...
        -1.0  1.6 -1.0];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs30,Xs3T,[-1 1; -1 1]*2.0,[140;140],5,.7,.85);
axis([-1 1 -1 1]*1.7 + [netSX netSX .4 .4]);
text(labX,labY,'\textbf{g}','Units','Normalized','fontsize',10,'fontweight','bold');
text(-0.04,labColY,'Displacement','Units','Normalized','fontsize',10);


%% h: Finite Network
nW = .02;
xSh = .3;
ySh = xSh/sqrt(3);
subplot(NRow,NCol,cellM{10}); cla;
Xu0 = [-1.0  0.0;...
        1.6 -0.90];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs30,Xs3T,[-1 1; -1 1]*2,[100;100],0,0,0,.15);
[Xu3,~] = construct_network(Xs30,Xs3T,Xu0,conn3,0,1);
Xu3 = Xu3(1:2,:);
visualize_network(Xs30,Xu3,conn3);
line_coordinates(Xs30(:,1:2),.35,.05,lw_d);
line_coordinates(Xs30(:,2:3),.35,.05,lw_d);
text(Xs30(1,1)-.70, mean(Xs30(2,1:2))+.15, '$\mathrm{d_1}$','fontsize',10);
text(Xs30(1,3)+.15, mean(Xs30(2,2:3))+.15, '$\mathrm{d_2}$','fontsize',10);
axis([-1 1 -1 1]*1.7 + [netSX netSX .4 .4]);
text(labX,labY,'\textbf{h}','Units','Normalized','fontsize',10,'fontweight','bold');


%% i: Simulate
subplot(NRow,NCol,cellM{11}); cla;
[XC3,~] = sim_motion(Xs30,Xu3,conn3,.005,700,[Xs30,Xu3],0);   % Simulate
dM0 = sqrt(sum(diff(Xs30,1,2).^2));
dMT = sqrt(sum(diff(Xs3T,1,2).^2));
d1 = sqrt(squeeze(sum((diff(XC3(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC3(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-dM0(1)) + abs(d2-dM0(2));    dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-2.16) + abs(d2-1.93);        dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-2.50) + abs(d2-2.25);        dInd3 = find(diff3==min(diff3),1);
diff4 = abs(d1-dMT(1)) + abs(d2-dMT(2));    dInd4 = find(diff4==min(diff4),1);
plot(d1,d2,'k-','linewidth',.5);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd1),d2(dInd1),'ko','markersize',2*ms);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd4),d2(dInd4),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd4),d2(dInd4),'ko','markersize',2*ms);
hold off;
visualize_network(XC3(:,1:3,dInd1)/8+[d1(dInd1)+.1 d2(dInd1)-.4]',...
                  XC3(:,4:5,dInd1)/8+[d1(dInd1)+.1 d2(dInd1)-.4]',conn3,netSc);
visualize_network(XC3(:,1:3,dInd2)/8+[d1(dInd2)+.15 d2(dInd2)-.35]',...
                  XC3(:,4:5,dInd2)/8+[d1(dInd2)+.15 d2(dInd2)-.35]',conn3,netSc);
visualize_network(XC3(:,1:3,dInd3)/8+[d1(dInd3)+.2 d2(dInd3)-.3]',...
                  XC3(:,4:5,dInd3)/8+[d1(dInd3)+.2 d2(dInd3)-.3]',conn3,netSc);
visualize_network(XC3(:,1:3,dInd4)/8+[d1(dInd4)+.2 d2(dInd4)-.3]',...
                  XC3(:,4:5,dInd4)/8+[d1(dInd4)+.2 d2(dInd4)-.3]',conn3,netSc);
set(gca,'visible',1,'XTick',[1.2 2.8],'YTick',[1.2 2.8],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);
axis([.9 3.1 .9 3.1]);
text(labX,labY+.2,'\textbf{i}','Units','Normalized','fontsize',10,'fontweight','bold');


%% j: Infinitesimal and Finite
subplot(NRow,NCol,cellM{13}); cla;
Xs40 = [-s/2 0 s/2;...
        -1/2 1 -1/2];
Xs4T = Xs40*1.7;
dXs4 = [-1.0  0.0  0.0;...
         0.0  0.0  0.0]/2;
visualize_conic_finite(Xs40,Xs4T,[-2 2;-2 2],[100;100],0,.7,0);
visualize_conic(Xs40,dXs4,[-2 2;-2 2],[100;100],0,1,0);
axis([-1 1 -1 1]*2 + [netSX netSX .27 .27]);
text(labX,labY,'\textbf{j}','Units','Normalized','fontsize',10,'fontweight','bold');
text(-0.31,labColY,'Velocity + Displacement','Units','Normalized','fontsize',10);


%% k: Construct Infinitesimal and Finite Network
subplot(NRow,NCol,cellM{14}); cla;
% Xu4 = [ 0.950 -1.690;...
%        -1.394 -0.080];
Xu4 = [-0.86 -0.86;...
       -1.45  1.47];
conn4 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs40,Xs4T,[-2 2;-2 2],[100;100],0,0,0,.15);
visualize_conic(Xs40,dXs4,[-2 2;-2 2],[100;100],0,0,0,.2);
visualize_network(Xs40,Xu4,conn4);
line_coordinates(Xs40(:,1:2),.35,.05,lw_d);
line_coordinates(Xs40(:,2:3),.35,.05,lw_d);
text(Xs40(1,1)-.75, mean(Xs40(2,1:2))+.1, '$\mathrm{d_1}$','fontsize',10);
text(Xs40(1,3)+.15, mean(Xs40(2,2:3))+.1, '$\mathrm{d_2}$','fontsize',10);
axis([-1 1 -1 1]*2 + [netSX netSX .27 .27]);
text(labX,labY,'\textbf{k}','Units','Normalized','fontsize',10,'fontweight','bold');


%% l: Simulate Infinitesimal and Finite
subplot(NRow,NCol,cellM{15}); cla;
dM0 = sqrt(sum(diff(Xs40,1,2).^2));
dMT = sqrt(sum(diff(Xs4T,1,2).^2));
[XC4,~] = sim_motion(Xs40,Xu4,conn4,.005,400,[Xs40,Xu4],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC4(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC4(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-dM0(1)) + abs(d2-dM0(2));    dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-2.27) + abs(d2-1.88);        dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-2.72) + abs(d2-2.26);        dInd3 = find(diff3==min(diff3),1);
diff4 = abs(d1-dMT(1)) + abs(d2-dMT(2));    dInd4 = find(diff4==min(diff4),1);
plot(d1,d2,'k-','linewidth',.5);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot([d1(dInd1)-.3,d1(dInd1)+.3],[d2(dInd1),d1(dInd1)],'-','color',cArrow);
plot(d1(dInd1),d2(dInd1),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd1),d2(dInd1),'ko','markersize',2*ms);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd4),d2(dInd4),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd4),d2(dInd4),'ko','markersize',2*ms);
hold off;
visualize_network(XC4(:,1:3,dInd1)/8+[d1(dInd1)-.2 d2(dInd1)+.25]',...
                  XC4(:,4:5,dInd1)/8+[d1(dInd1)-.2 d2(dInd1)+.25]',conn4,netSc);
visualize_network(XC4(:,1:3,dInd2)/8+[d1(dInd2)-.3 d2(dInd2)+.2]',...
                  XC4(:,4:5,dInd2)/8+[d1(dInd2)-.3 d2(dInd2)+.2]',conn4,netSc);
visualize_network(XC4(:,1:3,dInd3)/8+[d1(dInd3)-.3 d2(dInd3)+.1]',...
                  XC4(:,4:5,dInd3)/8+[d1(dInd3)-.3 d2(dInd3)+.1]',conn4,netSc);
visualize_network(XC4(:,1:3,dInd4)/8+[d1(dInd4)-.3 d2(dInd4)-.15]',...
                  XC4(:,4:5,dInd4)/8+[d1(dInd4)-.3 d2(dInd4)-.15]',conn4,netSc);
set(gca,'visible',1,'XTick',[1.2 2.8],'YTick',[1.2 2.8],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);
axis([.9 3.1 .9 3.1]);
text(labX,labY+.2,'\textbf{l}','Units','Normalized','fontsize',10,'fontweight','bold');


%% Legend
cMS = [cellM{[4; 8; 12; 16]}];
subplot(NRow,NCol,cMS(:)); cla;
set(gca,'visible',0);
ms = 3;
mhs = .8;
bw = 0.5;
lw = 1.2;                        % Line Width
bs = .4;
ws = .08;
C_SN = [255 100 100]/255;
C_UN = [100 100 255]/255;
C_SS = [100 100 255;...         % Color of Solution Space
        100 200 255]/255;    
C_E  = [0 0 0 .5];              % Color of Edge
LW_SA = 1;
LW_SS = 1;
tShX = .2;
tShY = ws/2-.04;
hold on;
% Specified Node
plot(1, bs, 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN);
plot(1, bs, 'ko', 'linewidth', bw, 'markersize', ms*2);
text(1+tShX, bs + tShY, 'designed node', 'fontsize', 10);
% Variable Node
plot(1, 0, 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN);
plot(1, 0, 'ko', 'linewidth', bw, 'markersize', ms*2);
text(1+tShX, tShY, 'variable node', 'fontsize', 10);
% Designed Motions
quiver(4, bs+ws/2, .6, 0, 0, 'filled', 'MaxHeadSize', mhs, 'linewidth', LW_SA, 'color', cArrow);
quiver(4, bs+ws/2, .6, 0, 0, 'filled', 'MaxHeadSize', mhs, 'linewidth', LW_SA/4, 'color', [1 1 1]);
quiver(4, bs-ws/2, .6, 0, 0, 'MaxHeadSize', mhs, 'linewidth', LW_SA*2/3, 'color', cArrow);
text(4.6+tShX, bs + tShY, 'designed motion', 'fontsize', 10);
% Variable Motions
quiver(4, ws/2, .6, 0, 0, 'filled', 'MaxHeadSize', mhs, 'linewidth', LW_SA, 'color', C_SS(1,:));
quiver(4, ws/2, .6, 0, 0, 'filled', 'MaxHeadSize', mhs, 'linewidth', LW_SA/4, 'color', [1 1 1]);
quiver(4, -ws/2, .6, 0, 0, 'MaxHeadSize', mhs, 'linewidth', LW_SA*2/3, 'color', C_SS(1,:));
text(4.6+tShX, tShY, 'variable motion', 'fontsize', 10);
% Instantaneous Solution Space
line([7.5 8], [ws ws], 'color', C_SS(1,:), 'linewidth', LW_SS);
line([7.5 8], [ws ws], 'color', [1 1 1], 'linewidth', LW_SS/4);
line([7.5 8], [0 0], 'color', C_SS(1,:), 'linewidth', LW_SS);
line([7.5 8], -[ws ws], 'color', C_SS(2,:), 'linewidth', LW_SS);
text(8.0+tShX, tShY, 'variable positions', 'fontsize', 10);
% Design
line([7.5 8], [bs bs], 'color', cArrow, 'linewidth', .4);
plot(7.75,bs,'o','markersize',ms/2,'linewidth',ms/2,'color',cArrow);
plot(7.75,bs,'ko','markersize',ms);
text(8+tShX, bs+tShY, 'designed distance', 'fontsize', 10);
hold off;
axis([0 10 -1*bs-ws 1.8*bs+ws]+[.7 .7 .2 .2]);


%% Size and Save Figure
fName = 'figure1';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.6 -2.4 23.7 17.1]*13.5/19;
fig.PaperSize = [19 13.5]*14/19;
saveas(fig, ['Figures/' fName], 'pdf');