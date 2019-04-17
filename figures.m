%% Figure 1: Designing a Single Module
% Prepare Space
clear; clc;

% Plot Preliminaries
cArrow = [255 100 100]/255;
cDist = [200 200 200]/255;
ms = 1.5;
lw = 1.5;
lw_d = .3;
labX = -.3;
labY = 0.94;
netSX = .2;
netSc = .3;
labRowX = -0.6;
labColY = 1.2;

% Subplot Indices
nRow = [8 8 8];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [3 3 3 3];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(1); clf;

% Global Annotation
annotation('arrow','HeadLength',8,'HeadWidth',12,'color',[.7 .7 .7],...
           'linewidth',4,'position',[.063 .73 0 -.07]);
annotation('arrow','HeadLength',8,'HeadWidth',12,'color',[.7 .7 .7],...
           'linewidth',4,'position',[.063 .435 0 -.07]);
annotation('line','linewidth',.2,'position',[.27 .1 0 .83],'color',[.9 .9 .9]);
annotation('line','linewidth',.2,'position',[.47 .1 0 .83],'color',[.9 .9 .9]);
annotation('line','linewidth',.2,'position',[.67 .1 0 .83],'color',[.9 .9 .9]);
annotation('line','linewidth',.2,'position',[.85 .1 0 .83],'color',[.9 .9 .9]);


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
axis([-3.5 3.5 -3.5 3.5] + [netSX netSX -.2 -.2]);
text(labX,labY,'a','Units','Normalized','fontsize',10,'fontweight','bold');
text(.05,labColY,'4-bar linkage','Units','Normalized','fontsize',10);
text(labRowX,.3,'Design','rotation',90,'Units','Normalized','fontsize',10);


% b: Define Distances
nW = .1;
subplot(NRow,NCol,cellM{2});
hold on;
line(Xs1(1,1:2)-.8,Xs1(2,1:2),'color','k','LineWidth',.5);
line(Xs1(1,1:2)-.8 + [nW -nW],[Xs1(2,1) Xs1(2,1)],'color','k','LineWidth',lw_d);
line(Xs1(1,1:2)-.8 + [nW -nW],[Xs1(2,2) Xs1(2,2)],'color','k','LineWidth',lw_d);
line(Xs1(1,3:4)+.8,Xs1(2,3:4),'color','k','LineWidth',lw_d);
line(Xs1(1,3:4)+.8 + [nW -nW],[Xs1(2,3) Xs1(2,3)],'color','k','LineWidth',lw_d);
line(Xs1(1,3:4)+.8 + [nW -nW],[Xs1(2,4) Xs1(2,4)],'color','k','LineWidth',lw_d);
visualize_network(Xs1,[],conn1);
text(Xs1(1,1)-2.2, mean(Xs1(2,1:2)), 'd_1', 'fontsize', 10);
text(Xs1(1,3)+1.2, mean(Xs1(2,1:2)), 'd_2', 'fontsize', 10);
hold off;
axis([-3.5 3.5 -3.5 3.5] + [netSX netSX -.2 -.2]);
text(labX,labY,'b','Units','Normalized','fontsize',10,'fontweight','bold');
text(labRowX,.2,'Construct','rotation',90,'Units','Normalized','fontsize',10);


% c: Motion Plot
subplot(NRow,NCol,cellM{3});
dX1 = [0 0 0 0; -1 1 1 -1];
[XCa,~] = sim_motion(Xs1,[],conn1,.01,140,dX1,0);   % Simulate
[XCb,~] = sim_motion(Xs1,[],conn1,.01,140,-dX1,0);  % Simulate
XC = cat(3,flip(XCa,3),XCb);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,3:4,:),[],2)).^2)));
diff1 = abs(d1-2) + abs(d2-2);       dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-1.143) + abs(d2-3.5); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-3.5) + abs(d2-1.143); dInd3 = find(diff3==min(diff3),1);
plot(d1,d2,'k-','linewidth',.4);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'ko','markersize',ms,'linewidth',lw);
hold off;
visualize_network(XC(:,:,dInd1)/8+[d1(dInd1)+.5 d2(dInd1)+.3]',[],conn1,netSc);
visualize_network(XC(:,:,dInd2)/8+[d1(dInd2)+.5 d2(dInd2)+0]',[],conn1,netSc);
visualize_network(XC(:,:,dInd3)/8+[d1(dInd3)+0 d2(dInd3)+.5]',[],conn1,netSc);
set(gca,'visible',1,'XTick',[1 4],'YTick',[1 4],'fontsize',10);
xlabel('d_1', 'Position', [2.5,.6]);
ylabel('d_2', 'Position', [.6 2.5]);
axis([1 4 1 4]);
text(labX,labY+.2,'c','Units','Normalized','fontsize',10,'fontweight','bold');
text(labRowX,.1,'Trajectory','rotation',90,'Units','Normalized','fontsize',10);


% d: Solution Space Example
subplot(NRow,NCol,cellM{4});
s = sqrt(3);
Xs2 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
dXs2 = [-s/2 0 s/2;...
        -1/2 1 -1/2]/2.5;
visualize_conic(Xs2,dXs2,[-2 2; -2 2],[100;100],3,1,1);
axis([-1 1 -1 1]*1.5 + [netSX netSX 0 0]);
text(labX,labY,'d','Units','Normalized','fontsize',10,'fontweight','bold');
text(.12,labColY,'Velocity','Units','Normalized','fontsize',10);


% e: Constructed Network
nW = .02;
xSh = .22;
ySh = xSh/sqrt(3);
subplot(NRow,NCol,cellM{5});
Xu0 = [-s/2 0;...
        1/2 -1];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,~] = construct_network(Xs2,dXs2,Xu0,conn2,0,0);
Xu2 = Xu2(1:2,:);
visualize_conic(Xs2,dXs2,[-2 2; -2 2],[100;100],0,0,0,.15);
visualize_network(Xs2,Xu2,conn2);
hold on;
line(Xs2(1,1:2)-xSh,Xs2(2,1:2)+ySh,'color','k','LineWidth',lw_d);
line(Xs2(1,1)-xSh+nW*[s -s], Xs2(2,1)+ySh+nW*[-1 1],'color','k','LineWidth',lw_d);
line(Xs2(1,2)-xSh+nW*[s -s], Xs2(2,2)+ySh+nW*[-1 1],'color','k','LineWidth',lw_d);
line(Xs2(1,2:3)+xSh,Xs2(2,2:3)+ySh,'color','k','LineWidth',lw_d);
line(Xs2(1,2)+xSh+nW*[s -s], Xs2(2,2)+ySh+nW*[1 -1],'color','k','LineWidth',lw_d);
line(Xs2(1,3)+xSh+nW*[s -s], Xs2(2,3)+ySh+nW*[1 -1],'color','k','LineWidth',lw_d);
text(Xs2(1,1)-.8, mean(Xs2(2,1:2))+.0, 'd_1','fontsize',10);
text(Xs2(1,3)+.2, mean(Xs2(2,2:3))+.0, 'd_2','fontsize',10);
hold off;
axis([-1 1 -1 1]*1.5 + [netSX netSX 0 0]);
text(labX,labY,'e','Units','Normalized','fontsize',10,'fontweight','bold');

% f: Simulate
subplot(NRow,NCol,cellM{6});
[XC2,~] = sim_motion(Xs2,Xu2,conn2,.005,920,[Xs2,Xu2],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC2(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC2(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-s) + abs(d2-s);       dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-1.96) + abs(d2-2.74); dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-1.24) + abs(d2-2.48); dInd3 = find(diff3==min(diff3),1);
plot(d1,d2,'k-','linewidth',.4);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot([d1(dInd1)-.2,d1(dInd1)+.2],[d2(dInd1)-.2,d1(dInd1)+.2],'-','color',cArrow);
plot(d1(dInd1),d2(dInd1),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd1),d2(dInd1),'ko','markersize',2*ms);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'ko','markersize',ms,'linewidth',lw);
hold off;
visualize_network(XC2(:,1:3,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',...
                  XC2(:,4:5,dInd1)/8+[d1(dInd1)+.2 d2(dInd1)-.2]',conn2,netSc);
visualize_network(XC2(:,1:3,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',...
                  XC2(:,4:5,dInd2)/8+[d1(dInd2)+.3 d2(dInd2)+0]',conn2,netSc);
visualize_network(XC2(:,1:3,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)+.3]',...
                  XC2(:,4:5,dInd3)/8+[d1(dInd3)-.1 d2(dInd3)+.3]',conn2,netSc);
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3],'fontsize',10);
xlabel('d_1', 'Position', [2,.6]);
ylabel('d_2', 'Position', [.6 2]);
axis([.9 3.1 .9 3.1]);
text(labX,labY+.2,'f','Units','Normalized','fontsize',10,'fontweight','bold');

% g: Finite Solution Space
subplot(NRow,NCol,cellM{7})
s = sqrt(3);
Xs30 = [-s/2 0 s/2;...
        -1/2 1 -1/2];
Xs3T = [-0.8  0.0  0.8;...
        -1.0  1.6 -1.0];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs30,Xs3T,[-1 1; -1 1]*2.0,[140;140],5,.7,.85);
axis([-1 1 -1 1]*1.7 + [netSX netSX .4 .4]);
text(labX,labY,'g','Units','Normalized','fontsize',10,'fontweight','bold');
text(.0,labColY,'Displacement','Units','Normalized','fontsize',10);

% h: Finite Network
subplot(NRow,NCol,cellM{8})
Xu0 = [-1.0  0.0;...
        1.6 -0.90];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs30,Xs3T,[-1 1; -1 1]*2,[100;100],0,0,0,.15);
[Xu3,~] = construct_network(Xs30,Xs3T,Xu0,conn3,0,1);
Xu3 = Xu3(1:2,:);
visualize_network(Xs30,Xu3,conn3);
axis([-1 1 -1 1]*1.7 + [netSX netSX .4 .4]);
text(labX,labY,'h','Units','Normalized','fontsize',10,'fontweight','bold');

% i: Simulate
subplot(NRow,NCol,cellM{9});
[XC3,~] = sim_motion(Xs30,Xu3,conn3,.005,700,[Xs30,Xu3],0);   % Simulate
dM0 = sqrt(sum(diff(Xs30,1,2).^2));
dMT = sqrt(sum(diff(Xs3T,1,2).^2));
d1 = sqrt(squeeze(sum((diff(XC3(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC3(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-dM0(1)) + abs(d2-dM0(2));    dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-2.27) + abs(d2-2.00);        dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-dMT(1)) + abs(d2-dMT(2));    dInd3 = find(diff3==min(diff3),1);
plot(d1,d2,'k-','linewidth',.4);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot(d1(dInd1),d2(dInd1),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd1),d2(dInd1),'ko','markersize',2*ms);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd3),d2(dInd3),'ko','markersize',2*ms);
hold off;
visualize_network(XC3(:,1:3,dInd1)/8+[d1(dInd1)+.1 d2(dInd1)-.4]',...
                  XC3(:,4:5,dInd1)/8+[d1(dInd1)+.1 d2(dInd1)-.4]',conn3,netSc);
visualize_network(XC3(:,1:3,dInd2)/8+[d1(dInd2)+.2 d2(dInd2)-.3]',...
                  XC3(:,4:5,dInd2)/8+[d1(dInd2)+.2 d2(dInd2)-.3]',conn3,netSc);
visualize_network(XC3(:,1:3,dInd3)/8+[d1(dInd3)+.2 d2(dInd3)-.3]',...
                  XC3(:,4:5,dInd3)/8+[d1(dInd3)+.2 d2(dInd3)-.3]',conn3,netSc);
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3],'fontsize',10);
xlabel('d_1', 'Position', [2,.6]);
ylabel('d_2', 'Position', [.6 2]);
axis([.9 3.1 .9 3.1]);
text(labX,labY+.2,'i','Units','Normalized','fontsize',10,'fontweight','bold');

% j: Infinitesimal and Finite
subplot(NRow,NCol,cellM{10});
Xs40 = [-s/2 0 s/2;...
        -1/2 1 -1/2];
Xs4T = Xs40*1.7;
dXs4 = [-1.0  0.0  0.0;...
         0.0  0.0  0.0]/2;
visualize_conic_finite(Xs40,Xs4T,[-2 2;-2 2],[100;100],0,.7,0);
visualize_conic(Xs40,dXs4,[-2 2;-2 2],[100;100],0,1,0);
axis([-1 1 -1 1]*2 + [netSX netSX .27 .27]);
text(labX,labY,'j','Units','Normalized','fontsize',10,'fontweight','bold');
text(.05,labColY,'Vel. + Disp.','Units','Normalized','fontsize',10);

% k: Construct Infinitesimal and Finite Network
subplot(NRow,NCol,cellM{11});
% Xu4 = [ 0.950 -1.690;...
%        -1.394 -0.080];
Xu4 = [-0.86 -0.86;...
       -1.45  1.47];
conn4 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
visualize_conic_finite(Xs40,Xs4T,[-2 2;-2 2],[100;100],0,0,0,.15);
visualize_conic(Xs40,dXs4,[-2 2;-2 2],[100;100],0,0,0,.2);
visualize_network(Xs40,Xu4,conn4);
axis([-1 1 -1 1]*2 + [netSX netSX .27 .27]);
text(labX,labY,'k','Units','Normalized','fontsize',10,'fontweight','bold');

% l: Simulate Infinitesimal and Finite
subplot(NRow,NCol,cellM{12});
dM0 = sqrt(sum(diff(Xs40,1,2).^2));
dMT = sqrt(sum(diff(Xs4T,1,2).^2));
[XC4,~] = sim_motion(Xs40,Xu4,conn4,.005,400,[Xs40,Xu4],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC4(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC4(:,2:3,:),[],2)).^2)));
diff1 = abs(d1-dM0(1)) + abs(d2-dM0(2));    dInd1 = find(diff1==min(diff1),1);
diff2 = abs(d1-2.62) + abs(d2-2.00);        dInd2 = find(diff2==min(diff2),1);
diff3 = abs(d1-dMT(1)) + abs(d2-dMT(2));    dInd3 = find(diff3==min(diff3),1);
plot(d1,d2,'k-','linewidth',.4);
hold on;
plot([1 4],[1 4], '--', 'color', [200 200 200]/255);
plot([d1(dInd1)-.3,d1(dInd1)+.3],[d2(dInd1),d1(dInd1)],'-','color',cArrow);
plot(d1(dInd1),d2(dInd1),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd1),d2(dInd1),'ko','markersize',2*ms);
plot(d1(dInd2),d2(dInd2),'ko','markersize',ms,'linewidth',lw);
plot(d1(dInd3),d2(dInd3),'o','markersize',ms,'linewidth',lw,'color',cArrow);
plot(d1(dInd3),d2(dInd3),'ko','markersize',2*ms);
hold off;
visualize_network(XC4(:,1:3,dInd1)/8+[d1(dInd1)-.2 d2(dInd1)+.3]',...
                  XC4(:,4:5,dInd1)/8+[d1(dInd1)-.2 d2(dInd1)+.3]',conn4,netSc);
visualize_network(XC4(:,1:3,dInd2)/8+[d1(dInd2)+.1 d2(dInd2)-.3]',...
                  XC4(:,4:5,dInd2)/8+[d1(dInd2)+.1 d2(dInd2)-.3]',conn4,netSc);
visualize_network(XC4(:,1:3,dInd3)/8+[d1(dInd3)-.3 d2(dInd3)-.15]',...
                  XC4(:,4:5,dInd3)/8+[d1(dInd3)-.3 d2(dInd3)-.15]',conn4,netSc);
set(gca,'visible',1,'XTick',[1 3],'YTick',[1 3],'fontsize',10);
xlabel('d_1', 'Position', [2,.6]);
ylabel('d_2', 'Position', [.6 2]);
axis([.9 3.1 .9 3.1]);
text(labX,labY+.2,'l','Units','Normalized','fontsize',10,'fontweight','bold');


% Legend
subplot(NRow,NCol,figM(:,end)); cla;
set(gca,'visible',0);
ms = 3;
bw = 0.5;
lw = 1.2;                        % Line Width
C_SN = [255 100 100]/255;
C_UN = [100 100 255]/255;
C_SS = [100 100 255;...         % Color of Solution Space
        100 200 255]/255;    
C_E  = [0 0 0 .5];              % Color of Edge
LW_SA = 1;
LW_SS = 1;
tShX = .9;
tShY = .008;
hold on;
% Specified Node
plot(.5, 6/5, 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN);
plot(.5, 6/5, 'ko', 'linewidth', bw, 'markersize', ms*2);
text(tShX, 6/5 + tShY, 'designed node', 'fontsize', 10);
% Variable Node
plot(.5, 5/5, 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN);
plot(.5, 5/5, 'ko', 'linewidth', bw, 'markersize', ms*2);
text(tShX, 5/5 + tShY, 'variable node', 'fontsize', 10);
% Designed Motions
quiver(.2, 4.1/5, .6, 0, 0, 'filled', 'linewidth', LW_SA, 'color', cArrow);
quiver(.2, 4.1/5, .6, 0, 0, 'filled', 'linewidth', LW_SA/4, 'color', [1 1 1]);
quiver(.2, 3.9/5, .6, 0, 0, 'linewidth', LW_SA*2/3, 'color', cArrow);
text(tShX, 4/5 + tShY, 'designed motion', 'fontsize', 10);
% Variable Motions
quiver(.2, 3.1/5, .6, 0, 0, 'filled', 'linewidth', LW_SA, 'color', C_SS(1,:));
quiver(.2, 3.1/5, .6, 0, 0, 'filled', 'linewidth', LW_SA/4, 'color', [1 1 1]);
quiver(.2, 2.9/5, .6, 0, 0, 'linewidth', LW_SA*2/3, 'color', C_SS(1,:));
text(tShX, 3/5 + tShY, 'variable motion', 'fontsize', 10);
% Instantaneous Solution Space
line([.2 .8], [2.2 2.2]/5, 'color', C_SS(1,:), 'linewidth', LW_SS);
line([.2 .8], [2.2 2.2]/5, 'color', [1 1 1], 'linewidth', LW_SS/4);
line([.2 .8], [2 2]/5, 'color', C_SS(1,:), 'linewidth', LW_SS);
line([.2 .8], [1.8 1.8]/5, 'color', C_SS(2,:), 'linewidth', LW_SS);
text(tShX, 2/5 + tShY, 'variable positions', 'fontsize', 10);
% Design
line([.2 .8], [1 1]/5, 'color', cArrow, 'linewidth', .4);
plot(.5,1/5,'o','markersize',ms/2,'linewidth',ms/2,'color',cArrow);
plot(.5,1/5,'ko','markersize',ms);
text(tShX, 1/5 + tShY, 'designed distance', 'fontsize', 10);
hold off;
axis([0 1 -1/5 8/5]);



% Size and Save Figure
fName = 'figure1';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.0 -.43 19.0 9.4];
fig.PaperSize = [19 8.9];
saveas(fig, ['Figures/' fName], 'pdf');



%% Figure 2
% Prepare Space
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


% a: Two Modules
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


% b: Coupled Modules
subplot(NRow,NCol,cellM{2}); cla;
[Xs1a,conn1a] = tesselate_network(Xs1,conn1,[4;0],[2;1]);
visualize_network(Xs1a,[],conn1a,pSc);
text(-2-tSx, tSy, 'd_1','fontsize',10);
text(2-tSx, tSy, 'd_2','fontsize',10);
text(6-tSx, tSy, 'd_3','fontsize',10);
axis([-5 5 -2 2]+[2 2 0 0]);
text(labX,labY,'b','Units','Normalized','fontsize',10,'fontweight','bold');


% c: Cobweb Plot
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


% d: 2 FP + Super Stability
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


% e: Tesselate
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


% f: Simulate
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


% g: Limit Cycle
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


% h: Tesselate Limit Cycle
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


% i: Simulate
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
pInd2 = 270;
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


% j: Chaos
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


% k: Tesselate Chaos
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


% l: Simulate
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


% Size and Save Figure
fName = 'figure2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters'; 
fig.PaperPosition = [-2.28 -0.3 24.4 8.9];
fig.PaperSize = [19 8.5];
saveas(fig, ['Figures/' fName], 'pdf');


%% Figure 3




%% Animate
fig = figure(4); clf;
D1 = sqrt(squeeze(sum(diff(XC,1,2).^2)));
D1 = (D1(1:max(conna(:,1)-1),:));
% XCP = permute(XCd, [2 1 3]);
XCP = XC;
fName = 'animation.gif';
nSV = 1;
dT = 0.03;
nV = [1:3];
nS = 2;
for i = 1:100:size(XC,3)
    cla;
    hold on;
    plot(d1,d2,'k-');
%     plot(d2,d1);
%     plot(d2,d1,'k-');
%     plot(d1,d2);
    plot([0 3],[0 3], '--', 'color', [200 200 200]/255);
    pInd3 = i;
    dP = [D1(:,pInd3)';D1(:,pInd3)']; dP = dP(:);
    dPa = dP(1:end-1); dPb = [0;dP(3:end)];
    line(dP(1:end-1),[0;dP(3:end)],'color',cTr3);
    visualize_network(XCP(:,1:max(conna(:,1)),i)/13 + [1.1;1.4],...
                      XCP(:,[(max(conna(:,1))+1):max(conna(:,2))],i)/13 + [1.1;1.4], conna,.8);
    axis([.8 2.5 .8 2.5]);
    hold off;
%     axis([min(min(min(XCP(1,:)))) max(max(max(XCP(1,:)))) min(min(min(XCP(2,:))))  max(max(max(XCP(2,:))))]);
    drawnow;

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




% subplot(NRow,NCol,cellM{12});










