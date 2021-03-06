% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');

% Subplot Indices
NCol = 20;
colBal = {[1 1 1], [1 1 1], [1 1 1], 1};
rowBal = [9 9 9 14]; 
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
axM = [-1 14 0 4];
labX = -0.02;
labY = 0.95;
labColY = 1.3;

fig = figure(3); clf;
set(gcf, 'Renderer', 'painters', 'Position', [10 10 24.0 9.5], 'Units', 'centimeters'); 
s = sqrt(3);
pSc = 0.5;
lw_d = .5;
lSh = .5;
nW = .1;

annotation('line','linewidth',.5,'position',[.380 .44 0 .4],'color',[.9 .9 .9]);
annotation('line','linewidth',.5,'position',[.658 .44 0 .4],'color',[.9 .9 .9]);
annotation('line','linewidth',.5,'position',[.135 .42 .765 0],'color',[.9 .9 .9]);


%% a: Slope 1 Module
subplot(NRow,NCol,cellM{1}); cla;

% Single Module
Xs1 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
dXs1 = [ 0.0  0.0 -0.0;...
         0.7 -0.7  0.7];
lSh = .6;
Xu1 = [[- lSh; 1] [ lSh; 1]];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Flipped Module
Xs1f = [Xs1(1,:); -Xs1(2,:)+.5];
dXs1f = [dXs1(1,:); -dXs1(2,:)];
Xu1f = [Xu1(1,:); -Xu1(2,:)+.5];

% Combined Module
Xs1c = [Xs1 Xs1f+[s/2;0]];
Xu1c = [Xu1 Xu1f+[s/2;0]];
conn1c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[3;1]);
dXs1G = zeros(size(Xs1a)); 
dXs1G(2,1:2:size(Xs1a,2)) = .7;
dXs1G(2,2:2:size(Xs1a,2)) = -.7;

% Show Module
construct_motion(Xs1+[1.1;0], dXs1, Xu1+[1.1;0], conn1, 1, 1, pSc);
construct_motion(Xs1a+[6.6;0],dXs1G,Xu1a+[6.6;0],conn1a,1, 1, pSc);
line_coordinates(Xs1(:,1:2)+[1.1;0],1.1,nW,lw_d);
line_coordinates(Xs1(:,2:3)+[1.1;0],1.1,nW,lw_d);
axis(axM - [1 1 1.3 1.3]);

% Text
text(.04, .6, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(.31, .6, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(labX,labY,'\textbf{a}','Units','Normalized','fontsize',10);
text(labX+.12,labY,'Module','Units','Normalized','fontsize',10);
text(labX+.6,labY,'Combined','Units','Normalized','fontsize',10);
text(.17,labColY,'\textbf{Uniform:} $\mathrm{\delta d_2/\delta d_1 = 1}$','Units','Normalized',...
     'fontsize',10);
drawnow;


%% b: Slope 1 Propagate Combined Motion
subplot(NRow,NCol,cellM{4}); cla;
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[10;1]);
[XC,fC] = sim_motion(Xs1a,Xu1a,conn1a,.05,75,[Xs1a Xu1a],0);
visualize_network(XC(:,1:size(Xs1a,2),1),...
                  XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),1),conn1a, pSc);
axis(axM+[2.5 2.5 -1.1 -1.1]);
text(labX,labY-.1,'\textbf{b}','Units','Normalized','fontsize',10);
text(labX+.4,labY-.1,'Motion','Units','Normalized','fontsize',10);
drawnow;

subplot(NRow,NCol,cellM{7}); cla;
visualize_network(XC(:,1:size(Xs1a,2),end),...
                  XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),end),conn1a, pSc);
axis(axM+[2.5 2.5 -2.3 -2.3]);
drawnow;


%% c: Slope -2 Module
subplot(NRow,NCol,cellM{2}); cla;

% Single Module
Xs2 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
dXs2 = [ 0.0  0.0  0.0;...
         -0.2  0.0  -0.3]*2.5;
Xu2 = [ 0.30  0.8;...
       -0.50  1.3];
xV = s/2;
Xu2 = [[0.0; -0.5], [xV; .5*xV+1]];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Flipped Module
Xs2f = [Xs2(1,:); -Xs2(2,:)+.5];
dXs2f = [dXs2(1,:); -dXs2(2,:)];
Xu2f = [Xu2(1,:); -Xu2(2,:)+.5];
[Uss,Uus,err] = construct_motion(Xs2f+[6;0], dXs2f, Xu2f+[6;0], conn2,0,0);

% Combined Module
Xs2c = [Xs2 Xs2f+[s/2;0]];
Xu2c = [Xu2 Xu2f+[s/2;0]];
conn2c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[3;1]);
dXs2G = zeros(size(Xs2a)); 
r = Uss(2,1)/Uss(2,3);
dXs2G(2,end) = 2.3;
dXs2G(2,end-1) = 0;
dXs2G(2,end-2) = dXs2G(2,end)*r^1;
dXs2G(2,end-3) = -dXs2G(2,end)*(r^2-r);
dXs2G(2,end-4) = dXs2G(2,end)*(r^3-r^2+r);
dXs2G(2,end-5) = -dXs2G(2,end)*(r^4-r^3+r^2-r);
dXs2G(2,end-6) = dXs2G(2,end)*(r^5-r^4+r^3-r^2+r);
dXs2G(2,end-7) = -dXs2G(2,end)*(r^6-r^5+r^6-r^5+r^2-r);
dXs2G(2,:) = (dXs2G(2,:)-dXs2G(2,end-5))*1.0;

% Show Module
construct_motion(Xs2+[1.1;0], dXs2, Xu2+[1.1;0], conn2, 1, 1, pSc);
construct_motion(Xs2a+[6.6;0],dXs2G,Xu2a+[6.6;0],conn2a,1,1, pSc);
line_coordinates(Xs2(:,1:2)+[1.1;0],.5,nW,lw_d);
line_coordinates(Xs2(:,2:3)+[1.1;0],.5,nW,lw_d);
axis(axM - [1 1 1.3 1.3]);

% Text
text(.04, .4, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(.31, .4, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(labX,labY,'\textbf{c}','Units','Normalized','fontsize',10);
text(labX+.12,labY,'Module','Units','Normalized','fontsize',10);
text(labX+.6,labY,'Combined','Units','Normalized','fontsize',10);
text(-.12,labColY,'\textbf{Partially Localized Right:} $\mathrm{\delta d_2/\delta d_1 = 1.5}$',...
     'Units','Normalized','fontsize',10);
drawnow;


%% d: Slope -2 Propagate Combined Motion
subplot(NRow,NCol,cellM{5}); cla;
[Xs2a,Xu2a,conn2a] = tesselate_network_old(Xs2c,Xu2c,conn2c,[s;0],[10;1]);
[XC,fC] = sim_motion(Xs2a,Xu2a,conn2a,.05, 50, [Xs2a Xu2a],0);
visualize_network(XC(:,1:size(Xs2a,2),1),...
                  XC(:,[1:size(Xu2a,2)]+size(Xs2a,2),1),conn2a, pSc);
axis(axM+[2.5 2.5 -1.1 -1.1]);
text(labX,labY-.1,'\textbf{d}','Units','Normalized','fontsize',10);
text(labX+.4,labY-.1,'Motion','Units','Normalized','fontsize',10);

subplot(NRow,NCol,cellM{8}); cla;
visualize_network(XC(:,1:size(Xs2a,2),end),...
                  XC(:,[1:size(Xu2a,2)]+size(Xs2a,2),end),conn2a, pSc);
axis(axM+[2.5 2.5 -2.3 -2.3]);
drawnow;


%% e: Slope 0 Module
subplot(NRow,NCol,cellM{3}); cla;

% Single Module
Xs3 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
dXs3 = [ 0.0  0.0 -0.0;...
         0.8  0.0 -0.0];
lSh = 0.3;
Xu3 = [[-lSh; -.5] [.25*s; .25]+[-1;s]*lSh/2];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Single Module Flipped
Xs3f = [Xs3(1,:); -Xs3(2,:)+.5];
dXs3f = [dXs3(1,:); -dXs3(2,:)];
Xu3f = [Xu3(1,:); -Xu3(2,:)+.5];

% Combined Module
Xs3c = [Xs3 Xs3f+[s/2;0]];
Xu3c = [Xu3 Xu3f+[s/2;0]];
conn3c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[3;1]);
dXs3G = zeros(size(Xs3a)); dXs3G(:,1) =  dXs3(:,1);

% Visualize
construct_motion(Xs3+[1.1;0], dXs3, Xu3+[1.1;0], conn3, 1, 1, pSc);
construct_motion(Xs3a+[6.6;0],dXs3G,Xu3a+[6.6;0],conn3a,1,1, pSc);
line_coordinates(Xs2(:,1:2)+[1.1;0],.8,nW,lw_d);
line_coordinates(Xs2(:,2:3)+[1.1;0],.8,nW,lw_d);

% Text
axis(axM - [1 1 1.3 1.3]);
text(.04, .55, '$\mathrm{d_1}$', 'Units', 'Normalized', 'fontsize',10);
text(.31, .55, '$\mathrm{d_2}$', 'Units', 'Normalized', 'fontsize',10);
text(labX,labY,'e','Units','Normalized','fontsize',10,'fontweight','bold');
text(labX+.12,labY,'Module','Units','Normalized','fontsize',10);
text(labX+.6,labY,'Combined','Units','Normalized','fontsize',10);
text(-.02,labColY,'\textbf{Fully Localized Left:} $\mathrm{\delta d_2/\delta d_1 = 0}$',...
     'Units','Normalized','fontsize',10);
drawnow;


%% f: Slope 0 Propagate Combined Motion
subplot(NRow,NCol,cellM{6}); cla;
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3c,Xu3c,conn3c,[s;0],[14;1]);
[XC,fC] = sim_motion(Xs3a,Xu3a,conn3a,.05,345,-[Xs3a Xu3a],0);
visualize_network(XC(:,1:size(Xs3a,2),1),...
                  XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),1),conn3a, pSc);
axis(axM+[4.5 4.5 -1.1 -1.1]);
text(labX,labY-.1,'\textbf{f}','Units','Normalized','fontsize',10,'fontweight','bold');
text(labX+.4,labY-.1,'Motion','Units','Normalized','fontsize',10);
drawnow;

subplot(NRow,NCol,cellM{9}); cla;
visualize_network(XC(:,1:size(Xs3a,2),end),...
                  XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),end),conn3a, pSc);
axis(axM+[4.5 4.5 -2.3 -2.3]);
drawnow;


%% g: Place Modules
subplot(NRow,NCol,cellM{10}); cla;

% Arrow Annotations
annotation('arrow','HeadLength',8,'HeadWidth',12,'color',[.7 .7 .7],...
           'linewidth',4,'position',[.285 .25 .05 0]);
annotation('arrow','HeadLength',8,'HeadWidth',12,'color',[.7 .7 .7],...
           'linewidth',4,'position',[.50 .25 .05 0]);
annotation('arrow','HeadLength',8,'HeadWidth',12,'color',[.7 .7 .7],...
           'linewidth',4,'position',[.71 .25 .05 0]);

% Construct Unit
Xs3 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs3f = [Xs3(1,:); -Xs3(2,:)+.5];
lSh = .3;
Xu3 = [[-lSh; -.5] [.25*s; .25]+[-1;s]*lSh/2];

% Xu3 = [-Xu3(1,:); Xu3(2,:)];
Xu3f = [Xu3(1,:); -Xu3(2,:)+.5];
Xs3c = [Xs3 Xs3f+[s/2;0]]; Xs3c = [Xs3c(1,:); -Xs3c(2,:)];
Xu3c = [Xu3 Xu3f+[s/2;0]]; Xu3c = [Xu3c(1,:); -Xu3c(2,:)];
conn1c = [1 7; 2 7; 3 7; 1 8; 2 8; 3 8; 4 9; 5 9; 6 9; 4 10; 5 10; 6 10];
nR = 5;
[Xs3a1,Xu3a1,conn3a1] = tesselate_network_old(Xs3c,Xu3c,conn1c,[s;0],[2*nR;1]);
[Xs3a2,Xu3a2,conn3a2] = tesselate_network_old(Xs3c,Xu3c,conn1c,[s;0],[nR;1]);
% Remove Offset
Xu3a1 = Xu3a1-Xs3a1(:,1); Xs3a1 = Xs3a1-Xs3a1(:,1);
Xu3a2 = Xu3a2-Xs3a2(:,1); Xs3a2 = Xs3a2-Xs3a2(:,1);

% Modified Unit: Branch 1
Xs3a2 = [Xs3a2(1,:); -Xs3a2(2,:)]; 
Xu3a2 = [Xu3a2(1,:); -Xu3a2(2,:)]; 
Rz = rotz(-60); Rz = Rz(1:2,1:2);
Xs3a2 = Rz*Xs3a2 + [nR*s/2;nR*s^2/2];
Xu3a2 = Rz*Xu3a2 + [nR*s/2;nR*s^2/2];

% Motions
dXs3a1 = zeros(size(Xs3a1)); dXs3a1(2,1) = -1;
dXs3a2 = zeros(size(Xs3a2)); dXs3a2(1,1) =  1;

% Visualize
construct_motion(Xs3a1,dXs3a1,Xu3a1,conn3a1,1,1,pSc/2);
construct_motion(Xs3a2+[0;1],dXs3a2,Xu3a2+[0;1],conn3a2,1,1,pSc/2);
axis([-.1 80 0 11.0]-[0 0 2 2]);
text(-.006,labY-.1,'\textbf{g}','Units','Normalized','fontsize',10);
text(.39,labY+.05,'\textbf{Mechanical AND Gate}','Units','Normalized','fontsize',10);
drawnow;


%% Combine and Simulate Modules
Xs3a = [Xs3a1 Xs3a2];
Xu3a = [Xu3a1 Xu3a2];
conn3a = conn3a1;
conn3a = [[conn3a(:,1), conn3a(:,2)+max(conn3a2(:,1))];...
          [conn3a2(:,1)+max(conn3a(:,1)), conn3a2(:,2)+max(conn3a(:,2))]];
[Xs3a,Xu3a,conn3a] = tesselate_network_old(Xs3a,Xu3a,conn3a,[1;1],[1;1]);

XsDot = [zeros(size(Xs3a)) zeros(size(Xu3a))];
XsDot(2,1) = -1;
[XC,fC] = sim_motion(Xs3a,Xu3a,conn3a,.01,3300,XsDot,0);


%% h: Visualize
pInd = 600;
rotV = rotz(0.8); rotV = rotV(1:2,1:2);
sh = [24; -.123];
visualize_network(rotV*XC(:,1:size(Xs3a,2),pInd)+sh,...
                  rotV*XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),pInd)+sh,conn3a, pSc/2);
text(.28,labY-.15,'\textbf{h}','Units','Normalized','fontsize',10);
drawnow;


%% i: More Collapsed
pInd = 1500;
XCD = -(XC(:,:,pInd)-XC(:,:,pInd-1));
XCD = 2*XCD / sqrt(sum(diag(XCD*XCD')));
R = rigidity(XC(:,:,pInd),conn3a);
[ua,sa,va] = svds(R,2,'smallest');
va = -2*va/sqrt(va'*va);
rotV = rotz(1.9); rotV = rotV(1:2,1:2);
sh = [45; -.609];
visualize_network(rotV*XC(:,1:size(Xs3a,2),pInd)+sh,...
                  rotV*XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),pInd)+sh,conn3a, pSc/2);
text(.56,labY-.15,'\textbf{i}','Units','Normalized','fontsize',10);
drawnow;


%% j: More Collapsed
pInd = 3300;
XCD = -(XC(:,:,pInd)-XC(:,:,pInd-1));
XCD = 2*XCD / sqrt(sum(diag(XCD*XCD')));
R = rigidity(XC(:,:,pInd),conn3a);
[ua,sa,va] = svds(R,2,'smallest');
va = -2*va/sqrt(va'*va);
rotV = rotz(0.52); rotV = rotV(1:2,1:2);
sh = [65; -0.712];
visualize_network(rotV*XC(:,1:size(Xs3a,2),pInd)+sh,...
                  rotV*XC(:,[1:size(Xu3a,2)]+size(Xs3a,2),pInd)+sh,conn3a, pSc/2);
text(.83,labY-.15,'\textbf{j}','Units','Normalized','fontsize',10);
drawnow;


%% Size and Save Figure
fName = 'figure3';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters'; 
fig.PaperPosition = [-2.9 -1.02 24.0 9.5];
fig.PaperSize = [19 8.2];
saveas(fig, ['Figures/' fName], 'pdf');


%% Animation
% fig = figure(4); clf;
% fName = 'animation.gif';
% dT = 0.03;
% dT = 0.1;
% nXSh = -2.5;
% nYSh = 2.3;
% nSV = 1;
% Xs1a = Xs2a;
% Xu1a = Xu2a;
% XC = XC2a;
% conn1a = conn2a;
% % XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);
% 
% D1 = sqrt(squeeze(sum(diff(XC(:,1:size(Xs1a,2),:),1,2).^2)));
% D1 = D1(1:end,:); 
% 
% for i = 11:40:1941
%     cla;
%     dP = [D1(:,i)';D1(:,i)']; dP = dP(:);
%     dPa = dP(1:end-1); dPb = [1;dP(3:end)];
%     plot(d1,d2,'k-');
%     hold on;
% %     [Us, Uu, ~] = construct_motion(XC(:,1:size(Xs1a,2),i),...
% %                                    XdotV,...
% %                                    XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),...
% %                                    conn1a, 0, 0);
% %     Us = 0.3*Us/sqrt(sum(sum(Us.^2))); 
% %     Uu = 0.3*Uu/sqrt(sum(sum(Uu.^2)));
% %     XdotV = Us;
% %     quiver(XC(1,:,i)/10 + nXSh, XC(2,:,i)/10 + nYSh,[Us(1,:) Uu(1,:)], [Us(2,:) Uu(2,:)],0);
%     plot([0 4],[0 4],'LineStyle','--','color',[0 0 0 .5]);
%     line(dP(1:end-1),[1;dP(3:end)]);
%     visualize_network(XC(:,1:size(Xs1a,2),i)/5 + [nXSh;nYSh],...
%                       XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i)/5 + [nXSh;nYSh],conn1a,.7);
%     axis([-4.0 3 1.6 3]);
% %     construct_motion(XC(:,1:size(Xsa,2),i), XC(:,1:size(Xsa,2),i+1)-XC(:,1:size(Xsa,2),i), XC(:,[1:size(Xua,2)]+size(Xsa,2),i), conn, 20, 20);
% %     axis([min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))]);
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


%% Animate Network Alone
% fig = figure(5); clf;
% fName = 'animation_net.gif';
% dT = 0.03;
% % dT = 0.1;
% nXSh = 1.6;
% nYSh = 2;
% nSV = 1;
% Xs1a = Xs3a;
% Xu1a = Xu3a;
% conn1a = conn3a;
% XdotV = XC(:,1:size(Xs1a,2),1); XdotV = XdotV-mean(XdotV,2);
% 
% for i = 1:20:size(XC,3)
%     cla;
% %     R = rigidity(XC(:,:,i),conn1a);
% %     [u,sa,va] = svds(R,1,'smallest');
%     visualize_network(XC(:,1:size(Xs1a,2),i),...
%                       XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i),conn1a,1);
%     hold on;
% %     quiver(XC(1,:,i),XC(2,:,i),-va(1:size(XC,2))',-va([1:size(XC,2)]+size(XC,2))',0);
%     hold off;
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


