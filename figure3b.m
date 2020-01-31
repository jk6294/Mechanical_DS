% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');

%% Figure dimensions
fig = figure(3); clf;
spA = [19 9];
set(gcf, 'renderer', 'painters',...
         'Position', [10 10 spA],...
         'Units', 'centimeters'); 
set(gcf, 'renderer', 'painters',...
         'Position', [10 10 spA],...
         'Units', 'centimeters'); 
spMarg = [0 -.02 .05 .1];
spa = [0.00 0.50 0.25 0.50] - spMarg;
spb = [0.25 0.50 0.24 0.50] - spMarg;
spc = [0.55 0.50 0.24 0.50] - spMarg;
spd = [0.79 0.50 0.25 0.50] - spMarg;
spe = [0.00 0.00 1/3  0.45] - spMarg;
spf = [1/3  0.00 1/3  0.45] - spMarg;
spg = [2/3  0.00 1/3  0.45] - spMarg;

% Plot item dimensions
LW = 1;
pSc = .7;

% Label dimension
labX =  0;
labY =  1.15;

% Colors
C_SN = [255 100 100]/255;       % Color of Specified Node
C_UN = [100 100 255]/255;       % Color of Unspecified Node
cTr1 = [126 200 255]/255;
cTr2 = [115 140 200]/255;
cTr3 = [015 082 186]/255;

% Other parameters
s = sqrt(3);
FS = 10;


%% a: Initial Design Visualization
subplot('position',spa); cla;
nxSh = s/2;
nySh = .3;
Xs = [-s/2 -nxSh  s/2;...
      -1/2  nySh -1/2];
Xf = [-s/2  nxSh  s/2;...
      -1/2  nySh -1/2];

% Sample conic
[Q,W,v0,err] = construct_conic(Xs,Xf,1);
P = [W([1:2],:) v0(1:2);...
     zeros(1,2) 1]^-1;
Q1 = P'*Q*P;
% Final
P2 = [W([1:2]+2,:) v0([1:2]+2);...
     zeros(1,2) 1]^-1;
Q2 = P2'*Q*P2;
[CB1 sInd] = sample_conic(Q1,[-1 1 -1 1] - [0 0 .5 .5],800,-90);
% Remap
S = [W v0] * P * [CB1; ones(1,size(CB1,2))];
CB2 = S([1:2]+2,:);

% Obtain conic coordinates
% thL1 = linspace(thSt,thSt+203,50);
% thL2 = linspace(thSt,thSt+153,50);


scatter(CB1(1,:),CB1(2,:),2,winter(size(CB1,2)));
hold on;
% plot([0 rM*cosd(thL1(1))], [0 rM*sind(thL1(1))], 'k-', 'linewidth', LW);
% plot([0 rM*cosd(thL1(end))], [0 rM*sind(thL1(end))], 'k:', 'linewidth', LW);
% plot([0 rM*cosd(thL2(end))], [0 rM*sind(thL2(end))], 'k:', 'linewidth', LW);
visualize_conic_finite(Xs,Xf,[-2 2; -2 2],[0;0],0,.7,0,1,1.5);
hold on;
% plot(rM*cosd(thL1(end)),rM*sind(thL1(end)),'o','markersize',3,'linewidth',3,'color',C_UN);
% plot(rM*cosd(thL1(end)),rM*sind(thL1(end)),'ko','markersize',6,'linewidth',.75);
% plot(rM*cosd(thL2(end)),rM*sind(thL2(end)),'o','markersize',3,'linewidth',3,'color',C_UN);
% plot(rM*cosd(thL2(end)),rM*sind(thL2(end)),'ko','markersize',6,'linewidth',.75);
% plot(.2*cosd(thL1), .2*sind(thL1), 'k-', 'linewidth', LW);
% plot(.35*cosd(thL2), .35*sind(thL2), 'k-', 'linewidth', LW);
line_coordinates(Xs(:,1:2),.3,.08,1);
line_coordinates(Xs(:,2:3),-.3,.08,1);
hold off;
axis([-1 1 -1 1]*1.5);
set(gca,'visible',0);
text(labX,labY,'\textbf{a} Select 2 fixed points as','Units','Normalized','fontsize',FS,'fontweight','bold');
text(labX,labY-.09,'initial and final positions','Units','Normalized','fontsize',FS,'fontweight','bold');
text(.1,.6,'$d_1$','Units','Normalized','fontsize',FS);
text(.6,.2,'$d_2$','Units','Normalized','fontsize',FS);
% text(.36,.4,'$\theta_1$','Units','Normalized','fontsize',FS);
% text(.5,.37,'$\theta_2$','Units','Normalized','fontsize',FS);


%% b: Phase Diagram% Evalulate slope at initial position
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Phase diagram
% Evalulate slope at initial position
Cm1 = sInd;
Cm2 = size(CB1,2)-sInd;
% Distances
d1s = sqrt((Xs(:,1)-Xs(:,2))' * (Xs(:,1)-Xs(:,2)));
d2s = sqrt((Xs(:,2)-Xs(:,3))' * (Xs(:,2)-Xs(:,3)));
d1f = sqrt((Xf(:,1)-Xf(:,2))' * (Xf(:,1)-Xf(:,2)));
d2f = sqrt((Xf(:,2)-Xf(:,3))' * (Xf(:,2)-Xf(:,3)));
% Change in distances
d1dot1 = zeros(Cm1,Cm2);
d2dot1 = zeros(Cm1,Cm2);
d1dot2 = zeros(Cm1,Cm2);
d2dot2 = zeros(Cm1,Cm2);

tic
fprintf([repmat('.',1,Cm1) '\n\n']);
parfor i = 1:Cm1
    fprintf('\b=\n');
    for j = 1:Cm2
        [Uss,Uus,err] = construct_motion(Xs,Xs,CB1(:,[i,sInd+j]),conn,0,0);
        d1dot1(i,j) = (Xs(:,1)-Xs(:,2))' * (Uss(:,1)-Uss(:,2));
        d2dot1(i,j) = (Xs(:,2)-Xs(:,3))' * (Uss(:,2)-Uss(:,3));
    end
end
fprintf([repmat('.',1,Cm1) '\n\n']);
parfor i = 1:Cm1
    fprintf('\b=\n');
    for j = 1:Cm2
        [Uss,Uus,err] = construct_motion(Xf,Xf,CB2(:,[i,sInd+j]),conn,0,0);
        d1dot2(i,j) = (Xf(:,1)-Xf(:,2))' * (Uss(:,1)-Uss(:,2));
        d2dot2(i,j) = (Xf(:,2)-Xf(:,3))' * (Uss(:,2)-Uss(:,3));
    end
end
toc
d1dot1 = d1dot1/d1s;
d2dot1 = d2dot1/d2s;
d1dot2 = d1dot2/d1f;
d2dot2 = d2dot2/d2f;
dRat1 = d2dot1./d1dot1;
dRat2 = d2dot2./d1dot2;
disp('done');


%% b-c: Draw Phase Diagrams
figure(3);
subplot('position',spb); cla;
dRat1c = dRat1.*dRat2;

% Phase diagram buffer
aW1 = round(Cm1/10); aW2 = round(Cm2/10); aW = [aW1 aW2];

% Phase diagram parameters
dRat1M = [abs(dRat1c) -1*ones(Cm1,aW2); -1*ones(aW1,aW2+Cm2)];
dRat1M(isnan(dRat1M)) = -1;

% Set up colormap and axis
colormap([1 1 1; parula(2^11)]);
col = winter(200);
col1 = flipud(col(1:floor(Cm1/(Cm1+Cm2)*size(col,1)),:));
col2 = flipud(col(ceil(Cm1/(Cm1+Cm2)*size(col,1)):size(col,1),:));


% %% Sim
% ind1 = floor(linspace(1,Cm2,20)) + sInd;
% ind2 = floor(linspace(1,Cm1,20));
% nFor = 300;
% nBac = 300;
% nTot = nFor + nBac;
% Zd1 = zeros(length(ind1), length(ind2), nTot);
% Zd2 = zeros(length(ind1), length(ind2), nTot);
% 
% for i = 1:length(ind1)
%     disp(i/length(ind1));
%     for j = 1:length(ind2)
%         Xu1 = CB1(:,[ind1(i) ind2(j)]);
%         
%         % Map simulations
%         [Xc1a,~] = sim_motion(Xs,Xu1,conn,.01,300,[Xs Xu1],0);    % Simulate
%         [Xc1b,~] = sim_motion(Xs,Xu1,conn,.01,300,-[Xs Xu1],0);   % Simulate
%         Xc1 = cat(3,flip(Xc1b,3), Xc1a);
%         
%         % Distances
%         Zd1(i,j,:) = sqrt(squeeze(sum((diff(Xc1(:,1:2,:),[],2)).^2)));
%         Zd2(i,j,:) = sqrt(squeeze(sum((diff(Xc1(:,2:3,:),[],2)).^2)));
%     end
% end
% 
% 
% %% Plot
% figb = figure(4); clf;
% for i = 1:length(ind1)
%     for j = 1:length(ind2)
%         subplot('Position',[1-(i)/length(ind1),1-(j)/length(ind2),...
%                             1/length(ind1), 1/length(ind2)]);
%         da = squeeze(Zd1(i,j,:));
%         db = squeeze(Zd2(i,j,:));
%         if(dRat1M(ind2(j),ind1(i)-sInd) < 1)
%             plot(da, db, 'linewidth', 1,'color',[0 0 1]);
%         else
%             plot(da, db, 'linewidth', 1,'color',[0 1 0]);
%         end
%         d = [da; db];
%         hold on;
%         plot([min(d) max(d)], [min(d) max(d)], '--', 'color', [1 1 1]*.7);
%         plot(d1s,d2s,'ko','markersize',2,'linewidth',2);
%         plot(d2s,d1s,'ko','markersize',2,'linewidth',2);
%         set(gca,'xtick',[],'ytick',[]);
%         axis([min(d) max(d) min(d) max(d)]);
%     end
%     drawnow;
% end
% 
% print('FillPageFigure','-dpdf','-bestfit')


%
% Select module
n1P = 480;
n2P = 237;
uInd1 = [n1P+sInd n2P];
% Unspecified node positions
Xu1 = CB1(:,uInd1);
Xu1f = CB2(:,uInd1);

[Uss1,~,~] = construct_motion(Xs,Xs,Xu1,conn,0,0);
[Uss2,~,~] = construct_motion(Xf,Xf,Xu1f,conn,0,0);

imagesc(tanh(dRat1M));
hold on;
plot(uInd1(1)-sInd,uInd1(2),'s','linewidth',1.8,'markersize',7,'color',cTr1);
scatter(linspace(0,Cm2,size(col2,1)), ones(1,size(col2,1))*Cm1+aW1*.8,1,col2,'s','linewidth',4);
scatter(ones(1,size(col1,1))*Cm2+aW2*.8,linspace(0,Cm1,size(col1,1)),1,col1,'s','linewidth',4);
hold off;
rot180 = rotz(180); rot180 = rot180(1:2,1:2);
visualize_network(rot180*Xs*40+[n1P;n2P],rot180*Xu1*40+[n1P;n2P],conn,...
                  'scale',pSc*.8,'scolor',C_SN,'ucolor',cTr1);
visualize_network(rot180*Xf*40+[n1P-100;n2P],rot180*Xu1f*40+[n1P-100;n2P],conn,...
                  'scale',pSc*.8,'scolor',C_SN,'ucolor',cTr1);
caxis([-.001 1]);
set(gca,'xtick',[],'ytick',[],'ydir','reverse','xdir','reverse','visible',0);
h = colorbar('location', 'east');
h.Position = [0.452 .57 .015 .35];
h.Ticks = [];
text(0,0,'0','Units','Normalized','fontsize',FS);
text(-0.13,0.96,'$2\pi$','Units','Normalized','fontsize',FS);
text(1.03,0,'$2\pi$','Units','Normalized','fontsize',FS);
text(.5,-.09,'$\theta_1$','Units','Normalized','fontsize',FS);
text(-.12,.5,'$\theta_2$','Units','Normalized','fontsize',FS);
text(1.2,.1,'$0$','Units','Normalize','fontsize',FS);
text(1.2,.96,'$1$','Units','Normalize','fontsize',FS);
text(1.22,.22,'$\tanh(|\delta d_2 / \delta d_1|)$','Units','Normalize','fontsize',FS,'rotation',90);
text(labX-.03,labY,'\textbf{b} Calculate slope $|\delta d_2/\delta d_2|$ at','Units','Normalized','fontsize',FS,'fontweight','bold');
text(labX-.01,labY-.09,'each variable node placement','Units','Normalized','fontsize',FS,'fontweight','bold');

%
subplot('position',spc); cla;
dRat1T = (abs(dRat1.*dRat2) > 1)*1;
dRat2M = [abs(dRat1T) -1*ones(Cm1,aW2); -1*ones(aW1,aW2+Cm2)];
dRat2M(isnan(dRat2M)) = -1;
imagesc(dRat2M);
hold on;
plot(uInd1(1)-sInd,uInd1(2),'s','linewidth',1.8,'markersize',7,'color',cTr1);
scatter(linspace(0,Cm2,size(col2,1)), ones(1,size(col2,1))*Cm1+aW1*.8,1,col2,'s','linewidth',4);
scatter(ones(1,size(col1,1))*Cm2+aW2*.8,linspace(0,Cm1,size(col1,1)),1,col1,'s','linewidth',4);
hold off;
caxis([-.001 1]);
text(.44,.89,'stable','Units','Normalized','fontsize',FS,'color',[1 1 1]);
text(.36,.2,'unstable','Units','Normalized','fontsize',FS);
text(0,0,'0','Units','Normalized','fontsize',FS);
text(-0.13,.96,'$2\pi$','Units','Normalized','fontsize',FS);
text(1.03,0,'$2\pi$','Units','Normalized','fontsize',FS);
text(.5,-.09,'$\theta_1$','Units','Normalized','fontsize',FS);
text(-.12,.5,'$\theta_2$','Units','Normalized','fontsize',FS);
text(labX-.1,labY,'\textbf{c} ~~~Unstable modules have','Units','Normalized','fontsize',FS,'fontweight','bold');
text(labX-.1,labY-.09,'~~~~~~~~~slope $|\delta d_2/\delta d_2| > 1$','Units','Normalized','fontsize',FS,'fontweight','bold');
set(gca,'xtick',[],'ytick',[],'ydir','reverse','xdir','reverse','visible',0);


% d: Maps
% Parameters
sRatd = spA(1)*spd(3)/(spA(2)*spd(4));
nSc = 2.4; nSh = 1.1;
netSc = .15; 

% Map simulations
[Xc1a,~] = sim_motion(Xs,Xu1,conn,.01,170,[Xs Xu1],0);    % Simulate
[Xc1b,~] = sim_motion(Xs,Xu1,conn,.01,80,-[Xs Xu1],0);   % Simulate
Xc1 = cat(3,flip(Xc1b,3), Xc1a);

% Map distances
d11 = sqrt(squeeze(sum((diff(Xc1(:,1:2,:),[],2)).^2)));
d21 = sqrt(squeeze(sum((diff(Xc1(:,2:3,:),[],2)).^2)));
d1all = [d11 d21]';
d1all = sum(abs(d1all - flipud(d1all)));
fpInd = find(d1all == min(d1all));

aSh = [1; 1];
subplot('position',spd); cla;
plot([1 3],[1 3], '--', 'color', [200 200 200]/255);
visualize_network(Xc1(:,1:3,fpInd)/4+aSh,Xc1(:,4:5,fpInd)/4+aSh,conn);
hold on;
plot(d11,d21,'-','linewidth',2,'color',cTr1);
plot(d1s,d2s,'k*','linewidth',1.5,'markersize',9);
plot(d1f,d2f,'k*','linewidth',1.5,'markersize',9);
hold off;
drawnow;
axis([min([d11;d21]) max([d11;d21]) min([d11;d21]) max([d11;d21])]);

% Build chain
Xsc = [Xs [s/2; 0.3]];
Xuc = [Xu1 [Xu1f(1,:);-Xu1f(2,:)-.2]];
connc = [1 5; 2 5; 3 5; 1 6; 2 6; 3 6; 2 7; 3 7; 4 7; 2 8; 3 8; 4 8];
[Xsa1,Xua1,conna1] = tesselate_network_old(Xsc,Xuc,connc,[s;0],[5;1]);
subplot('position',spe); cla;
visualize_network(Xsa1,Xua1,conna1);


%%


% % Simulate
[XC1,fC1] = sim_motion(Xsa1,Xua1,conna1,.01,1000,[Xsa1 Xua1],0);
% [XC2,fC2] = sim_motion(Xsa2,Xua2,conna2,.02,700,[Xsa2 Xua2],0);
% [XC3,fC3] = sim_motion(Xsa3,Xua3,conna3,.01,1000,[Xsa3 Xua3],0);
disp(max(fC1));
% 
% %
nS1 = size(Xsa1,2);
% nS2 = size(Xsa2,2);
% nS3 = size(Xsa3,2);
% 
% 
% %% e-g: Visualize
% sRate = spA(1)*spe(3)/(spA(2)*spe(4));
% 
subplot('position',spe); cla;
dc = squeeze(sqrt(sum(diff(XC1,1,2).^2)));
dc = dc([1:10],:);

visualize_network(XC1(:,1:nS1,1)*netSc+[.86;1.2],...
                  XC1(:,nS1+1:end,1)*netSc+[.86;1.2],conna1,...
                  'scale',pSc,'scolor',C_SN,'ucolor',cTr1);
% visualize_network(XC1(:,1:nS1,end)*netSc+[.7;0.4],...
%                   XC1(:,nS1+1:end,end)*netSc+[.7;0.4],conna1,...
%                   'scale',pSc,'scolor',C_SN,'ucolor',cTr1);
% axis([[0 1]*sRate-.08 [0 1]+.05]*1.6);
% text(labX,labY,'\textbf{e} ~~~Stable networks change shape','Units','Normalized','fontsize',FS,'fontweight','bold');
% text(labX,labY-.11,'~~~~~~~~~~~~~~from the $d_1$ end','Units','Normalized','fontsize',FS,'fontweight','bold');
% 
% subplot('position',spf); cla;
% delXC2 = XC2(:,:,2) - XC2(:,:,1);
% [Us2,Uu2,err] = construct_motion(XC2(:,1:nS1,1),delXC2(:,1:nS1),...
%                                  XC2(:,nS1+1:end,1),conna2,0,0);
% quiver(XC2(1,:,1)*netSc+.5,XC2(2,:,1)*netSc+1.2,...
%        [Us2(1,:) Uu2(1,:)], [Us2(2,:) Uu2(2,:)],'linewidth',1,'color','g');
% visualize_network(XC2(:,1:nS1,1)*netSc+[.5;1.2],...
%                   XC2(:,nS1+1:end,1)*netSc+[.5;1.2],conna2,...
%                   'scale',pSc,'scolor',C_SN,'ucolor',cTr2);
% visualize_network(XC2(:,1:nS1,end)*netSc+[.5;0.4],...
%                   XC2(:,nS1+1:end,end)*netSc+[.5;0.4],conna2,...
%                   'scale',pSc,'scolor',C_SN,'ucolor',cTr2);
% axis([[0 1]*sRate-.1 [0 1]+.05]*1.6);
% text(labX-.1,labY,'\textbf{f} ~~~Marginally stable networks change','Units','Normalized','fontsize',FS,'fontweight','bold');
% text(labX-.1,labY-.11,'~~~~~~shape through the whole network','Units','Normalized','fontsize',FS,'fontweight','bold');
% 
% subplot('position',spg); cla;
% delXC3 = XC3(:,:,2) - XC3(:,:,1);
% [Us3,Uu3,err] = construct_motion(XC3(:,1:nS1,1),delXC3(:,1:nS1),...
%                                  XC3(:,nS1+1:end,1),conna3,0,0);
% quiver(XC3(1,:,1)*netSc+.23,XC3(2,:,1)*netSc+1.2,...
%        [Us3(1,:) Uu3(1,:)], [Us3(2,:) Uu3(2,:)],'linewidth',1,'color','g');
% visualize_network(XC3(:,1:nS1,1)*netSc+[.23;1.2],...
%                   XC3(:,nS1+1:end,1)*netSc+[.23;1.2],conna2,...
%                   'scale',pSc,'scolor',C_SN,'ucolor',cTr3);
% visualize_network(XC3(:,1:nS1,end)*netSc+[.4;0.4],...
%                   XC3(:,nS1+1:end,end)*netSc+[.4;0.4],conna2,...
%                   'scale',pSc,'scolor',C_SN,'ucolor',cTr3);
% axis([[0 1]*sRate-.16 [0 1]+.05]*1.6);
% text(labX,labY,'\textbf{g} ~~~Unstable networks change shape','Units','Normalized','fontsize',FS,'fontweight','bold');
% text(labX,labY-.11,'~~~~~~~~~~~~~~~~from the $d_2$ end','Units','Normalized','fontsize',FS,'fontweight','bold');


figure(9); clf;
animate_network(XC1,conna1,'distance',[d11 d21]','nu',size(Xua1,2),...
                'nframe',100,'figwidth',10);
% Animation

                        
%% Save
% fName = 'figure3b';
% set(gcf, 'Renderer', 'painters');
% fig.PaperPositionMode = 'manual';
% fig.PaperUnits = 'centimeters';
% fig.PaperPosition = [0 0 19 9];
% fig.PaperSize = [19 9];
% saveas(fig, ['Figures/' fName], 'pdf');

