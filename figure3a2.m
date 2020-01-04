% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');

%% Figure dimensions
fig = figure(3); clf;
spA = [19 9.5];
set(gcf, 'renderer', 'painters',...
         'Position', [80 12 spA],...
         'Units', 'centimeters'); 
spMarg = [-.015 -.03 .05 .1];
spa = [0.00 0.50 0.25 0.50] - spMarg;
spb = [0.00 0.00 0.25 0.50] - spMarg;
spc = [0.24 0.00 0.50 1.00] - spMarg;
spd = [0.75 0.50 0.25 0.50] - spMarg;
spe = [0.00 0.00 1/3  0.45] - spMarg;
spf = [1/3  0.00 1/3  0.45] - spMarg;
spg = [2/3  0.00 1/3  0.45] - spMarg;

% Plot item dimensions
LW = 1;
pSc = .56;

% Label dimension
labX =  -.07;
labY =  1.13;

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
Xs = [-s/2 0 s/2;...
      -1/2 1 -1/2];
rM = 1.7;
Xf = rM*Xs;

% Sample conic
[Q,W,v0,err] = construct_conic(Xs,Xf,1);
P = [W([1:2],:) v0(1:2);...
     zeros(1,2) 1]^-1;
Q = P'*Q*P;
nP = 360*2;
% thSt = 103;
thSt = 90;
CB1 = sample_conic(Q,[-2 2 -2 2],nP,thSt);
Cm1 = size(CB1,2);

% Obtain conic coordinates
thL1 = linspace(thSt,thSt+203,50);
thL2 = linspace(thSt,thSt+153,50);


scatter(CB1(1,:),CB1(2,:),2,cool(Cm1));
hold on;
plot([0 rM*cosd(thL1(1))], [0 rM*sind(thL1(1))], 'k-', 'linewidth', LW);
plot([0 rM*cosd(thL1(end))], [0 rM*sind(thL1(end))], 'k:', 'linewidth', LW);
plot([0 rM*cosd(thL2(end))], [0 rM*sind(thL2(end))], 'k:', 'linewidth', LW);
visualize_conic_finite(Xs,Xf,[-2 2; -2 2],[0;0],0,.7,0,1,1.5);
hold on;
plot(rM*cosd(thL1(end)),rM*sind(thL1(end)),'o','markersize',3,'linewidth',3,'color',C_UN);
plot(rM*cosd(thL1(end)),rM*sind(thL1(end)),'ko','markersize',6,'linewidth',.75);
plot(rM*cosd(thL2(end)),rM*sind(thL2(end)),'o','markersize',3,'linewidth',3,'color',C_UN);
plot(rM*cosd(thL2(end)),rM*sind(thL2(end)),'ko','markersize',6,'linewidth',.75);
plot(.2*cosd(thL1), .2*sind(thL1), 'k-', 'linewidth', LW);
plot(.35*cosd(thL2), .35*sind(thL2), 'k-', 'linewidth', LW);
line_coordinates(Xs(:,1:2),.35,.08,1);
line_coordinates(Xs(:,2:3),.35,.08,1);
hold off;
axis([-1 1 -1 1]*2);
set(gca,'visible',0);
text(labX,labY,'\textbf{a} Select 2 fixed points as','Units','Normalized','fontsize',FS,'fontweight','bold');
text(labX,labY-.09,'initial and final positions','Units','Normalized','fontsize',FS,'fontweight','bold');
text(.24,.67,'$d_1$','Units','Normalized','fontsize',FS);
text(.69,.67,'$d_2$','Units','Normalized','fontsize',FS);
text(.35,.43,'$\theta_1$','Units','Normalized','fontsize',FS);
text(.465,.37,'$\theta_2$','Units','Normalized','fontsize',FS);


%% b: Fixed Points
% Distances
d1s = sqrt((Xs(:,1)-Xs(:,3))' * (Xs(:,1)-Xs(:,3)));
d2s = sqrt((Xs(:,2)-Xs(:,3))' * (Xs(:,2)-Xs(:,3)));
d1f = sqrt((Xf(:,1)-Xf(:,3))' * (Xf(:,1)-Xf(:,3)));
d2f = sqrt((Xf(:,2)-Xf(:,3))' * (Xf(:,2)-Xf(:,3)));
p = [d1s 2.0 2.5 d1f;...
     d2s 1.8 2.2 d2f];
xq = linspace(d1s-1, d1f+1, 100);
yq = spline(p(1,:), p(2,:), xq);
xqs = [d1s d1s + .001];
yqs = spline(p(1,:), p(2,:), xqs);
sl = diff(yqs)/diff(xqs);
subplot('position',spb); cla;
plot([0 4], [0 4], '--', 'color', [200 200 200]/255, 'linewidth',1);
hold on;
plot([d1s d1f], [d2s d2f], 'k*', 'markersize', 8, 'linewidth', 1);
plot([0 4], [0 4]*sl + (1-sl)*d1s, 'color', 'b', 'linewidth', 1);
plot(xq,yq,'k-','linewidth',1);
hold off;
axis([1.2 3.5 1.2 3.5]);
set(gca,'xtick',[],'ytick',[]);
% Axis and labels
text(labX,labY,'\textbf{b}','Units','Normalized','fontsize',FS,'fontweight','bold');
text(1.03,0,'$d_1$','Units','Normalized','fontsize',FS);
text(labX,1.05,'$d_2$','Units','Normalized','fontsize',FS);
% Markers
text(.15,.35,'$D_1^*$','Units','Normalized','fontsize',FS);
text(.65,.86,'$D_2^*$','Units','Normalized','fontsize',FS);
text(.21,.12,'slope $=\delta d_2 / \delta d_1$','Units','Normalized',...
             'fontsize',FS,'color','b');


%% c: Phase Diagram Evalulate slope at initial position
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Change in distances
d1dot1 = zeros(Cm1);
d2dot1 = zeros(Cm1);

UssM1 = zeros([size(Xs),Cm1,Cm1]);
th = linspace(0,360,Cm1+1); th = th(1:end-1);
tic
fprintf([repmat('.',1,Cm1) '\n\n']);
parfor i = 1:Cm1
    fprintf('\b=\n');
    for j = 1:Cm1
        if i ~= j
            [Uss,Uus,err] = construct_motion(Xs,Xs,CB1(:,[i,j]),conn,0,0);
            UssM1(:,:,i,j) = Uss;
            d1dot1(i,j) = (Xs(:,1)-Xs(:,2))' * (Uss(:,1)-Uss(:,2));
            d2dot1(i,j) = (Xs(:,2)-Xs(:,3))' * (Uss(:,2)-Uss(:,3));
        end
    end
end
toc
d1dot1 = d1dot1/d1s;
d2dot1 = d2dot1/d2s;
disp('done');


%% c: Draw Phase Diagrams
nScb = 20; 
nSh1b = [0;-50];
nSh1c = [27; -42];
nSh1d = [0; -50];
subplot('position',spc); cla;
dRat1 = d2dot1./d1dot1;
aW = round(Cm1/20);
dRat1M = [nan(aW,aW+Cm1); nan(Cm1,aW) abs(dRat1)];
colormap([1 1 1; parula(2^11)]);
lx = round(linspace(27,Cm1/2-50,6)); ly = Cm1 - lx;
lx = lx+aW+1; ly = ly+aW+1;
l0x = [63  121 180 219 268 285 309];
l0y = [320 348 397 458 603 668 752];
l1x = [319 360 421 490 532 580 639 700];
l1y = [131 221 264 304 377 439 468 484];
l0xq = linspace(min(l0x), max(l0x), 100);
l0yq = pchip(l0x,l0y,l0xq);
C_B = [0.2422    0.1504    0.6603];

imagesc(tanh((dRat1M)),'alphadata',.7);
hold on;
contour(dRat1M,[1 1],'r','linewidth',1.5);
plot(l0xq,l0yq,'-','linewidth',1.5,'color',C_B);
for i = 1:length(lx)
    uIndi = [lx(i) ly(i)];
    Xui = CB1(:,uIndi-aW);
    plot(uIndi(1),uIndi(2),'s','linewidth',1.5,'markersize',5,'color',cTr1);
    visualize_network( Xs*nScb+uIndi'+nSh1b,...
                      Xui*nScb+uIndi'+nSh1b,conn,pSc,C_SN,cTr1);
    hold on;
end
for i = 1:length(l0x)
    uIndi = [l0x(i) l0y(i)];
    Xui = CB1(:,uIndi-aW);
    plot(uIndi(1),uIndi(2),'s','linewidth',1.5,'markersize',5,'color',cTr3);
    visualize_network( Xs*nScb+uIndi'+nSh1c,...
                      Xui*nScb+uIndi'+nSh1c,conn,pSc,C_SN,cTr3);
    hold on;
end
for i = 1:length(l1x)
    uIndi = [l1x(i) l1y(i)];
    Xui = CB1(:,uIndi-aW);
    plot(uIndi(1),uIndi(2),'s','linewidth',1.5,'markersize',5,'color',cTr2);
    visualize_network( Xs*nScb+uIndi'+nSh1d,...
                      Xui*nScb+uIndi'+nSh1d,conn,pSc,C_SN,cTr2);
    hold on;
end
scatter(linspace(0,Cm1-.4*aW,100)+1.25*aW, ones(1,100)*16,1,cool(100),'s',...
        'linewidth',3);
scatter(ones(1,100)*16,linspace(0,Cm1-.4*aW,100)+1.25*aW,1,cool(100),'s',...
        'linewidth',3);
hold off;
caxis([-.001 1]);
set(gca,'xtick',[],'ytick',[],'ydir','normal','visible',0);
h = colorbar('location', 'east','position',[0.715 .1 .015 .805]);
delete(findall(gcf,'type','annotation'));
a = annotation('textbox',h.Position,'FitBoxToText','off','FaceAlpha',0.3,...
    'EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]);
h.Ticks = tanh(1); h.TickLabels = {};

% Axis Ticks
text(0.01,0.015,'0','Units','Normalized','fontsize',FS);
text(-0.04,0.98,'$2\pi$','Units','Normalized','fontsize',FS);
text(1.02,0.015,'$2\pi$','Units','Normalized','fontsize',FS);
text(.5,-0.01,'$\theta_1$','Units','Normalized','fontsize',FS);
text(-.04,.5,'$\theta_2$','Units','Normalized','fontsize',FS);
% Colorbar Ticks
text(1.07,.075,'$0$','Units','Normalize','fontsize',FS);
text(1.07,.96,'$1$','Units','Normalize','fontsize',FS);
text(1.085,.71,'$\tanh(1)$','Units','Normalize','fontsize',FS,...
               'rotation',90,'color','r');
annotation('line',[.715 .731], [.714 .714], 'Units','Normalize',...
           'color','r','linewidth',1.5);
annotation('line',[.715 .731], [.1 .1], 'Units','Normalize',...
           'color',C_B,'linewidth',1.5);
text(1.085,.38,'$\tanh(|\delta d_2 / \delta d_1|)$','Units','Normalize','fontsize',FS,'rotation',90);
% Description
text(labX+.03,1.055,'\textbf{c} Calculate slope $|\delta d_2/\delta d_2|$ at each variable node placement','Units','Normalized','fontsize',FS,'fontweight','bold');
text(labX-.01,1.02,'','Units','Normalized','fontsize',FS,'fontweight','bold');
text(.3,.38,'Unstable','Units','Normalize','fontsize',FS);
text(.59,.65,'Stable','Units','Normalize','fontsize',FS);


%% d: Motion of Chains
Xu1 = CB1(:,[97  335]-aW);
Xu2 = CB1(:,[380 675]-aW);
Xu3 = CB1(:,[l0x(5) l0y(5)]-aW);

% Map
% figure(1); clf;
% [XCp,fCp] = sim_motion(Xs,Xu2,conn,.02,400,[Xs Xu2],0);
% d1p = squeeze(sqrt(sum(diff(XCp(:,[1,2],:),1,2).^2)));
% d2p = squeeze(sqrt(sum(diff(XCp(:,[2,3],:),1,2).^2)));
% 
% plot(d1p,d2p);
% hold on;
% plot(d2p,d1p);
% plot([1 4], [1 4]);
% hold off;


%%

% Chain
[Xsa1,Xua1,conna1] = network_chain_x(Xs(:,1:2),Xu1,ones(1,10));
[Xsa2,Xua2,conna2] = network_chain_x(Xs(:,1:2),Xu2,ones(1,10));
[Xsa3,Xua3,conna3] = network_chain_x(Xs(:,1:2),Xu3,ones(1,10));

% Simulate
[XC1,fC1] = sim_motion(Xsa1,Xua1,conna1,.02,400,[Xsa1 Xua1],0);
[XC2,fC2] = sim_motion(Xsa2,Xua2,conna2,.02,400,[Xsa2 Xua2],0);
[XC3,fC3] = sim_motion(Xsa3,Xua3,conna3,.02,400,[Xsa3 Xua3],0);

% Sizes
nS1 = size(Xsa1,2);
nS2 = size(Xsa2,2);
nS3 = size(Xsa3,2);


%% d: Plot
netSc = .1;
netSh1 = [.6; 1.4];
netSh2 = [.2; 0.9];
sRate = spA(1)*spd(3)/(spA(2)*spd(4));

% Average difference
xcd1 = mean(XC1(:,9:nS1,end) - XC1(:,9:nS1,1),2);
xcd2 = mean(XC2(:,1:4,end) - XC2(:,1:4,1),2);

subplot('position',spd); cla;
delXC1 = XC1(:,:,2) - XC1(:,:,1);
[Us1,Uu1,err] = construct_motion(XC1(:,1:nS1,1),delXC1(:,1:nS1),...
                                 XC1(:,nS1+1:end,1),conna1,0,0);
visualize_network(XC1(:,1:nS1,end)*netSc+netSh1-xcd1*netSc,...
                  XC1(:,nS1+1:end,end)*netSc+netSh1-xcd1*netSc,conna1,pSc,C_SN,cTr1,.3);
hold on;
quiver(XC1(1,:,1)*netSc+netSh1(1),XC1(2,:,1)*netSc+netSh1(2),...
       [Us1(1,:) Uu1(1,:)], [Us1(2,:) Uu1(2,:)],'linewidth',1,'color','g');
hold off;
visualize_network(XC1(:,1:nS1,1)*netSc+netSh1,...
                  XC1(:,nS1+1:end,1)*netSc+netSh1,conna1,pSc,C_SN,cTr1);
              
delXC2 = XC2(:,:,2) - XC2(:,:,1);
[Us2,Uu2,err] = construct_motion(XC2(:,1:nS1,1),delXC2(:,1:nS1),...
                                 XC2(:,nS1+1:end,1),conna1,0,0);
visualize_network(XC2(:,1:nS1,end)*netSc+netSh2-xcd2*netSc,...
                  XC2(:,nS1+1:end,end)*netSc+netSh2-xcd2*netSc,conna1,pSc,C_SN,cTr1,.3);
hold on;
quiver(XC2(1,:,1)*netSc+netSh2(1),XC2(2,:,1)*netSc+netSh2(2),...
       [Us2(1,:) Uu2(1,:)], [Us2(2,:) Uu2(2,:)],'linewidth',1,'color','g');
hold off;
visualize_network(XC2(:,1:nS1,1)*netSc+netSh2,...
                  XC2(:,nS1+1:end,1)*netSc+netSh2,conna1,pSc,C_SN,cTr1);
axis([[0 1]*sRate-.08 [0 1]+.05]*1.6);


%% Save
fName = 'figure3a2';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 19 9.5];
fig.PaperSize = [19 9.5];
saveas(fig, ['Figures/' fName], 'pdf');



% 
% 
% %% e-g: Visualize
% sRate = spA(1)*spe(3)/(spA(2)*spe(4));
% 
% subplot('position',spe); cla;
% delXC1 = XC1(:,:,2) - XC1(:,:,1);
% [Us1,Uu1,err] = construct_motion(XC1(:,1:nS1,1),delXC1(:,1:nS1),...
%                                  XC1(:,nS1+1:end,1),conna1,0,0);
% quiver(XC1(1,:,1)*netSc+.86,XC1(2,:,1)*netSc+1.2,...
%        [Us1(1,:) Uu1(1,:)], [Us1(2,:) Uu1(2,:)],'linewidth',1,'color','g');
% visualize_network(XC1(:,1:nS1,1)*netSc+[.86;1.2],...
%                   XC1(:,nS1+1:end,1)*netSc+[.86;1.2],conna1,pSc,C_SN,cTr1);
% visualize_network(XC1(:,1:nS1,end)*netSc+[.7;0.4],...
%                   XC1(:,nS1+1:end,end)*netSc+[.7;0.4],conna1,pSc,C_SN,cTr1);
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
%                   XC2(:,nS1+1:end,1)*netSc+[.5;1.2],conna2,pSc,C_SN,cTr2);
% visualize_network(XC2(:,1:nS1,end)*netSc+[.5;0.4],...
%                   XC2(:,nS1+1:end,end)*netSc+[.5;0.4],conna2,pSc,C_SN,cTr2);
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
%                   XC3(:,nS1+1:end,1)*netSc+[.23;1.2],conna2,pSc,C_SN,cTr3);
% visualize_network(XC3(:,1:nS1,end)*netSc+[.4;0.4],...
%                   XC3(:,nS1+1:end,end)*netSc+[.4;0.4],conna2,pSc,C_SN,cTr3);
% axis([[0 1]*sRate-.16 [0 1]+.05]*1.6);
% text(labX,labY,'\textbf{g} ~~~Unstable networks change shape','Units','Normalized','fontsize',FS,'fontweight','bold');
% text(labX,labY-.11,'~~~~~~~~~~~~~~~~from the $d_2$ end','Units','Normalized','fontsize',FS,'fontweight','bold');
% 
%                         


