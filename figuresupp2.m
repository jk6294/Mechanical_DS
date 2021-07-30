% Specifics on Iterated Maps
clear; clc;
set(groot,'defaulttextinterpreter','latex');

fig = figure(2); clf;
flw = [19.8 7.8];
set(gcf, 'Renderer', 'painters', 'Position', [0 15 flw], 'Units', 'centimeters'); 
set(gcf, 'Renderer', 'painters', 'Position', [0 15 flw], 'Units', 'centimeters'); 

pSc = .8;
cBTr = [100 150 255]/255;
cBOp = [000 000 255]/255;
cB   = [100 150 255]/255;
cRTr = [255 255 255]/255;
cROp = [000 000 000]/255;
cR   = [255 255 255]/255;


%% a: Module Combination
subplot(1,2,1); cla;
pLim = [-1 17];
lSh = .5;
nW = .1;
lw = .7;


% Single Modules
s = sqrt(3);
Xs1 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
Xs1p = [Xs1(1,:);-Xs1(2,:)+.5];
Xu1 = [-0.86 -0.86;...
       -1.45  1.47];
Xu1p = [Xu1(1,:);-Xu1(2,:)+.5];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Combined Modules
Xs1c = [Xs1 [s;1]];
Xu1c = [Xu1 Xu1p+[s/2;0]];
conn1c = [1 5; 2 5; 3 5; 1 6; 2 6; 3 6; 2 7; 3 7; 4 7; 2 8; 3 8; 4 8];

% Combined Modules 2
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[2;1]);
[Xs1aa,Xu1aa,conn1aa] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[3;1]);

% 2 modules
% M1
Sh = [0;14];
visualize_network(Xs1+Sh,Xu1+Sh,conn1,'scale',pSc,'scolor',[cRTr;cROp;cROp],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1(1,nP))+Sh(1)-0.9,mean(Xs1(2,nP))+Sh(2)+2.2,'$\mathrm{l_1}$','fontsize',10);
nP = 2:3;
line_coordinates(Xs1(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1(1,nP))+Sh(1)-0.3,mean(Xs1(2,nP))+Sh(2)+2.2,'$\mathrm{l_2}$','fontsize',10);
% M2
Sh = [1.5*s;14];
visualize_network(Xs1p+Sh,Xu1p+Sh,conn1,'scale',pSc,'scolor',[cROp;cROp;cRTr],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1p(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1p(1,nP))+Sh(1)-0.2,mean(Xs1p(2,nP))+Sh(2)-1.7,'$\mathrm{l_2^\prime}$','fontsize',10);
nP = 2:3;
line_coordinates(Xs1p(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1p(1,nP))+Sh(1)-0.0,mean(Xs1p(2,nP))+Sh(2)-1.7,'$\mathrm{l_3}$','fontsize',10);
% M3
Sh = [6*s;14];
visualize_network(Xs1c+Sh,Xu1c+Sh,conn1c,'scale',pSc,'scolor',[cRTr;cROp;cROp;cRTr],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.9,mean(Xs1c(2,nP))+Sh(2)+2.2,'$\mathrm{l_1}$','fontsize',10);
text(7.5,14.2,'$\boldmath{\rightarrow}$','fontsize',14);
nP = 3:4;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.0,mean(Xs1c(2,nP))+Sh(2)-1.7,'$\mathrm{l_3}$','fontsize',10);


% 4 modules
% M1
Sh = [0;8];
visualize_network(Xs1c+Sh,Xu1c+Sh,conn1c,'scale',pSc,...
                  'scolor',[cRTr;cRTr;cROp;cROp],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.9,mean(Xs1c(2,nP))+Sh(2)+2.2,'$\mathrm{l_1}$','fontsize',10);
nP = 3:4;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.0,mean(Xs1c(2,nP))+Sh(2)-1.7,'$\mathrm{l_3}$','fontsize',10);
% M2
Sh = [2*s;8];
visualize_network(Xs1c+Sh,Xu1c+Sh,conn1c,'scale',pSc,...
                  'scolor',[cROp;cROp;cRTr;cRTr],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.9,mean(Xs1c(2,nP))+Sh(2)+2.2,'$\mathrm{l_3^\prime}$','fontsize',10);
nP = 3:4;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.0,mean(Xs1c(2,nP))+Sh(2)-1.7,'$\mathrm{l_5}$','fontsize',10);
% M3
Sh = [6*s;8];
visualize_network(Xs1a+Sh,Xu1a+Sh,conn1a,'scale',pSc,...
                  'scolor',[cRTr;cRTr;cROp;cROp;cRTr;cRTr],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1a(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1a(1,nP))+Sh(1)-0.9,mean(Xs1a(2,nP))+Sh(2)+2.2,'$\mathrm{l_1}$','fontsize',10);
text(7.5,8.2,'$\boldmath{\rightarrow}$','fontsize',14);
nP = 5:6;
line_coordinates(Xs1a(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1a(1,nP))+Sh(1)-0.0,mean(Xs1a(2,nP))+Sh(2)-1.7,'$\mathrm{l_5}$','fontsize',10);


% 6 modules
% M1
Sh = [0;2];
visualize_network(Xs1a+Sh,Xu1a+Sh,conn1a,'scale',pSc,...
                  'scolor',[cRTr;cRTr;cRTr;cRTr;cROp;cROp],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1a(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1a(1,nP))+Sh(1)-0.9,mean(Xs1a(2,nP))+Sh(2)+2.2,'$\mathrm{l_1}$','fontsize',10);
nP = 5:6;
line_coordinates(Xs1a(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1a(1,nP))+Sh(1)-0.0,mean(Xs1a(2,nP))+Sh(2)-1.7,'$\mathrm{l_5}$','fontsize',10);
% M2
Sh = [3*s;2];
visualize_network(Xs1c+Sh,Xu1c+Sh,conn1c,'scale',pSc,...
                  'scolor',[cROp;cROp;cRTr;cRTr],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.9,mean(Xs1c(2,nP))+Sh(2)+2.2,'$\mathrm{l_5^\prime}$','fontsize',10);
nP = 3:4;
line_coordinates(Xs1c(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1c(1,nP))+Sh(1)-0.0,mean(Xs1c(2,nP))+Sh(2)-1.7,'$\mathrm{l_7}$','fontsize',10);
% M3
Sh = [6*s;2];
visualize_network(Xs1aa+Sh,Xu1aa+Sh,conn1aa,'scale',pSc,...
                  'scolor',[cRTr;cRTr;cRTr;cRTr;cROp;cROp;cRTr;cRTr],'ucolor',cBTr);
nP = 1:2;
line_coordinates(Xs1aa(:,nP)+Sh,'lSh',lSh,'nW',nW,'lw',lw);
text(mean(Xs1aa(1,nP))+Sh(1)-0.9,mean(Xs1aa(2,nP))+Sh(2)+2.2,'$\mathrm{l_1}$','fontsize',10);
text(7.5,2.2,'$\boldmath{\rightarrow}$','fontsize',14);
nP = 7:8;
line_coordinates(Xs1aa(:,nP)+Sh,'lSh',-lSh,'nW',nW,'lw',lw);
text(mean(Xs1aa(1,nP))+Sh(1)-0.0,mean(Xs1aa(2,nP))+Sh(2)-1.7,'$\mathrm{l_7}$','fontsize',10);
           


% Plot Parameters
axis([pLim-.5,pLim]);
text(-0.05,0.98,'\textbf{a}','units','normalized','fontsize',10);
text(-0.05,0.65,'\textbf{b}','units','normalized','fontsize',10);
text(-0.05,0.31,'\textbf{c}','units','normalized','fontsize',10);



%% b: Cobweb
subplot(1,2,2); cla;
s = sqrt(3);
pLim = [1.6 3.1];
pSh = 0;
netSC = .7;

Xs1 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
Xu1 = [-0.86 -0.86;...
       -1.45  1.47];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Simulate: Single module for map
[XC,fC] = sim_motion10(Xs1,Xu1,conn1,.02,95,[Xs1 Xu1],0);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
% Simulate: Combined network for propagation
[XC,fC] = sim_motion10(Xs1aa,Xu1aa,conn1aa,.02,200,[Xs1aa Xu1aa],0);
D1 = sqrt(squeeze(sum(diff(XC,1,2).^2)));
D1 = D1(1:size(Xs1aa,2)-1,:);
pInd = 200;
dP = [D1(:,pInd)';D1(:,pInd)']; dP = dP(:);
dPa = dP(1:end-1); dPb = [1;dP(3:end)];


% Plot
arrHL = 5;
arrHW = 6;
caF = .8;

plot(d1,d2,'k-','linewidth',2);
hold on;
plot(pLim,pLim, '--', 'color', [200 200 200]/255);
line(dP(1:end-1),[1;dP(3:end)],'color',cBTr,'linewidth',1);
for i = 1:length(dP)-2
    ah = annotation('arrow','HeadLength',arrHL,'HeadWidth',arrHW,'color',cBTr);
    set(ah,'parent',gca);
    set(ah,'position',[dPa(i) dPb(i) diff(dPa(i:i+1))*caF diff(dPb(i:i+1))*caF]);
end
cRShT = ones(size(Xs1aa,2)-1,3).*linspace(0,1,size(Xs1aa,2)-1)'; cRSh(:,1)=1;
plot(dP(1),1.6,'s','linewidth',4,'markersize',4,'color',cRShT(1,:));
plot(dP(1),1.6,'ks','linewidth',1,'markersize',8);
for i = 2:2:length(dP)-1
    plot(dP(i),dP(i+1),'s','linewidth',4,'markersize',4,'color',cRShT(i/2+1,:));
    plot(dP(i),dP(i+1),'ks','linewidth',1,'markersize',8);
end
Sh = [2;2.8];
visualize_network(XC(:,1:size(Xs1aa,2),pInd)/10+Sh,...
                  XC(:,[1:size(Xu1aa,2)]+size(Xs1aa,2),pInd)/10+Sh,...
                  conn1aa,'scolor',cR,'ucolor',cB);
line_coordinates(XC(:,1:2,pInd)/10+Sh,'lSh',.07,'nW',.02);
line_coordinates(XC(:,7:8,pInd)/10+Sh,'lSh',-.07,'nW',.02);
text(dP(1)-.1,1.7,'$\mathrm{l_1}$','fontsize',10);
for i = 2:length(dP)/2
    text(dP(2*i-2)-.1,dP(2*i-1)+.1,['$\mathrm{l_', num2str(i), '}$'],'fontsize',10);
end
hold on;
plot(1.75,2.86,'s','linewidth',4,'markersize',4,'color',cRShT(1,:));
plot(1.75,2.86,'ks','linewidth',1,'markersize',8);
plot(2.6,2.79,'s','linewidth',4,'markersize',4,'color',cRShT(7,:));
plot(2.6,2.79,'ks','linewidth',1,'markersize',8);
hold off;
axis([pLim+pSh pLim-pSh]);
set(gca,'visible',1,'XTick',[],'YTick',[],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);

% Labels
text(-0.05,0.98,'\textbf{d}','units','normalized','fontsize',10);
text(2.3,1.8,'$\mathrm{l_{k+1} = f(l_k)}$','fontsize',10);
text(1.9,2.3,'$\mathrm{l_{k+1} = l_k}$','fontsize',10,'color',[200 200 200]/255);
xlabel('$\mathrm{l_k}$','fontsize',10);
ylabel('$\mathrm{l_{k+1}}$','fontsize',10);


%% Size and Save Figure
fName = 'figure2supp';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.5 -.2 flw];
fig.PaperSize = [16.5 7.2];
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl'); 