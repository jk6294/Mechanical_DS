% Specifics on Iterated Maps
clear; clc;
set(groot,'defaulttextinterpreter','latex');

fig = figure(1); clf;
flw = [19.8 7.8];
set(gcf, 'Renderer', 'painters', 'Position', [-25 25 flw], 'Units', 'centimeters'); 


%% a: Module: Full Motion
subplot(1,2,1);
s = sqrt(3);
pLim = [1.2 4.2];
netSC = .7;

Xs1 = [-s/2 0 s/2;...
       -1/2 1 -1/2];
Xu1 = [-0.86 -0.86;...
       -1.45  1.47];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

[XC,fC] = sim_motion(Xs1,Xu1,conn1,.02,430,[Xs1 Xu1],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
plot(d1,d2,'k-','linewidth',1);
hold on;
plot(pLim,pLim, '--', 'color', [200 200 200]/255);
pIndL = [[0:22:100]+1 135 156 180 210 340 360 390];
vI = [ 0.00  0.10  0.20  0.40  0.60  0.10 -0.10 -0.20 -0.20 -0.30 -0.30 -0.25;...
      -0.25 -0.30 -0.25 -0.10  0.00  0.25  0.30  0.30  0.00  0.00  0.00  0.10];
for i = 1:length(pIndL)
    hold on;
    pInd = pIndL(i);
    vSh = [d1(pInd);d2(pInd)] + vI(:,i);
    plot(d1(pInd),d2(pInd),'ko','linewidth',3,'markersize',3);
    visualize_network(XC(:,1:3,pInd)/8+vSh,XC(:,4:5,pInd)/8+vSh,conn1,netSC);
end
hold on;
pInd = pIndL(5);
vSh = [d1(pInd);d2(pInd)] + vI(:,5);
line_coordinates(XC(:,1:2,pInd)/8+vSh,.1,.03,1);
line_coordinates(XC(:,2:3,pInd)/8+vSh,.1,.03,1);
hold off;
axis([pLim pLim]);
set(gca,'visible',1,'XTick',[],'YTick',[],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);

% Labels
text(vSh(1)-.4,vSh(2)+.15,'$\mathrm{d_1}$','fontsize',10);
text(vSh(1)+.2,vSh(2)+.15,'$\mathrm{d_2}$','fontsize',10);
text(vSh(1)+0,vSh(2)+1,'$\mathrm{d_2 = d_1}$','fontsize',10,'color',[200 200 200]/255);
text(vSh(1)-.5,vSh(2)+.7,'$\mathrm{g(d_1,d_2)=0}$','fontsize',10);
text(-0.05,0.98,'\textbf{a}','units','normalized','fontsize',10);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);


%% b: Module: Partial Motion
subplot(1,2,2); cla;
pLim = [1.2 4.2];

% Section 1
pR = [1:91];
plot(d1(pR),d2(pR),'-','linewidth',1,'color',[126 200 255]/255);
hold on;
% Section 2
pR = [93:200];
plot(d1(pR),d2(pR),'-','linewidth',1,'color',[115 140 200]/255);
plot(pLim,pLim, '--', 'color', [200 200 200]/255);
% Section 3
pR = [205:295];
plot(d1(pR),d2(pR),'-','linewidth',1,'color',[015 082 186]/255);
plot(pLim,pLim, '--', 'color', [200 200 200]/255);
% Section 4
pR = [300:405];
plot(d1(pR),d2(pR),'-','linewidth',1,'color',[000 000 000]/255);
plot(pLim,pLim, '--', 'color', [200 200 200]/255);
hold off;
axis([pLim pLim]);
set(gca,'visible',1,'XTick',[],'YTick',[],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);
                
% Labels
text(2,3.8,'$\mathrm{d_2 = f_2(d_1)}$','fontsize',10,'color',[115 140 200]/255);
text(3,2.3,'$\mathrm{d_2 = f_1(d_1)}$','fontsize',10,'color',[126 200 255]/255);
text(1.3,2.8,'$\mathrm{d_2 = f_3(d_1)}$','fontsize',10,'color',[015 082 186]/255);
text(1.3,2.5,'$\mathrm{d_2 = f_4(d_1)}$','fontsize',10,'color',[000 000 000]/255);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);
text(-0.05,0.98,'\textbf{a}','units','normalized','fontsize',10);




%% Size and Save Figure
fName = 'figure1supp';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.5 .2 flw];
fig.PaperSize = [16.5 7.8];
saveas(fig, ['Figures/' fName], 'pdf');