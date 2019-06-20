% Specifics on Iterated Maps
clear; clc;
set(groot,'defaulttextinterpreter','latex');

fig = figure(4); clf;
flw = [19.8 7.8];
set(gcf, 'Renderer', 'painters', 'Position', [0 15 flw], 'Units', 'centimeters'); 


%% a: Module: Map
subplot(1,2,1); cla;
s = sqrt(3);
pLim = [0.6 2.5];
netSC = .75;

Xu1 = [0 1/2;...
       0 s/2];
Xs1 = [1 2;...
       0 0];
conn1 = [1 3; 2 3; 1 4; 2 4];

[XC,fC] = sim_motion(Xs1,Xu1,conn1,.01,202,[Xs1 Xu1],0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,3:4,:),[],2)).^2)));
plot(d1,d2,'k-','linewidth',1);
hold on;
plot(pLim,pLim, '--', 'color', [200 200 200]/255);
pIndL = [.5 15 29 44 59 69 80]*2;
vI = [-0.2 -0.34 -0.34 -0.25 -0.15  0.00  0.03;...
      -0.15  0.05  0.10  0.15  0.20  0.20  0.15];
for i = 1:length(pIndL)
    hold on;
    pInd = pIndL(i);
    vSh = [d1(pInd);d2(pInd)] + vI(:,i);
    visualize_network(XC(:,1:2,pInd)/8+vSh,XC(:,3:4,pInd)/8+vSh,conn1,netSC,...
                      [255 100 100]/255, [255 100 100]/255);
end
hold on;
for i = [2 3 4 5 7]
    pInd = pIndL(i);
    plot(d1(pInd),d2(pInd),'ko','linewidth',3,'markersize',3);
end
for i = [1 6]
    pInd = pIndL(i);
    plot(d1(pInd),d2(pInd),'k*','linewidth',1.5,'markersize',9);
end
pInd = 1;
vSh = [d1(pInd);d2(pInd)] + vI(:,1);
line_coordinates(XC(:,1:2,pInd)/8+vSh,-.07,.02,1);
line_coordinates(XC(:,3:4,pInd)/8+vSh,.07,.02,1);
hold off;
axis([pLim pLim]);
set(gca,'visible',1,'XTick',[],'YTick',[],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);

% Labels
text(d1(pIndL(1))+.1,d2(pIndL(1))-.01,'$\mathrm{D_1^*}$','fontsize',10);
text(d1(pIndL(6))-.05,d2(pIndL(6))-.15,'$\mathrm{D_2^*}$','fontsize',10);
text(vSh(1)-.15,vSh(2)+.15,'$\mathrm{d_2}$','fontsize',10);
text(vSh(1)+.15,vSh(2)-.14,'$\mathrm{d_1}$','fontsize',10);
text(vSh(1)+.52,vSh(2)+.4,'$\mathrm{d_2 = d_1}$','fontsize',10,'color',[200 200 200]/255);
text(vSh(1)+1.1,vSh(2)+.55,'$\mathrm{d_2 = f(d_1)}$','fontsize',10);
text(-0.05,0.98,'\textbf{a}','units','normalized','fontsize',10);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);


%% b: Planar: Map
subplot(1,2,2); cla;
s = sqrt(3);
pLim = [0.6 5.4];
netSC = .75;

Xs2 = [ 0  0  0  1.2  1.2  3  5.2  5.2;...
        0 -2  2 -0.4  0.4  0 -2.0  2.0];
Xs2 = [-Xs2(1,:)+5.2;Xs2(2,:)];
conn2 = [1 2; 1 3; 1 4; 1 5; 2 4; 3 5; 4 6; 4 7; 5 6; 5 8; 6 7; 6 8];

[XC,fC] = sim_motion(Xs2,[],conn2,.05,81,Xs2,0);   % Simulate
d1 = sqrt(squeeze(sum((diff(XC(:,7:8,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
plot(d1,d2,'k-','linewidth',1);
hold on;
plot(pLim,pLim, '--', 'color', [200 200 200]/255);
pIndL = [1 29 59 81];
vI = [ 0.05 -0.30  0.10  0.30;...
      -0.65 -0.70 -0.40  0.00];
pIndsTr = [1 2 4; 1 3 5; 4 6 7; 5 6 8];
for i = 1:length(pIndL)
    pInd = pIndL(i);
    vSh = [d1(pInd);d2(pInd)] + vI(:,i);
    hold on;
    for j = 1:size(pIndsTr,1)
        patch(XC(1,pIndsTr(j,:),pInd)/6+vSh(1),...
              XC(2,pIndsTr(j,:),pInd)/6+vSh(2),...
              [160 160 255]/255);
    end
    visualize_network(XC(:,:,pInd)/6+vSh,[],conn2,netSC,...
                      [255 100 100]/255, [255 100 100]/255);
end
hold on;
for i = [2 3]
    pInd = pIndL(i);
    plot(d1(pInd),d2(pInd),'ko','linewidth',3,'markersize',3);
end
for i = [1 4]
    pInd = pIndL(i);
    plot(d1(pInd),d2(pInd),'k*','linewidth',1.5,'markersize',9);
end
pInd = 1;
vSh = [d1(pInd);d2(pInd)] + vI(:,1);
line_coordinates(XC(:,2:3,pInd)/6+vSh,-.15,.05,1);
line_coordinates(XC(:,7:8,pInd)/6+vSh,.15,.05,1);
hold off;
axis([pLim pLim]);
set(gca,'visible',1,'XTick',[],'YTick',[],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);

% Labels
text(d1(pIndL(1))+.2,d2(pIndL(1))-.02,'$\mathrm{D_1^*}$','fontsize',10);
text(d1(pIndL(4))-.5,d2(pIndL(4))-.00,'$\mathrm{D_2^*}$','fontsize',10);
text(vSh(1)-.43,vSh(2),'$\mathrm{d_1}$','fontsize',10);
text(vSh(1)+1.08,vSh(2),'$\mathrm{d_2}$','fontsize',10);
text(4.00,5.10,'$\mathrm{d_2 = d_1}$','fontsize',10,'color',[200 200 200]/255);
text(1.00,3.5,'$\mathrm{d_2 = f(d_1)}$','fontsize',10);
text(-0.05,0.98,'\textbf{b}','units','normalized','fontsize',10);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);


%% Size and Save Figure
fName = 'figure4supp';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.5 -.2 flw];
fig.PaperSize = [16.5 7.2];
saveas(fig, ['Figures/' fName], 'pdf');