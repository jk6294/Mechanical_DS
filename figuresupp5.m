% Specifics on Iterated Maps
clear; clc;
set(groot,'defaulttextinterpreter','latex');

fig = figure(5); clf;
flw = [19.8 7.8];
set(gcf, 'Renderer', 'painters', 'Position', [0 15 flw], 'Units', 'centimeters'); 


%% a: Module: Construction
subplot(1,2,1); cla;
axis([0 4.5 0 4.5] - [0 0 1 1]*.5);

% Plot Parameters
lw_o = 2;
lw_c = 1;

% Square
s2 = sqrt(2);
L = 2;
c = L*(sqrt(2)/(2+2*sqrt(2)));
ca = L/(1+sqrt(2));
cSh = [2 2]';

line(L*s2/2*[0 0]+cSh(1), [-L*s2/2 c]+cSh(2),...
     'linewidth',lw_c,'color','k','LineStyle','--');
hold on;
line(L*s2/2*[0 0]+cSh(1), [c L*s2/2]+cSh(2),...
     'linewidth',lw_c,'color','k','LineStyle',':');
line(L*s2/2*[-1 0]+cSh(1), [0 c]+cSh(2),...
     'linewidth',lw_c,'color','k','LineStyle',':');
line(L*s2/2*[1 0]+cSh(1), [0 c]+cSh(2),...
     'linewidth',lw_c,'color','k','LineStyle',':');
 
line(L*s2/2*[-1 0 1 0 -1]+cSh(1), L*s2/2*[0 1 0 -1 0]+cSh(2),...
     'linewidth',lw_o,'color','k');
line_coordinates(L*s2/2*[0 0; -1 1]+cSh,L*.85,L/30,lw_c);
line_coordinates(L*s2/2*[-1 1; 0 0]+cSh,L*.85,L/30,lw_c);
line_coordinates([0 0; L*s2/2 c]+cSh,L*.1,L/30,lw_c);
line_coordinates(L*s2/2*[0 1; -1 0]+cSh,-L*.1,L/30,lw_c);
hold on;
plot(cSh(1), c+cSh(2),'ko', 'linewidth',1,'markersize',5);
hold off;
text(0.0,2,'$\mathrm{d_1}$','fontsize',10);
text(2.0,3.9,'$\mathrm{d_2}$','fontsize',10);
text(2.3,3.3,'$\mathrm{\frac{L}{1+\sqrt{2}}}$','fontsize',10);
text(3.0,1,'$\mathrm{L}$','fontsize',10);
text(-0.05,0.98,'\textbf{a}','units','normalized','fontsize',10);
set(gca,'visible',0);

legend('Mountain','Valley','location','southeast');


%% b: Planar: Map
subplot(1,2,2); cla;
% Plot Parameters
pLim = [0.2 1.75];
netSC = .75;
% Network Parameters
s2 = sqrt(2);
L = 1;
c = L*(sqrt(2)/(2+2*sqrt(2)));
zSh = pLim(1)-mean(pLim);

Xs2 = L*s2/2*[-1 0 1 0 0; 0 1 0 -1 c; 0 0 0 0 0];
dXs2 = [0 0 0 0 0; 0 0 0 0 0; 1 -1 1 -1 0];
dXs2b = [0 0 0 0 0; 0 0 0 0 0; 1 -1 1 -1 -1];
conn2 = [1 2; 1 4; 2 3; 3 4; 1 5; 2 5; 3 5; 4 5];

[XCa,fCa] = sim_motion3D(Xs2,[],conn2,.001,1700,-dXs2,0);   % Simulate
[XCb,fCb] = sim_motion3D(Xs2,[],conn2,.02,50,dXs2b,0);    % Simulate
XC = cat(3,flip(XCa,3),XCb);
pIndsTr = [1 2 5; 2 3 5; 3 4 5; 1 4 5];

d1 = sqrt(squeeze(sum((diff(XC(:,[2,4],:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,[1 3],:),[],2)).^2)));
dInd = (abs(d1-0.58575)+abs(d2-0.58575));
dInd = find(dInd == min(dInd));
pIndL = [dInd 600 950 size(XCa,3) size(XC,3)];
nSc = 1/3;

plot3(d1,d2,zSh*ones(length(d1),1),'k-','linewidth',1);
for i = 1:length(pIndL)
    pInd = pIndL(i);
    vSh = [d1(pInd);d2(pInd);-zSh*.5];
    hold on;
    for j = 1:size(pIndsTr,1)
        patch(XC(1,pIndsTr(j,:),pInd)*nSc+vSh(1),...
              XC(2,pIndsTr(j,:),pInd)*nSc+vSh(2),...
              XC(3,pIndsTr(j,:),pInd)*nSc+vSh(3),...
              [160 220 255]/255);
    end
    visualize_network(XC(:,:,pInd)*nSc+vSh,[],conn2,.3,...
                      [255 100 100]/255, [255 100 100]/255);

end
hold on;
for i = [1 4]
    pInd = pIndL(i);
    vSh = [d1(pInd);d2(pInd);0];
    plot3(d1(pInd)*[1 1],d2(pInd)*[1 1],zSh*[1,-.5],'k--','linewidth',1);
    plot3(d1(pInd),d2(pInd),zSh,'k*','markersize',10,'linewidth',2);
end
for i = setdiff(1:length(pIndL),[1 4])
    pInd = pIndL(i);
    vSh = [d1(pInd);d2(pInd);0];
    plot3(d1(pInd)*[1 1],d2(pInd)*[1 1],zSh*[1,-.5],'k--','linewidth',1);
    plot3(d1(pInd),d2(pInd),zSh,'ko','markersize',3,'linewidth',3);
end

plot3(pLim,pLim,zSh*ones(2,1),'--','color',[200 200 200]/255);
hold off;
view(45,30);
axis([pLim,pLim,pLim-mean(pLim)]);
camlight(40, 30); lighting gouraud; material([.4 1 0]);

text(-0.05,0.98,'\textbf{b}','units','normalized','fontsize',10);
text(d1(size(XCa,3))-.10,d2(size(XCa,3))+.16,zSh,'$\mathrm{D_1^*}$','fontsize',10);
text(d1(dInd)+.2,d2(dInd)-.14,zSh,'$\mathrm{D_2^*}$','fontsize',10);

set(gca,'visible',1,'XTick',[],'YTick',[],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);
xlabel('$\mathrm{d_1}$','fontsize',10);
ylabel('$\mathrm{d_2}$','fontsize',10);


%% Size and Save Figure
fName = 'figure5supp';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.5 -.2 flw];
fig.PaperSize = [16.5 7.2];
saveas(fig, ['Figures/' fName], 'pdf');