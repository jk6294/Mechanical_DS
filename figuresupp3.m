% Specifics on Iterated Maps
clear; clc;
set(groot,'defaulttextinterpreter','latex');

fig = figure(3); clf;
flw = [19.8 7.8];
set(gcf, 'Renderer', 'painters', 'Position', [0 15 flw], 'Units', 'centimeters'); 
set(gcf, 'Renderer', 'painters', 'Position', [0 15 flw], 'Units', 'centimeters'); 

pSc = .8;
cBTr = [200 200 255]/255;
cBOp = [000 000 255]/255;
cB   = [100 100 255]/255;
cRTr = [255 200 200]/255;
cROp = [255 000 000]/255;
cR   = [255 100 100]/255;


%% a: Cobweb
subplot(1,2,1); cla;
s = sqrt(3);
pLim = [.3 1.9];
pSh = 0;
netSC = .7;

% Modules
load chaos_network.mat;
Xs1 = Xsc5(:,1:3,1);
Xu1 = Xuc5(:,1:2,1);
Xs1p = [Xs1(1,:); -Xs1(2,:)];
Xu1p = [Xu1(1,:); -Xu1(2,:)];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Map and map limits
d2 = maps(Xs1,Xu1);
d2f = matlabFunction(d2);
dkMax = double(fzero(matlabFunction(diff(d2)),1));
dMax = double(subs(d2,dkMax));
dMin = double(subs(d2,dMax));
r = [dMin + (dMax-dMin)*.001, dMax - (dMax-dMin)*.001];

% Lyapunov
nV = 50;
cRShT = linspace(.5,1,nV)';
cRShT = [cRShT ones(nV,1)*.5 flipud(cRShT)];
tV = linspace(0,.01,nV)+dMin+.008;
lyap_exp(d2,tV,[0 10],10,cRShT);

% Plot
hold on;
plot([0 3], [0 3], '--', 'color', [200 200 200]/255);
plot(linspace(dMin,dMax+.01,2000),d2f(linspace(dMin,dMax+.01,2000)),'k-','linewidth',1);
plot(mean(tV),pLim(1),'ks','markersize',8,'linewidth',1);
text(mean(tV)+.05,pLim(1)+.05,'$\mathrm{l_1}$','fontsize',10);
plot(r,pLim(2)-.1*[1 1],'k-','linewidth',1);
plot(r(1)*[1 1],(pLim(2)-.1)+[-.03 .03],'k-','linewidth',1);
plot(r(2)*[1 1],(pLim(2)-.1)+[-.03 .03],'k-','linewidth',1);
hold off;

axis([pLim+pSh pLim-pSh]);
set(gca,'visible',1,'XTick',[],'YTick',[],'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',[]);
drawnow;

% Labels
text(-0.05,0.98,'\textbf{a}','units','normalized','fontsize',10);
text(1,pLim(2)-.05,'r','fontsize',10);
text(1.5,.45,'$\mathrm{l_{k+1} = f(l_k)}$','fontsize',10);
text(1.52,1.6,'$\mathrm{l_{k+1} = l_k}$','fontsize',10,'color',[200 200 200]/255,'rotation',45);
xlabel('$\mathrm{l_k}$','fontsize',10);
ylabel('$\mathrm{l_{k+1}}$','fontsize',10);


%% b: Lyapunov Exponents
subplot(1,2,2); cla;

% Exponent Calculation
nV = 1000;
tV = linspace(r(1),r(2),nV);
d2 = maps(Xs1,Xu1);
l = lyap_exp(d2,tV,[10000,50000],1,[]);
plot(tV,l, 'ko','markersize',1,'linewidth',1);
axis([r mean(l)+std(l)*[-8 8]]);
set(gca,'visible',1,'XTick',[],'YTick',round(mean(l)+std(l)*[-4 0 4],3),'fontsize',10,...
                    'XTickLabel',[],'YTickLabel',round(mean(l)+std(l)*[-4 0 4],3));
text(-0.05,0.98,'\textbf{b}','units','normalized','fontsize',10);
xlabel('r');
ylabel('$\lambda$');
drawnow;


%% Size and Save Figure
fName = 'figure3supp';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [-1.5 -.2 flw];
fig.PaperSize = [16.5 7.2];
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl'); 