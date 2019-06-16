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
pLim = [.6 2];
pSh = 0;
netSC = .7;

% Modules
Xs1 = [-s/2  0.0  s/2;...
        -0.5  0.0 -0.5]*1.5;
Xu1 = [ 0.15 -0.22;...
       -0.50 -0.90]*1.5;
Xs1p = [Xs1(1,:); -Xs1(2,:)];
Xu1p = [Xu1(1,:); -Xu1(2,:)];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

% Combined
Xs1c = [-s/2  0  s/2  s;...
        -1/2  0 -1/2  0]*1.5;
Xu1c = [Xu1 Xu1p+[s;-1]*1.5/2];
conn1c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs1aa,Xu1aa,conn1aa] = tesselate_network_old(Xs1c,Xu1c,conn1c,[1.5*s;0],[4;1]);

% Map
nV = 50;
cRShT = linspace(.5,1,nV)';
cRShT = [cRShT ones(nV,1)*.5 flipud(cRShT)];
d2 = maps(Xs1,Xu1);
tV = linspace(1,1.005,nV)-.241;
lyap_exp(d2,tV,[0 11],11,cRShT);
% Find Map Range
dkMax = double(fzero(matlabFunction(diff(d2)),1));
dMax = double(subs(d2,dkMax));
dMin = double(subs(d2,dMax));
r = [dMin, dMax];
% Plot
hold on;
plot([0 3], [0 3], '--', 'color', [200 200 200]/255);
plot(mean(tV),pLim(1),'ks','markersize',8,'linewidth',1);
text(mean(tV)+.05,pLim(1)+.05,'$\mathrm{d_1}$','fontsize',10);
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
text(mean(pLim),pLim(2)-.05,'r','fontsize',10);
text(1.46,.67,'$\mathrm{d_{k+1} = f(d_k)}$','fontsize',10);
text(.9,.86,'$\mathrm{d_{k+1} = d_k}$','fontsize',10,'color',[200 200 200]/255);
xlabel('$\mathrm{d_k}$','fontsize',10);
ylabel('$\mathrm{d_{k+1}}$','fontsize',10);


%% b: Lyapunov Exponents
subplot(1,2,2); cla;

% Exponent Calculation
nV = 1000;
tV = linspace(r(1),r(2),nV);
Xs1 = [-s/2  0.0  s/2;...
        -0.5  0.0 -0.5]*1.5;
Xu1 = [ 0.15 -0.22;...
       -0.50 -0.90]*1.5;
d2 = maps(Xs1,Xu1);
l = lyap_exp(d2,tV,[5000,20000],1,[]);
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