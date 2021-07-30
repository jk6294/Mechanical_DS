% Figure 4: Stability phase diagram
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(4); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 9.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  3.75  9.50];...
        [ 4.50  0.00  9.50  9.50];...
        [14.25  0.00  4.75  9.50]];
% Fontsize
FS = 10;
% Distance visualization parameters
nW = .08;
lw = .5;
gr = 0.8;
% Scaling parameters
scc = .25;
    
% Adjust Position
for i = 1:3
    if(i==2)
        subp(i,:) = subp(i,:) + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))]*2.5;
    else
        subp(i,:) = subp(i,:) + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
    end
end
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');

% Default Name-Value Pairs
NVTitle = {'Units','centimeters','fontsize',FS};
NVTextH = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','right'};
NVTextR = {'Units','Normalized','fontsize',FS};


%% Color
% Coefficients
nT = 720;
cR = .299;
cG = .587;
cB = .114;
br = .6;

% Gradient along red/blue, sweep red
DLin = linspace(0,1,nT);
RLin = .75*ones(1,nT);
BLin = linspace(0,1,nT);
GLin = sqrt((br^2 - cB*BLin.^2 - cR*RLin.^2)/cG);
CP = [RLin; GLin; BLin]';
CP([1 end],:)


%% a: Initial Design Visualization
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;

s = sqrt(3);
Xs = [-s/2 0 s/2;...
      -1/2 1 -1/2];
R = [-2 2; -2 2];
rM = 1.5;
Xf = rM*Xs;
sc = .10;

% Legend
xst = .00; yst = .94;
ysh = .04;
visualize_network([xst;yst],[],[1 1],'msize',5);
visualize_network([xst;yst-ysh],[],[1 1],'msize',5,'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);
visualize_network([xst+.16;yst],[],[1 1],'msize',5,'scolor',interp1(DLin,CP,tanh(0)));
visualize_network([xst+.20;yst],[],[1 1],'msize',5,'scolor',interp1(DLin,CP,tanh(1)));
visualize_network([xst+.24;yst],[],[1 1],'msize',5,'scolor',interp1(DLin,CP,1));
line([0 .07]+.15,[yst yst]-ysh,'linewidth',1,'color',[1 1 1]*.7);
text(xst+.07,yst,'start',NVTextR{:});
text(xst+.07,yst-ysh,'end',NVTextR{:});
text(1,yst,'added',NVTextRA{:});
text(1,yst-ysh,'solution',NVTextRA{:});

% Obtain conic coordinates
thL1 = linspace(90,-60,nT);
thL2 = linspace(90,-120,nT);
Xu = rM*[cosd([thL1(end) thL2(end)]);sind([thL1(end) thL2(end)])];

% Row 1
sh = [.18;.63];
text(.5,.82,'add nodes at $\theta_1$ and $\theta_2$',NVTextH{:});
visualize_conic_finite(Xs*sc+sh,Xf*sc+sh,R*sc+sh,'ucolori',[1 1 1]*.7,'overlay',.99);
hold on;
% line markers
plot([0 0]*sc+sh(1), [0 rM/2.5]*sc+sh(2), 'k-', 'linewidth', .5);
plot([0 Xu(1,1)]*sc+sh(1), [0 Xu(2,1)]*sc+sh(2), 'k:', 'linewidth', 1);
plot([0 Xu(1,2)]*sc+sh(1), [0 Xu(2,2)]*sc+sh(2), 'k:', 'linewidth', 1);
% angle markers
visualize_network(Xu*sc+sh,[],[1 1],'msize',5,'scolor',interp1(DLin,CP,tanh(1)));
% angle measures
plot(.35*cosd(thL1)*sc+sh(1), .35*sind(thL1)*sc+sh(2), 'k-', 'linewidth', .5);
plot(.20*cosd(thL2)*sc+sh(1), .2*sind(thL2)*sc+sh(2), 'k-', 'linewidth', .5);
% Network
visualize_network(Xs*sc+sh,[],[1 1],'msize',5);
visualize_network(Xf*sc+sh,[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8,'msize',5);
hold off;
% line_coordinates(Xs(:,1:2)*sc+sh,'lSh',.03,'nw',.007,'style','-','lw',.5);
% line_coordinates(Xs(:,2:3)*sc+sh,'lSh',.03,'nw',.007,'style','-','lw',.5);

% Text
text(labX,subp(pInd,4)+labY,'\textbf{a} \hspace{1.5mm}\textit{find solution space}',NVTitle{:});
text(.59,.595,'$\theta_1$',NVTextR{:});
text(.45,.58,'$\theta_2$',NVTextR{:});


% Obtain conic coordinates
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
sh3 = [.18;.14];
% Small motion
[XC,fC] = sim_motion(Xs,Xu,conn,.01,30,[Xs Xu],0);
XC = XC - XC(:,2,:) + XC(:,2,1);

visualize_network(XC(:,1:3,end)*sc+sh3,XC(:,4:5,end)*sc+sh3,conn,'msize',5,'ucolor',interp1(DLin,CP,tanh(1)),...
                  'nalpha',.3,'lalpha',.3);
visualize_network(XC(:,1:3,1)*sc+sh3,XC(:,4:5,1)*sc+sh3,conn,'msize',5,'ucolor',interp1(DLin,CP,tanh(1)));
% initial distance marker
line_coordinates(XC(:,1:2,1)*sc+sh3,'lSh',.03,'nw',.007,'style','-','lw',.5);
line_coordinates(XC(:,2:3,1)*sc+sh3,'lSh',.03,'nw',.007,'style','-','lw',.5);
% final distance marker
line_coordinates(XC(:,1:2,end)*sc+sh3,'lSh',.03,'nw',.007,'style','-','lw',.5,'color',[1 1 1]*.3);
line_coordinates(XC(:,2:3,end)*sc+sh3,'lSh',.03,'nw',.007,'style','-','lw',.5,'color',[1 1 1]*.3);
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);

% Text
text(labX,subp(pInd,4)-5.8,'\textbf{b} \hspace{1.5mm}\textit{construct network}',NVTitle{:});
text(.21,.2,'$\delta d_1$',NVTextH{:});
text(.76,.2,'$\delta d_2$',NVTextH{:});
text(.5,.31,'slope = $\delta d_2 / \delta d_1$',NVTextH{:});


%% c: Phase Diagram Evalulate slope at initial position
% angle list
th1 = linspace(90,-270,nT);
th2 = linspace(90,-270,nT);
Xua = zeros(2,2,nT,nT);

% Change in distances
d1dot1 = zeros(nT);
d2dot1 = zeros(nT);

UssM1 = zeros([size(Xs),nT,nT]);
th = linspace(0,360,nT+1); th = th(1:end-1);
tic
fprintf([repmat('.',1,nT) '\n\n']);
parfor i = 1:nT
    fprintf('\b=\n');
    for j = 1:nT
        if i ~= j
            Xu = rM*[cosd([th1(i) th2(j)]);sind([th1(i) th2(j)])];
            Xua(:,:,i,j) = Xu;
            [Uss,Uus,~] = construct_motion(Xs,Xs,Xu,conn,0,0);
            UssM1(:,:,i,j) = Uss;
            d1dot1(i,j) = (Xs(:,1)-Xs(:,2))' * (Uss(:,1)-Uss(:,2));
            d2dot1(i,j) = (Xs(:,2)-Xs(:,3))' * (Uss(:,2)-Uss(:,3));
        end
    end
end
toc
disp('done');


%% c: Draw Phase Diagrams
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;

nScb = 20; 
nSh1d = [0; -50];

dRat1 = d2dot1./d1dot1;
dRat1M = abs(dRat1);
colormap([1 1 1; CP]);

dRat2 = dRat1.*(dRat1M<.5);
dRat2(dRat2(:)==0) = nan;
dRat3 = dRat1.*(dRat1M>.5);
dRat3(dRat3(:)==0) = nan;

imagesc(tanh((dRat1M)),'alphadata',.6);
hold on;

% Metastable
C1 = contour(dRat1M,[1 1],'color',interp1(DLin,CP,tanh(1)),'linewidth',2);
% Superstable
C0 = contour(dRat2,[0 0],'color',interp1(DLin,CP,tanh(0)),'linewidth',1,'linestyle',':');
% Boundaries
line([1 1 nT nT 1 1], [0 nT nT 1 1 nT],'linewidth',2,'color','k');


% Trace slope 1 line 1
lx = round(linspace(27,nT/2-50,6)); lx = [lx; nT - lx+1];
sInd11 = round(linspace(27,nT/2-50,100)); sInd11 = [sInd11; nT - sInd11+1];
nSh1b = [0;-50];
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    plot(lx(1,i),lx(2,i),'ko','linewidth',.7,'markersize',5);
    visualize_network( Xs*nScb+lx(:,i)+nSh1b,...
                      Xui*nScb+lx(:,i)+nSh1b,conn,'ucolor',interp1(DLin,CP,tanh(1)),'msize',3);
end

% Trace slope 1 line 2
nP = 55; nJ = 5;
sInd12 = zeros(2,nP+1);
sInd12(:,1:2) = [[285;474] [292;477]];
[~,I] = min(sum(abs(sInd12(:,1) - C1)));
sInd12(:,1) = C1(:,I);
for i = 2:nP
    sloPr = sInd12(:,i)-sInd12(:,i-1);
    [~,I] = min(sum(abs(sInd12(:,i)+sloPr*nJ/norm(sloPr) - C1)));
    sInd12(:,i+1) = C1(:,I);
end
lx = round(sInd12(:,floor(linspace(7,nP-7,4))));
nSh1b = [linspace(35,45,length(lx));...
         linspace(-25,-20,length(lx))];
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    plot(lx(1,i),lx(2,i),'ks','linewidth',.7,'markersize',5);
    visualize_network( Xs*nScb+lx(:,i)+nSh1b(:,i),...
                      Xui*nScb+lx(:,i)+nSh1b(:,i),conn,'ucolor',interp1(DLin,CP,tanh(1)),'msize',3);
end

% Trace slope 1 line 3
nP = 40; nJ = 5;
sInd13 = zeros(2,nP+1);
sInd13(:,1:2) = [[234;526] [239;545]];
[~,I] = min(sum(abs(sInd13(:,1) - C1)));
sInd13(:,1) = C1(:,I);
for i = 2:nP
    sloPr = sInd13(:,i)-sInd13(:,i-1);
    [~,I] = min(sum(abs(sInd13(:,i)+sloPr*nJ/norm(sloPr) - C1)));
    sInd13(:,i+1) = C1(:,I);
end
lx = round(sInd13(:,floor(linspace(3,nP-7,3))));
nSh1b = [linspace(50,35,length(lx));...
         linspace(-15,-35,length(lx))];
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    plot(lx(1,i),lx(2,i),'kh','linewidth',.7,'markersize',5);
    visualize_network( Xs*nScb+lx(:,i)+nSh1b(:,i),...
                      Xui*nScb+lx(:,i)+nSh1b(:,i),conn,'ucolor',interp1(DLin,CP,tanh(1)),'msize',3);
end


% Trace slope 0 line 1
nP = 90; nJ = 5;
sInd01 = zeros(2,nP+1);
sInd01(:,1:2) = [[467;57] [469;67]];
[~,I] = min(sum(abs(sInd01(:,1) - C0)));
sInd01(:,1) = C0(:,I);
for i = 2:nP
    sloPr = sInd01(:,i)-sInd01(:,i-1);
    [~,I] = min(sum(abs(sInd01(:,i)+sloPr*nJ/norm(sloPr) - C0)));
    sInd01(:,i+1) = C0(:,I);
end
lx = round(sInd01(:,floor(linspace(3,nP-2,7))));
nSh1b = [linspace(-50,-20,length(lx));...
         linspace(0,45^2,length(lx)).^(1/2)];
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    plot(lx(1,i),lx(2,i),'k+','linewidth',.7,'markersize',5);
    visualize_network( Xs*nScb+lx(:,i)+nSh1b(:,i),...
                      Xui*nScb+lx(:,i)+nSh1b(:,i),conn,'ucolor',interp1(DLin,CP,tanh(0)),'msize',3);
end

% Trace slope 0 line 2
nP = 40; nJ = 5;
sInd02 = zeros(2,nP+1);
sInd02(:,1:2) = [[262;47] [262;52]];
[~,I] = min(sum(abs(sInd02(:,1) - C0)));
sInd02(:,1) = C0(:,I);
for i = 2:nP
    sloPr = sInd02(:,i)-sInd02(:,i-1);
    [~,I] = min(sum(abs(sInd02(:,i)+sloPr*nJ/norm(sloPr) - C0)));
    sInd02(:,i+1) = C0(:,I);
end
lx = round(sInd02(:,floor(linspace(3,nP-2,4))));
nSh1b = [linspace(40,40,length(lx));...
         linspace(-20,-5,length(lx))];
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    plot(lx(1,i),lx(2,i),'kx','linewidth',.7,'markersize',5);
    visualize_network( Xs*nScb+lx(:,i)+nSh1b(:,i),...
                      Xui*nScb+lx(:,i)+nSh1b(:,i),conn,'ucolor',interp1(DLin,CP,tanh(0)),'msize',3);
end


% Trace slope 0 line 3
nP = 57; nJ = 5;
sInd03 = zeros(2,nP+1);
sInd03(:,1:2) = [[700;44] [700;49]];
[~,I] = min(sum(abs(sInd03(:,1) - C0)));
sInd03(:,1) = C0(:,I);
for i = 2:nP
    sloPr = sInd03(:,i)-sInd03(:,i-1);
    [~,I] = min(sum(abs(sInd03(:,i)+sloPr*nJ/norm(sloPr) - C0)));
    sInd03(:,i+1) = C0(:,I);
end
lx = round(sInd03(:,floor(linspace(3,nP-2,4))));
nSh1b = [linspace(-45,-45,length(lx));...
         linspace(-20,-5,length(lx))];
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    plot(lx(1,i),lx(2,i),'k*','linewidth',.7,'markersize',5);
    visualize_network( Xs*nScb+lx(:,i)+nSh1b(:,i),...
                      Xui*nScb+lx(:,i)+nSh1b(:,i),conn,'ucolor',interp1(DLin,CP,tanh(0)),'msize',3);
end

hold off;
caxis([-.001 1]);
set(gca,'xtick',[],'ytick',[],'ydir','normal','visible',0);
h = colorbar('location', 'east','position',[0.21 .1 .021 .805]);
delete(findall(gcf,'type','annotation'));
a = annotation('textbox',h.Position,'FitBoxToText','off','FaceAlpha',0.4,...
    'EdgeColor',[1 1 1]*0,'BackgroundColor',[1 1 1]);
h.Ticks = tanh(1); h.TickLabels = {};
annotation(gcf,'textarrow',[1 1]*.2603, [1 1]*.5,'String','$|\delta d_2 / \delta d_1|$',...
    'HeadStyle', 'none', 'LineStyle', 'none',NVTextH{:}, 'TextRotation',90,...
    'interpreter','latex','color','k');

% Axis Ticks
text(-.025,-.025,'0',NVTextR{:});
text(-0.055,0.99,'$2\pi$',NVTextR{:});
text(1.0,-.025,'$2\pi$',NVTextR{:});
text(.5,-0.036,'$\theta_1$',NVTextR{:});
text(-.05,.5,'$\theta_2$',NVTextR{:});
% Colorbar Ticks
text(-.094,.03,'$0$',NVTextH{:},'color',interp1(DLin,CP,tanh(0)));
text(-.140,.74,'$1$',NVTextH{:},'color',interp1(DLin,CP,tanh(1)));
text(-.094,.97,'$\infty$',NVTextH{:},'color',interp1(DLin,CP,1));
annotation('line',[0 .0195]+.2108, [.714 .714], 'Units','Normalize',...
           'color',interp1(DLin,CP,tanh(1)),'linewidth',2);
annotation('line',[0 .0195]+.2108, [0 0]+.105, 'Units','Normalize',...
           'color',interp1(DLin,CP,tanh(0)),'linewidth',2);
% Description
text(labX-.4,subp(pInd,4)+.33,'\textbf{c} \hspace{4mm}\textit{calculate slope at each added node position:} $|\delta d_2/\delta d_1|$',NVTitle{:});
text(.45,.41,'stable',NVTextR{:});
text(.525,.5,'unstable',NVTextR{:});


%% d: Networks
% Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
sc = .055;
sha = [.0;0];
sh1 = [.09;.79];
sh2 = [.09;.56];
sh3 = [.09;.33];
sh4 = [.09;.10];
sh5 = [.39;.79];
sh6 = [.39;.56];
sh7 = [.39;.33];
sh8 = [.39;.10];
R = [-2 2; -2 2];
CR = [.9 .7 .5];
lSh = .01162;

% Metastable
% Symmetric 1
visualize_conic_finite(Xs*sc+sh1,Xf*sc+sh1,R*sc+sh1,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd11(:,floor(linspace(1,50,5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    if(i < length(lx))
        visualize_network(Xs*sc+sh1+(i-1)*sha,Xui*sc+sh1+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(1)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        visualize_network(Xs*sc+sh1+(i-1)*sha,Xui*sc+sh1+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(1)),'msize',5);
    end
end
plot(sh1(1)-.12,sh1(2)+.109,'ko','linewidth',.7,'markersize',5,'clipping',0);
% Symmetric 2
visualize_conic_finite(Xs*sc+sh2,Xf*sc+sh2,R*sc+sh2,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd11(:,floor(linspace(51,100,5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    if(i < length(lx))
        visualize_network(Xs*sc+sh2+(i-1)*sha,Xui*sc+sh2+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(1)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        visualize_network(Xs*sc+sh2+(i-1)*sha,Xui*sc+sh2+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(1)),'msize',5);
    end
end
plot(sh2(1)-.12,sh2(2)+.109,'ko','linewidth',.7,'markersize',5,'clipping',0);
% Asymmetric 1
visualize_conic_finite(Xs*sc+sh3,Xf*sc+sh3,R*sc+sh3,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd12(:,floor(linspace(1,28,5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    if(i < length(lx))
        visualize_network(Xs*sc+sh3+(i-1)*sha,Xui*sc+sh3+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(1)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        visualize_network(Xs*sc+sh3+(i-1)*sha,Xui*sc+sh3+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(1)),'msize',5);
    end
end
plot(sh3(1)-.12,sh3(2)+.109,'ks','linewidth',.7,'markersize',5,'clipping',0);
% Asymmetric 2
visualize_conic_finite(Xs*sc+sh4,Xf*sc+sh4,R*sc+sh4,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd13(:,floor(linspace(1,size(sInd13,2),5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    if(i < length(lx))
        visualize_network(Xs*sc+sh4+(i-1)*sha,Xui*sc+sh4+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(1)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        visualize_network(Xs*sc+sh4+(i-1)*sha,Xui*sc+sh4+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(1)),'msize',5);
    end
end
plot(sh4(1)-.12,sh4(2)+.109,'kh','linewidth',.7,'markersize',5,'clipping',0);


% Superstable
% 2 colinear 1
visualize_conic_finite(Xs*sc+sh5,Xf*sc+sh5,R*sc+sh5,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd01(:,floor(linspace(size(sInd01,2)/1.7,size(sInd01,2)/1.1,5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    XfM = mean(Xui,2); XfD = diff(Xui,1,2);
    if(i < length(lx))
        line_coordinates(([-XfD XfD]*.65+XfM)*sc+sh5,'lSh',-lSh,...
                         'style','-','lw',.5,'color',interp1(DLin,CP,tanh(0)),'lalpha',.1*i);
        visualize_network(Xs*sc+sh5+(i-1)*sha,Xui*sc+sh5+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(0)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        line_coordinates(([-XfD XfD]*.65+XfM)*sc+sh5,'lSh',-lSh,...
                         'style','-','lw',.5,'color',interp1(DLin,CP,tanh(0)));
        visualize_network(Xs*sc+sh5+(i-1)*sha,Xui*sc+sh5+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(0)),'msize',5);
    end
end
visualize_network(Xs(:,1)*sc+sh5+(i-1)*sha,[],[1 1],'msize',5);
plot(sh5(1)-.1,sh5(2)+.112,'k+','linewidth',.7,'markersize',5,'clipping',0);
% 2 colinear 2
visualize_conic_finite(Xs*sc+sh6,Xf*sc+sh6,R*sc+sh6,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd01(:,floor(linspace(15,size(sInd01,2)/1.7,5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    XfM = mean(Xui,2); XfD = diff(Xui,1,2);
    if(i < length(lx))
        line_coordinates(([-XfD XfD]*.65+XfM)*sc+sh6,'lSh',-lSh,...
                         'style','-','lw',.5,'color',interp1(DLin,CP,tanh(0)),'lalpha',.1*i);
        visualize_network(Xs*sc+sh6+(i-1)*sha,Xui*sc+sh6+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(0)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        line_coordinates(([-XfD XfD]*.65+XfM)*sc+sh6,'lSh',-lSh,...
                         'style','-','lw',.5,'color',interp1(DLin,CP,tanh(0)));
        visualize_network(Xs*sc+sh6+(i-1)*sha,Xui*sc+sh6+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(0)),'msize',5);
    end
end
visualize_network(Xs(:,1)*sc+sh6+(i-1)*sha,[],[1 1],'msize',5);
plot(sh6(1)-.1,sh6(2)+.112,'k+','linewidth',.7,'markersize',5,'clipping',0);
% 1 colinear 1
visualize_conic_finite(Xs*sc+sh7,Xf*sc+sh7,R*sc+sh7,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd02(:,floor(linspace(1,size(sInd02,2)/2,5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    if(i < length(lx))
        visualize_network(Xs*sc+sh7+(i-1)*sha,Xui*sc+sh7+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(0)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        XfM = mean(Xs(:,2:3),2); XfD = diff(Xs(:,2:3),1,2);
        line_coordinates(([-XfD XfD]*1+XfM)*sc+sh7,'lSh',-lSh,...
                         'style','-','lw',.5,'color',interp1(DLin,CP,tanh(0)));
        visualize_network(Xs*sc+sh7+(i-1)*sha,Xui*sc+sh7+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(0)),'msize',5);
    end
end
visualize_network(Xs(:,2:3)*sc+sh7+(i-1)*sha,[],[1 1],'msize',5);
plot(sh7(1)-.1,sh7(2)+.112,'kx','linewidth',.7,'markersize',5,'clipping',0);
% 1 colinear 2
visualize_conic_finite(Xs*sc+sh8,Xf*sc+sh8,R*sc+sh8,'ucolori',[1 1 1]*.8,'overlay',.99);
lx = round(sInd03(:,floor(linspace(5,size(sInd03,2)/2,5))));
for i = 1:length(lx)
    Xui = Xua(:,:,lx(1,i),lx(2,i));
    if(i < length(lx))
        visualize_network(Xs*sc+sh8+(i-1)*sha,Xui*sc+sh8+(i-1)*sha,[1 1],...
            'ucolor',interp1(DLin,CP,tanh(0)),'nalpha',.1*i,'lalpha',.1*i,'msize',5);
    else
        XfM = mean(Xs(:,2:3),2); XfD = diff(Xs(:,2:3),1,2);
        line_coordinates(([-XfD XfD]*1+XfM)*sc+sh8,'lSh',-lSh,...
                         'style','-','lw',.5,'color',interp1(DLin,CP,tanh(0)));
        visualize_network(Xs*sc+sh8+(i-1)*sha,Xui*sc+sh8+(i-1)*sha,conn,...
            'ucolor',interp1(DLin,CP,tanh(0)),'msize',5);
    end
end
visualize_network(Xs(:,2:3)*sc+sh8+(i-1)*sha,[],[1 1],'msize',5);
plot(sh8(1)-.1,sh8(2)+.112,'k*','linewidth',.7,'markersize',5,'clipping',0);



axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);

% Text
text(labX-.25,subp(pInd,4)+labY,'\textbf{d} \hspace{1.5mm}\textit{metastable}',NVTitle{:});
text(.16,.955,'$|\delta d_2 / \delta d_1|=1$',NVTextH{:});
text(-.02,.90,'symmetric 1',NVTextR{:});
text(-.02,.67,'symmetric 2',NVTextR{:});
text(-.02,.44,'asymmetric 1',NVTextR{:});
text(-.02,.21,'asymmetric 2',NVTextR{:});
text(labX+2.5,subp(pInd,4)+labY,'\textbf{e} \hspace{1.5mm}\textit{superstable}',NVTitle{:});
text(.79,.955,'$|\delta d_2 / \delta d_1|=0$',NVTextH{:});
text(.65,.90,'2 colinear',NVTextR{:});
text(.65,.67,'2 colinear',NVTextR{:});
text(.65,.44,'1 colinear',NVTextR{:});
text(.65,.21,'1 colinear',NVTextR{:});


%% Save
fName = 'figure4b';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 19 9.5];
fig.PaperSize = [19 9.5];
saveas(fig, ['Figures/' fName], 'pdf');


