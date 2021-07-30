% Figure 4: Stability phase diagram
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
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
% Adjust Position
for i = 1:3
    if(i==2)
        subp(i,:) = subp(i,:) + [ fMarg(1) fMarg(3)...
                                 -sum(fMarg(1:2)) -sum(fMarg(3:4))]*2.5;
    else
        subp(i,:) = subp(i,:) + [ fMarg(1) fMarg(3)...
                                 -sum(fMarg(1:2)) -sum(fMarg(3:4))];
    end
end
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');


%% Figure Parameters
FS = 10;            % Fontsize
gr = 0.9;           % Gray
al = 0.2;           % Alpha
    
% Default Name-Value Pairs
NVTitle  = {'Units','centimeters','fontsize',FS};
NVTextH  = {'Units','Normalized','fontsize',FS,...
           'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,...
            'HorizontalAlignment','right'};
NVTextR  = {'Units','Normalized','fontsize',FS};

% Color
nT = 2000;
o = [1 1 1];
% Gradient: Slope. Interpolate between 2 colors
DLin = linspace(0,1,nT);
C3a = [197 066 085]/255;
C3b = [031 172 204]/255;
CP3 = interp1([0 1],[C3a;C3b],linspace(0,1,nT));


%% Construct Network
s = sqrt(3);
Xs = [-s/2 0 s/2;...                    % Initial position
      -1/2 1 -1/2];
rM = 1.5;
Xf = rM*Xs;                             % Final position
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];  % Network connectivity


%% Create MATLAB Function for Calculating Slope
% Define variables and constants
dx1 = [-1/2;-s/2];
dx2 = [0;0];
syms('th1s','real')
syms('th2s','real');
x4 = rM*[cos(th1s);sin(th1s)];
x5 = rM*[cos(th2s);sin(th2s)];

% Node 4 motion
M4 = (Xs(:,1:2)-x4)';
b4 = sum((Xs(:,1:2) - x4).*[dx1 dx2])';
dx4 = M4\b4;
% Node 5 motion
M5 = (Xs(:,1:2)-x5)';
b5 = sum((Xs(:,1:2) - x5).*[dx1 dx2])';
dx5 = M5\b5;
% Node 3 motion
M3 = (Xs(:,3)-[x4 x5])';
b3 = sum((Xs(:,3)-[x4 x5]) .* [dx4 dx5])';
dx3 = M3\b3;
% Slope
d2sym = simplify([1/2 -s/2]*dx3);
d2symf = matlabFunction(d2sym);

% Sample phase space
nS = 720;                                           % Number of samples
th = linspace(pi/2,-pi*3/2,nS);                     % Angle range
Xua = @(t1,t2) rM*[cos([t1 t2]);sin([t1 t2])];      % Added node positions
[T1,T2] = meshgrid(th,th);
dSl = arrayfun(d2symf,T1,T2);                       % Slope at all angles
dSl(isinf(dSl)) = nan;
dSl(1:1+nS:end)=nan;


%% a: Initial Design Visualization
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .10;                       % Scale network
sh = [.18;.63];                 % Shift network
R = [-2 2; -2 2];               % Range of conic visualization
mSs = 5;                        % Node size: start nodes
mSf = 3;                        % Node size: end nodes

% Draw legend
xst = .00; yst = .94;
ysh = .04;
visualize_network([xst;yst],[],[1 1],'msize',mSs);
visualize_network([xst;yst-ysh],[],[1 1],'msize',mSf,'scolor',o);
visualize_network([xst+.16;yst],[],[1 1],'msize',mSs,...
                  'scolor',interp1(DLin,CP3,tanh(0)));
visualize_network([xst+.20;yst],[],[1 1],'msize',mSs,...
                  'scolor',interp1(DLin,CP3,tanh(1)));
visualize_network([xst+.24;yst],[],[1 1],'msize',mSs,...
                  'scolor',interp1(DLin,CP3,1));
line([0 .07]+.152,[yst yst]-ysh,'linewidth',1,'color',o*.7);

% Obtain conic coordinates
thL1 = linspace(90,-60,nS);
thL2 = linspace(90,-120,nS);
Xu = rM*[cosd([thL1(end) thL2(end)]);sind([thL1(end) thL2(end)])];

% Draw solution space
visualize_conic_finite(Xs*sc+sh,Xf*sc+sh,R*sc+sh,'ucolori',o*.7,...
                       'overlay',.99);
% Draw line markers
plot([0 0]*sc+sh(1), [0 rM/2.5]*sc+sh(2), 'k-', 'linewidth', .5);
plot([0 Xu(1,1)]*sc+sh(1), [0 Xu(2,1)]*sc+sh(2), 'k:', 'linewidth', 1);
plot([0 Xu(1,2)]*sc+sh(1), [0 Xu(2,2)]*sc+sh(2), 'k:', 'linewidth', 1);
% Draw added nodes
visualize_network(Xu*sc+sh,[],[1 1],'msize',mSs,'scolor',...
                  interp1(DLin,CP3,tanh(1)));
% Draw angle measures
plot(.35*cosd(thL1)*sc+sh(1),0.35*sind(thL1)*sc+sh(2),'k-','linewidth',.5);
plot(.20*cosd(thL2)*sc+sh(1),0.20*sind(thL2)*sc+sh(2),'k-','linewidth',.5);
% Draw start and end nodes
visualize_network(Xs*sc+sh,[],[1 1],'msize',mSs);
visualize_network(Xf*sc+sh,[],[1 1],'scolor',o,'msize',mSf);

% Text
texta = '\textbf{a} \hspace{1.5mm}find solution space';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(xst+.07,yst,'start',NVTextR{:});
text(xst+.07,yst-ysh,'end',NVTextR{:});
text(1,yst,'added',NVTextRA{:});
text(1,yst-ysh,'solution',NVTextRA{:});
text(.5,.82,'add nodes at $\theta_4$ and $\theta_5$',NVTextH{:});
text(.60,.65,'$\theta_4$',NVTextR{:});
text(.45,.58,'$\theta_5$',NVTextR{:});


% b: Motion of Constructed Network
pInd = 1;
subplot('position',subpN(pInd,:)); hold on;

% Parameters
sc = .10;                           % Scale network
sh = [.18;.14];                    % Shift network
mSs = 5;                            % Node size: start nodes
mSf = 3;                            % Node size: end nodes

% Simulate network: small motion
[XC,fC] = sim_motion10(Xs,Xu,conn,.01,30,[Xs Xu],0);
XC = XC - XC(:,2,:) + XC(:,2,1);    % Correct offset

% Draw original and slightly moved network
visualize_network(XC(:,1:3,end)*sc+sh,XC(:,4:5,end)*sc+sh,conn,...
                  'msize',mSs,'scolor',o*gr.^.2,'ucolor',...
                  interp1(DLin,CP3,tanh(1)).^.2,'lalpha',.3);
visualize_network(XC(:,1:3,1)*sc+sh,XC(:,4:5,1)*sc+sh,conn,'msize',...
                  mSs,'ucolor',interp1(DLin,CP3,tanh(1)));
% final distance marker
line_coordinates(XC(:,1:2,end)*sc+sh,'lSh',.03,'nw',.007,'style','-',...
                 'lw',.5,'color',o*gr);
line_coordinates(XC(:,2:3,end)*sc+sh,'lSh',.03,'nw',.007,'style','-',...
                 'lw',.5,'color',o*gr);
% initial distance marker
line_coordinates(XC(:,1:2,1)*sc+sh,'lSh',.03,'nw',.007,'style','-',...
                 'lw',.5);
line_coordinates(XC(:,2:3,1)*sc+sh,'lSh',.03,'nw',.007,'style','-',...
                 'lw',.5);

% Text
textb = '\textbf{b} \hspace{1.5mm}construct network';
text(labX,subp(pInd,4)-5.8,textb,NVTitle{:});
text(.21,.2,'$d l_1$',NVTextH{:});
text(.76,.2,'$d l_2$',NVTextH{:});
text(.5,.31,'slope = $d l_2 / d l_1$',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% c: Draw Phase Diagrams
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = 20; 

% Draw phase diagram
colormap([1 1 1; CP3]);
imagesc(tanh(abs(dSl)),'alphadata',.6);

% Draw metastable points
C1 = contour(abs(dSl),[1 1],'color',interp1(DLin,CP3,tanh(1)),...
             'linewidth',2);
         
% Draw superstable points
dV2 = dSl.*(abs(dSl)<.5);
dV2(dV2(:)==0) = nan;
C0 = contour(dV2,[0 0],'color',interp1(DLin,CP3,tanh(0)),...
             'linewidth',1,'linestyle',':');
         
% Draw figure border
line([1 1 nS nS 1 1], [0 nS nS 1 1 nS],'linewidth',2,'color','k');

% Trace curve and draw networks: metastable
nP1 = [78 55 40];                               % Number of trace points
sI1 = cell(1,3);                                % Store trace points
nE1 = [6 4 3];                                  % Number of examples
rE1  = [1 7 3; nP1-[7 7 7]];                    % Range of examples
sh1 = cell(1,3);                                % Network shift
sL1 = ['ko';'ks';'kh'];                         % Marker symbol list
sh1{1} = [linspace(  0,  0,nE1(1));...
          linspace(-50,-50,nE1(1))];
sh1{2} = [linspace( 35, 45,nE1(2));...
          linspace(-25,-20,nE1(2))];
sh1{3} = [linspace( 50, 35,nE1(3));...
          linspace(-15,-35,nE1(3))];
sII1 = cat(3,[ 27  34 ; 694 687],...            % Initial points to trace
             [285 292 ; 474 477],...
             [234 239 ; 526 545]);
for i = 1:3
    % Trace curve
    [~,I] = min(sum(abs(sII1(:,1,i) - C1)));
    sIp = [C1(:,I) sII1(:,2,i) zeros(2,nP1(i)-1)];
    for j = 2:nP1(i)
        sloP = sIp(:,j)-sIp(:,j-1);
        [~,I] = min(sum(abs(sIp(:,j)+sloP*5/norm(sloP) - C1)));
        sIp(:,j+1) = C1(:,I);
    end
    sI1{i} = sIp;           % Store traced curve
    
    % Draw examples
    lx = round(sIp(:,floor(linspace(rE1(1,i),rE1(2,i),nE1(i)))));
    for j = 1:nE1(i)
        Xui = Xua(th(lx(1,j)),th(lx(2,j)));
        plot(lx(1,j),lx(2,j),sL1(i,:),'linewidth',.7,'markersize',5);
        visualize_network( Xs*sc+lx(:,j)+sh1{i}(:,j),...
                          Xui*sc+lx(:,j)+sh1{i}(:,j),conn,...
                          'ucolor',interp1(DLin,CP3,tanh(1)),'msize',3);
    end
end

% Trace curve and draw networks: superstable
nP0 = [90 40 57];                               % Number of trace points
sI0 = cell(1,3);                                % Store trace points
nE0 = [7 4 4];                                  % Number of examples
rE0  = [3 3 3; nP0-[2 2 2]];                    % Range of examples
sh0 = cell(1,3);                                % Network shift
sL0 = ['k+';'kx';'k*'];                         % Marker symbol list
sh0{1} = [linspace(-50,-20,nE0(1));...
          linspace(  0, 45^2,nE0(1)).^(1/2)];
sh0{2} = [linspace( 40, 40,nE0(2));...
          linspace(-20, -5,nE0(2))];
sh0{3} = [linspace(-45,-45,nE0(3));...
          linspace(-20, -5,nE0(3))];
sII0 = cat(3,[467 469 ; 057 067],...            % Initial points to trace
             [262 262 ; 047 052],...
             [700 700 ; 044 049]);
for i = 1:3
    % Trace curve
    [~,I] = min(sum(abs(sII0(:,1,i) - C0)));
    sIp = [C0(:,I) sII0(:,2,i) zeros(2,nP0(i)-1)];
    for j = 2:nP0(i)
        sloP = sIp(:,j)-sIp(:,j-1);
        [~,I] = min(sum(abs(sIp(:,j)+sloP*5/norm(sloP) - C0)));
        sIp(:,j+1) = C0(:,I);
    end
    sI0{i} = sIp;           % Store traced curve
    
    % Draw examples
    lx = round(sIp(:,floor(linspace(rE0(1,i),rE0(2,i),nE0(i)))));
    for j = 1:nE0(i)
        Xui = Xua(th(lx(1,j)),th(lx(2,j)));
        plot(lx(1,j),lx(2,j),sL0(i,:),'linewidth',.7,'markersize',5);
        visualize_network( Xs*sc+lx(:,j)+sh0{i}(:,j),...
                          Xui*sc+lx(:,j)+sh0{i}(:,j),conn,...
                          'ucolor',interp1(DLin,CP3,0),'msize',3);
    end
end

% Colorbar
delete(findall(gcf,'type','annotation'));
caxis([-.001 1]);
h = colorbar('location', 'east','position',[0.21 .1 .021 .805]);
a = annotation('textbox',h.Position,'FitBoxToText','off','FaceAlpha',...
               0.4,'EdgeColor',o*0,'BackgroundColor',o);
h.Ticks = tanh(1); h.TickLabels = {};

% Text
textc = '\textbf{c}\hspace{5.5mm}calculate slope at each added node position: $|d l_2/d l_1|$';
text(labX-.4,subp(pInd,4)+.31,textc,NVTitle{:});
text(.45,.41,'stable',NVTextR{:});
text(.525,.5,'unstable',NVTextR{:});
text(-.025,-.025,'0',NVTextR{:});
text(-0.055,0.99,'$2\pi$',NVTextR{:});
text(1.0,-.025,'$2\pi$',NVTextR{:});
text(.5,-0.036,'$\theta_4$',NVTextR{:});
text(-.045,.5,'$\theta_5$',NVTextR{:});
text(-.094,.03,'$0$',NVTextH{:},'color',interp1(DLin,CP3,tanh(0)));
text(-.140,.74,'$1$',NVTextH{:},'color',interp1(DLin,CP3,tanh(1)));
text(-.094,.97,'$\infty$',NVTextH{:},'color',interp1(DLin,CP3,1));
annotation('line',[0 .0195]+.2108, [.714 .714], 'Units','Normalize',...
           'color',interp1(DLin,CP3,tanh(1)),'linewidth',2);
annotation('line',[0 .0195]+.2108, [0 0]+.105, 'Units','Normalize',...
           'color',interp1(DLin,CP3,tanh(0)),'linewidth',2);
annotation(gcf,'textarrow',[1 1]*.2569, [1 1]*.38,'String',...
    '$|d l_2 / d l_1|$','HeadStyle','none','LineStyle','none',...
    NVTextH{:}, 'TextRotation',90,'interpreter','latex','color','k');
annotation(gcf,'textarrow',[1 1]*.25, [1 1]*.64,'String','stable',...
    'HeadStyle','none','LineStyle','none',NVTextH{:},'TextRotation',90,...
    'interpreter','latex','color','k');
annotation(gcf,'textarrow',[1 1]*.2603, [1 1]*.81,'String','unstable',...
    'HeadStyle','none','LineStyle','none',NVTextH{:},'TextRotation',90,...
    'interpreter','latex','color','k');

% Axes
axis([0 nS 0 nS]);
set(gca,'xtick',[],'ytick',[],'ydir','normal','visible',0);
drawnow;


%% d: Networks
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .055;                              % Scale network
shx = [.09 .39];                        % Shift: x
shy = [.79 .56 .33 .10];                % Shift: y
pR1 = [ 3 35; 40 70; 1 28; 1 40];       % Plot range: metastable
pR0 = [15 50; 55 85; 1 20; 5 30];       % Plot range: superstable
R = [-2 2; -2 2];                       % Range of solution space plot
lSh = .01162;                           % Distance bar shift
sL1 = ['ko';'ko';'ks';'kh'];            % Marker symbol list: metastable
sL0 = ['k+';'k+';'kx';'k*'];            % Marker symbol list: superstable

% Draw examples units: metastable
for i = 1:4
    sh = [shx(1);shy(i)];
    visualize_conic_finite(Xs*sc+sh,Xf*sc+sh,R*sc+sh,'ucolori',...
                           o*.8,'overlay',.99);
    lx = round(sI1{max(i-1,1)}(:,floor(linspace(pR1(i,1),pR1(i,2),5))));
    for j= 1:length(lx)
        Xui = Xua(th(lx(1,j)),th(lx(2,j)));
        if(j < length(lx))
            visualize_network(Xs*sc+sh,Xui*sc+sh,[1 1],...
                'ucolor',interp1(DLin,CP3,tanh(1)),'nalpha',.1*j,...
                'lalpha',.1*j,'msize',5);
        else
            visualize_network(Xs*sc+sh,Xui*sc+sh,conn,...
                'ucolor',interp1(DLin,CP3,tanh(1)),'msize',5);
        end
    end
    plot(sh(1)-.12,sh(2)+.109,sL1(i,:),'linewidth',.7,'markersize',5,...
         'clipping',0);
end

% Draw examples units: super
for i = 1:4
    sh = [shx(2);shy(i)];
    visualize_conic_finite(Xs*sc+sh,Xf*sc+sh,R*sc+sh,'ucolori',...
                           o*.8,'overlay',.99);
    lx = round(sI0{max(i-1,1)}(:,floor(linspace(pR0(i,1),pR0(i,2),5))));
    for j= 1:length(lx)
        Xui = Xua(th(lx(1,j)),th(lx(2,j)));
        if(i<3)
            XfM = mean(Xui,2); XfD = diff(Xui,1,2);
        else
            XfM = mean(Xs(:,2:3),2); XfD = diff(Xs(:,2:3),1,2)*1.6;
        end
        if(j < length(lx))
            line_coordinates(([-XfD XfD]*.65+XfM)*sc+sh,'lSh',-lSh,...
                             'style','-','lw',.5,'color',...
                             interp1(DLin,CP3,tanh(0)),'lalpha',.1*i);
            visualize_network(Xs*sc+sh,Xui*sc+sh,[1 1],...
                              'ucolor',interp1(DLin,CP3,tanh(0)),...
                              'nalpha',.1*j,'lalpha',.1*j,'msize',5);
        else
            line_coordinates(([-XfD XfD]*.65+XfM)*sc+sh,'lSh',-lSh,...
                             'style','-','lw',.5,'color',...
                             interp1(DLin,CP3,tanh(0)));
            visualize_network(Xs*sc+sh,Xui*sc+sh,conn,...
                'ucolor',interp1(DLin,CP3,tanh(0)),'msize',5);
        end
    end
    plot(sh(1)-.1,sh(2)+.112,sL0(i,:),'linewidth',.7,'markersize',5,...
         'clipping',0);
end

% Text
textd = '\textbf{d} \hspace{1.5mm}metastable';
text(labX-.25,subp(pInd,4)+labY,textd,NVTitle{:});
text(.16,.955,'$|d l_2 / d l_1|=1$',NVTextH{:});
text(-.02,.90,'symmetric 1',NVTextR{:});
text(-.02,.67,'symmetric 2',NVTextR{:});
text(-.02,.44,'asymmetric 1',NVTextR{:});
text(-.02,.21,'asymmetric 2',NVTextR{:});
texte = '\textbf{e} \hspace{1.5mm}superstable';
text(labX+2.5,subp(pInd,4)+labY,texte,NVTitle{:});
text(.79,.955,'$|d l_2 / d l_1|=0$',NVTextH{:});
text(.65,.90,'2 colinear',NVTextR{:});
text(.65,.67,'2 colinear',NVTextR{:});
text(.65,.44,'1 colinear',NVTextR{:});
text(.65,.21,'1 colinear',NVTextR{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% Save
fName = 'figure4c';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');