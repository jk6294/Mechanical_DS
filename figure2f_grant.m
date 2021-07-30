% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(2); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 7.5];
% Margins in cm, [l,r,d,u]
fMarg = [.4 .0 .4 .0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 2.75  4.75 4.75];...
        [ 0.00 0.00  4.75 2.75];....
        [ 5.00 0.00 14.00 7.50]];
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
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
al = 0.3;           % Alpha

% Default Name-Value Pairs
NVTitle  = {'Units','centimeters','fontsize',FS};
NVTextH  = {'Units','Normalized','fontsize',FS,...
            'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,...
            'HorizontalAlignment','right'};
NVTextR  = {'Units','Normalized','fontsize',FS};

% Color
nT = 1000;
o = [1 1 1];
% Gradient: Distance. Interpolate between 3 colors
DLin = linspace(sqrt(2)-.01,2+.01,nT);
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));


%% Define module
% Parameters
L = 1;
Rz = rotz(45); Rz = Rz(1:2,1:2);

% Single module
Xs = [0  0  0  2*L;...
      L -L  0  0];
conn = [1 3; 1 4; 2 3; 2 4];

% Double module
Xsc = [0  0  0  2*L  L  L;...
       L -L  0  0    0 -2*L];
connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];

% Tesselate module
[Xsa, conna] = tesselate_network(Rz*Xsc, connc, [sqrt(2);0], [4;1]);


%% Simulate networks
% Initial conditions
Xs0 = [-1 -1 0 0; 0 0 0 0];
Xsc0 = [0 0 0 0 0  0;...
        0 0 0 0 -1 0];
Xsa0 = 0*Xsa; Xsa0(1,end-1) = 1;

% Simulate networks
[XC,fC] = sim_motion(Xs,[],conn,.001,756,Xs0,0);        % 1 unit
[XCc,fCc] = sim_motion(Xsc,[],connc,.001,1122,Xsc0,0);  % 2 units
[XCa,fCa] = sim_motion(Xsa,[],conna,.001,4379,Xsa0,0);  % Chain

% Simulation errors
disp(['mean simulation error: '...
     [num2str(mean(fC)) '  ' num2str(mean(fCc))] '  ' num2str(mean(fCc))]);

 
% Translate networks to align with origin
for i = 1:size(XC,3)
    XC(:,:,i) = XC(:,:,i) - XC(:,4,i);
end
for i = 1:size(XCa,3)
    XCa(:,:,i) = XCa(:,:,i) - XCa(:,1,i);
end

% Calculate distances between nodes
D = squeeze(sqrt(sum(diff(XC,1,2).^2)));                % 1 unit
D = D([1,3],:);
Dc = squeeze(sqrt(sum(diff(XCc,1,2).^2)));              % 2 units
Dc = Dc(1:2:end,:);
Da = squeeze(sqrt(sum(diff(XCa,1,2).^2)));              % CHain
Da = Da(1:2:end,:);

% Colors
CPI = interp1(DLin,CP1,D(1,:));
CPJ = interp1(DLin,CP1,D(2,:));


%% a: single module conformational motion
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .075;          % Scale drawing
aX0 = 1.27;         % Map axes start
aXF = 2.07;         % Map axes end

% Draw map ticks
line(sqrt([2 2]), [0 .04]+aX0, 'color', interp1(DLin,CP1,sqrt(2)),...
     'linewidth', .5);
line([2 2], [0 .04]+aX0, 'color', interp1(DLin,CP1,2),...
     'linewidth', .5);
line([0 .04]+aX0, sqrt([2 2]), 'color', interp1(DLin,CP1,sqrt(2)),...
     'linewidth', .5);
line([0 .04]+aX0, [2 2], 'color', interp1(DLin,CP1,2),...
     'linewidth', .5);

% Draw distance colormap
scatter(D(1,:),ones(1,size(D,2))*aX0+.012,20,interp1(DLin,CP1,D(1,:)),...
        'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.012,D(2,:),20,interp1(DLin,CP1,D(2,:)),...
        'filled','marker','s');
 
% Draw map axes
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
    
% Draw map and y=x lines
plot([1.3 2.0],[1.3 2.0], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);

% Draw spline
sCo11 = [0.70 0.50 1.00;...
         0.00 0.50 1.00] .* [.09;-.17] + [1.87;2.03];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5);

% Draw example networks
plInd = [1 230 480 756];
plSh = [ .2  .2 .22 .3;...
        -.1 -.1 -.09 -.04];
lSh = -.01;
for i = 1:length(plInd)
    pI = plInd(i);
    XP = XC(:,:,pI)*sc + [D(1,pI);D(2,pI)]+plSh(:,i);
    line_coordinates(XP(:,1:2), 'lSh',lSh, 'lw',1, 'color',CPI(pI,:));
    line_coordinates(XP(:,3:4), 'lSh',-lSh, 'lw',1, 'color',CPJ(pI,:));
    visualize_network(XP,[],conn);
end
scatter(D(1,plInd),D(2,plInd),40,[interp1(DLin,CP1,D(1,plInd(1)));...
       zeros(2,3);interp1(DLin,CP1,D(1,plInd(4)))],...
       'filled','marker','s','linewidth',.5);

% Draw large template network
XP = Xs*sc*2 + [D(1,pI);D(2,pI)]+[.08;.57];
visualize_network(XP,[],conn,'lalpha',.3,'scolor',o*gr^.3);
line_coordinates(XP(:,1:2), 'lSh',-.035, 'lw',.5,'style','-','nw',.01);
line_coordinates(XP(:,3:4), 'lSh',0.035, 'lw',.5,'style','-','nw',.01);

% Text
texta = '\textbf{a}\hspace{.2cm}each unit $k$ maps $l_k$ to $l_{k+1}$';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(.011,.15,'$\sqrt{2}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(2)));
text(.011,.7,'$2$',NVTextRA{:},'color',interp1(DLin,CP1,2));
text(.83,0.02,'$l_k$',NVTextR{:});
text(-.087,.84,'$l_{k+1}$',NVTextR{:});
text(.23,.30,'$l_{k+1} = l_k$',NVTextR{:},'rotation',45,'color',o*gr);
text(.401,.8,'$l_{k+1} = f(l_k)$',NVTextR{:});
text(.14,0.71,'$l_k$',NVTextR{:});

% Axes
axis([1.24 2.3 1.24 2.3]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% b: single module with defined distances
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .2;                % Scale network
sh1 = [.60; .38];       % Shift network
sh2 = [1.75; .38];
lSh = .1; nW = .02;     % Distance bar

% Draw colored markers
scatter(-.12,.72,40,interp1(DLin,CP1,D(1,plInd(1))),...
        'filled','marker','s','linewidth',.5,'clipping',0);
scatter(1.02,.72,40,interp1(DLin,CP1,D(1,plInd(end))),...
        'filled','marker','s','linewidth',.5,'clipping',0);
    
% Draw module 1 and motion
pI1 = 300;
visualize_network(XC(:,:,pI1)*sc+sh1,[],conn,'lalpha',al,...
                  'scolor',o*sqrt(gr));
visualize_network(XC(:,:,1)*sc+sh1,[],conn);
line_coordinates(XC(:,1:2,pI1)*sc+sh1,'lSh',-lSh+.03,'nW',nW,...
                 'color',CPI(pI1,:).^.2);
line_coordinates(XC(:,1:2,1)*sc+sh1,'lSh',-lSh-.03,'nW',nW,...
                 'color', CPI(1,:));
line_coordinates(XC(:,3:4,pI1)*sc+sh1,'lSh',lSh+.12,'nW',nW,...
                 'color', CPJ(pI1,:).^.2);
line_coordinates(XC(:,3:4,1)*sc+sh1,'lSh',lSh+.175,'nW',nW,...
                 'color', CPJ(1,:));

% Draw module 2 and motion
pI2 = size(XC,3)-300;
visualize_network(XC(:,:,pI2)*sc+sh2,[],conn,'lalpha',al,...
                  'scolor',o*sqrt(gr));
visualize_network(XC(:,:,end)*sc+sh2,[],conn);
line_coordinates(XC(:,1:2,pI2)*sc+sh2,'lSh',-lSh+.02,'nw',nW,...
                 'color',CPI(pI2,:).^.2);
line_coordinates(XC(:,1:2,end)*sc+sh2,'lSh',-lSh-.02,'nw',nW,...
                 'color',CPI(end,:));
line_coordinates(XC(:,3:4,pI2)*sc+sh2,'lSh',lSh+.1,'nw',nW,...
                 'color',CPJ(pI2,:).^.2);
line_coordinates(XC(:,3:4,end)*sc+sh2,'lSh',lSh+.155,'nw',nW,...
                 'color',CPJ(end,:));

% Text
% textb = '\textbf{b}\hspace{.2cm}stability of unit geometry';
% text(labX,subp(pInd,4)+labY+.07,textb,NVTitle{:});
text(-0.075, 0.38,'$d l_k$',NVTextR{:});
text( 0.140, 0.75,'$d l_{k+1}$',NVTextR{:});
text( 0.540, 0.38,'$d l_k$',NVTextR{:});
text( 0.785, 0.75,'$d l_{k+1}$',NVTextR{:});
text( 0.138, 0.06,'unstable',NVTextH{:});
text( 0.138,-0.09,'$|d l_k| < |d l_{k+1}|$',NVTextH{:});
text( 0.755, 0.06,'stable',NVTextH{:});
text( 0.755,-0.09,'$|d l_k| > |d l_{k+1}|$',NVTextH{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c-e: module combination
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = 1/15;                                      % Scale network
sh = [[.1;.845] [.1;.47] [.1;.08]];             % Shift network: initial
shi = [.25;0];                                  % Shift network: indexed
scX = [linspace(0,shi(1),15); zeros(1,15)];     % Dots connecting nodes
lSh = .007;                                     % Distance bar offset
grL = [gr*ones(15,2) linspace(gr,1,15)'];       % Fading of dots

% Draw modules
pSh = [1 3600 size(XCa,3)];                     % Index of lattice to draw
sPr = {'$l_1$', '$l_2$', '$l_3$'};              % Text of length names
for i = 1:3                                     % Iterate across modules
    for j = 1:3                                 % Iterate across lattices
        XCp = XCa(:,(1:6)+2*(j-1),pSh(i))*sc+sh(:,i) + shi*(j-1);
        % Drwa dotted lines connecting module nodes
        scatter(scX(1,:)+XCp(1,3),scX(2,:)+XCp(2,3),.2,o.*grL(:,j));
        scatter(scX(1,:)+XCp(1,4),scX(2,:)+XCp(2,4),.2,o.*grL(:,j));
        % Draw module connection animation with transparent modules
        if(i == 1)
            visualize_network(XCp(:,3:6),[],conn,'lalpha',(1-gr)/2,...
                              'scolor',o*gr^.2);
            visualize_network(XCp(:,3:6)+[.04;0],[],conn,...
                              'lalpha',(1-gr),'scolor',o*gr^.3);
            if(j==1)
                line_coordinates(XCp(:,3:4), 'lSh',lSh, 'lw',1,...
                                 'color', interp1(DLin,CP1,Da(j,pSh(i))));
            end
        end
        % Draw length measurements
        line_coordinates(XCp(:,1:2), 'lSh',(-1)^j*lSh, 'lw',1,...
                         'color', interp1(DLin,CP1,Da(j,pSh(i))));
        % Draw networks
        visualize_network(XCp(:,1:4),[],conn);
        % Text
        text((XCp(1,2))/sRat(pInd)-.045-.005*(-1)^j,...
              XCp(2,2)-.03*(-1)^(j),sPr{j},NVTextH{:},...
              'color', interp1(DLin,CP1,Da(j,pSh(i))));
    end
end
text(.085,.91,'$f(l_1)$',NVTextH{:},'color',interp1(DLin,CP1,Da(2,1)));

% Draw lattice
sh2 = [-.04 -.02 0];                % Shift: Lattice
sh3 = [1.07 1.07 1.07;...
        0.09 0.09 0.09];
for i = 1:3                         % Iterate across lattices
    XCp = XCa(:,:,pSh(i))*sc+sh(:,i) + sh3(:,i);
    for j = 1:size(Da,1)            % Iterate across distances
        % Draw distance measurement
        line_coordinates(XCp(:,(-1:0)+2*j),'lSh',(-1)^j*lSh,'lw',1,...
                         'color',interp1(DLin,CP1,Da(j,pSh(i))));
        % x-axis text
        if(j < 3 || j == 9)
            text(XCp(1,-1+2*j)+.03,sh(2,i)-.11+sh2(i),num2str(j),...
                 'fontsize',FS,'horizontalalignment','center');
        elseif(j==4)
            text(XCp(1,-1+2*j)+.03,sh(2,i)-.11+sh2(i),'$\cdots$',...
                 'fontsize',FS,'horizontalalignment','center');
        end
    end
    % Draw distance squares on plot
    scatter(XCp(1,1:2:end)+.03,(sh(2,i)-.188+sh2(i))+Da(:,pSh(i))'/11,20,...
            interp1(DLin,CP1,Da(:,pSh(i))),'filled','marker','s','linewidth',2);
    % Draw network
    visualize_network(XCp,[],conna);
    % Draw tick marks
    axV = [0 .6 0 .14] + [1.16 1.16 sh(2,i)-.08+sh2(i) sh(2,i)-.08+sh2(i)];
    line([0 .01]+axV(1), [0 0]+sh(2,i)-.188+sh2(i)+2/11,...
         'color',interp1(DLin,CP1,2));
    line([0 .01]+axV(1), [0 0]+sh(2,i)-.188+sh2(i)+sqrt(2)/11,...
         'color',interp1(DLin,CP1,sqrt(2)));
    line([1;1].*(XCp(1,1:2:end)+.03),...
         [0;.01].*ones(1,size(XCp(1,1:2:end),2))+sh(2,i)-.08+sh2(i),...
         'color','k');
    % Draw plot axes
    line(axV([1 1 2]),axV([4 3 3]),'color','k','linewidth',.5);
    % Text
    text(1.15,sh(2,i)-.19+sh2(i)+2/11,'$2$','fontsize',FS,...
         'horizontalalignment','right','color',interp1(DLin,CP1,2));
    text(1.15,sh(2,i)-.19+sh2(i)+sqrt(2)/11,'$\sqrt{2}$','fontsize',FS,...
         'horizontalalignment','right','color',interp1(DLin,CP1,sqrt(2)));
    text(1.15,sh(2,i)-.19+sh2(i)+2.7/11,'$l_k$','fontsize',FS,...
         'horizontalalignment','right');
    text(1.78,sh(2,i)-.08+sh2(i),'$k$','fontsize',FS);
end

% Text
textc = ['\textbf{b}\hspace{.2cm}fixed point: if $l_1 = 2$, ',...
         'then $l_{k+1} = f(2) = 2$'];
textd = ['\textbf{c}\hspace{.2cm}unstable: if $l_1 < 2$, ',...
         'then $l_{k+1}$ moves away from 2'];
texte = ['\textbf{e}\hspace{.2cm}fixed point: if $l_1 = \sqrt{2}$, ',...
         'then $l_{k+1} = f\left(\sqrt{2}\right) = \sqrt{2}$'];
text(0,.974,textc,NVTextR{:});
text(0,.6,textd,NVTextR{:});
text(0,.23,texte,NVTextR{:});
text(.11,.71,'glue together',NVTextR{:},'color',o*gr^2);
text(1,.975,'periodic',NVTextRA{:});
text(1,.925,'lattice',NVTextRA{:});
text(1,.23,'periodic',NVTextRA{:});
text(1,.18,'lattice',NVTextRA{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'figure2fg';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 -2.67 fSize];
fig.PaperSize = [19 4.85];
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');