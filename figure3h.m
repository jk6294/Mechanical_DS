% Figure 3: Designing large finite conformatoinal geometries
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% figure dimensions
fig = figure(3); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 10];
% Margins in cm, [l,r,d,u]
fMarg = [.0 .4 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  5.00 10.00];...
        [ 5.30  5.00  5.00  5.00];...
        [ 5.30  0.00  5.25  5.25];...
        [10.50  6.50  4.75  3.50];...
        [10.50  3.15  4.75  3.50];...
        [10.50  0.00  4.75  3.00];...
        [15.25  0.00  3.75 10.00]];
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


%% figure parameters
FS = 10;            % Fontsize
gr = 0.9;           % Gray
dgr = 1;           % Dark gray

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
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));
% Gradient: Distance c. Interpolate lightness of one color
C2a = [133 37 83]/255;
C2b = C2a.^.1;
CP2 = interp1([0 1],[C2a;C2b],linspace(0,1,nT));


%% define module
% Parameters
sq = sqrt(3)/2;
L = 1;                                          % Initial length
LF = 1.7;                                       % Final length

% Module: increase c
Xs2 = [-sq   0.00   sq;...
       -0.5  1.00  -0.5]*L;
Xf2 = Xs2 * LF - [0;LF-L];
ds2 = sqrt(sum((Xf2(:,1)-Xf2(:,2)).^2));        % Final distance

% module: decrease c
Xs1 = [-sq 0  sq;...
      -.5 1 -.5]*L;
th = 15;
Xf1 = ds2*[[-sind(th);-cosd(th)] [0;0] [sind(th);-cosd(th)]] + [0;L];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

Xu2 = [-sq*L -sq*L ;...
       -sqrt((L*LF)^2 - (sq*L)^2) sqrt((L*LF)^2 - (sq*L)^2)];
Xu10 = [sq 0; 1.76 -1.1];
[Xu1,~] = construct_network(Xs1, Xf1, Xu10, conn, 0, 1);
Xu1 = Xu1(1:2,:);


%% simulate networks
% Run simulations
[XC2,fC2] = sim_motion10(Xs2,Xu2,conn,.01,181,[Xs2 Xu2],0);
D2 = squeeze(sqrt(sum(diff(XC2(:,[1 2 3 1],:),1,2).^2)));
[XC1,fC1] = sim_motion10(Xs1,Xu1,conn,.01,176,[Xs1 Xu1],0);
D1 = squeeze(sqrt(sum(diff(XC1(:,[1 2 3 1],:),1,2).^2)));

% Simulation errors
disp(['mean simulation error: '...
     [num2str(mean(fC1)) '  ' num2str(mean(fC2))]]);

% Correct tilt
for i = 1:size(XC1,3)
    Rz = rotz(atan2d(diff(XC1(2,[1 3],i)),diff(XC1(1,[1 3],i)))); 
    Rz = Rz(1:2,1:2);
    XC1(:,:,i) = Rz'*XC1(:,:,i);
end
for i = 1:size(XC2,3)
    Rz = rotz(atan2d(diff(XC2(2,[1 3],i)),diff(XC2(1,[1 3],i)))); 
    Rz = Rz(1:2,1:2);
    XC2(:,:,i) = Rz'*XC2(:,:,i);
end

% Number of solution spaces to sample
DLin = linspace(2*sq-.01,ds2+.01,nT);           % Distance range, d1,d2
DLin2 = linspace(1.52,3.97,nT);                 % Distance range, c
CPaI = interp1(DLin,CP1,D1(1,:));
CPaJ = interp1(DLin,CP1,D1(2,:));
CPaK = interp1(DLin2,CP2,D1(3,:));
CPbI = interp1(DLin,CP1,D2(1,:));
CPbJ = interp1(DLin,CP1,D2(2,:));
CPbK = interp1(DLin2,CP2,D2(3,:));


%% Tesselate quadrifolium
% Parameters
a = 16.2;                                   % Quadrifolium radius
f = @(tho)[a*sin(2*tho) .* cos(tho);...     % Quadrifolium equation
           a*sin(2*tho) .* sin(tho)];

% Trace quadrifolium
nv = 1000;
thv = linspace(0,2.1*pi,nv) + pi/4;         % Sampled angles
xv = f(thv);                                % Sampled points

% Initialize 
xP = zeros(2,113);                          % Triangle tesselation points
xP(:,1) = f(thv(1))*1.084543;               % Initial triangle point
thop = thv(1);

% Tesselate quadrifolium with isoscelese triangles
options = optimset('TolFun',1e-15,'TolX',1e-15);
for i = 1:size(xP,2)
    [thop, fV] = fzero(@(t) (sqrt(sum((f(t)-xP(:,i)).^2))-ds2/2),...
                             thop+[.01,.2],options);
    xP(:,i+1) = xP(:,i) + 2*(f(thop)-xP(:,i));
end

% Bin triangle distances
Dc = sqrt(sum((xP(:,1:end-2)-xP(:,3:end)).^2));     % Triangle distances: c
[C,~,IC] = unique(round(Dc,3),'stable');            % Bin to 3rd decimal
Xuc = zeros(2,2,length(C));                         % Unique node positions

% Construct units for each unique bin
Xu0a = [-3.5*sq 3.5*sq; 2.2 2.2];                   % Initial conditions
Xu0b = [-sq+.1 -sq+.1; -2 2.5];
Xu0c = [-sq*2 sq*2;  -3.3 -3.3];
for i = 1:length(C)
    % Different initial conditions based on distance c
    if(C(i)<2.2)
        Xu0 = Xu0a;
    elseif(C(i)>4)
        Xu0 = Xu0c;
    else
        Xu0 = Xu0b;
    end
    Xfp = [[-C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)] [0;0]...
           [C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)]] + [0;Xf2(2,2)];
    Xup = construct_network(Xs1,Xfp,Xu0,conn,0,1);
    Xuc(:,:,i) = Xup(1:2,:);
end

% Construct chain for quadrifolium network at the start position
[Xscc,Xucc,conncc,CSSc] = network_chain_x(Xs1(:,1:2),Xuc,IC',...
                                          interp1(DLin2,CP2,Dc));


%% a: design finite motion
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .05;                   % Scale drawing
sh1 = [.10 .10 .10;...      % Shift drawing
       .67 .39 .09];
sh2 = [.37 .37 .37;...
       .67 .39 .09];
R = [-1 1; -1 1.3]*1.7;     % Solution space ranges
lSh = .003;                 % Distance bar
mSs = 5;                    % Size: start node
mSf = 3;                    % Size: end node

% Define coordinates of rainbow node
thB = linspace(0,2*pi,360);
rB = .00176*mSf/subpN(pInd,4);
CPB = interp1(linspace(0,2*pi,nT),CP2,[thB fliplr(thB)]);

% Row 1: draw start and end positions, and distances d1,d2
line_coordinates(Xs1(:,[1,2])*sc+sh1(:,1),'lSh',lSh,'color',CPaI(1,:));
line_coordinates(Xs1(:,[2,3])*sc+sh1(:,1),'lSh',lSh,'color',CPaJ(1,:));
line_coordinates(Xf1(:,[1,2])*sc+sh1(:,1),'lSh',lSh,'color',CPaI(end,:));
line_coordinates(Xf1(:,[2,3])*sc+sh1(:,1),'lSh',lSh,'color',CPaJ(end,:));
visualize_network(Xs1*sc+sh1(:,1),[],[1 1],'msize',mSs);
visualize_network(Xf1*sc+sh1(:,1),[],[1 1],'msize',mSf,...
                  'scolor',o*dgr,'bcolor','k');
line_coordinates(Xs2(:,[1,2])*sc+sh2(:,1),'lSh',lSh,'color',CPaI(1,:));
line_coordinates(Xs2(:,[2,3])*sc+sh2(:,1),'lSh',lSh,'color',CPaJ(1,:));
line_coordinates(Xf2(:,[1,2])*sc+sh2(:,1),'lSh',-lSh,'color',CPaI(end,:));
line_coordinates(Xf2(:,[2,3])*sc+sh2(:,1),'lSh',-lSh,'color',CPaJ(end,:));
visualize_network(Xs2*sc+sh2(:,1),[],[1 1],'msize',mSs);
visualize_network(Xf2*sc+sh2(:,1),[],[1 1],'msize',mSf,...
                  'scolor',o*dgr,'bcolor','k');

% Row 2: draw solution spaces and distances c
visualize_conic_finite(Xs1*sc+sh1(:,2),Xf1*sc+sh1(:,2),R*sc+sh1(:,2),...
                       'ucolori',CPaK(end,:),'ucolorf',o,'overlay',.99);
line_coordinates(Xf1(:,[3,1])*sc+sh1(:,2),'lSh',lSh,'color',CPaK(end,:));
visualize_network(Xs1*sc+sh1(:,2),[],[1 1],'msize',mSs);
visualize_network(Xf1*sc+sh1(:,2),[],[1 1],'msize',mSf,...
                  'scolor',o*dgr,'bcolor',CPaK(end,:));
visualize_conic_finite(Xs2*sc+sh2(:,2),Xf2*sc+sh2(:,2),R*sc+sh2(:,2),...
                       'ucolori',CPbK(end,:),'ucolorf',o,'overlay',.99);
line_coordinates(Xf2(:,[3,1])*sc+sh2(:,2),'lSh',lSh,'color',CPbK(end,:));
visualize_network(Xs2*sc+sh2(:,2),[],[1 1],'msize',mSs);
visualize_network(Xf2*sc+sh2(:,2),[],[1 1],'msize',mSf,...
                  'scolor',o*dgr,'bcolor',CPbK(end,:));

% Row 3: draw constructed networks by adding extra nodes on solution space
visualize_conic_finite(Xs1*sc+sh1(:,3),Xf1*sc+sh1(:,3),R*sc+sh1(:,3),...
                       'ucolori',CPaK(end,:),'ucolorf',o,...
                       'overlay',.99,'lw',.5);
visualize_conic_finite(Xs2*sc+sh2(:,3),Xf2*sc+sh2(:,3),R*sc+sh2(:,3),...
                       'ucolori',CPbK(end,:),'ucolorf',o,...
                       'overlay',.99,'lw',.5);
XCP = XC1(:,:,1)*sc+sh1(:,3);
visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'ucolor',CPaK(end,:),...
                  'msize',5);
XCP = XC2(:,:,1)*sc+sh2(:,3);
visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'ucolor',CPbK(end,:),...
                  'msize',5);

% Legend
yPos = .92; nPo = 19;
visualize_network([.053;yPos+.005],[],[1 1],'msize',mSs);
visualize_network([.15;yPos+.005],[],[1 1],'scolor',o*dgr,'msize',mSf);
scatter(rB*cos(thB)+.17,rB*sin(thB)+yPos+.005,.25,CPB(1:2:end,:),...
        'filled','linewidth',.001)
scatter(linspace(.0,.17,nPo)+.018, ones(1,nPo)*yPos-.06, 1,...
        CP1(floor(linspace(1,size(CP1,1),nPo)),:),'filled', 'marker', 's');
scatter(linspace(.0,.17,nPo)+.275, ones(1,nPo)*yPos+.005, 1,...
        CP2(floor(linspace(1,size(CP2,1),nPo)),:),'filled', 'marker', 's');
scatter(linspace(.0,.17,size(CP2,1))+.275, ones(1,size(CP2,1))*yPos-.06,...
        1, CP2, 'filled', 'marker', 's');
text(.025,yPos+.035,'start',NVTextR{:});
text(.27,yPos+.035,'end',NVTextR{:});
text(.75,yPos+.035,'length $c$',NVTextH{:});
text(.025,yPos-.035,'length $l_1,l_2$',NVTextR{:});
text(.75,yPos-.035,'solution at $c$',NVTextH{:});

% Text
texta = '\textbf{a}\hspace{0.45cm}design unit geometry';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(0,sh1(2,1)+.13,'i.\hspace{1mm}fix start and end positions',NVTextR{:});
text(.22,sh1(2,1)+.09,'example 1',NVTextH{:});
text(.77,sh1(2,1)+.09,'example 2',NVTextH{:});
text(0,sh1(2,2)+.13,'ii. ~solve for solution space',NVTextR{:});
text(0,sh1(2,3)+.13,'iii. ~add nodes on solution',NVTextR{:});
% Node coordinates
text(Xs1(1,1)*sc+sh1(1,1)-.013,Xs1(2,1)*sc+sh1(2,1),'$1$',...
     'fontsize',FS*4/4,'horizontalalignment','right','color',o*gr^4);
text(Xs1(1,2)*sc+sh1(1,1)+.015,Xs1(2,2)*sc+sh1(2,1),'$2$',...
     'fontsize',FS*4/4,'horizontalalignment','left','color',o*gr^4);
text(Xs1(1,3)*sc+sh1(1,1)+.015,Xs1(2,3)*sc+sh1(2,1),'$3$',...
     'fontsize',FS*4/4,'horizontalalignment','left','color',o*gr^4);
% Edge lengths
text(.10,sh1(2,1)+.025,'$l_1$',NVTextR{:});
text(.32,sh1(2,1)+.025,'$l_2$',NVTextR{:});
text(.66,sh1(2,1)+.02,'$l^\bullet$',NVTextR{:},'color',CPaI(1,:));
text(.85,sh1(2,1)+.02,'$l^\bullet$',NVTextR{:},'color',CPaJ(1,:));
text(.68,sh1(2,1)-.065,'$l^\circ$',NVTextR{:},'color',CPaI(end,:));
text(.83,sh1(2,1)-.065,'$l^\circ$',NVTextR{:},'color',CPaJ(end,:));
text(.21,.27,'$c^\circ$',NVTextR{:},'color',CPaK(end,:));
text(.75,.27,'$c^\circ$',NVTextR{:},'color',CPbK(end,:));

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: conformational motion
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Draw map and line y=x
aX0 = 1.2; aXF = 3.25;                  % Axes min and max locations
plot(D1(1,:),D1(2,:),'k-','linewidth',1);
plot(D2(1,:),D2(2,:),'k-','linewidth',1);
plot([aX0 aXF],[aX0 aXF], '-', 'color', o*gr, 'linewidth',.5);

% Draw axes and ticks
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line([1 1]*D1(1,1), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([1 1]*D1(1,end), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [1 1]*D1(1,1), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [1 1]*D1(1,end), 'color', 'k', 'linewidth', .5);

% Draw examples
nE = 4;                 % Number of examples
sc = .13;               % Scale
pIa = floor(linspace(1,size(XC1,3),nE));
sha = [linspace(-.25,-.2,nE);...
       linspace(.15,.25,nE)];
pIb = floor(linspace(1,size(XC2,3),nE));
shb = [linspace(.23,.27,nE);...
       linspace(-.26,-.1,nE)];
for i = 1:nE
    % unit a
    pa = pIa(i);
    scatter(D1(1,pa),D1(2,pa),40,'k','filled','linewidth',.5,'marker','s');
    XCP = XC1(:,:,pa)*sc+D1(1:2,pa)+sha(:,i);
    if(i == 4)
        line_coordinates(XCP(:,[3 1]),'lSh',.02,'color',CPaK(pa,:));
    end
    visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'ucolor',CPaK(end,:));
    % unit b
    pb = pIb(i);
    scatter(D2(1,pb),D2(2,pb),40,'k','filled','linewidth',.5,'marker','s');
    XCP = XC2(:,:,pb)*sc+D2(1:2,pb)+shb(:,i);
    if(i == 4)
        line_coordinates(XCP(:,[3 1]),'lSh',.02,'color',CPbK(pb,:));
    end
    visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'ucolor',CPbK(end,:));
end

% Text
textb = '\textbf{b}\hspace{0.3cm}conformational motion';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
text(.89,0.045,'$l_1$',NVTextR{:});
text(.03,.91,'$l_2$',NVTextR{:});
text(.08,.15,'start',NVTextR{:});
text(.78,.87,'end',NVTextR{:});
text(D1(1,1),aX0-(aXF-aX0)*.06,num2str(D1(1,1),3),...
     'Horizontalalignment','center','fontsize',FS);
text(D1(1,end),aX0-(aXF-aX0)*.07,num2str(D1(1,end),3),...
     'Horizontalalignment','center','fontsize',FS);

% Axes
axis([1.1 3.6 1.1 3.6]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: Range of solution spaces
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
nSS = 15;                       % Number of solutoin spaces to sample
sc = .13;                       % Scale networks
sh = [0.45;0.42];               % Shift networks
thR = linspace(15,42.2683,nSS); % Angles of solution spaces to draw
R2 = [-1 1; -1 1]*3.2;          % Range of solution space to draw

% Draw solution space
for i = 1:nSS
    Xfp = ds2*[[-sind(thR(i));-cosd(thR(i))] [0;0],...
               [sind(thR(i));-cosd(thR(i))]] + [0;Xf2(2,2)];
    visualize_conic_finite(Xs1*sc+sh(:,1),Xfp*sc+sh(:,1),R2*sc+sh(:,1),...
                           'ucolori',CP2(ceil(i/nSS*size(CP2,1)),:),...
                           'ucolorf',o,'overlay',.99);
end

% Draw start and end node positions
for i = 1:nSS
    Xfp = ds2*[[-sind(thR(i));-cosd(thR(i))] [0;0],...
               [sind(thR(i));-cosd(thR(i))]] + [0;Xf2(2,2)];
    visualize_network(Xs1*sc+sh(:,1),[],[1 1],'msize',mSs);
    visualize_network(Xfp*sc+sh(:,1),[],[1 1],'scolor',o*dgr,'msize',mSf,...
                      'bcolor',interp1(DLin2,CP2,diff(Xfp(1,[1 3]))));
end

% Draw rainbow node
rB = .00176*mSf/subpN(pInd,4);
scatter(rB*cos(thB)+0*sc+sh(1,1),rB*sin(thB)+1*sc+sh(2,1),...
        .25,CPB(1:2:end,:),'filled','linewidth',.001)

% Text
textc = '\textbf{c}\hspace{.35cm}range of solution spaces';
text(labX,subp(pInd,4)-.75,textc,NVTitle{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% d: quadrifolium trace
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .03;               % Scale
shl = [.7; .55];        % Shift

% Draw quadrifolium curve
plot(xv(1,:)*sc+shl(1),xv(2,:)*sc+shl(2),'k-','linewidth',.5);

% Text
textd = '\textbf{d}\hspace{.54cm}trace desired curve';
text(labX,subp(pInd,4)+labY,textd,NVTitle{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% d: quadrifolium tesselation
pInd = 5;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .03;               % Scale trace
sc2 = .085;             % Scale example networks
sc3 = .4;               % Scale unit triangle symbol
sh = [.7; .48];         % Shift Networks
sh1a = [.14;.54];
sh1b = [.14;.12];
sh1c = [.07;.65];
sh2a = [1.295;.54];
sh2b = [1.295;.12];
sh2c = [.765;.65];

% Draw quadrifolium trace
plot(xv(1,:)*sc+sh(1),xv(2,:)*sc+sh(2),'-','linewidth',.5,'color',o*gr);
% Draw all tesselated isoscelese triangles
for i = 1:size(xP,2)-2
    line_coordinates(xP(:,[0 2]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin2,CP2,Dc(i)),'style','-');
    line_coordinates(xP(:,[0 1]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin,CP1,ds2),'style','-');
end
% Fill in a subset of triangles
for i = 3:16
    fill(xP(1,[0:2 0]+i)*sc+sh(1),xP(2,[0:2 0]+i)*sc+sh(2),...
         interp1(DLin2,CP2,Dc(i)),'linewidth',.5);
end
for i = 3:17
    line_coordinates(xP(:,[0 2]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin2,CP2,Dc(i)),'style','-');
    line_coordinates(xP(:,[0 1]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin,CP1,ds2),'style','-');
end

% Simulate Examples
[Xsm1,Xum1,connm1,CSScm1] = network_chain_x(Xs1(:,1:2),Xuc,IC(3),...
                                            interp1(DLin2,CP2,Dc));
[Xsm2,Xum2,connm2,CSScm2] = network_chain_x(Xs1(:,1:2),Xuc,IC(4),...
                                            interp1(DLin2,CP2,Dc));
[XCm1,fCm1] = sim_motion10(Xsm1,Xum1,connm1,.01,199,[Xsm1,Xum1],0);
[XCm2,fCm2] = sim_motion10(Xsm2,Xum2,connm2,.01,178,[Xsm2,Xum2],0);

% Simulation errors
disp(['mean simulation error: '...
     [num2str(mean(fCm1)) '  ' num2str(mean(fCm2))]]);
 
% Distance calculation
DM1 = sqrt(squeeze(sum(diff(XCm1(:,[1 2 3 1],:),1,2).^2)));
DM2 = sqrt(squeeze(sum(diff(XCm2(:,[1 2 3 1],:),1,2).^2)));

% Draw unit 1
Rz1 = rotz(3); Rz1 = Rz1(1:2,1:2);
Xpl1a = Rz1*XCm1(:,:,end)*sc2+sh1a; 
Xpl1b = XCm1(:,:,1)*sc2+sh1b;

% Draw shape
fill(Xpl1a(1,[1 2 3 1])*sc3+sh1c(1),Xpl1a(2,[1 2 3 1])*sc3+sh1c(2),...
     interp1(DLin2,CP2,DM1(3,end)),'linewidth',.01,'edgecolor','w');
line_coordinates(Xpl1a(:,[1 2])*sc3+sh1c,'color',...
                 interp1(DLin,CP1,DM1(1,end)),'lSh',0,'lw',.5,'style','-');
line_coordinates(Xpl1a(:,[2 3])*sc3+sh1c,'color',...
                 interp1(DLin,CP1,DM1(2,end)),'lSh',0,'lw',.5,'style','-');
line_coordinates(Xpl1a(:,[1 3])*sc3+sh1c,'color',...
                 interp1(DLin2,CP2,DM1(3,end)),'lSh',0,'lw',.5,'style','-');

% Draw network 1
fill(Xpl1a(1,[1 2 3 1]),Xpl1a(2,[1 2 3 1]),...
     interp1(DLin2,CP2,DM1(3,end)),'linewidth',.01,'edgecolor','w');
line_coordinates(Xpl1a(:,[1 2]),'color',interp1(DLin,CP1,DM1(1,end)),'lSh',0);
line_coordinates(Xpl1a(:,[2 3]),'color',interp1(DLin,CP1,DM1(2,end)),'lSh',0);
line_coordinates(Xpl1a(:,[1 3]),'color',interp1(DLin2,CP2,DM1(3,end)),'lSh',0);
visualize_network(Xpl1a(:,1:3),Rz1*Xpl1a(:,4:5),connm1,...
                  'ucolor',interp1(DLin2,CP2,DM1(3,end)));
visualize_network(Xpl1b(:,1:3),Xpl1b(:,4:5),connm1,...
                  'ucolor',interp1(DLin2,CP2,DM1(3,end)));

% Draw unit 2
Rz2 = rotz(2); Rz2 = Rz2(1:2,1:2);
Xpl2a = Rz2*XCm2(:,:,end)*sc2+sh2a; 
Xpl2b = XCm2(:,:,1)*sc2+sh2b;

% Draw shape 2
fill(Xpl2a(1,[1 2 3 1])*sc3+sh2c(1),Xpl2a(2,[1 2 3 1])*sc3+sh2c(2),...
     interp1(DLin2,CP2,DM2(3,end)),'linewidth',.01,'edgecolor','w');
line_coordinates(Xpl2a(:,[1 2])*sc3+sh2c,'color',...
                 interp1(DLin,CP1,DM2(1,end)),'lSh',0,'lw',.5,'style','-');
line_coordinates(Xpl2a(:,[2 3])*sc3+sh2c,'color',...
                 interp1(DLin,CP1,DM2(2,end)),'lSh',0,'lw',.5,'style','-');
line_coordinates(Xpl2a(:,[1 3])*sc3+sh2c,'color',... 
                 interp1(DLin2,CP2,DM2(3,end)),'lSh',0,'lw',.5,'style','-');

% Draw network 2
fill(Xpl2a(1,[1 2 3 1]),Xpl2a(2,[1 2 3 1]),...
     interp1(DLin2,CP2,DM2(3,end)),'linewidth',.01,'edgecolor','w');
line_coordinates(Xpl2a(:,[1 2]),'color',interp1(DLin,CP1,DM2(1,end)),'lSh',0);
line_coordinates(Xpl2a(:,[2 3]),'color',interp1(DLin,CP1,DM2(2,end)),'lSh',0);
line_coordinates(Xpl2a(:,[1 3]),'color',interp1(DLin2,CP2,DM2(3,end)),'lSh',0);
visualize_network(Xpl2a(:,1:3),Xpl2a(:,4:5),connm1,...
                  'ucolor',interp1(DLin2,CP2,DM2(3,end)));
visualize_network(Xpl2b(:,1:3),Xpl2b(:,4:5),connm1,...
                  'ucolor',interp1(DLin2,CP2,DM2(3,end)));

% Text
texte2 = 'isoceles triangle tessellation';
text(labX,subp(pInd,4)+labY+.3,texte2,NVTitle{:});
text(.08,.98,'unit',NVTextH{:});
text(.08,.755,'end',NVTextH{:});
text(.08,.32,'start',NVTextH{:});
text(.92,.98,'unit',NVTextH{:});
text(.92,.755,'end',NVTextH{:});
text(.92,.32,'start',NVTextH{:});

% Axis
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% d: combine individual units
pInd = 6;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .1;                    % Scale networks
sh1 = [0.49;0.61];          % Shift networks: offset
shI = [.4;-.0];             % Shift networks: iteratively
nP = 10;                    % Number of points in dots connecting nodes
nSh = 4.2;                  % Network shift amount

% Draw example units
XsP = Xscc(:,:,1)*sc+sh1;
XuP = Xucc(:,:,1)*sc+sh1;
for i = 3:5
    if(i<5)
        scatter(linspace(0,shI(1),nP)+XsP(1,i+1)+shI(1)*(i-nSh),...
                ones(1,nP)*XsP(2,i+1),.3,o*.8);
        scatter(linspace(0,shI(1),nP)+XsP(1,i+2)+shI(1)*(i-nSh),...
                ones(1,nP)*XsP(2,i+2),.3,o*.8);
    else
        scatter(linspace(0,shI(1),nP)+XsP(1,i+1)+shI(1)*(i-nSh),...
                ones(1,nP)*XsP(2,i+1),.3,o.*(linspace(.8,1,nP)'));
        scatter(linspace(0,shI(1),nP)+XsP(1,i+2)+shI(1)*(i-nSh),...
                ones(1,nP)*XsP(2,i+2),.3,o.*(linspace(.8,1,nP)'));
    end
    visualize_network(XsP(:,(0:2)+i)+shI*(i-nSh),...
                      XuP(:,(-1:0)+2*i)+shI*(i-nSh),conn,...
                      'ucolor',interp1(DLin2,CP2,Dc(i)));
end

% Draw example lattice
for i = 3:16
    visualize_network(XsP(:,(0:2)+i)+shI*(3-nSh)+[0;-.45],...
                      XuP(:,(-1:0)+2*i)+shI*(3-nSh)+[0;-.45],conn,...
                      'ucolor',interp1(DLin2,CP2,Dc(i)));
end

% Draw spline arrows
sCo11 = [0.60 0.00 0.00 1.00;...
         1.00 0.70 0.30 0.00] .* [.1;.58] + [-.03;.65];
sCo12 = [0.87 1.00 0.20 0.00;...
         1.00 0.20 0.10 0.00] .* [1;.42] + [.7;.8];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd));
plot_spline(sCo12,'head',1,'headpos',1,'ratio',sRat(pInd));

% Text
texte3 = 'combine units in tessellation';
text(labX,subp(pInd,4)+labY-.15,texte3,NVTitle{:});
text(.82,.2,'$\dots$', NVTextR{:});

% Axis
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0);
drawnow;


%% e: Conformational motion
pInd = 7;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .0098;                 % Scale
sh1 = [0.345;0.015];        % Shift transition example networks
sh2 = [0.28; 0.24];         % Shift final configuration network

% Simulation data
load quadrifolium.mat;
% Uncomment the below 4 lines to run the simulation
% XCc0 = zeros(2,size([Xscc Xucc],2));
% XCc0(1,1) = -1;
% [XCcc,fCcc] = sim_motion10(Xscc,Xucc,conncc,1.5,1157,XCc0,0);
% disp(['mean simulation error: ' num2str(mean(fCcc))]);

% Correct rotation and offset in simulations
XCcc = XCcc - XCcc(:,size(Xscc,2),:);
% XCcc(1,:,:) = -XCcc(1,:,:);
for i = 1:size(XCcc,3)
    Rz = rotz(atan2d(diff(XCcc(2,[-2 0]+size(Xscc,2),i)),...
              diff(XCcc(1,[-2 0]+size(Xscc,2),i))));
    Rz2 = rotz(270);
    XCcc(:,:,i) = Rz2(1:2,1:2)*Rz(1:2,1:2)'*XCcc(:,:,i);
end

%% Draw example networks through the progression along conformational change
pI = [1 100 520];
for i = 1:length(pI)
    la = .1 + .15*i;
    visualize_network(XCcc(:,unique(conncc(:,1)),pI(i))*sc+sh1,...
              XCcc(:,unique(conncc(:,2)),pI(i))*sc+sh1,...
              conncc + [0 max(conncc(:,1))-size(Xscc,2)],...
              'lalpha',la,'msize',2,'ucolor',...
              CSSc(1:length(unique(conncc(:,2))),:).^(.15*i));
end

% Draw final configuration
Rz = rotz(45); Rz = Rz(1:2,1:2);
visualize_network(Rz*XCcc(:,1:size(Xscc,2),end)*sc+sh2,...
                  Rz*XCcc(:,(1:size(Xucc,2))+size(Xscc,2),end)*sc+sh2,...
                  conncc,'msize',2.5,'ucolor',CSSc);

% Text
texte = '\textbf{e}\hspace{1.5mm}conformational motion';
text(labX-.1,subp(pInd,4)+labY,texte,NVTitle{:});
text(.8,.95,'start',NVTextH{:},'color',o*gr^2);
text(.65,.71,'middle',NVTextH{:},'color',o*gr^4);
text(.5,.31,'end',NVTextH{:});
text(.29,.79,'$l_1$',NVTextH{:},'color',o*gr^4);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);


%% Save
% fName = 'figure3h';
% set(gcf, 'Renderer', 'painters');
% fig.PaperPositionMode = 'manual';
% fig.PaperUnits = 'centimeters';
% fig.PaperPosition = [0 0 fSize];
% fig.PaperSize = fSize;
% saveas(fig, ['Figures/' fName], 'pdf');
% set(gcf, 'Renderer', 'opengl');
    
