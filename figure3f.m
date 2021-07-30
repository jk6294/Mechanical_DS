% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(3); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 9.5];
% Margins in cm, [l,r,d,u]
fMarg = [.0 .4 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 0.00 4.50 9.50];...
        [ 4.50 5.00 5.00 4.50];...
        [ 4.50 0.00 5.00 5.00];...
        [ 9.50 4.75 4.75 4.75];...
        [ 9.50 0.00 4.75 4.75];...
        [14.25 6.00 4.75 3.50];...
        [14.25 0.00 4.75 6.00]];
% Fontsize
FS = 10;
% Distance visualization parameters
lSh = .3;
nW = .08;
lw = .5;
gr = 0.8;
% Scaling parameters
scc = .25;
    
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.17;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');

% Default Name-Value Pairs
NVTitle = {'Units','centimeters','fontsize',FS};
NVTextH = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','center'};
NVTextRA = {'Units','Normalized','fontsize',FS,'HorizontalAlignment','right'};
NVTextR = {'Units','Normalized','fontsize',FS};


%% Define module
sh = sqrt(3)/2;
L = 1;
LF = 1.7;

% module: increase c
Xs2 = [-sh 0  sh;...
      -.5 1 -.5]*L;
Xf2 = Xs2 * LF;
ds2 = sqrt(sum((Xf2(:,1)-Xf2(:,2)).^2));

% module: decrease c
Xs1 = [-sh 0  sh;...
      -.5 1 -.5]*L;
th = 15;
Xf1 = ds2*[[-sind(th);-cosd(th)] [0;0] [sind(th);-cosd(th)]] + [0;Xf2(2,2)];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

Xu2 = [-sh*L -sh*L ;...
       -sqrt((L*LF)^2 - (sh*L)^2) sqrt((L*LF)^2 - (sh*L)^2)];
Xu10 = [sh 0; 1.76 -1.1];
[Xu1,~] = construct_network(Xs1, Xf1, Xu10, conn, 0, 1);
Xu1 = Xu1(1:2,:);

% Simulation
[XC2,fC2] = sim_motion(Xs2,Xu2,conn,.001,1808,[Xs2 Xu2],0);
D2 = squeeze(sqrt(sum(diff(XC2(:,[1 2 3 1],:),1,2).^2)));
[XC1,fC1] = sim_motion(Xs1,Xu1,conn,.001,1758,[Xs1 Xu1],0);
D1 = squeeze(sqrt(sum(diff(XC1(:,[1 2 3 1],:),1,2).^2)));

% correct tilt
for i = 1:size(XC1,3)
    Rz = rotz(atan2d(diff(XC1(2,[1 3],i)),diff(XC1(1,[1 3],i)))); 
    Rz = Rz(1:2,1:2);
    XC1(:,:,i) = Rz'*XC1(:,:,i);
%     XC1(:,:,i) = XC1(:,:,i) - mean(XC1(:,:,i),2);
end
for i = 1:size(XC2,3)
    Rz = rotz(atan2d(diff(XC2(2,[1 3],i)),diff(XC2(1,[1 3],i)))); 
    Rz = Rz(1:2,1:2);
    XC2(:,:,i) = Rz'*XC2(:,:,i);
%     XC2(:,:,i) = XC2(:,:,i) - mean(XC2(:,1:3,i),2);
end

% color
nSS = 15;

% Coefficients
cR = .299;
cG = .587;
cB = .114;
br = 1;

% Gradient along red/blue, sweep red
DLin = linspace(min([D1(:);D2(:)])-.01,max([D1(:);D2(:)])+.01,size(D1,2));
RLin = linspace(0,1,size(D1,2));
GLin = .1 * ones(1,size(D1,2));
BLin = sqrt(br^2 - cG*GLin.^2 - cR*RLin.^2);
CP = [RLin; GLin; BLin]';

% Gradient along blue/green, sweep blue
br = .8;
BLin = linspace(0,1,nSS).^(.5);
RLin = 0 * ones(1,nSS);
GLin = sqrt(br^2 - cB*BLin.^2 - cR*RLin.^2);
CS = flipud([RLin; GLin; BLin]');

% CP = flipud(winter(size(XC,3)));
CPaI = interp1(DLin,CP,D1(1,:));
CPaJ = interp1(DLin,CP,D1(2,:));
CPaK = interp1(DLin,CP,D1(3,:));
CPbI = interp1(DLin,CP,D2(1,:));
CPbJ = interp1(DLin,CP,D2(2,:));
CPbK = interp1(DLin,CP,D2(3,:));


%% Combined network
[Xsa,Xua,conna,CSS] = network_chain_x(Xs1(:,1:2),cat(3,Xu2,Xu1),[1 2 1 2 1 2 1 2],CS([15,1],:));
% Simulate
[XCa,fCa] = sim_motion(Xsa,Xua,conna,.005,2096,[Xsa Xua],0);
% Distances
Dc = squeeze(sqrt(sum(diff(XCa(:,[1:2:size(Xsa,2)],:),1,2).^2)));
Dc = [Dc; squeeze(sqrt(sum(diff(XCa(:,[2:2:size(Xsa,2)],:),1,2).^2)))];
Da = squeeze(sqrt(sum(diff(XCa(:,1:size(Xsa,2),:),1,2).^2)));
[~,I] = min(sum(abs(Da-ds2)));


%% Animate
% figure(4); clf; 
% animate_network(XCa,conna,'ucolor',CSS,'nu',size(Xua,2));


%% a: design finite motion
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Plot parameters
sc = .05; 
sh1 = [.10 .10 .10;...
       .62 .35 .09];
sh2 = [.35 .35 .35;...
       .62 .35 .09];
C_I = [1 1 1]*.92;
C_F = [1 1 1]*.6;
C_SS = [170 150 255]/255;
R = [-1 1; -1 1.3]*1.7;
lSh = .004;

% Title
texta = '\textbf{a}\hspace{0.2cm}design unit geometry';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});

% Text
text(.0,.75,'fix start and end positions',NVTextR{:});
text(.0,.47,'distances between nodes',NVTextR{:});
text(.0,.22,'solve for solution space',NVTextR{:});

% Row 1
visualize_network(Xs1*sc+sh1(:,1),[],[1 1]);
visualize_network(Xf1*sc+sh1(:,1),[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);
visualize_network(Xs2*sc+sh2(:,1),[],[1 1]);
visualize_network(Xf2*sc+sh2(:,1),[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);

% Row 2
line_coordinates(Xs1(:,[1,2])*sc+sh1(:,2),'lSh',lSh,'color',CPaI(1,:));
line_coordinates(Xs1(:,[2,3])*sc+sh1(:,2),'lSh',lSh,'color',CPaJ(1,:));
line_coordinates(Xs1(:,[3,1])*sc+sh1(:,2),'lSh',lSh,'color',CPaK(1,:));
line_coordinates(Xf1(:,[1,2])*sc+sh1(:,2),'lSh',lSh,'color',CPaI(end,:));
line_coordinates(Xf1(:,[2,3])*sc+sh1(:,2),'lSh',lSh,'color',CPaJ(end,:));
line_coordinates(Xf1(:,[3,1])*sc+sh1(:,2),'lSh',lSh,'color',CPaK(end,:));
visualize_network(Xs1*sc+sh1(:,2),[],[1 1]);
visualize_network(Xf1*sc+sh1(:,2),[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);
line_coordinates(Xs2(:,[1,2])*sc+sh2(:,2),'lSh',lSh,'color',CPbI(1,:));
line_coordinates(Xs2(:,[2,3])*sc+sh2(:,2),'lSh',lSh,'color',CPbJ(1,:));
line_coordinates(Xs2(:,[3,1])*sc+sh2(:,2),'lSh',lSh,'color',CPbK(1,:));
line_coordinates(Xf2(:,[1,2])*sc+sh2(:,2),'lSh',lSh,'color',CPbI(end,:));
line_coordinates(Xf2(:,[2,3])*sc+sh2(:,2),'lSh',lSh,'color',CPbJ(end,:));
line_coordinates(Xf2(:,[3,1])*sc+sh2(:,2),'lSh',lSh,'color',CPbK(end,:));
visualize_network(Xs2*sc+sh2(:,2),[],[1 1]);
visualize_network(Xf2*sc+sh2(:,2),[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);

% Row 3
visualize_conic_finite(Xs1*sc+sh1(:,3),Xf1*sc+sh1(:,3),R*sc+sh1(:,3),...
                       'ucolori',CS(1,:),'ucolorf',[1 1 1],'overlay',.99);
visualize_network(Xs1*sc+sh1(:,3),[],[1 1]);
visualize_network(Xf1*sc+sh1(:,3),[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);

visualize_conic_finite(Xs2*sc+sh2(:,3),Xf2*sc+sh2(:,3),R*sc+sh2(:,3),...
                       'ucolori',CS(15,:),'ucolorf',[1 1 1],'overlay',.99);
visualize_network(Xs2*sc+sh2(:,3),[],[1 1]);
visualize_network(Xf2*sc+sh2(:,3),[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);

% Legend
yPos = .83;
visualize_network([.05;yPos+.12],[],[1 1],'scolor',C_I);
visualize_network([.05;yPos+.08],[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);
scatter(linspace(0,.1,size(CP,1)), ones(1,size(CP,1))*yPos+.04, 1, interp1(1:nSS,CS,linspace(1,nSS,size(CP,1))), 'marker', 's');
scatter(linspace(0,.1,size(CP,1)), ones(1,size(CP,1))*yPos, 1, CP, 'marker', 's');
text(.33,yPos+.12,'start position',NVTextR{:});
text(.33,yPos+.08,'end position',NVTextR{:});
text(.33,yPos+.04,'solution space',NVTextR{:});
text(.33,yPos,'distance',NVTextR{:});

% Text
% text(.3,.77,'$d_1^0$',NVTextH{:},'color',CPaI(1,:));
% text(.7,.77,'$d_2^0$',NVTextH{:},'color',CPaJ(1,:));
% text(.39,.86,'$d_1^*$',NVTextH{:},'color',CPaI(end,:));
% text(.61,.86,'$d_2^*$',NVTextH{:},'color',CPaJ(end,:));
% text(.5,.688,'$c^0$',NVTextH{:},'color',CPaK(1,:));
% text(.5,.64,'$c^*$',NVTextH{:},'color',CPaK(end,:));

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: Units
pInd = 2;
subplot('position',subpN(pInd,:)); cla;

% Title
textb = '\textbf{b}\hspace{0.2cm}construct units';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Text
text(.0,.89,'add nodes on solution space',NVTextR{:});
text(.0,.41,'connect nodes with edges',NVTextR{:});

% plot
sc = .09;
sh1 = [0.25 0.25; .65 .18];
sh2 = [0.80 0.80; .65 .18];

% solutoin spaces
visualize_conic_finite(Xs1*sc+sh1(:,1),Xf1*sc+sh1(:,1),R*sc+sh1(:,1),...
                       'ucolori',CS(1,:),'ucolorf',[1 1 1],'overlay',.2);
visualize_conic_finite(Xs2*sc+sh2(:,1),Xf2*sc+sh2(:,1),R*sc+sh2(:,1),...
                       'ucolori',CS(end,:),'ucolorf',[1 1 1],'overlay',.2);
% Overlay networks
XCP = XC1(:,:,1)*sc+sh1(:,1);
visualize_network(XCP(:,1:3),XCP(:,4:5),[1 1],'scolor',C_I,'ucolor',CS(1,:));
XCP = XC2(:,:,1)*sc+sh2(:,1);
visualize_network(XCP(:,1:3),XCP(:,4:5),[1 1],'scolor',C_I,'ucolor',CS(end,:));
% Constructed networks
XCP = XC1(:,:,1)*sc+sh1(:,2);
visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'scolor',C_I,'ucolor',CS(1,:));
XCP = XC2(:,:,1)*sc+sh2(:,2);
visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'scolor',C_I,'ucolor',CS(end,:));

% Axis
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0);


%% conformational motion
% Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
hold on;

% Title
text(labX,subp(pInd,4)+labY,'unit conformational motion',NVTitle{:});

% Map
CPa = winter(size(XC1,3));
CPb = winter(size(XC2,3));
hold on;
plot(D1(1,:),D1(2,:),'k-','linewidth',1);
plot(D2(1,:),D2(2,:),'k-','linewidth',1);
% scatter(Da(1,:),Da(2,:),1,CPa,'o','linewidth',.1);
% scatter(Db(1,:),Db(2,:),1,CPb,'o','linewidth',.1);

% Axes
aX0 = 1.2; aXF = 3.3;
plot([aX0 aXF],[aX0 aXF], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line([1 1]*D1(1,1), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([1 1]*D1(1,end), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [1 1]*D1(1,1), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [1 1]*D1(1,end), 'color', 'k', 'linewidth', .5);
text(.91,0.045,'$d_1$',NVTextR{:});
text(.03,.93,'$d_2$',NVTextR{:});
% text(.56,.44,'$d_2 = f_a(d_1)$',NVTextH{:},'rotation',45);
% text(.45,.55,'$d_2 = f_b(d_1)$',NVTextH{:},'rotation',45);
text(D1(1,1),aX0-(aXF-aX0)*.05,num2str(D1(1,1),3),'Horizontalalignment','center','fontsize',FS);
text(D1(1,end),aX0-(aXF-aX0)*.05,num2str(D1(1,end),3),'Horizontalalignment','center','fontsize',FS);
axis([1.1 3.6 1.1 3.6]);

% Examples
nE = 4;
plInda = floor(linspace(1,size(XC1,3),nE));
plSha = [linspace(-.25,-.2,nE);...
         linspace(.15,.25,nE)];
plIndb = floor(linspace(1,size(XC2,3),nE));
plShb = [linspace(.23,.27,nE);...
         linspace(-.26,-.1,nE)];
sc3 = .13;
for i = 1:nE
    % unit a
    pIa = plInda(i);
    plot(D1(1,pIa),D1(2,pIa),'ks','linewidth',3,'markersize',3);
    XCP = XC1(:,:,pIa)*sc3+D1(1:2,pIa)+plSha(:,i);
    if(i==1 || i == 4)
        line_coordinates(XCP(:,1:2),'lSh',.02, 'nW',.0, 'lw',1, 'color',CPaI(pIa,:));
        line_coordinates(XCP(:,2:3),'lSh',.02, 'nW',.0, 'lw',1, 'color',CPaJ(pIa,:));
    end
    visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'scolor',C_I,'ucolor',CS(1,:));
    
    % unit b
    pIb = plIndb(i);
    plot(D2(1,pIb),D2(2,pIb),'ks','linewidth',3,'markersize',3);
    XCP = XC2(:,:,pIb)*sc3+D2(1:2,pIb)+plShb(:,i);
    if(i==1 || i == 4)
        line_coordinates(XCP(:,1:2),'lSh',.02, 'nW',.0, 'lw',1, 'color',CPbI(pIb,:));
        line_coordinates(XCP(:,2:3),'lSh',.02, 'nW',.0, 'lw',1, 'color',CPbJ(pIb,:));
    end
    visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'scolor',C_I,'ucolor',CS(end,:));
end

% Axis limits
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: Combine to form curvature
pInd = 4;
subplot('position',subpN(pInd,:)); cla;

% Title
textc = '\textbf{c}\hspace{.2cm}combine to form curve';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});

% Combine
sc = .07;
sh1 = [.15;.81];
sh2 = [.15;.51];
sh3 = [.2;.18];
XP = XCa(:,:,1)*sc+sh1;
shI = [.25;0];
nP = 10;
hold on;
for i = 1:3
    scatter(linspace(0,shI(1),nP)+XP(1,i+1)+shI(1)*(i-1),ones(1,nP)*XP(2,i+1),.2,[1 1 1]*.8);
    scatter(linspace(0,shI(1),nP)+XP(1,i+2)+shI(1)*(i-1),ones(1,nP)*XP(2,i+2),.2,[1 1 1]*.8);
    visualize_network(XP(:,[0:2]+i)+shI*(i-1),XP(:,(-1:0)+size(Xsa,2)+2*i)+shI*(i-1),conn,...
                      'ucolor',CS(1+14*(mod(i,2)),:));
end

% initial network
for i = 1:4
    line_coordinates(XCa(:,[-1,1]+2*i,1)*sc+sh2,'lSh',-.01,'color',interp1(DLin,CP,Dc(i,1)));
    line_coordinates(XCa(:,[0,2]+2*i,1)*sc+sh2,'lSh',.01,'color',interp1(DLin,CP,Dc(i+4,1)));
end
visualize_network(XCa(:,1:size(Xsa,2),1)*sc+sh2,XCa(:,(size(Xsa,2)+1):end,1)*sc+sh2,conna,'ucolor',CSS);

% final network
for i = 1:4
    line_coordinates(XCa(:,[-1,1]+2*i,end)*sc+sh3,'lSh',-.01,'color',interp1(DLin,CP,Dc(i,end)));
    line_coordinates(XCa(:,[0,2]+2*i,end)*sc+sh3,'lSh',.01,'color',interp1(DLin,CP,Dc(i+4,end)));
end
visualize_network(XCa(:,1:size(Xsa,2),end)*sc+sh3,XCa(:,(size(Xsa,2)+1):end,end)*sc+sh3,conna,'ucolor',CSS);

% text
text(1,.52,'initial',NVTextRA{:});
text(1,.44,'position',NVTextRA{:});
text(1,.13,'final',NVTextRA{:});
text(1,.05,'position',NVTextRA{:});

% Axis
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0);


%% d: Range of solution spaces
pInd = 5;
subplot('position',subpN(pInd,:)); cla;

% Title
textd = '\textbf{d}\hspace{.2cm}range of solution spaces';
text(labX,subp(pInd,4)+labY,textd,NVTitle{:});

sc = .18;
sh1 = [0.5;0.5];
nM = 1;
thR = linspace(15,30,nSS*nM);
R2 = [-1 1; -1 1]*2.9;
for i = 1:nSS*nM
    Xfp = ds2*[[-sind(thR(i));-cosd(thR(i))] [0;0] [sind(thR(i));-cosd(thR(i))]] + [0;Xf2(2,2)];
    visualize_conic_finite(Xs1*sc+sh1(:,1),Xfp*sc+sh1(:,1),R2*sc+sh1(:,1),...
                           'ucolori',CS(ceil(i/nM),:),'ucolorf',[1 1 1],'overlay',.99);
    visualize_network(Xs1*sc+sh1(:,1),[],[1 1]);
    visualize_network(Xfp*sc+sh1(:,1),[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);
    drawnow;
end

% Axis
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0);


%% e: Tracing
pInd = 6;
subplot('position',subpN(pInd,:)); cla;


% Title
texte = '\textbf{e}\hspace{.2cm}tesselate desired curve';
text(labX,subp(pInd,4)+labY,texte,NVTitle{:});

sc = .02;
shl = [.7; .5];

% Hypotrochoid
Rv = 5;
rv = 3;
dv = 5;
nv = 1000000;
nCh = ceil(nv/20);
thv = linspace(0,2.09*pi*lcm(rv,Rv)/Rv,nv) + 2*pi*lcm(rv,Rv)/(10*Rv);
thp = [1 3 5 7 9] * 2*pi*lcm(rv,Rv)/(10*Rv);
f = @(tho)[(Rv-rv)*cos(tho) + dv*cos((Rv-rv)*tho/rv);...
           (Rv-rv)*sin(tho) - dv*sin((Rv-rv)*tho/rv)];
xvo = f(thv);
xvpo = f(thp);
nISh = 22;
scC = 2.7;
xv = xvo*scC;
xvp = xvpo*scC;
plot(xv(1,:)*sc+shl(1),xv(2,:)*sc+shl(2),'k-','linewidth',1);
hold on;

% Trace
nM = 1;
xP = zeros(2,nM+1);
% xP(:,1) = xv(:,1) - xvp(:,1)*.1378671673714;
% xP(:,1) = xv(:,1) - xvp(:,1)*.125577249384;
xP(:,1) = xv(:,1) - xvp(:,1)*.1487034107;
cI = 1;
i = 1;

% Optimization
options = optimset('TolFun',1e-15,'TolX',1e-15);
while(cI < nv - nCh)
    % Compute distance initial condition
    dP = sqrt(sum((xP(:,i)-xv(:,(1:nCh)+cI+10)).^2));
    dA = abs(dP-ds2/2);
    cIP = find(diff(sign(diff(dA)))==2,1)+11;
    cI = cI + cIP;

    if(isempty(cI))
        break;
    end
    % Find better minimum with optimization
    [thop, fV] = fminsearch(@(tho) abs(sqrt(sum((f(tho)*scC-xP(:,i)).^2))-ds2/2),...
                            thv(cI),options);
    xP(:,i+1) = xP(:,i) + 2*(f(thop)*scC-xP(:,i));
    i = i+1;
end

fvp = f(thp);
line([zeros(1,5); fvp(1,:)]+shl(1), [zeros(1,5); fvp(2,:)]+shl(2), 'color','k');


% Distances
cVec = sqrt(sum((xP(:,1:end-2) - xP(:,3:end)).^2));
cVecL = linspace(min(cVec)-.01,max(cVec)+.01,size(CSS,1));
nSk = 11;
for i = 1:size(xP,2)-2
    fill(xP(1,[0:2 0]+i)*sc+shl(1),xP(2,[0:2 0]+i)*sc+shl(2),...
         interp1(cVecL,CSS,cVec(i)),'linewidth',.5);
    if(mod(i,nSk) == 1)
        plot(xP(1,i)*sc+shl(1),xP(2,i)*sc+shl(2),'ks','markersize',5,'linewidth',5);
    end
end
% line([xP(1,1:end-1);xP(1,2:end)]*sc+sh(1), [xP(2,1:end-1);xP(2,2:end)]*sc+sh(2),...
%      'linewidth',.5,'color',CP(end,:));

% Text
text(0,.05,'hypotrochoid',NVTextR{:});

% Axis
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0);

disp(cVec(1:2*nSk))


%% Trifolium
pInd = 6;
% subplot('position',subpN(pInd,:)); cla;
figure(5); clf;


% Title
texte = '\textbf{e}\hspace{.2cm}tesselate desired curve';
text(labX,subp(pInd,4)+labY,texte,NVTitle{:});

sc = .03;
shl = [.7; .5];

% Hypotrochoid
Rv = 5;
rv = 3;
dv = 5;
nv = 1000000;
nCh = ceil(nv/20);
thv = linspace(0,2.1*pi,nv) + pi/4;
thp = [1/4 3/4 5/4 7/4] * pi;
a = 6;
f = @(tho)[a*sin(2*tho) .* cos(tho);...
           a*sin(2*tho) .* sin(tho)];
xvo = f(thv);
xvpo = f(thp);
nISh = 22;
scC = 2.7;
xv = xvo*scC;
xvp = xvpo*scC;
plot(xv(1,:)*sc+shl(1),xv(2,:)*sc+shl(2),'k-','linewidth',1);
hold on;

% Trace
nM = 1;
xP = zeros(2,nM+1);
% xP(:,1) = xv(:,1) - xvp(:,1)*.1378671673714;
% xP(:,1) = xv(:,1) - xvp(:,1)*.125577249384;
xP(:,1) = xv(:,1) + xvp(:,1)*.084543;
cI = 1;
i = 1;

% Optimization
options = optimset('TolFun',1e-15,'TolX',1e-15);
while(cI < nv - nCh)
    % Compute distance initial condition
    dP = sqrt(sum((xP(:,i)-xv(:,(1:nCh)+cI+10)).^2));
    dA = abs(dP-ds2/2);
    cIP = find(diff(sign(diff(dA)))==2,1)+11;
%     plot(xv(1,(1:nCh)+cI)*sc+shl(1),xv(2,(1:nCh)+cI)*sc+shl(2),'-','linewidth',1);
    cI = cI + cIP;

    if(isempty(cI))
        break;
    end
    % Find better minimum with optimization
    [thop, fV] = fminsearch(@(tho) abs(sqrt(sum((f(tho)*scC-xP(:,i)).^2))-ds2/2),...
                            thv(cI),options);
    xP(:,i+1) = xP(:,i) + 2*(f(thop)*scC-xP(:,i));
    i = i+1;
end

fvp = f(thp);
line([zeros(1,length(thp)); fvp(1,:)]+shl(1), [zeros(1,length(thp)); fvp(2,:)]+shl(2), 'color','k');
plot(xP(1,:)*sc+shl(1),xP(2,:)*sc+shl(2), 'ko', 'linewidth',4,'markersize',4);

% Distances
cVec = sqrt(sum((xP(:,1:end-2) - xP(:,3:end)).^2));
cVecL = linspace(min(cVec)-.01,max(cVec)+.01,size(CSS,1));
nSk = 28;
for i = 1:size(xP,2)-2
    fill(xP(1,[0:2 0]+i)*sc+shl(1),xP(2,[0:2 0]+i)*sc+shl(2),...
         interp1(cVecL,CSS,cVec(i)),'linewidth',.5);
    if(mod(i,nSk) == 1)
        plot(xP(1,i)*sc+shl(1),xP(2,i)*sc+shl(2),'ks','markersize',5,'linewidth',5);
    end
end
% line([xP(1,1:end-1);xP(1,2:end)]*sc+sh(1), [xP(2,1:end-1);xP(2,2:end)]*sc+sh(2),...
%      'linewidth',.5,'color',CP(end,:));

% Text
text(0,.05,'hypotrochoid',NVTextR{:});

% Axis
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0);

disp(cVec(1:2*nSk))


%% f: Solve for solution spaces
% pInd = 7;
% subplot('position',subpN(pInd,:)); cla;
% figure(4); clf;

[C,IA,IC] = unique(round(cVec,3),'stable');

Xuc = zeros(2,2,length(C));
nI = 1;
sc = .1;
shl = [.4;.5];
Xu0a = [-3.5*sh 3.5*sh; 2.2 2.2];
Xu0b = [-sh+.1 -sh+.1; -2 2.5];
% Xu0c = [-sh*2 sh*2;  -3.1 -3.1];
Xu0c = [-sh*2 sh*2;  -3.3 -3.3];
for i = 1:length(C)
    cla;
%     if(C(i)<1)
    if(C(i)<2.2)
        Xu0 = Xu0a;
    elseif(C(i)>4)
        Xu0 = Xu0c;
    else
        Xu0 = Xu0b;
    end
    for j = 1:nI
        Xfp = [[-C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)] [0;0] [C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)]] + [0;Xf2(2,2)];
        Xup = construct_network(Xs1,Xfp,Xu0,conn,0,1);
    end
    Xuc(:,:,i) = Xup(1:2,:);
    visualize_conic_finite(Xs1*sc+shl,Xfp*sc+shl,2*R*sc+shl);
    visualize_network(Xs1*sc+shl,Xuc(:,:,i)*sc+shl,conn,'ucolor',[0 0 0]);
    
    ax = [0 sRat(pInd) 0 1];
    axis(ax);
    drawnow;
end


%% Test solutions
figure(5); clf;
for i = 1:length(C)
    cla;
    X0 = zeros(2,5); X0(:,1) = -diff(Xs1(:,1:2),1,2);
    [XCp,fCp] = sim_motion(Xs1,Xuc(:,:,i),conn,.01,400,X0,0);
    % Distances
    Dp = squeeze(sqrt(sum(diff(XCp(:,[1 2 3 1],:),1,2).^2)));
    hold on;
%     plot([1.5 ds2], [1.5 ds2]);
    plot(Dp(1,:),Dp(3,:),'color','k');
    plot(ds2,C(i), 'ks','linewidth',5,'markersize',5);
    plot(ds2,ds2, 'ks','linewidth',5,'markersize',5);
    plot(Dp(1,:),Dp(2,:),'color',CP(75*i,:));
    hold off;
    axis([min(C)-.1 max(C)+.1 min(C)-.1 max(C)+.1]);
    drawnow;
    pause;
end


%% Assemble
% Concatenate
[Xscc,Xucc,conncc,CSSc] = network_chain_x(Xs1(:,1:2),Xuc,IC',interp1(cVecL,CSS,cVec));
[XCcc,fCcc] = sim_motion8(Xscc,Xucc,conncc,.5,3467,[Xscc Xucc],0);

% ax = [0 sRat(pInd) 0 1];
% axis(ax);

%% Plot
figure(4); clf;
animate_network(XCcc(:,:,1:3467),conncc,'nu',size(Xucc,2),'figwidth',20,'nframe',400,'save',1,'scolor',[1 1 1]*.85);


%% Save
fName = 'figure3f';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');
    
