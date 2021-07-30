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
subp = [[ 0.00 0.00 4.00 9.50];...
        [ 4.00 5.50 5.50 4.00];...
        [ 4.00 0.00 5.50 5.50];...
        [11.40 0.00 7.60 7.60]];
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
NVTextR = {'Units','Normalized','fontsize',FS};


%% Define module
sh = sqrt(3)/2;
L = 1;
LF = 1.7;
Xs2 = [-sh 0  sh;...
      -.5 1 -.5]*L;
Xf2 = Xs2 * LF;
ds2 = sqrt(sum((Xf2(:,1)-Xf2(:,2)).^2));

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
% figure(4); clf; plot(D2(1,:),D2(2,:)); hold on; plot(D2(2,:),D2(1,:)); hold off;


%% Color
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


%% a: design finite motion
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Plot parameters
sc = .08; 
sh1 = [.2;.76];
sh2 = [.2;.35];
C_I = [1 1 1]*.92;
C_F = [1 1 1]*.6;
C_SS = [170 150 255]/255;
R = [-1 1; -1 1]*2;
lSh = .006;

% Title
texta = '\textbf{a}\hspace{0.2cm}set node positions';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(.5,.95,'decrease distance $c$',NVTextH{:});
text(.5,.53,'increase distance $c$',NVTextH{:});
% Network 1
visualize_conic_finite(Xs1*sc+sh1,Xf1*sc+sh1,R*sc+sh1,...
                       'ucolori',CS(1,:),'ucolorf',[1 1 1],'overlay',.99);
line_coordinates(Xs1(:,[1,2])*sc+sh1,'lSh',lSh,'color',CPaI(1,:));
line_coordinates(Xs1(:,[2,3])*sc+sh1,'lSh',lSh,'color',CPaJ(1,:));
line_coordinates(Xs1(:,[3,1])*sc+sh1,'lSh',lSh,'color',CPaK(1,:));
line_coordinates(Xf1(:,[1,2])*sc+sh1,'lSh',lSh,'color',CPaI(end,:));
line_coordinates(Xf1(:,[2,3])*sc+sh1,'lSh',lSh,'color',CPaJ(end,:));
line_coordinates(Xf1(:,[3,1])*sc+sh1,'lSh',lSh,'color',CPaK(end,:));
visualize_network(Xs1*sc+sh1,[],[1 1]);
visualize_network(Xf1*sc+sh1,[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);
% Network 2
visualize_conic_finite(Xs2*sc+sh2,Xf2*sc+sh2,R*sc+sh2,...
                       'ucolori',CS(15,:),'ucolorf',[1 1 1],'overlay',.99);
line_coordinates(Xs2(:,[1,2])*sc+sh2,'lSh',lSh,'color',CPbI(1,:));
line_coordinates(Xs2(:,[2,3])*sc+sh2,'lSh',lSh,'color',CPbJ(1,:));
line_coordinates(Xs2(:,[3,1])*sc+sh2,'lSh',lSh,'color',CPbK(1,:));
line_coordinates(Xf2(:,[1,2])*sc+sh2,'lSh',lSh,'color',CPbI(end,:));
line_coordinates(Xf2(:,[2,3])*sc+sh2,'lSh',lSh,'color',CPbJ(end,:));
line_coordinates(Xf2(:,[3,1])*sc+sh2,'lSh',lSh,'color',CPbK(end,:));
visualize_network(Xs2*sc+sh2,[],[1 1]);
visualize_network(Xf2*sc+sh2,[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);

% Legend
yPos = 0;
visualize_network([.05;yPos+.12],[],[1 1],'scolor',C_I);
visualize_network([.05;yPos+.08],[],[1 1],'scolor',[1 1 1]*.0,'bcolor',[1 1 1]*.8);
scatter(linspace(0,.1,size(CP,1)), ones(1,size(CP,1))*yPos+.04, 1, interp1(1:nSS,CS,linspace(1,nSS,size(CP,1))), 'marker', 's');
scatter(linspace(0,.1,size(CP,1)), ones(1,size(CP,1))*yPos, 1, CP, 'marker', 's');
text(.33,yPos+.12,'initial position',NVTextR{:});
text(.33,yPos+.08,'final position',NVTextR{:});
text(.33,yPos+.04,'solution space',NVTextR{:});
text(.33,yPos,'distance',NVTextR{:});

% Text
text(.3,.77,'$d_1^0$',NVTextH{:},'color',CPaI(1,:));
text(.7,.77,'$d_2^0$',NVTextH{:},'color',CPaJ(1,:));
text(.39,.86,'$d_1^*$',NVTextH{:},'color',CPaI(end,:));
text(.61,.86,'$d_2^*$',NVTextH{:},'color',CPaJ(end,:));
text(.5,.688,'$c^0$',NVTextH{:},'color',CPaK(1,:));
text(.5,.64,'$c^*$',NVTextH{:},'color',CPaK(end,:));

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: Units
pInd = 2;
subplot('position',subpN(pInd,:)); cla;

% Title
textb = '\textbf{b}\hspace{0.2cm}add 2 nodes on solution space';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
text(.5,.45,'connect nodes with edges',NVTextH{:});


% plot
sc = .08;
sh1 = [0.4  0.4 ; .76 .2];
sh2 = [1.05 1.05; .76 .2];

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


ax = [0 sRat(pInd) 0 1];
axis(ax);


set(gca,'visible',0);



%% c: design infinitesimal
% Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
hold on;

% Title
textc = '\textbf{c}\hspace{.2cm}unit conformational motion';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});

% Map
CPa = winter(size(XC1,3));
CPb = winter(size(XC2,3));
hold on;
plot(D1(1,:),D1(2,:),'k-','linewidth',1);
plot(D2(1,:),D2(2,:),'k-','linewidth',1);
% scatter(Da(1,:),Da(2,:),1,CPa,'o','linewidth',.1);
% scatter(Db(1,:),Db(2,:),1,CPb,'o','linewidth',.1);

% Axes
aX0 = 1.3; aXF = 3.25;
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
axis([1.2 3.5 1.2 3.5]);

% Examples
nE = 4;
plInda = floor(linspace(1,size(XC1,3),nE));
plSha = [linspace(-.25,-.2,nE);...
         linspace(.15,.25,nE)];
plIndb = floor(linspace(1,size(XC2,3),nE));
plShb = [linspace(.2,.25,nE);...
         linspace(-.21,-.1,nE)];
sc3 = .1;
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


%% label
% % Plot
% pInd = 3;
% subplot('position',subpN(pInd,:)); cla;
% delete(findall(gca,'type','annotation'))
% hold on;
% 
% % Set final positions
% 
% 
% % Axis limits
% set(gca,'visible',0,'xtick',[],'ytick',[]);
% drawnow;




%% a: iterated map
% % Plot
% pInd = 1;
% subplot('position',subpN(pInd,:)); cla;
% delete(findall(gca,'type','annotation'))
% hold on;
% 
% % Title
% texta2 = '\textbf{a}\hspace{.2cm}conformational motion of a module';
% text(labX,subp(pInd,4)+labY,texta2,NVTitle{:});
% 
% % Map
% CP = winter(size(XC,3));
% aX0 = 0.8; aXF = 2.8; aXD = aXF - aX0;
% hold on;
% plot([aX0 aXF],[aX0 aXF], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
% scatter(D(1,:),D(2,:),1,CP,'o','linewidth',.1);
% 
% % Axes
% arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
% arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
% line([1 1]*D(1,1), [0 aXD/30]+aX0, 'color', 'k', 'linewidth', .5);
% line([1 1]*D(1,fpI), [0 aXD/30]+aX0, 'color', 'k', 'linewidth', .5);
% line([1 1]*D(1,fpL), [0 aXD/30]+aX0, 'color', 'k', 'linewidth', .5);
% line([0 aXD/30]+aX0, [1 1]*D(1,1), 'color', 'k', 'linewidth', .5);
% line([0 aXD/30]+aX0, [1 1]*D(1,fpI), 'color', 'k', 'linewidth', .5);
% line([0 aXD/30]+aX0, [1 1]*D(1,fpL), 'color', 'k', 'linewidth', .5);
% text(.1,.005,'$D_a$',NVTextR{:});
% text(.41,.005,'$D^*$',NVTextR{:});
% text(.61,.005,'$D_b$',NVTextR{:});
% text(.895,0.045,'$d_1$',NVTextR{:});
% text(.02,.92,'$d_2$',NVTextR{:});
% axis([aX0 aXF aX0 aXF] + [-.05 .15 -.05 .15]*aXD);
% 
% % Examples
% scc = .22;
% plInd = [1 fpI fpL];
% plSh = [ 0.40  0.35  0.20;...
%          0.20  0.20  0.20];
% for i = 1:length(plInd)
%     pI = plInd(i);
%     XCP = XC(:,:,pI)*scc + [D(1,pI);D(2,pI)]+plSh(:,i);
%     plot(D(1,pI),D(2,pI),'s','color',CP(pI,:),'linewidth',4,'markersize',4);
%     visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'lcolor',CP(pI,:),'scolor',[1 1 1]*.7);
%     line_coordinates(XCP(:,1:2), lSh, nW, lw, 'color', [1 1 1]*.5);
%     line_coordinates(XCP(:,2:3), lSh, nW, lw, 'color', [1 1 1]*.5);
%     if(i == 2)
%         visualize_network(XCP(:,1:3)-.8,XCP(:,4:5)-.8,conn,'scolor',[1 1 1]*.7);
%         line_coordinates(XCP(:,1:2)-.8, lSh, nW, lw, 'color', [1 1 1]*.5);
%         line_coordinates(XCP(:,2:3)-.8, lSh, nW, lw, 'color', [1 1 1]*.5);
%     end
% end
% 
% % Lables
% text(.07,.77,'$D_a$',NVTextR{:});
% text(.31,.84,'$D_b$',NVTextR{:});
% text(.445,.63,'$D^*$',NVTextR{:});
% text(.68,.63,'$D^*$',NVTextR{:});
% text(.65,.33,'$D_b$',NVTextR{:});
% text(.87,.25,'$D_a$',NVTextR{:});
% text(.11,.30,'$d_1$',NVTextR{:});
% text(.35,.30,'$d_2$',NVTextR{:});
% set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
% drawnow;

                        
%% Save
fName = 'figure3e';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');

