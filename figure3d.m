% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(3); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 7.6];
% Margins in cm, [l,r,d,u]
fMarg = [.0 .4 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 0.00 3.80 7.60];...
        [ 3.80 0.00 7.60 7.60];...
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
Xs = [-sh 0  sh;...
      -.5 1 -.5]*L;
Xf = Xs * LF;
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

Xua = [-sh*L -sh*L ;...
       -sqrt((L*LF)^2 - (sh*L)^2) sqrt((L*LF)^2 - (sh*L)^2)];
Xub = [0    LF/sqrt(2);...
       -LF  LF/sqrt(2)];
   
% Simulation
[XCa,fCa] = sim_motion(Xs,Xua,conn,.001,1808,[Xs Xua],0);
[XCb,fCb] = sim_motion(Xs,Xub,conn,.001,1691,[Xs Xub],0);
Da = squeeze(sqrt(sum(diff(XCa(:,[1 2 3 1],:),1,2).^2)));
Db = squeeze(sqrt(sum(diff(XCb(:,[1 2 3 1],:),1,2).^2)));
% figure(4); clf; plot(Da(1,:),Da(2,:)); hold on; plot(Db(1,:),Db(2,:)); plot([1.7 3],[1.7 3]); hold off;


%% a: design finite motion
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
texta = '\textbf{a}\hspace{0.25cm}set node positions';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});

% Plot parameters
sc = .1; 
sh = [.24;.75];
sh2 = [.24;.75];
C_I = [1 1 1]*.92;
C_F = [1 1 1]*.6;
C_SS = [170 150 255]/255;

% Draw network
visualize_conic_finite(Xs*sc+sh2,Xf*sc+sh2,[-5 5;-5 5]*sc+sh2,...
        'scolori',C_I,'scolorf',C_F,'ucolorf',[1 1 1],'ucolori',C_SS);
line_coordinates(Xs(:,[1 2])*sc+sh2, lSh*sc, nW*sc, lw, 'color', [1 1 1]*.7);
line_coordinates(Xs(:,[2 3])*sc+sh2, lSh*sc, nW*sc, lw, 'color', [1 1 1]*.7);
line_coordinates(Xs(:,[1 3])*sc+sh2,-lSh*sc, nW*sc, lw, 'color', [1 1 1]*.7);
% line_coordinates(Xf(:,[1 2])*sc+sh2, lSh*sc, nW*sc, lw, 'color', [1 1 1]*.7);
% line_coordinates(Xf(:,[2 3])*sc+sh2, lSh*sc, nW*sc, lw, 'color', [1 1 1]*.7);
% line_coordinates(Xf(:,[1 3])*sc+sh2,-lSh*sc, nW*sc, lw, 'color', [1 1 1]*.7);
text(.26,.81,'$d_1$',NVTextR{:});
text(.67,.81,'$d_2$',NVTextR{:});
text(.48,.70,'$c$',NVTextR{:});

% Legend
visualize_network([.07;.52],[],[1 1],'scolor',C_I);
visualize_network([.07;.46],[],[1 1],'scolor',C_F);
line([.02 .12],[.4 .4],'linewidth',1,'color',C_SS);
text(.3,.52,'initial position',NVTextR{:});
text(.3,.46,'final position',NVTextR{:});
text(.3,.4,'solution space',NVTextR{:});

% Units
sc2 = .5;
sha = [.14;.085]; shb = [.38;.085];
visualize_conic_finite(Xs*sc*sc2+sha,Xf*sc*sc2+sha,[-5 5;-5 5]*sc*sc2+sha,...
        'ucolorf',[1 1 1],'ucolori',C_SS,'overlay',.2);
visualize_network(Xs*sc*sc2+sha,Xua*sc*sc2+sha,conn,'scolor',C_I,'ucolor',C_SS);
visualize_conic_finite(Xs*sc*sc2+shb,Xf*sc*sc2+shb,[-5 5;-5 5]*sc*sc2+shb,...
        'ucolorf',[1 1 1],'ucolori',C_SS,'overlay',.2);
visualize_network(Xs*sc*sc2+shb,Xub*sc*sc2+shb,conn,'scolor',C_I,'ucolor',C_SS);
textb = '\textbf{b}\hspace{0.25cm}place 2 more nodes';
text(labX,subp(pInd,4)-5.2,textb,NVTitle{:});
text(.3,.21,'unit a',NVTextH{:});
text(.81,.21,'unit b',NVTextH{:});

% Axis limits
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% c: design infinitesimal
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textc = '~\textbf{c}\hspace{.15cm}both units reach initial and final positions';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});

% Map
CPa = winter(size(XCa,3));
CPb = winter(size(XCb,3));
hold on;
scatter(Da(1,:),Da(2,:),1,CPa,'o','linewidth',.1);
scatter(Db(1,:),Db(2,:),1,CPb,'o','linewidth',.1);

% Axes
aX0 = 1.3; aXF = 3.25;
plot([aX0 aXF],[aX0 aXF], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line([1 1]*Da(1,1), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([1 1]*Da(1,end), [0 .03]+aX0, 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [1 1]*Da(1,1), 'color', 'k', 'linewidth', .5);
line([0 .03]+aX0, [1 1]*Da(1,end), 'color', 'k', 'linewidth', .5);
text(.91,0.045,'$d_1$',NVTextR{:});
text(.03,.93,'$d_2$',NVTextR{:});
text(.56,.44,'$d_2 = f_a(d_1)$',NVTextH{:},'rotation',45);
text(.45,.55,'$d_2 = f_b(d_1)$',NVTextH{:},'rotation',45);
text(Da(1,1),aX0-(aXF-aX0)*.05,num2str(Da(1,1),3),'Horizontalalignment','center','fontsize',FS);
text(Da(1,end),aX0-(aXF-aX0)*.05,num2str(Da(1,end),3),'Horizontalalignment','center','fontsize',FS);
axis([1.2 3.5 1.2 3.5]);

% Examples
nE = 4;
plInda = floor(linspace(1,size(XCa,3),nE));
plSha = [linspace(.2,.25,nE);...
         linspace(-.21,-.1,nE)];
plIndb = floor(linspace(1,size(XCb,3),nE));
plShb = [linspace(-.25,-.2,nE);...
         linspace(.15,.25,nE)];
sc3 = .1;
for i = 1:nE
    % unit a
    pIa = plInda(i);
    plot(Da(1,pIa),Da(2,pIa),'s','color',CPa(pIa,:),'linewidth',3,'markersize',3);
    if(i<nE)
        XCP = XCa(:,:,pIa)*sc3+Da(1:2,plInda(i+1))+plSha(:,i+1);
        visualize_network(XCP(:,1:3),XCP(:,4:5),conn,...
            'scolor',C_I,'ucolor',C_SS,'lcolor',CPa(pIa,:),'nalpha',.2,'lalpha',.2);
    end
    XCP = XCa(:,:,pIa)*sc3+Da(1:2,pIa)+plSha(:,i);
    visualize_network(XCP(:,1:3),XCP(:,4:5),conn,...
        'scolor',C_I,'ucolor',C_SS,'lcolor',CPa(pIa,:));
    % unit b
    pIb = plIndb(i);
    plot(Db(1,pIb),Db(2,pIb),'s','color',CPb(pIb,:),'linewidth',3,'markersize',3);
    if(i<nE)
        XCP = XCb(:,:,pIb)*sc3+Db(1:2,plIndb(i+1))+plShb(:,i+1);
        visualize_network(XCP(:,1:3),XCP(:,4:5),conn,...
            'scolor',C_I,'ucolor',C_SS,'lcolor',CPb(pIb,:),'nalpha',.2,'lalpha',.2);
    end
    XCP = XCb(:,:,pIb)*sc3+Db(1:2,pIb)+plShb(:,i);
    visualize_network(XCP(:,1:3),XCP(:,4:5),conn,...
        'scolor',C_I,'ucolor',C_SS,'lcolor',CPb(pIb,:));
end

% Axis limits
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% label
% Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Set final positions


% Axis limits
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;




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
fName = 'figure3d';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');

