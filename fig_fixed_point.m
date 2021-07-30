% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(3); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 7.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 5.00  2.50 2.50];...
        [ 0.00 2.50  2.50 2.50];....
        [ 0.00 0.00  2.50 2.50];...
        [ 2.75 5.00  2.50 2.50];...
        [ 2.75 2.50  2.50 2.50];....
        [ 2.75 0.00  2.50 2.50];...
        [ 5.75 2.50  5.00 5.00];...
        [ 5.75 0.00  3.00 4.50];....
        [ 5.75 0.00  2.50 2.50]];
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','centimeters');


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
NVTextHu  = {'fontsize',FS,'HorizontalAlignment','center'};
NVTextRAu = {'fontsize',FS,'HorizontalAlignment','right'};
NVTextRu  = {'fontsize',FS};

% Color
nT = 1000;
o = [1 1 1];
% Gradient: Distance. Interpolate between 3 colors
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));
C2a = [133 37 83]/255;


%% Analysis
% Define random module
Xsr = [-1.0  0.2  0.5;...
       -0.5  0.0 -0.9];
Xur = [-0.5  0.5;...
       -0.3 -0.2];

% Parameters
L1 = 1;
L2 = 1.5;
Rz1 = rotz(-30); Rz1 = Rz1(1:2,1:2);
Rz2 = rotz(-17); Rz2 = Rz2(1:2,1:2);
os = [0;-1]*L1;
of = [0;-1]*L2;

% Single module
Xs = [Rz1*os [0;0] Rz1'*os];
Xf = [Rz2*of [0;0] Rz2'*of];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
ax = [-1 1 -1 1]*1.4 + [0 0 -1 -1]*.5;

% Solution space
R = [ax(1:2); ax(3:4)] + [1 -1]*.4;
Xui = [.5 0;...
       0.33 -1.15];
[Xu,~] = construct_network(Xs,Xf,Xui,conn,0,1);

% Simulate
[XC,fC] = sim_motion10(Xs,Xu(1:2,:),conn,.01,70,[Xs Xu(1:2,:)],0);
D = squeeze(sqrt(sum(diff(XC(:,[1 2 3],:),1,2).^2)));
DLin = linspace(L1,L2,nT);
CPa1 = interp1(DLin,CP1,D(1,:),'linear','extrap');
CPa2 = interp1(DLin,CP1,D(2,:),'linear','extrap');
% [~,mI] = min(sum(abs(D-L2)));


%% a: Random module
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));
% Plot
visualize_network(Xsr,[],[1 1]);
line_coordinates(Xsr(:,1:2), 'lSh',.22, 'lw',.5,'style','-','nw',.05);
line_coordinates(Xsr(:,2:3), 'lSh',.22, 'lw',.5,'style','-','nw',.05);
% Text
texta = '\textbf{a}\hspace{.2cm}require 2 unique lengths';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(Xsr(1,1)-.3,Xsr(2,1),'$a$',NVTextHu{:},'color',o*gr^4);
text(Xsr(1,2)+.1,Xsr(2,2)+.35,'$b$',NVTextHu{:},'color',o*gr^4);
text(Xsr(1,3)+.25,Xsr(2,3)-.15,'$c$',NVTextHu{:},'color',o*gr^4);
text(mean(Xsr(1,1:2))-.1,mean(Xsr(2,1:2))+.55,'$l_k$',NVTextHu{:});
text(mean(Xsr(1,2:3))+.35,mean(Xsr(2,2:3))+.1,'$l_{k+1}$',NVTextRu{:});
% Axes
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% b: Start
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
% Plot
line_coordinates(Xs(:,[1,2]),'lSh',.05,'color',CPa1(1,:));
line_coordinates(Xs(:,[2,3]),'lSh',.05,'color',CPa2(1,:));
visualize_network(Xs,[],[1 1]);
% Text
texta = '\textbf{b}\hspace{.2cm}fixed point: $l_k = l_{k+1}$';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(0,.45,'$l_k = l_{k+1} = l^\bullet$',NVTextHu{:});
% text(mean(Xs(1,2:3))+.75,mean(Xs(2,2:3))+.2,'$l_2 = l^\bullet$',NVTextHu{:});
% Axes
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% c: Start
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
% Plot
line_coordinates(Xf(:,[1,2]),'lSh',.05,'color',CPa1(end,:));
line_coordinates(Xf(:,[2,3]),'lSh',.05,'color',CPa2(end,:));
visualize_network(Xf,[],[1 1]);
% Text
texta = '\textbf{c}\hspace{.75cm}end';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(0,.45,'$l_k = l_{k+1} = l^\circ$',NVTextHu{:});
% Axes
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% d: Added
pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
% Plot
visualize_network(Xs,Xur,conn,'ucolor',C2a);
% Text
texta = '\textbf{d}\hspace{.2cm}4 motions';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(Xs(1,1)-.3,Xs(2,1)-.1,'$a$',NVTextHu{:},'color',o*gr^4);
text(Xs(1,2),Xs(2,2)+.35,'$b$',NVTextHu{:},'color',o*gr^4);
text(Xs(1,3)+.3,Xs(2,3)-.1,'$c$',NVTextHu{:},'color',o*gr^4);
text(Xur(1,1)-.3,Xur(2,1)+.1,'$d$',NVTextHu{:},'color',C2a);
text(Xur(1,2)+.3,Xur(2,2)+.1,'$e$',NVTextHu{:},'color',C2a);
% Axes
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% e: Start
pInd = 5;
subplot('position',subpN(pInd,:)); cla; hold on;
% Plot
visualize_conic_finite(Xs,Xf,R,'ucolori',C2a,'ucolorf',o,'overlay',.99);
visualize_network(Xs,Xu(1:2,:),conn,'ucolor',C2a);
% Text
texta = '\textbf{e}\hspace{.2cm}start space';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
% Axes
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% f: Start
pInd = 6;
subplot('position',subpN(pInd,:)); cla; hold on;
% Plot
visualize_conic_finite(Xf,Xs,R,'ucolori',C2a,'ucolorf',o,'overlay',.99);
visualize_network(Xf,Xu(3:4,:),conn,'ucolor',C2a);
% Text
texta = '\textbf{f}\hspace{.2cm}end space';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
% Axes
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% g: Start
pInd = 7;
subplot('position',subpN(pInd,:)); cla; hold on;
aX0 = .77; aXF = 1.65;                  % Axes min and max locations
nE = 4;
sc = .13;
pI = floor(linspace(1,size(XC,3),nE));
sha = [linspace(-.1,-.03,nE);...
       linspace(.2^5,.20^5,nE).^(1/5)] + D(:,pI);
% Plot
scatter(D(1,:),ones(1,size(D,2))*aX0+.012,20,interp1(DLin,CP1,D(1,:)),...
        'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.012,D(2,:),20,interp1(DLin,CP1,D(1,:)),...
        'filled','marker','s');
line([1 1]*L1, [0;.03]+aX0, 'color', CP1(1,:));
line([1 1]*L2, [0;.03]+aX0, 'color', CP1(end,:));
line([0;.03]+aX0, [1 1]*L1, 'color', CP1(1,:));
line([0;.03]+aX0, [1 1]*L2, 'color', CP1(end,:));
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
hold on;
plot([aX0 aXF],[aX0 aXF], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);
for i = 1:nE
    scatter(D(1,pI(i)),D(2,pI(i)),40,'k','filled','linewidth',.5,'marker','s');
    visualize_network(XC(:,1:3,pI(i))*sc+sha(:,i),...
                      XC(:,4:5,pI(i))*sc+sha(:,i),conn,'ucolor',C2a);
end
scatter(D(1,pI([1 4])),D(2,pI([1 4])),80,CP1([1 end],:),...
       'filled','marker','s','linewidth',.5);
scatter(D(1,pI([1 4])),D(2,pI([1 4])),40,ones(2,3),'marker','s','linewidth',.1);
hold off;

% Text
texta = '\textbf{g}\hspace{.2cm}conformational motion';
text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
text(aXF+.03,aX0,'$l_k$',NVTextRu{:});
text(aX0,aXF+.05,'$l_{k+1}$',NVTextRu{:});
text(aX0-.06,L1,'$l^{\bullet}$',NVTextRu{:},'color',CP1(1,:));
text(aX0-.06,L2,'$l^{\circ}$',NVTextRu{:},'color',CP1(end,:));
text(D(1),D(1)-.08,'$l_{k+1} = f(l_k)$',NVTextRu{:});
text(mean([aXF,aX0]),mean([aXF,aX0])-.08,'$l_k = l_{k+1}$',NVTextRu{:},...
     'rotation',45,'color',o*gr^2);
% Axes
axis([.75 1.8 .75 1.8]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%%
% 
% % Double module
% Xsc = [0  0  0  2*L  L  L;...
%        L -L  0  0    0 -2*L];
% connc = [1 3; 1 4; 2 3; 2 4; 3 5; 3 6; 4 5; 4 6];
% 
% % Tesselate module
% [Xsa, conna] = tesselate_network(Rz1*Xsc, connc, [sqrt(2);0], [4;1]);
% 
% 
% %% Simulate networks
% % Initial conditions
% Xs0 = [-1 -1  0  0;...
%         0  0  0  0];
% Xsc0 = [0  0  0  0  0  0;...
%         0  0  0  0 -1  0];
% Xsa0 = 0*Xsa; Xsa0(1,end-1) = 1;
% 
% % Simulate networks
% [XC ,fC ] = sim_motion10(Xs, [],conn, .01,77 ,Xs0 ,0);        % 1 unit
% [XCa,fCa] = sim_motion10(Xsa,[],conna,.01,438,Xsa0,0);  % Chain
% 
% % Simulation errors
% disp(['mean simulation error: '...
%      [num2str(mean(fC)) '  ' num2str(mean(fCa))]]);
% 
%  
% % Translate networks to align with origin
% for i = 1:size(XC,3)
%     XC(:,:,i) = XC(:,:,i) - XC(:,4,i);
% end
% for i = 1:size(XCa,3)
%     XCa(:,:,i) = XCa(:,:,i) - XCa(:,1,i);
% end
% 
% % Calculate distances between nodes
% D = squeeze(sqrt(sum(diff(XC,1,2).^2)));                % 1 unit
% D = D([1,3],:);
% Da = squeeze(sqrt(sum(diff(XCa,1,2).^2)));              % CHain
% Da = Da(1:2:end,:);
% 
% % Colors
% CPI = interp1(DLin,CP1,D(1,:));
% CPJ = interp1(DLin,CP1,D(2,:));
% 
% 
% %% a: single module conformational motion
% pInd = 1;
% subplot('position',subpN(pInd,:)); cla; hold on;
% delete(findall(gca,'type','annotation'));
% 
% % Parameters
% sc = .075;          % Scale drawing
% aX0 = 1.27;         % Map axes start
% aXF = 2.07;         % Map axes end
% 
% % Draw map ticks
% line(sqrt([2 2]), [0 .04]+aX0, 'color', interp1(DLin,CP1,sqrt(2)),...
%      'linewidth', .5);
% line([2 2], [0 .04]+aX0, 'color', interp1(DLin,CP1,2),...
%      'linewidth', .5);
% line([0 .04]+aX0, sqrt([2 2]), 'color', interp1(DLin,CP1,sqrt(2)),...
%      'linewidth', .5);
% line([0 .04]+aX0, [2 2], 'color', interp1(DLin,CP1,2),...
%      'linewidth', .5);
% 
% % Draw distance colormap
% scatter(D(1,:),ones(1,size(D,2))*aX0+.012,20,interp1(DLin,CP1,D(1,:)),...
%         'filled','marker','s');
% scatter(ones(1,size(D,2))*aX0+.012,D(2,:),20,interp1(DLin,CP1,D(2,:)),...
%         'filled','marker','s');
%  
% % Draw map axes
% line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
% arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
% arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
%     
% % Draw map and y=x lines
% plot([1.3 2.0],[1.3 2.0], '-', 'color', o*gr, 'linewidth',.5);
% plot(D(1,:),D(2,:),'k-','linewidth',1);
% 
% % Draw spline
% sCo11 = [0.70 0.50 1.00;...
%          0.00 0.50 1.00] .* [.09;-.17] + [1.87;2.03];
% plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
%             'linewidth',.5,'headwidth',5,'headlength',5);
% 
% % Draw example networks
% plInd = ceil([1 23 48 size(XC,3)]);
% plSh = [ .21  .2 .22 .3;...
%         -.1 -.1 -.09 -.04];
% lSh = -.01;
% for i = 1:length(plInd)
%     pI = plInd(i);
%     XP = XC(:,:,pI)*sc + [D(1,pI);D(2,pI)]+plSh(:,i);
%     line_coordinates(XP(:,1:2), 'lSh',lSh, 'lw',1, 'color',CPI(pI,:));
%     line_coordinates(XP(:,3:4), 'lSh',-lSh, 'lw',1, 'color',CPJ(pI,:));
%     visualize_network(XP,[],conn);
% end
% scatter(D(1,plInd(2:3)),D(2,plInd(2:3)),40,zeros(2,3),...
%        'filled','marker','s','linewidth',.5);
% scatter(D(1,plInd([1 4])),D(2,plInd([1 4])),80,...
%        [interp1(DLin,CP1,D(1,plInd(1)));interp1(DLin,CP1,D(1,plInd(4)))],...
%        'filled','marker','s','linewidth',.5);
% scatter(D(1,plInd([1 4])),D(2,plInd([1 4])),40,ones(2,3),'marker','s','linewidth',.1);
%    
% % Draw large template network
% XP = Xs*sc*2 + [D(1,pI);D(2,pI)]+[.08;.57];
% visualize_network(XP,[],conn,'lalpha',.3,'scolor',o*gr^.3);
% line_coordinates(XP(:,1:2), 'lSh',-.035, 'lw',.5,'style','-','nw',.01);
% line_coordinates(XP(:,3:4), 'lSh',0.035, 'lw',.5,'style','-','nw',.01);
% 
% % Text
% texta = '\textbf{a}\hspace{.2cm}each unit $k$ maps $l_k$ to $l_{k+1}$';
% text(labX,subp(pInd,4)+labY,texta,NVTitle{:});
% text(.011,.15,'$\sqrt{2}$',NVTextRA{:},'color',interp1(DLin,CP1,sqrt(2)));
% text(.011,.7,'$2$',NVTextRA{:},'color',interp1(DLin,CP1,2));
% text(.83,0.02,'$l_k$',NVTextR{:});
% text(-.087,.84,'$l_{k+1}$',NVTextR{:});
% text(.23,.30,'$l_{k+1} = l_k$',NVTextR{:},'rotation',45,'color',o*gr);
% text(.401,.8,'$l_{k+1} = f(l_k)$',NVTextR{:});
% text(.14,0.71,'$l_k$',NVTextR{:});
% 
% % Axes
% axis([1.24 2.3 1.24 2.3]);
% set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
% drawnow;
% 
% 
% %% b: single module with defined distances
% pInd = 2;
% subplot('position',subpN(pInd,:)); cla; hold on;
% delete(findall(gca,'type','annotation'));
% 
% % Parameters
% sc = .2;                % Scale network
% sh1 = [.50; .42];       % Shift network
% sh2 = [1.65; .42];
% lSh = .1; nW = .02;     % Distance bar
% 
% % Draw colored markers
% scatter(-.10,.75,80,interp1(DLin,CP1,D(1,plInd(1))),...
%         'filled','marker','s','linewidth',.5,'clipping',0);
% scatter(-.10,.75,40,[1 1 1],'marker','s','linewidth',.1,'clipping',0);
% scatter(1.02,.75,80,interp1(DLin,CP1,D(1,plInd(end))),...
%         'filled','marker','s','linewidth',.5,'clipping',0);
% scatter(1.02,.75,40,[1 1 1],'marker','s','linewidth',.1,'clipping',0);
%     
% % Draw module 1 and motion
% pI1 = 30;
% line_coordinates(XC(:,1:2,1)*sc+sh1,'lSh',-0.02,'color',CPI(1,:));
% line_coordinates(XC(:,3:4,1)*sc+sh1,'lSh',0.02,'color', CPJ(1,:));
% line([-.07 -.09 -.09 -.07]+XC(1,1,1)*sc+sh1(1), squeeze(XC(2,1,[1 1 pI1 pI1]))*sc+sh1(2),'color','k','linewidth',.4);
% line([-.07 -.09 -.09 -.07]+XC(1,2,1)*sc+sh1(1), squeeze(XC(2,2,[1 1 pI1 pI1]))*sc+sh1(2),'color','k','linewidth',.4);
% line(squeeze(XC(1,3,[1 1 pI1 pI1]))*sc+sh1(1), [.48 .50 .50 .48]+XC(2,2,1)*sc+sh1(2),'color','k','linewidth',.4);
% visualize_network(XC(:,:,pI1)*sc+sh1,[],conn,'lalpha',al,...
%                   'scolor',o*sqrt(gr));
% visualize_network(XC(:,:,1)*sc+sh1,[],conn);
% 
% 
% % Draw module 2 and motion
% pI2 = size(XC,3)-30;
% line_coordinates(XC(:,1:2,end)*sc+sh2,'lSh',-0.02,'color',CPI(end,:));
% line_coordinates(XC(:,3:4,end)*sc+sh2,'lSh',0.02,'color', CPJ(end,:));
% line([-.07 -.09 -.09 -.07]+XC(1,1,end)*sc+sh2(1), squeeze(XC(2,1,[end end pI2 pI2]))*sc+sh2(2),'color','k','linewidth',.4);
% line([-.07 -.09 -.09 -.07]+XC(1,2,end)*sc+sh2(1), squeeze(XC(2,2,[end end pI2 pI2]))*sc+sh2(2),'color','k','linewidth',.4);
% line(squeeze(XC(1,3,[end end pI2 pI2]))*sc+sh2(1), [.38 .40 .40 .38]+XC(2,2,end)*sc+sh2(2),'color','k','linewidth',.4);
% visualize_network(XC(:,:,pI2)*sc+sh2,[],conn,'lalpha',al,...
%                   'scolor',o*sqrt(gr));
% visualize_network(XC(:,:,end)*sc+sh2,[],conn);
% 
% % Text
% textb = '\textbf{b}\hspace{.2cm}stability of unit geometry';
% text(labX,subp(pInd,4)+labY,textb,NVTitle{:});
% text(-0.075, 0.42,'$d l_k$',NVTextR{:});
% text( 0.110, 0.66,'$d l_{k+1}$',NVTextR{:});
% text( 0.540, 0.42,'$d l_k$',NVTextR{:});
% text( 0.76, 0.66,'$d l_{k+1}$',NVTextR{:});
% text( 0.138, 0.06,'unstable',NVTextH{:});
% text( 0.138,-0.09,'$|d l_k| < |d l_{k+1}|$',NVTextH{:});
% text( 0.755, 0.06,'stable',NVTextH{:});
% text( 0.755,-0.09,'$|d l_k| > |d l_{k+1}|$',NVTextH{:});
% 
% % Axes
% axis([0 sRat(pInd) 0 1]);
% set(gca,'visible',0,'xtick',[],'ytick',[]);
% drawnow;
% 
% 
% %% c-e: module combination
% pInd = 3;
% subplot('position',subpN(pInd,:)); cla; hold on;
% delete(findall(gca,'type','annotation'));
% 
% % Parameters
% sc = 1/15;                                      % Scale network
% sh = [[.05;.88] [.05;.51] [.05;.12]];             % Shift network: initial
% shi = [.35;0];                                  % Shift network: indexed
% scX = [linspace(0,shi(1),15); zeros(1,15)];     % Dots connecting nodes
% lSh = .007;                                     % Distance bar offset
% grL = [gr*ones(15,2) linspace(gr,1,15)'];       % Fading of dots
% 
% % Draw modules
% pSh = [1 360 size(XCa,3)];                     % Index of lattice to draw
% sPr = {'$l_1$', '$l_2=f(l_1)$', '$l_3=f(l_2)$'};              % Text of length names
% for i = 1:3                                     % Iterate across modules
%     for j = 1:3                                 % Iterate across lattices
%         XCp = XCa(:,(1:6)+2*(j-1),pSh(i))*sc+sh(:,i) + shi*(j-1);
%         % Drwa dotted lines connecting module nodes
%         scatter(scX(1,:)+XCp(1,3),scX(2,:)+XCp(2,3),.2,o.*grL(:,j));
%         scatter(scX(1,:)+XCp(1,4),scX(2,:)+XCp(2,4),.2,o.*grL(:,j));
%         % Draw module connection animation with transparent modules
%         if(i == 1 & j == 1)
%             visualize_network(XCp(:,3:6),[],conn,'lalpha',(1-gr)/2,'scolor',o*gr^.2);
%             visualize_network(XCp(:,3:6)+[.03;0],[],conn,'lalpha',(1-gr),'scolor',o*gr^.3);
%             line_coordinates(XCp(:,3:4), 'lSh',lSh, 'lw',1,'color', interp1(DLin,CP1,Da(j,pSh(i))));
%         end
%         if(j == 3)
%             arrow([1.05 1.15], [1 1]*mean(XCp(2,3:4)), sRat(pInd),'linewidth',2,'headwidth',10,'headlength',8,'lp',.8,'color',[1 1 1]*.95);
%         end
%         % Draw length measurements
%         line_coordinates(XCp(:,1:2), 'lSh',(-1)^j*lSh, 'lw',1,...
%                          'color', interp1(DLin,CP1,Da(j,pSh(i))));
%         % Draw networks
%         visualize_network(XCp(:,1:4),[],conn);
%         % Text
%         text((XCp(1,2))/sRat(pInd)-.035-.002*(-1)^j,...
%               XCp(2,2)-.03*(-1)^(j)-.005,sPr{j},NVTextRA{:},...
%               'color', interp1(DLin,CP1,Da(j,pSh(i))));
%     end
% end
% 
% % Draw lattice
% sh2 = [-.125 -.11 -.05];                % Shift: Lattice
% sh3 = [1.22 1.22 1.22;...
%         0.00 0.00 0.00];
% shy = .188;
% for i = 1:3                         % Iterate across lattices
%     XCp = XCa(:,:,pSh(i))*sc+sh(:,i) + sh3(:,i);
%     for j = 1:size(Da,1)            % Iterate across distances
%         % Draw distance measurement
%         line_coordinates(XCp(:,(-1:0)+2*j),'lSh',(-1)^j*lSh,'lw',1,...
%                          'color',interp1(DLin,CP1,Da(j,pSh(i))));
%         % x-axis text
%         if(j < 3 || j == 9)
%             text(XCp(1,-1+2*j)+.03,sh(2,i)-.11+sh2(i),num2str(j),...
%                  'fontsize',FS,'horizontalalignment','center');
%         elseif(j==4)
%             text(XCp(1,-1+2*j)+.03,sh(2,i)-.11+sh2(i),'$\cdots$',...
%                  'fontsize',FS,'horizontalalignment','center');
%         end
%     end
%     % Draw distance squares on plot
%     scatter(XCp(1,1:2:end)+.03,(sh(2,i)-shy+sh2(i))+Da(:,pSh(i))'/11,20,...
%             interp1(DLin,CP1,Da(:,pSh(i))),'filled','marker','s','linewidth',2);
%     % Draw network
%     visualize_network(XCp,[],conna);
%     % Draw tick marks
%     axV = [0 .63 0 .14] + [1.23 1.23 sh(2,i)-.08+sh2(i) sh(2,i)-.08+sh2(i)];
%     line([0 .01]+axV(1), [0 0]+sh(2,i)-shy+sh2(i)+2/11,...
%          'color',interp1(DLin,CP1,2));
%     line([0 .01]+axV(1), [0 0]+sh(2,i)-shy+sh2(i)+sqrt(2)/11,...
%          'color',interp1(DLin,CP1,sqrt(2)));
%     line([1;1].*(XCp(1,1:2:end)+.03),...
%          [0;.01].*ones(1,size(XCp(1,1:2:end),2))+sh(2,i)-.08+sh2(i),...
%          'color','k','clipping',0);
%     % Draw plot axes
%     line(axV([1 1 2]),axV([4 3 3]),'color','k','linewidth',.5,'clipping',0);
%     % Text
%     text(1.22,sh(2,i)-.19+sh2(i)+2/11,'$2$','fontsize',FS,...
%          'horizontalalignment','right','color',interp1(DLin,CP1,2));
%     text(1.22,sh(2,i)-.19+sh2(i)+sqrt(2)/11,'$\sqrt{2}$','fontsize',FS,...
%          'horizontalalignment','right','color',interp1(DLin,CP1,sqrt(2)));
%     text(1.22,sh(2,i)-.19+sh2(i)+2.7/11,'$l_k$','fontsize',FS,...
%          'horizontalalignment','right');
%     text(1.87,sh(2,i)-.08+sh2(i),'$k$','fontsize',FS);
% end
% 
% % Text
% textc = ['\textbf{c}\hspace{.2cm}fixed point: if $l_1 = 2$, ',...
%          'then $l_{k+1} = f(2) = 2$'];
% textd = ['\textbf{d}\hspace{.2cm}unstable: if $l_1 < 2$, ',...
%          'then $l_{k+1}$ moves away from 2'];
% texte = ['\textbf{e}\hspace{.2cm}fixed point: if $l_1 = \sqrt{2}$, ',...
%          'then $l_{k+1} = f\left(\sqrt{2}\right) = \sqrt{2}$'];
% text(0,.974,textc,NVTextR{:});
% text(0,.6,textd,NVTextR{:});
% text(0,.24,texte,NVTextR{:});
% text(.035,.71,'glue together',NVTextR{:},'color',o*gr^2);
% text(.8,.974,'periodic lattice',NVTextH{:});
% text(.8,.24,'periodic lattice',NVTextH{:});
% 
% % Axes
% axis([0 sRat(pInd) 0 1]);
% set(gca,'visible',0,'xtick',[],'ytick',[]);
% drawnow;


%% Size and Save Figure
fName = 'fig_fixed_point';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');