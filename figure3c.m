% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;
set(groot,'defaulttextinterpreter','latex');


%% Figure dimensions
fig = figure(3); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 6.5];
% Margins in cm, [l,r,d,u]
fMarg = [.0 .4 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 0.00 6.50 6.50];...
        [ 6.10 0.00 4.75 6.50];...
        [10.45 0.00 4.75 6.50];...
        [14.80 0.00 4.75 6.50]];
% Fontsize
FS = 10;
% Distance visualization parameters
lSh = .1;
nW = .03;
lw = .5;
gr = 0.8;
% Scaling parameters
sc = .2;
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
L = 1; s = sqrt(3);
Xs = [-L -L  L;...
       0  L  0];
Xu = [ .4 0;...
       .4 -.8]*L;
XuF = [ .4 0;...
       -.4 -.8];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];

Xsc = [-L -L L L;...
        0  L 0 L];
% Limit cycle module
Xuc = [Xu [XuF(1,:); -XuF(2,:)+L]];
connc = [1 5; 2 5; 3 5; 1 6; 2 6; 3 6; 2 7; 3 7; 4 7; 2 8; 3 8; 4 8];
% Lattice
[Xsa,Xua,conna] = tesselate_network_old(Xsc,Xuc,connc,[2*L;0],[5,1]);

% Simulate
[XC,fC] = sim_motion(Xs,Xu,conn,.001,1645,[Xs Xu],0);
[XCc,fCc] = sim_motion(Xsc,Xuc,connc,.001,300,[Xsc Xuc],0);
[XCa,fCa] = sim_motion(Xsa,Xua,conna,.01,496,[Xsa Xua],1);

% Correct offset in XC
for i = 1:size(XC,3)
    XC(:,:,i) = XC(:,:,i) - mean(XC(:,[1 3],i),2);
    Rz = rotz(atan2d(diff(XC(2,[1 3],i)),diff(XC(1,[1 3],i)))); 
    Rz = Rz(1:2,1:2);
    XC(:,:,i) = Rz'*XC(:,:,i);
end
% Correct offset in XCa
for i = 1:size(XCa,3)
    XCa(:,:,i) = XCa(:,:,i) - mean(XCa(:,:,i),2);
    sl = diff(XCa(2,[1 11],i)) / diff(XCa(1,[1 11],i));
    Rz = rotz(atan2d(sl,1)); 
    Rz = Rz(1:2,1:2);
    XCa(:,:,i) = Rz'*XCa(:,:,i);
end

D = squeeze(sqrt(sum(diff(XC,1,2).^2))); D = D(1:2,:);
Da = squeeze(sqrt(sum(diff(XCa,1,2).^2))); Da = Da(1:11,:);
disp(max(fC));

% Find fixed point
[~,fpI] = min(abs(diff(D(:,1:1000),1,1)));
[~,fpL] = min(sum(abs(D(:,1001:end)-D([2,1],1))));
fpL = fpL + 1000;


%% a: iterated map
% Plot
pInd = 1;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
texta2 = '\textbf{a}\hspace{.2cm}conformational motion of a module';
text(labX,subp(pInd,4)+labY,texta2,NVTitle{:});

% Map
CP = winter(size(XC,3));
aX0 = 0.8; aXF = 2.8; aXD = aXF - aX0;
hold on;
plot([aX0 aXF],[aX0 aXF], '-', 'color', [1 1 1]*gr, 'linewidth',.5);
scatter(D(1,:),D(2,:),1,CP,'o','linewidth',.1);

% Axes
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
line([1 1]*D(1,1), [0 aXD/30]+aX0, 'color', 'k', 'linewidth', .5);
line([1 1]*D(1,fpI), [0 aXD/30]+aX0, 'color', 'k', 'linewidth', .5);
line([1 1]*D(1,fpL), [0 aXD/30]+aX0, 'color', 'k', 'linewidth', .5);
line([0 aXD/30]+aX0, [1 1]*D(1,1), 'color', 'k', 'linewidth', .5);
line([0 aXD/30]+aX0, [1 1]*D(1,fpI), 'color', 'k', 'linewidth', .5);
line([0 aXD/30]+aX0, [1 1]*D(1,fpL), 'color', 'k', 'linewidth', .5);
text(.1,.005,'$D_a$',NVTextR{:});
text(.41,.005,'$D^*$',NVTextR{:});
text(.61,.005,'$D_b$',NVTextR{:});
text(.895,0.045,'$d_1$',NVTextR{:});
text(.02,.92,'$d_2$',NVTextR{:});
axis([aX0 aXF aX0 aXF] + [-.05 .15 -.05 .15]*aXD);

% Examples
scc = .22;
plInd = [1 fpI fpL];
plSh = [ 0.40  0.35  0.20;...
         0.20  0.20  0.20];
for i = 1:length(plInd)
    pI = plInd(i);
    XCP = XC(:,:,pI)*scc + [D(1,pI);D(2,pI)]+plSh(:,i);
    plot(D(1,pI),D(2,pI),'s','color',CP(pI,:),'linewidth',4,'markersize',4);
    visualize_network(XCP(:,1:3),XCP(:,4:5),conn,'lcolor',CP(pI,:),'scolor',[1 1 1]*.7);
    line_coordinates(XCP(:,1:2), lSh, nW, lw, 'color', [1 1 1]*.5);
    line_coordinates(XCP(:,2:3), lSh, nW, lw, 'color', [1 1 1]*.5);
    if(i == 2)
        visualize_network(XCP(:,1:3)-.8,XCP(:,4:5)-.8,conn,'scolor',[1 1 1]*.7);
        line_coordinates(XCP(:,1:2)-.8, lSh, nW, lw, 'color', [1 1 1]*.5);
        line_coordinates(XCP(:,2:3)-.8, lSh, nW, lw, 'color', [1 1 1]*.5);
    end
end

% Lables
text(.07,.77,'$D_a$',NVTextR{:});
text(.31,.84,'$D_b$',NVTextR{:});
text(.445,.63,'$D^*$',NVTextR{:});
text(.68,.63,'$D^*$',NVTextR{:});
text(.65,.33,'$D_b$',NVTextR{:});
text(.87,.25,'$D_a$',NVTextR{:});
text(.11,.30,'$d_1$',NVTextR{:});
text(.35,.30,'$d_2$',NVTextR{:});
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% b: fixed point lattice
% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textb = '~\textbf{b}\hspace{.2cm}fixed point: $D^* ~\overrightarrow{{_f~}}~ D^*$';
text(labX,subp(pInd,4)+labY,textb,NVTitle{:});

% Draw networks
sc = .055;
Xfp1 = XC(:,:,fpI)*sc; Xfp2 = [Xfp1(1,:); -Xfp1(2,:) + Xfp1(2,2)];
sh = [.22; .83]; shx = [.25; 0];
line([Xfp1(1,[3,2]);Xfp2(1,[2,1])+shx(1)]+sh(1),...
     [Xfp1(2,[3,2]);Xfp2(2,[2,1])+shx(2)]+sh(2),...
     'linestyle',':','color',[1 1 1]*gr,'linewidth',1);
visualize_network(Xfp1(:,1:3)+sh,Xfp1(:,4:5)+sh,conn,'scolor',[1 1 1]*.7,'lcolor',CP(fpI,:));
visualize_network(Xfp2(:,1:3)+sh+shx,Xfp2(:,4:5)+sh+shx,conn,'scolor',[1 1 1]*.7,'lcolor',CP(fpI,:));
line_coordinates(Xfp1(:,1:2)+sh, lSh*.4, nW*.4, lw, 'color', [1 1 1]*.5);
line_coordinates(Xfp2(:,2:3)+sh+shx, -lSh*.4, nW*.4, lw, 'color', [1 1 1]*.5);
text(.13,.91,'$D^*$',NVTextR{:});
text(.45,.865,'$D^*$',NVTextR{:});
text(.75,.82,'$D^*$',NVTextR{:});

% Draw network lattice
sh = [.02; .51]; shx2 = [Xfp1(1,3);0];
for i = 1:5
    visualize_network(Xfp1(:,1:3)+sh+(2*i-1)*shx2,Xfp1(:,4:5)+sh+(2*i-1)*shx2,conn,'scolor',[1 1 1]*.7,'lcolor',CP(fpI,:));
    visualize_network(Xfp2(:,1:3)+sh+2*i*shx2,Xfp2(:,4:5)+sh+2*i*shx2,conn,'scolor',[1 1 1]*.7,'lcolor',CP(fpI,:));
end
text(.5,.67,'combine to form lattice',NVTextH{:});

% Draw network
sc = .09; sh = [.36;.12]; fpI1 = fpI+50; fpI2 = fpI-50;
visualize_network(XC(:,1:3,fpI1)*sc+sh,XC(:,4:5,fpI1)*sc+sh,conn,'nalpha',.4,'lalpha',.4,'lcolor',CP(fpI1,:),'scolor',[1 1 1]*.7);
visualize_network(XC(:,1:3,fpI2)*sc+sh,XC(:,4:5,fpI2)*sc+sh,conn,'lcolor',CP(fpI2,:),'scolor',[1 1 1]*.7);
line_coordinates(XC(:,1:2,fpI1)*sc+sh, lSh*.5, nW*.5, lw, 'color', [1 1 1]*.8);
line_coordinates(XC(:,1:2,fpI2)*sc+sh, lSh*.5, nW*.5, lw, 'color', [1 1 1]*.3);
line_coordinates(XC(:,2:3,fpI1)*sc+sh, lSh*.5, nW*.5, lw, 'color', [1 1 1]*.8);
line_coordinates(XC(:,2:3,fpI2)*sc+sh, lSh*.5, nW*.5, lw, 'color', [1 1 1]*.3);
text(.25,.24,'$\delta d_1$',NVTextR{:});
text(.64,.24,'$\delta d_2$',NVTextR{:});
text(.5,0,'$|\delta d_1| < |\delta d_2|$',NVTextH{:});
text(.5,.33,'unstable module',NVTextH{:});

% Axes
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% c: limit cycle lattice
% Plot
pInd = 3;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textc = '\textbf{c}\hspace{.2cm}2-cycle: $D_a ~\overrightarrow{{_f~}}~ D_b ~\overrightarrow{{_f~}}~ D_a$';
text(labX,subp(pInd,4)+labY,textc,NVTitle{:});

% Draw networks
sc = .055;
Xfp1 = XC(:,:,1)*sc; Xfp2 = [XC(1,:,fpL); -XC(2,:,fpL) + XC(2,2,fpL)]*sc;
sh = [.22; .84]; shx = [.25; 0];
line([Xfp1(1,[3,2]);Xfp2(1,[2,1])+shx(1)]+sh(1),...
     [Xfp1(2,[3,2]);Xfp2(2,[2,1])+shx(2)]+sh(2),...
     'linestyle',':','color',[1 1 1]*gr,'linewidth',1);
visualize_network(Xfp1(:,1:3)+sh,Xfp1(:,4:5)+sh,conn,'scolor',[1 1 1]*.7,'lcolor',CP(1,:));
visualize_network(Xfp2(:,1:3)+sh+shx,Xfp2(:,4:5)+sh+shx,conn,'scolor',[1 1 1]*.7,'lcolor',CP(fpL,:));
line_coordinates(Xfp1(:,1:2)+sh, lSh*.4, nW*.4, lw, 'color', [1 1 1]*.5);
line_coordinates(Xfp2(:,2:3)+sh+shx, -lSh*.4, nW*.4, lw, 'color', [1 1 1]*.5);
text(.05,.865,'$D_a$',NVTextR{:});
text(.81,.865,'$D_a$',NVTextR{:});

% Draw network lattice
sh = [.02; .52]; shx2 = [Xfp1(1,3);0];
for i = 1:5
    visualize_network(Xfp1(:,1:3)+sh+2*i*shx2,Xfp1(:,4:5)+sh+2*i*shx2,conn,'scolor',[1 1 1]*.7,'lcolor',CP(1,:));
    visualize_network(Xfp2(:,1:3)+sh+2*i*shx2,Xfp2(:,4:5)+sh+2*i*shx2,conn,'scolor',[1 1 1]*.7,'lcolor',CP(fpL,:));
end
text(.5,.67,'combine to form lattice',NVTextH{:});

% Draw network
sc = .07; sh = [.36;.14]; fpI1 = 1; fpI2 = 300;
visualize_network(XCc(:,1:4,fpI2)*sc+sh,XCc(:,5:8,fpI2)*sc+sh,connc,'nalpha',.4,'lalpha',.4,'lcolor',CP([fpI1 fpL],:),'scolor',[1 1 1]*.7);
visualize_network(XCc(:,1:4,fpI1)*sc+sh,XCc(:,5:8,fpI1)*sc+sh,connc,'lcolor',CP([fpI1 fpL],:),'scolor',[1 1 1]*.7);
line_coordinates(XCc(:,1:2,fpI1)*sc+sh, lSh*.5, nW*.5, lw, 'color', [1 1 1]*.8);
line_coordinates(XCc(:,1:2,fpI2)*sc+sh, lSh*.5, nW*.5, lw, 'color', [1 1 1]*.3);
line_coordinates(XCc(:,3:4,fpI1)*sc+sh, -lSh*.5, nW*.5, lw, 'color', [1 1 1]*.8);
line_coordinates(XCc(:,3:4,fpI2)*sc+sh, -lSh*.5, nW*.5, lw, 'color', [1 1 1]*.3);
text(.2,.18,'$\delta d_1$',NVTextR{:});
text(.7,.17,'$\delta d_3$',NVTextR{:});
text(.5,0,'$|\delta d_1| > |\delta d_3|$',NVTextH{:});
text(.5,.33,'stable modules',NVTextH{:});

% Axes
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% d: conformational motion
% Plot
pInd = 4;
subplot('position',subpN(pInd,:)); cla;
delete(findall(gca,'type','annotation'))
hold on;

% Title
textd = '\textbf{d}\hspace{.2cm}conformational motion';
text(labX,subp(pInd,4)+labY,textd,NVTitle{:});

% Draw networks
sc = .055;
plInds = [size(XCa,3) 300 1];
pSh = [[.345;.72] [.345;.45] [.345;.2]];
for i = 1:length(plInds)
    pI = plInds(i);
    [~,cInds] = min(abs(Da(:,pI)' - D(1,:)'));
    cInds = cInds(1:end-1);
    XCP = XCa(:,:,pI)*sc + pSh(:,i);
    visualize_network(XCP(:,1:12),XCP(:,13:end),conna,'lcolor',CP(cInds,:),'scolor',[1 1 1]*.7);
end

% Text
text(.5,.91,'unstable: $|\delta d_1| \ll |\delta d_k|$',NVTextH{:});
text(.5,.85,'motion begins at $d_k$',NVTextH{:});
text(.5,.06,'stable: $|\delta d_1| \gg |\delta d_k|$',NVTextH{:});
text(.5,.0,'motion begins at $d_1$',NVTextH{:});
text(0,.74,'$d_1$',NVTextR{:});
text(.89,.7,'$d_k$',NVTextR{:});

% Axes
ax = [0 sRat(pInd) 0 1];
axis(ax);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;

                        
%% Save
fName = 'figure3c';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');

