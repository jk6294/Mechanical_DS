% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
fig = figure(9); clf;
params_fig;
suppfig_energy_quad_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [8.6 16.6];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 8.00 8.60 8.60];...
        [ 0.00 2.00 8.60 6.00];...
        [ 0.00 0.00 8.60 2.00]];
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


%% a: quadrifolium tesselation
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .035;              % Scale trace
sh = [.5; .48];         % Shift Networks
plRange = (1:28)+20;

% Draw quadrifolium trace
plot(xv(1,:)*sc+sh(1),xv(2,:)*sc+sh(2),'-','linewidth',.5,'color',[0 0 0]);
% Draw all tesselated isoscelese triangles
for i = 1:size(xP,2)-2
    line_coordinates(xP(:,[0 2]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin2,CP2,Dc(i)),'style','-');
    line_coordinates(xP(:,[0 1]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin,CP1,ds2),'style','-');
end
% Fill in a subset of triangles
for i = plRange(1:end-2)
    fill(xP(1,[0:2 0]+i)*sc+sh(1),xP(2,[0:2 0]+i)*sc+sh(2),...
         interp1(DLin2,CP2,Dc(i)),'linewidth',.5);
end
for i = plRange
    line_coordinates(xP(:,[0 2]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin2,CP2,Dc(i)),'style','-');
    line_coordinates(xP(:,[0 1]+i)*sc+sh,'lw',.5,'lSh',0,...
        'color',interp1(DLin,CP1,ds2),'style','-');
end
% Draw nodes
visualize_network(xP(:,plRange)*sc+sh,[],[1 1]);

% Text
text(labX,subp(pInd,4)+labY,'\textbf{a}',NVTitle{:});
text(subp(pInd,3)/2,subp(pInd,4)+labY,'tesselate curve',NVTitleH{:});

% Axis
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% b: Draw units
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc2 = .08;             % Scale example networks
sh1a = [.14;.52];
sh1b = [.14;.92];
sh1c = [.0;.38];
sha = [0.4;0.0];
Rz = [1 0; 0 -1];

% Draw units
for i = 1:4
    plI1 = IC(plRange(i));
    % Orient and shift units
    Xspr = Rz^(i-1)*[Xs1 Xuc(:,:,plI1)];
    Xspr = Xspr - [Xspr(1,2); Xspr(2,1+mod(i,2))];
    Xfpr = Rz^(i-1)*[Xfp(:,:,plI1) Xucf(:,:,plI1)];
    Xfpr = Xfpr - [Xfpr(1,2); Xfpr(2,1+mod(i,2))];
    Xspr = Xspr*sc2+sh1a+sha*(i-1); 
    Xfpr = Xfpr*sc2+sh1b+sha*(i-1);
    
    % Connection between units
    line_coordinates([[0;0] sha-[ds1*sc2/2;0]]+Xspr(:,2), 'lSh',0,...
                     'lw',.5,'color',o*gr,'style','--');
    line_coordinates([[0;0] sha-[ds1*sc2/2;0]]+Xspr(:,3), 'lSh',0,...
                     'lw',.5,'color',o*gr,'style','--');
    if(i<4)
        text(sha(1)/2.0+Xspr(1,2), Xspr(2,2)-(-1)^i*.04,'join',...
             NVTexth{:},'color',o*gr);
        text(sha(1)/2.0+Xspr(1,2), Xspr(2,2)-(-1)^i*.035-.07,'2 and 1',...
             NVTexth{:},'color',o*gr);
        text(sha(1)/2.0+Xspr(1,2), Xspr(2,3)+(-1)^i*.04,'join',...
             NVTexth{:},'color',o*gr);
        text(sha(1)/2.0+Xspr(1,2), Xspr(2,3)+(-1)^i*.035-.07,'3 and 2',...
             NVTexth{:},'color',o*gr);
    end

    % Draw units in end configuration
    fill(Xfpr(1,[1 2 3 1]),Xfpr(2,[1 2 3 1]),...
         interp1(DLin2,CP2,C(plI1)),'linewidth',.01,'edgecolor','w');
    line_coordinates(Xfpr(:,[1 2]),'color',interp1(DLin,CP1,ds2),'lSh',0);
    line_coordinates(Xfpr(:,[2 3]),'color',interp1(DLin,CP1,ds2),'lSh',0);
    line_coordinates(Xfpr(:,[1 3]),'color',interp1(DLin2,CP2,C(plI1)),'lSh',0);
    visualize_network(Xfpr(:,1:3),Xfpr(:,4:5),conn,...
                      'ucolor',interp1(DLin2,CP2,ds2),'ms',ms);

    % Draw units in start configuration
    visualize_network(Xspr(:,1:3),Xspr(:,4:5),conn,...
                      'ucolor',interp1(DLin2,CP2,ds2),'ms',ms);
                  
    % Node labels
    for j = 1:3
        text(Xspr(1,j),Xspr(2,j)-.004,num2str(j),NVTexth{:},'fontsize',FS2);
        text(Xfpr(1,j),Xfpr(2,j)-.004,num2str(j),NVTexth{:},'fontsize',FS2);
    end
end

% Draw splines
sC1 = [0.86 0.60 0.20; 1.29 1.13 0.97];
sC2 = [0.95 0.58; 1.27 0.97];
sC3 = [0.97 0.93 0.93; 1.19 1.10 0.97];
sC4 = [1.06 1.01 1.20; 1.20 1.06 0.97];
plot_spline(sC1,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
plot_spline(sC2,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
plot_spline(sC3,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
plot_spline(sC4,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);

% Axis
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0);
drawnow;


%% c: Quadrifolium folding sequence
pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;


% Parameters
sc = .05;              % Scale drawing
sh = [4.9;0.15];       % Shift drawing

% Draw networks
pI =  [001 050 085 130 180 590];
pIL = [200 200 350 350 350 size(conncc,1)];
for i = 1:length(pI)
    la = .1 + .1*i;
    visualize_network(XCcc(:,unique(conncc(1:pIL(i),1)),pI(i))*sc+sh,...
                      XCcc(:,unique(conncc(1:pIL(i),2)),pI(i))*sc+sh,...
                      conncc(1:pIL(i),:) + [0 max(conncc(1:pIL(i)))-size(Xscc,2)],...
                      'lalpha',la,'msize',2,'ucolor',CSSc(1:length(unique(conncc(1:pIL(i),2))),:).^(.05*i));
end

% Start conformation
% visualize_network(XCcc(:,1:size(Xscc,2),1)*sc+sh,...
%                   XCcc(:,(1:size(Xucc,2))+size(Xscc,2),1)*sc+sh,...
%                   conncc,'msize',2,'ucolor',CSSc);

% End conformation
sh = [5.0; 0.5];
Rz = rotz(47.5); Rz = Rz(1:2,1:2);
visualize_network(Rz*XCcc(:,1:size(Xscc,2),end)*sc+sh,...
                  Rz*XCcc(:,(1:size(Xucc,2))+size(Xscc,2),end)*sc+sh,...
                  conncc,'msize',2,'ucolor',CSSc);

% Text
% textc = '...to design folding sequences?';
% text(labX,subp(pInd,4)+labY-.1,'\textbf{c}',NVTitle{:});
% text(fSize(1)-fMarg(1),subp(pInd,4)+labY-.1,textc,NVTitleR{:});

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'suppfig_energy_quad';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');