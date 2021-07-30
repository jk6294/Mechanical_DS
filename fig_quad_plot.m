% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
fig = figure(7); clf;
params_fig;
fig_map_design_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [17.8 7.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 0.00  7.50  7.50];...
        [ 7.50 0.00  7.50  7.50];...
        [15.00 0.00  3.00  7.50]];
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


%% a: module combination
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));


% Parameters
nSS = 15;                       % Number of solutoin spaces to sample
sc = .1;                       % Scale networks
sh = [0.3;0.85];               % Shift networks
thR = linspace(15,42.2683,nSS); % Angles of solution spaces to draw
R2 = [-1 1; -1 1]*3.2;          % Range of solution space to draw

% Draw solution space
for i = 1:nSS
    Xfp = sq*2*Ls*[[-sind(thR(i));-cosd(thR(i))] [0;0],...
               [sind(thR(i));-cosd(thR(i))]] + [0;Xf2(2,2)];
    visualize_conic_finite(Xs*sc+sh(:,1),Xfp*sc+sh(:,1),R2*sc+sh(:,1),...
                           'ucolori',CP2(ceil(i/nSS*size(CP2,1)),:),...
                           'ucolorf',o,'overlay',.99);
end

% Draw start and end node positions
for i = 1:nSS
    Xfp = ds2*[[-sind(thR(i));-cosd(thR(i))] [0;0],...
               [sind(thR(i));-cosd(thR(i))]];
    visualize_network(Xs*sc+sh(:,1),[],[1 1],'msize',5);
    visualize_network(Xfp*sc+sh(:,1),[],[1 1],'scolor',o*gr.^.3,'msize',3,...
                      'bcolor',interp1(DLin2,CP2,diff(Xfp(1,[1 3]))));
end


% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',1,'xtick',[],'ytick',[]);
drawnow;


%% b: single module conformational motion
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .1;
aX0 = 1.3;         % Map axes start
aXF = 3.15;         % Map axes end
sh1 = [1.4;2.2];
sh2 = [1.7;2.6];
sh3 = [2.0;3.0];
% Colors
CPI = interp1(DLin,CP1,D(1,:));
CPJ = interp1(DLin,CP1,D(2,:));
aMarg = .5;

% Draw map ticks
line([1 1]*Ls, [0 .06]+aX0, 'color', interp1(DLin,CP1,Ls),'linewidth', 1);
line([1 1]*Le, [0 .06]+aX0, 'color', interp1(DLin,CP1,Le),'linewidth', 1);
line([0 .06]+aX0, [1 1]*Ls, 'color', interp1(DLin,CP1,Ls),'linewidth', 1);
line([0 .06]+aX0, [1 1]*Le, 'color', interp1(DLin,CP1,Le),'linewidth', 1);

% Draw distance colormap
scatter(D(1,:),ones(1,size(D,2))*aX0+.015,15,interp1(DLin,CP1,D(1,:)),'filled','marker','s');
scatter(ones(1,size(D,2))*aX0+.015,D(2,:),15,interp1(DLin,CP1,D(2,:)),'filled','marker','s');

% Draw map axes
line([aX0 aX0 aXF],[aXF aX0 aX0],'color','k','linewidth',.5);
arrow([aX0 aXF], [aX0 aX0], sRat(pInd));
arrow([aX0 aX0], [aX0 aXF], sRat(pInd));
    
% Draw map and y=x lines
plot([aX0 aXF],[aX0 aXF], '-', 'color', o*gr, 'linewidth',.5);
plot(D(1,:),D(2,:),'k-','linewidth',1);

% Cobweb
mInter = 2.8;
[~,mInd] = min(abs(Da(1,:) - mInter));
CPa = interp1(DLin,CP1,Da(:,1));
CPb = interp1(DLin,CP1,Da(:,mInd));
CPc = interp1(DLin,CP1,Da(:,end));
cobweb(Da(:,1),aX0,'mdecay',0,'msize',2,'lcolor',CPa,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',1);
cobweb(Da(:,mInd),aX0,'mdecay',0,'msize',2,'lcolor',CPb,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',1);
cobweb(Da(:,end),aX0,'mdecay',0,'msize',2,'lcolor',CPc,'ncolor',o*0,...
       'arrowind',[1 4],'mtype','o','linestyle','--','linewidth',1);
   
% Scale and shift networks for drawing
XCapa = XCa(:,:,1) * sc + sh1;
XCapb = XCa(:,:,mInd) * sc + sh2;
XCapc = XCa(:,:,end) * sc + sh3;


% Network 1
for i = 1:size(Da,1)
    line_coordinates(XCapa(:,(0:1)+i), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPa(i,:));
end
visualize_network(XCapa(:,1:ns),XCapa(:,ns+1:end),conna,'ucolor',CP2a);
% Network 2
for i = 1:size(Da,1)
    line_coordinates(XCapb(:,(0:1)+i), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPb(i,:));
end
visualize_network(XCapb(:,1:ns),XCapb(:,ns+1:end),conna,'ucolor',CP2a);
% Network 3
for i = 1:size(Da,1)
    line_coordinates(XCapc(:,(0:1)+i), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPc(i,:));
end
visualize_network(XCapc(:,1:ns),XCapc(:,ns+1:end),conna,'ucolor',CP2a);


% Text
text(labX,subp(pInd,4)+labY,'\textbf{b}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'shape change',NVTitleH{:});
% % Axis legend
text(Ls,aX0-.08,['$l^\bullet=$' sprintf('%0.2f',Ls)],NVTexth{:},'color',interp1(DLin,CP1,Ls));
text(Le,aX0-.08,['$l^\circ=$' sprintf('%0.2f',Le)],NVTexth{:},'color',interp1(DLin,CP1,Le));
text(aX0-.05,Ls,'$l^\bullet$',NVTexth{:},'color',interp1(DLin,CP1,Ls));
text(aX0-.05,Le,'$l^\circ$',NVTexth{:},'color',interp1(DLin,CP1,Le));
% % axis labels
text(.89,.03,'$l_k$',NVTextL{:});
text(0,.9,'$l_{k+1}$',NVTextL{:});

% Axes
axis([D(1,1)-aMarg D(1,end)+aMarg D(1,1)-aMarg D(1,end)+aMarg]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% Size and Save Figure
fName = 'fig_quada';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');