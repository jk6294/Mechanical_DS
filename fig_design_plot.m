% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
fig = figure(7); clf;
params_fig;
fig_design_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [17.8 7.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .0 .2 .0];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  4.00  3.50  3.50];...
        [ 3.50  4.00  3.50  3.50];...
        [ 7.00  4.00  3.50  3.50];...
        [ 0.00  0.00  3.50  3.50];...
        [ 3.50  0.00  3.50  3.50];...
        [ 7.00  0.00  3.50  3.50];...
        [10.50  0.00  7.50  7.50]];
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


%% a-c: unit, start, and end
% Parameters
sc = .23;
sh = [.4;.65];
lSh1 = .015;
lSh2 = .03;
nw = .015;
lw = .5;
CP20 = interp1(DLin2,CP2,sqrt(sum(diff(Xs(:,[1,3]),1,2).^2)));
CP2a = interp1(DLin2,CP2,sqrt(sum(diff(Xf(:,[1,3]),1,2).^2)));

pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
% Lengths
line_coordinates(Xs(:,1:2)*sc+sh,'lSh',lSh1,'color',C1a);
line_coordinates(Xs(:,2:3)*sc+sh,'lSh',lSh1,'color',C1a);
line_coordinates(Xs(:,[1 3])*sc+sh,'lSh',-lSh1,'color',CP20);
line_coordinates(Xf(:,1:2)*sc+sh,'lSh',lSh2,'color',C1c);
line_coordinates(Xf(:,2:3)*sc+sh,'lSh',lSh2,'color',C1c);
line_coordinates(Xf(:,[1 3])*sc+sh,'lSh',-lSh2,'color',CP2a);
% Draw network
visualize_network(Xf*sc+sh,[],[1 1],'ms',ms,'scolor',o);
visualize_network(Xs*sc+sh,[],[1 1],'ms',ms,'scolor',o*gr^3);
% Text
text(labX,subp(pInd,4)+labY,'\textbf{a}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'unit $k$',NVTitleH{:});
% Node nabels
for i = 1:3
    text(Xf(1,i)*sc+sh(1),Xf(2,i)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
    text(Xs(1,i)*sc+sh(1),Xs(2,i)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o);
end
% Node coordinates
text(Xs(1,1)*sc+sh(1)+.05,Xs(2,1)*sc+sh(2)-.09,'$(x_1^\bullet,y_1^\bullet)$',NVTexth{:},...
     'fontsize', FS2, 'color', o*gr^4);
text(Xs(1,2)*sc+sh(1)-.12,Xs(2,2)*sc+sh(2)+.09,'$(x_2^\bullet,y_2^\bullet)$',NVTexth{:},...
     'fontsize', FS2, 'color', o*gr^4);
text(Xs(1,3)*sc+sh(1)-.05,Xs(2,3)*sc+sh(2)-.09,'$(x_3^\bullet,y_3^\bullet)$',NVTexth{:},...
     'fontsize', FS2, 'color', o*gr^4);
text(Xf(1,1)*sc+sh(1)+.00,Xf(2,1)*sc+sh(2)-.09,'$(x_1^\circ,y_1^\circ)$',NVTexth{:},...
     'fontsize', FS2, 'color', o*gr^2);
text(Xf(1,2)*sc+sh(1)+.12,Xf(2,2)*sc+sh(2)+.09,'$(x_2^\circ,y_2^\circ)$',NVTexth{:},...
     'fontsize', FS2, 'color', o*gr^2);
text(Xf(1,3)*sc+sh(1)-.00,Xf(2,3)*sc+sh(2)-.09,'$(x_3^\circ,y_3^\circ)$',NVTexth{:},...
     'fontsize', FS2, 'color', o*gr^2);
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;
% Lengths
% line_coordinates(Xs(:,1:2)*sc+sh,'lSh',lSh1,'color',C1a);
% line_coordinates(Xs(:,2:3)*sc+sh,'lSh',lSh1,'color',C1a);
% Draw network
R = [-1 1; -2 1.3]*1.7;     % Solution space ranges
visualize_conic_finite(Xs*sc+sh,Xf*sc+sh,R*sc+sh,'ucolori',CP2a,...
                       'ucolorf',o,'overlay',.99,'lw',.5);
visualize_network(Xs*sc+sh,Xu*sc+sh,conn,'ucolor',CP2a,'ms',ms);
% Node labels
for i = 1:5
    if(i<=3)
        text(XC(1,i,1)*sc+sh(1),XC(2,i,1)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
             'fontsize', FS2, 'color', o*gr^4);
    else
        text(XC(1,i,1)*sc+sh(1),XC(2,i,1)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o);
    end
end
% Edge labels
shLx = [-1  10  8  -1  9   14]*.01;
shLy = [ 0  4  -14  2 -10 -4]*.01;
for i = 1:size(conn,1)
    text(mean(XC(1,conn(i,:),1))*sc+sh(1)+shLx(i),...
         mean(XC(2,conn(i,:),1))*sc+sh(2)+shLy(i),...
         sprintf('%0.1f',L(i)),NVTextr{:},'fontsize',8,'color',o*gr); 
end

% Text
text(labX,subp(pInd,4)+labY,'\textbf{b}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'add extra nodes',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,2)+labY,'fix rigid rod lengths',...
     NVTitleH{:},'color',o*gr);
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


pInd = 3;
subplot('position',subpN(pInd,:)); cla; hold on;
% line_coordinates(Xf(:,1:2)*sc+sh,'lSh',lSh1,'lw',1,'color',C1c);
% line_coordinates(Xf(:,2:3)*sc+sh,'lSh',lSh1,'lw',1,'color',C1c);
visualize_network(XC(:,1:3,end)*sc+sh,XC(:,4:5,end)*sc+sh,conn,'ucolor',CP2a,'ms',ms);
% Node labels
for i = 1:5
    if(i<=3)
        text(XC(1,i,end)*sc+sh(1),XC(2,i,end)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
             'fontsize', FS2, 'color', o*gr^4);
    else
        text(XC(1,i,end)*sc+sh(1),XC(2,i,end)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o);
    end
end
text(labX,subp(pInd,4)+labY,'\textbf{c}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'end position',NVTitleH{:});
% Edge labels
shLx = [-3 -3  5  5  10  6]*.01;
shLy = [ 4  4  8 -5 -15 -5]*.01;
for i = 1:size(conn,1)
    text(mean(XC(1,conn(i,:),end))*sc+sh(1)+shLx(i),...
         mean(XC(2,conn(i,:),end))*sc+sh(2)+shLy(i),...
         sprintf('%0.1f',L(i)),NVTextr{:},'fontsize',8,'color',o*gr); 
end
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% d-f: unit, start, and end
% Parameters
sc = .23;
sh = [.4;.65];
lSh1 = .015;
lSh2 = .03;
nw = .015;
lw = .5;
CP2b = interp1(DLin2,CP2,sqrt(sum(diff(Xf2(:,[1,3]),1,2).^2)));

pInd = 4;
subplot('position',subpN(pInd,:)); cla; hold on;
% Lengths
line_coordinates(Xs2(:,1:2)*sc+sh,'lSh',lSh1,'color',C1a);
line_coordinates(Xs2(:,2:3)*sc+sh,'lSh',lSh1,'color',C1a);
line_coordinates(Xs2(:,[1 3])*sc+sh,'lSh',-lSh1,'color',C2a);
line_coordinates(Xf2(:,1:2)*sc+sh,'lSh',lSh2,'color',C1c);
line_coordinates(Xf2(:,2:3)*sc+sh,'lSh',lSh2,'color',C1c);
line_coordinates(Xf2(:,[1 3])*sc+sh,'lSh',-lSh2,'color',C2b);
% Draw network
visualize_network(Xf2*sc+sh,[],[1 1],'ms',ms,'scolor',o);
visualize_network(Xs2*sc+sh,[],[1 1],'ms',ms,'scolor',o*gr^4);
% Text
text(labX,subp(pInd,4)+labY,'\textbf{a}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'unit $k$',NVTitleH{:});
% Node nabels
for i = 1:3
    text(Xf2(1,i)*sc+sh(1),Xf2(2,i)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^3);
    text(Xs2(1,i)*sc+sh(1),Xs2(2,i)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o);
end
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


pInd = 5;
subplot('position',subpN(pInd,:)); cla; hold on;
% Lengths
% line_coordinates(Xs2(:,1:2)*sc+sh,'lSh',lSh1,'color',C1a);
% line_coordinates(Xs2(:,2:3)*sc+sh,'lSh',lSh1,'color',C1a);
% Draw network
R = [-1 1; -2 1.3]*1.7;     % Solution space ranges
visualize_conic_finite(Xs2*sc+sh,Xf2*sc+sh,R*sc+sh,'ucolori',CP2b,...
                       'ucolorf',o,'overlay',.99,'lw',.5);
visualize_network(Xs2*sc+sh,Xu2*sc+sh,conn,'ucolor',CP2b,'ms',ms);
% Node labels
for i = 1:5
    if(i<=3)
        text(XC2(1,i,1)*sc+sh(1),XC2(2,i,1)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
             'fontsize', FS2, 'color', o*gr^4);
    else
        text(XC2(1,i,1)*sc+sh(1),XC2(2,i,1)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o);
    end
end
% Edge labels
shLx = [-8 -3 12  0 12 14]*.01;
shLy = [-7  1 -3 -4 -1 -4]*.01;
for i = 1:size(conn,1)
    text(mean(XC2(1,conn(i,:),1))*sc+sh(1)+shLx(i),...
         mean(XC2(2,conn(i,:),1))*sc+sh(2)+shLy(i),...
         sprintf('%0.1f',L2(i)),NVTextr{:},'fontsize',8,'color',o*gr); 
end

% Text
text(labX,subp(pInd,4)+labY,'\textbf{b}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'add extra nodes',NVTitleH{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,2)+labY,'fix rigid rod lengths',...
     NVTitleH{:},'color',o*gr);
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


pInd = 6;
subplot('position',subpN(pInd,:)); cla; hold on;
% line_coordinates(Xf2(:,1:2)*sc+sh,'lSh',lSh1,'lw',1,'color',C1c);
% line_coordinates(Xf2(:,2:3)*sc+sh,'lSh',lSh1,'lw',1,'color',C1c);
visualize_network(XC2(:,1:3,end)*sc+sh,XC2(:,4:5,end)*sc+sh,conn,'ucolor',CP2b,'ms',ms);
% Node labels
for i = 1:5
    if(i<=3)
        text(XC2(1,i,end)*sc+sh(1),XC2(2,i,end)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
             'fontsize', FS2, 'color', o*gr^4);
    else
        text(XC2(1,i,end)*sc+sh(1),XC2(2,i,end)*sc+sh(2)-.005, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o);
    end
end
% Edge labels
shLx = [11 11 12  9 -1  9]*.01;
shLy = [15  3 -3 -6 -1  4]*.01;
for i = 1:size(conn,1)
    text(mean(XC2(1,conn(i,:),end))*sc+sh(1)+shLx(i),...
         mean(XC2(2,conn(i,:),end))*sc+sh(2)+shLy(i),...
         sprintf('%0.1f',L2(i)),NVTextr{:},'fontsize',8,'color',o*gr); 
end
text(labX,subp(pInd,4)+labY,'\textbf{c}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'end position',NVTitleH{:});
% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% g: conformational motion
pInd = 7;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .11;
aX0 = 1.3;         % Map axes start
aXF = 3.15;         % Map axes end
% Colors
CPI = interp1(DLin,CP1,D(1,:)); 
CPJ = interp1(DLin,CP1,D(2,:));
CPI2 = interp1(DLin,CP1,D2(1,:)); 
CPJ2 = interp1(DLin,CP1,D2(2,:));

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
plot(D2(1,:),D2(2,:),'k-','linewidth',1);

% Draw spline
% sCo11 = [0.00 1.00; 0.20 0.90] .* [.13;-.6] + [2.45;2.8];
% plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
%             'linewidth',.5,'headwidth',5,'headlength',5);

% Draw example networks
plInd = ceil([1 26 47 70 size(XC,3)]);
plSh = [ .22  .22  .22 .22 .22;...
        -.09 -.07 -.00 .05 .05];
plSh2 = [-.18 -.22 -.16 -.14 -.16;...
          .19  .20  .34  .35  .35];
lSh = -.015;
for i = 1:length(plInd)
    pI = plInd(i);
    XP = XC(:,:,pI)*sc + [D(1,pI);D(2,pI)]+plSh(:,i);
    XP2 = XC2(:,:,pI)*sc + [D2(1,pI);D2(2,pI)]+plSh2(:,i);
    % Draw unit 1
    line_coordinates(XP(:,1:2), 'lSh',-lSh, 'lw',1, 'color',CPI(pI,:));
    line_coordinates(XP(:,2:3), 'lSh',-lSh, 'lw',1, 'color',CPJ(pI,:));
    visualize_network(XP(:,1:3),XP(:,4:5),conn,'ucolor',CP2a);
    % Draw unit 2
    line_coordinates(XP2(:,1:2), 'lSh',-lSh, 'lw',1, 'color',CPI2(pI,:));
    line_coordinates(XP2(:,2:3), 'lSh',-lSh, 'lw',1, 'color',CPJ2(pI,:));
    visualize_network(XP2(:,1:3),XP2(:,4:5),conn,'ucolor',CP2b);
end
scatter(D(1,plInd(2:end-1)),D(2,plInd(2:end-1)),40,zeros(length(plInd)-2,3),...
       'filled','marker','s','linewidth',.5);
scatter(D2(1,plInd(2:end-1)),D2(2,plInd(2:end-1)),40,zeros(length(plInd)-2,3),...
       'filled','marker','s','linewidth',.5);
scatter([Ls Le],[Ls Le],80,...
       [interp1(DLin,CP1,D(1,plInd(1)));interp1(DLin,CP1,D(1,plInd(end)))],...
       'filled','marker','s','linewidth',.5);
scatter([Ls Le],[Ls Le],40,ones(2,3),'marker','s','linewidth',.1);

% Draw large template network
% XP = XC(:,:,1)*sc*2 + [D(1,1);D(2,1)]+[0.2;1.1];
% visualize_network(XP(:,1:3),XP(:,4:5),conn,'lalpha',.3,'scolor',o*gr^.3,...
%                   'ucolor',CP2a);
% line_coordinates(XP(:,1:2), 'lSh',.07, 'lw',.5,'style','-','nw',.01);
% line_coordinates(XP(:,2:3), 'lSh',.07, 'lw',.5,'style','-','nw',.01);

% Text
text(labX,subp(pInd,4)+labY,'\textbf{g}',NVTitle{:});
text(subp(pInd,3)/2-fMarg(1),subp(pInd,4)+labY,'shape change',NVTitleH{:});
% % Axis legend
text(Ls,aX0-.08,['$l^\bullet=$' sprintf('%0.2f',Ls)],NVTexth{:},'color',interp1(DLin,CP1,Ls));
text(Le,aX0-.08,['$l^\circ=$' sprintf('%0.2f',Le)],NVTexth{:},'color',interp1(DLin,CP1,Le));
text(aX0-.05,Ls,'$l^\bullet$',NVTexth{:},'color',interp1(DLin,CP1,Ls));
text(aX0-.05,Le,'$l^\circ$',NVTexth{:},'color',interp1(DLin,CP1,Le));
% % axis labels
% text(.89,.03,'$l_k$',NVTextL{:});
% text(0,.9,'$l_{k+1}$',NVTextL{:});
% % template labels
% text(XP(1,2),XP(2,2)+.25,'unit $k$',NVTexth{:},'color',o*gr);
% text(.18,0.7,'$l_k$',NVTextL{:});
% text(.39,.7,'$l_{k+1} = f(l_k)$',NVTextL{:});
% text(XP(1,1)-.05,XP(2,1)-.05,'1',NVTexth{:},'color',o*gr);
% text(XP(1,2),XP(2,2)+.08,'2',NVTexth{:},'color',o*gr);
% text(XP(1,3)+.05,XP(2,3)-.05,'3',NVTexth{:},'color',o*gr);

% Axes
aMarg = .5;
axis([D(1,1)-aMarg D(1,end)+aMarg D(1,1)-aMarg D(1,end)+aMarg]);
set(gca,'visible',0,'xtick',[],'ytick',[],'box',0);
drawnow;


%% Size and Save Figure
fName = 'fig_design2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');