% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
fig = figure(8); clf;
params_fig;
fig_map_design_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [17.8 7.5];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00 0.00  9.50  7.50];...
        [11.00 0.00  7.50  7.50]];
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
sc = .055;          % Scale drawing

sho1 = [.00;.76];
sho2 = [.00;.42];
sho3 = [.00;.02];
sh1 = [1.025;0];
sh2 = [0.96;0];
sh3 = [.775;0];
lSh = -.000;


% Parameters
plInd = 30;
XCa1 = XCa(:,:,plInd)*sc+sho1;
XCa2 = XCa(:,:,plInd)*sc+sho2;
XCa3 = XCa(:,:,plInd)*sc+sho3;
CPap = interp1(DLin,CP1,Da(:,plInd));
CP2a = interp1(DLin2,CP2,sqrt(sum(diff(Xf(:,[1,3]),1,2).^2)));


% Row 1
nPl = 1;
for i = 1:nPl+1
    line_coordinates(XCa1(:,[0 1]+i), 'lSh',-lSh, 'lw',1,...
                     'color',CPap(i,:));
end
line_coordinates(XCa1(:,(0:1)+nPl+1)+sh1, 'lSh',lSh, 'lw',1,'color',CPap(nPl+1,:));
line_coordinates(XCa1(:,(0:1)+nPl+2)+sh1, 'lSh',lSh, 'lw',1,'color',CPap(nPl+2,:));
line_coordinates([[0;0] sh1]+XCa1(:,nPl+1), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
line_coordinates([[0;0] sh1]+XCa1(:,nPl+2), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
% arrow([.15 -.15]+.66, [1 1]*.785, sRat(pInd),'color',o*gr);
visualize_network(XCa1(:,1:nPl+2),XCa1(:,[1:2*nPl]+ns),...
                  conna(1:nPl*6,:)-[0 ns-nPl-2],'ucolor',CP2a,'ms',ms);
visualize_network(XCa1(:,[1:3]+nPl)+sh1,XCa1(:,[1:2]+ns+nPl*2)+sh1,...
                  conn,'ucolor',CP2a,'ms',ms);
% Node labels
for i = 1:3
    text(XCa1(1,i,1),XCa1(2,i,1)-.003, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
    text(XCa1(1,i+1,1)+sh1(1),XCa1(2,i+1,1)+sh1(2)-.003, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
end
% Transition 1
sCo11 = [0.02 -0.00 0.02; 0.73 0.65 0.55];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
% Unit 1
yv = mean(XCa1(2,2:3));
text(XCa1(1,2),XCa1(2,1)+.19,'unit 1',NVTexth{:});
text(XCa1(1,1)-.02,yv,'$l_1$',NVTextlv{:},'color',CPap(1,:));
text(XCa1(1,3)-.01,yv,'$f(l_1)$',NVTextlv{:},'color',CPap(2,:));
% Unit 2
text(XCa1(1,3)+sh1(1),XCa1(2,1)+.19,'unit 2',NVTexth{:});
text(XCa1(1,2)+sh1(1)-.03,yv,'$l_2$',NVTextlv{:},'color',CPap(2,:));
text(XCa1(1,4)+sh1(1)-.02,yv-.02,'$f(l_2)$',NVTextlv{:},'color',CPap(3,:));
% Combine
text(.65,yv,'set $f(l_1) = l_2$',NVTexth{:},'color',CPap(2,:));
text(.65,.74,'join nodes 3 and 2',NVTexth{:},'color',o*gr);
text(.65,.925,'join nodes 2 and 1',NVTexth{:},'color',o*gr);


% Row 2
nPl = 2;
for i = 1:nPl+1
    line_coordinates(XCa2(:,[0 1]+i), 'lSh',-lSh, 'lw',1,...
                     'color',CPap(i,:));
end
line_coordinates(XCa2(:,(0:1)+nPl+1)+sh2, 'lSh',-lSh, 'lw',1,'color',CPap(nPl+1,:));
line_coordinates(XCa2(:,(0:1)+nPl+2)+sh2, 'lSh',-lSh, 'lw',1,'color',CPap(nPl+2,:));
line_coordinates([[0;0] sh2]+XCa2(:,nPl+1), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
line_coordinates([[0;0] sh2]+XCa2(:,nPl+2), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
% arrow([.15 -.15]+.66, [1 1]*.44, sRat(pInd),'color',o*gr);
visualize_network(XCa2(:,1:nPl+2),XCa2(:,[1:2*nPl]+ns),conna(1:nPl*6,:)-[0 ns-nPl-2],...
                  'scolor',o*gr^.2,'ucolor',CP2a.^.1,'lalpha',.1,'ms',ms);
visualize_network(XCa2(:,[1:3]+nPl-1),XCa2(:,[1:2]+ns+nPl*2-2),conn,'ucolor',CP2a,'ms',ms);
visualize_network(XCa2(:,[1:3]+nPl)+sh2,XCa2(:,[1:2]+ns+nPl*2)+sh2,conn,'ucolor',CP2a,'ms',ms);
% Node labels
for i = 1:3
    text(XCa2(1,i+1,1),XCa2(2,i+1,1)-.003, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
    text(XCa2(1,i+2,1)+sh2(1),XCa2(2,i+2,1)+sh2(2)-.003, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
end
% Transition 2
sCo11 = [0.02 0.00 0.02; [0.74 0.7 0.66]-.35];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
sCo11 = [0.02 0.00 0.02; [0.74 0.7 0.66]-.51];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
% Unit 2
yv = mean(XCa2(2,3:4));
text(XCa2(1,3),XCa2(2,1)+.19,'unit 2',NVTexth{:});
text(XCa2(1,4)-.02,yv,'$f(l_2)$',NVTextlv{:},'color',CPap(3,:));
% Unit 3
text(XCa2(1,4)+sh2(1),XCa2(2,1)+.19,'unit 3',NVTexth{:});
text(XCa2(1,3)+sh2(1)-.05,yv,'$l_3$',NVTextlv{:},'color',CPap(3,:));
text(XCa2(1,5)+sh2(1)+.00,yv+.02,'$f(l_3)$',NVTextlv{:},'color',CPap(4,:));
% Combine
text(.68,yv,'set $f(l_2) = l_3$',NVTexth{:},'color',CPap(3,:));
text(.68,.562,'join nodes 3 and 2',NVTexth{:},'color',o*gr);
text(.68,.4,'join nodes 2 and 1',NVTexth{:},'color',o*gr);


% Row 3
nPl = 6;
for i = 1:nPl+1
    line_coordinates(XCa3(:,[0 1]+i), 'lSh',-lSh, 'lw',1,...
                     'color',CPap(i,:));
end
line_coordinates(XCa3(:,(0:1)+nPl+1)+sh3, 'lSh',-lSh, 'lw',1,'color',CPap(nPl+1,:));
line_coordinates(XCa3(:,(0:1)+nPl+2)+sh3, 'lSh',-lSh, 'lw',1,'color',CPap(nPl+2,:));
line_coordinates([[0;0] sh3]+XCa3(:,nPl+1), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
line_coordinates([[0;0] sh3]+XCa3(:,nPl+2), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
% arrow([.15 -.15]+.66, [1 1]*.05, sRat(pInd),'color',o*gr);
visualize_network(XCa3(:,1:nPl+2),XCa3(:,[1:2*nPl]+ns),conna(1:nPl*6,:)-[0 ns-nPl-2],...
                  'scolor',o*gr^.2,'ucolor',CP2a.^.1,'lalpha',.1,'ms',ms);
visualize_network(XCa3(:,[1:3]+nPl-1),XCa3(:,[1:2]+ns+nPl*2-2),conn,'ucolor',CP2a,'ms',ms);
visualize_network(XCa3(:,[1:3]+nPl)+sh3,XCa3(:,[1:2]+ns+nPl*2)+sh3,conn,'ucolor',CP2a,'ms',ms);
% Node labels
for i = 1:3
    text(XCa3(1,i+nPl-1,1),XCa3(2,i+nPl-1,1)-.003, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
    text(XCa3(1,i+nPl,1)+sh3(1),XCa3(2,i+nPl,1)+sh3(2)-.003, num2str(i),NVTexth{:},...
         'fontsize', FS2, 'color', o*gr^4);
end
% Unit k
yv = mean(XCa3(2,7:8))-.005;
text(XCa3(1,nPl+1),XCa3(2,nPl)+.09,'unit $k$',NVTexth{:});
text(XCa3(1,8)-.02,yv,'$f(l_k)$',NVTextlv{:},'color',CPap(nPl+1,:));
% Unit k+1
text(XCa3(1,nPl+2)+sh3(1),XCa3(2,nPl)+.09,'unit $k+1$',NVTexth{:});
text(XCa3(1,nPl+1)+sh3(1)-.1,yv,'$l_{k+1}$',NVTextlv{:},'color',CPap(nPl+1,:));
text(XCa3(1,nPl+3)+sh3(1)-.01,yv+.02,'$f(l_{k+1})$',NVTextlv{:},'color',CPap(nPl+1,:));
% Combine
text(.01,.27,'$\vdots$',NVTexth{:},'color',o*gr);
text(.79,yv,'set $f(l_k) = l_{k+1}$',NVTexth{:},'color',CPap(nPl+2,:));
text(.79,.16,'join nodes 3 and 2',NVTexth{:},'color',o*gr);
text(.79,.02,'join nodes 2 and 1',NVTexth{:},'color',o*gr);


% Title
text(labX,subp(pInd,4)+labY,'\textbf{a}',NVTitle{:});
text(subp(pInd,3)/2 + subp(pInd,2),subp(pInd,4)+labY,'combine units by joining nodes',NVTitleH{:});


% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% b: single module conformational motion
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .1;
aX0 = 1.3;         % Map axes start
aXF = 3.15;         % Map axes end
sh1 = [1.4;2.2];
sh2 = [1.5;2.6];
sh3 = [1.6;3.0];
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
mInter = 2.93;
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
visualize_network(XCapa(:,1:ns),XCapa(:,ns+1:end),conna,...
                  'scolor',o*gr.^.1,'ucolor',CP2a.^.1,'lalpha',.1);
% Network 2
for i = 1:size(Da,1)
    line_coordinates(XCapb(:,(0:1)+i), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPb(i,:));
end
visualize_network(XCapb(:,1:ns),XCapb(:,ns+1:end),conna,...
                  'scolor',o*gr.^.1,'ucolor',CP2a.^.1,'lalpha',.1);
% Network 3
for i = 1:size(Da,1)
    line_coordinates(XCapc(:,(0:1)+i), 'lSh',(-1)^i*lSh, 'lw',1,'color',CPc(i,:));
end
visualize_network(XCapc(:,1:ns),XCapc(:,ns+1:end),conna,...
                  'scolor',o*gr.^.1,'ucolor',CP2a.^.1,'lalpha',.1);


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
fName = 'fig_map_design2';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');