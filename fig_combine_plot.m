% Figure 2: Networks form crystals at attractors
%% Prepare Space
clear; clc;
fig = figure(4); clf;
params_fig;
fig_combine_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [8.6 7.1];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .4 .0];
% Subplot position in cm [x,y,w,h]
subp = [0.00 0.00  10.00 7.50];
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


%% odule combination
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;
delete(findall(gca,'type','annotation'));

% Parameters
sc = .07;          % Scale drawing
sho1 = [.00;.83];
sho2 = [.00;.53];
sho3 = [.00;.06];
sh1 = [.89;0];
sh2 = [.84;0];
sh3 = [.533;0];

% Parameters
plInd = 335;
XCa1 = XCa(:,:,plInd)*sc+sho1;
XCa2 = XCa(:,:,plInd)*sc+sho2;
XCa3 = XCa(:,:,plInd)*sc+sho3;
CPap = interp1(DLin,CP1,Da(:,plInd));

% Draw
lSh = -.007;
line_coordinates(XCa1(:,1:2), 'lSh',lSh, 'lw',1,'color',CPap(1,:));
line_coordinates(XCa1(:,3:4), 'lSh',-lSh, 'lw',1,'color',CPap(2,:));
line_coordinates(XCa1(:,3:4)+sh1, 'lSh',-lSh, 'lw',1,'color',CPap(2,:));
line_coordinates(XCa1(:,5:6)+sh1, 'lSh',lSh, 'lw',1,'color',CPap(3,:));
line_coordinates([[0;0] sh1]+XCa1(:,3), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
line_coordinates([[0;0] sh1]+XCa1(:,4), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
visualize_network(XCa1(:,1:4),[],conn);
visualize_network(XCa1(:,3:6)+sh1,[],conn);
% arrow([.15 -.15]+.84, [1 1]*.81, sRat(pInd),'color',o*gr);

% Transition 1
sCo11 = [0.02 -0.01 0.02; 0.72 0.66 0.6];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
        
% Second units
% Draw
line_coordinates(XCa2(:,1:2), 'lSh',lSh, 'lw',1,'color',CPap(1,:));
line_coordinates(XCa2(:,3:4), 'lSh',-lSh, 'lw',1,'color',CPap(2,:));
line_coordinates(XCa2(:,5:6), 'lSh',-lSh, 'lw',1,'color',CPap(3,:));
line_coordinates(XCa2(:,5:6)+sh2, 'lSh',lSh, 'lw',1,'color',CPap(3,:));
line_coordinates(XCa2(:,7:8)+sh2, 'lSh',lSh, 'lw',1,'color',CPap(4,:));
line_coordinates([[0;0] sh2]+XCa2(:,5), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
line_coordinates([[0;0] sh2]+XCa2(:,6), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
visualize_network(XCa2(:,1:6),[],conna(1:8,:));
visualize_network(XCa2(:,5:8)+sh2,[],conn);
% arrow([.15 -.15]+.84, [1 1]*.445, sRat(pInd),'color',o*gr);

% Transition 2
sCo11 = [0.02 -0.01 0.02; [0.75 0.7 0.65]-.3];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
sCo11 = [0.02 -0.01 0.02; [0.75 0.7 0.65]-.55];
plot_spline(sCo11,'head',1,'headpos',1,'ratio',sRat(pInd),...
            'linewidth',.5,'headwidth',5,'headlength',5,'color',o*gr);
        
% Third units
for i = 1:8
    line_coordinates(XCa3(:,[1 2]+2*(i-1)), 'lSh',-(-1)^i*lSh, 'lw',1,...
                     'color',CPap(i,:));
end
line_coordinates(XCa3(:,15:16)+sh3, 'lSh',-lSh, 'lw',1,'color',CPap(8,:));
line_coordinates(XCa3(:,17:18)+sh3, 'lSh',lSh, 'lw',1,'color',CPap(9,:));
line_coordinates([[0;0] sh3]+XCa3(:,15), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
line_coordinates([[0;0] sh3]+XCa3(:,16), 'lSh',0, 'lw',.5,'color',o*gr,'style','--');
% arrow([.15 -.15]+.84, [1 1]*.05, sRat(pInd),'color',o*gr);
visualize_network(XCa3(:,1:16),[],conna(1:28,:));
visualize_network(XCa3(:,15:18)+sh3,[],conn);


% Text
% Label 1 unit 1
tsh = .05;
text(XCa1(1,1),XCa1(2,1)+.09,'unit 1',NVTextlv{:});
text(XCa1(1,1),XCa1(2,1)-tsh-.02,'$l_1$',NVTextlv{:},'color',CPap(1,:));
text(XCa1(1,4),XCa1(2,4)-tsh,'$f(l_1)$',NVTextlv{:},'color',CPap(2,:));
text(XCa1(1,3)+.035,XCa1(2,3)+.015,'3',NVTexth{:},'color',o*gr,'fontsize',8);
text(XCa1(1,4)+.015,XCa1(2,4)+.025,'4',NVTexth{:},'color',o*gr,'fontsize',8);
% Label 1 unit 2
text(XCa1(1,3)+sh1(1),XCa1(2,1)+.09,'unit 2',NVTextlv{:});
text(XCa1(1,3)+sh1(1)-.03,XCa1(2,4)-tsh,'$l_2$',NVTextlv{:},'color',CPap(2,:));
text(XCa1(1,6)+sh1(1)-.01,XCa1(2,6)+tsh+.02,'$f(l_2)$',NVTextlv{:},'color',CPap(3,:));
text(XCa1(1,3)+sh1(1),XCa1(2,3)-0.035,'1',NVTexth{:},'color',o*gr,'fontsize',8);
text(XCa1(1,4)+0.035+sh1(1),XCa1(2,4),'2',NVTexth{:},'color',o*gr,'fontsize',8);
% Label 1 combine nodes
text(.35,XCa1(2,4)-tsh,'set $f(l_1) = l_2$',NVTextL{:},'color',CPap(2,:));
text(.02,.67,'join nodes marked by dashed lines',NVTextL{:},'color',o*gr);

% Label 2 unit 2
text(XCa2(1,3),XCa2(2,1)+.09,'unit 2',NVTextlv{:});
text(XCa2(1,6)+.01,XCa2(2,6)+tsh-.01,'$f(l_2)$',NVTextlv{:},'color',CPap(3,:));
text(XCa2(1,5)+.015,XCa2(2,5)-.035,'3',NVTexth{:},'color',o*gr,'fontsize',8);
text(XCa2(1,6)+.015,XCa2(2,6)-.035,'4',NVTexth{:},'color',o*gr,'fontsize',8);
% Label 2 unit 3
text(XCa2(1,5)+sh2(1),XCa2(2,1)+.09,'unit 3',NVTextlv{:});
text(XCa2(1,5)+sh2(1)-.02,XCa2(2,6)+tsh-.01,'$l_3$',NVTextlv{:},'color',CPap(3,:));
text(XCa2(1,8)+sh2(1)-.02,XCa2(2,8)-tsh-.02,'$f(l_3)$',NVTextlv{:},'color',CPap(4,:));
text(XCa2(1,5)+sh2(1),XCa2(2,5)+0.035,'1',NVTexth{:},'color',o*gr,'fontsize',8);
text(XCa2(1,6)+0.035+sh2(1),XCa2(2,6),'2',NVTexth{:},'color',o*gr,'fontsize',8);
% Label 2 combine nodes
text(.40,XCa2(2,6)+tsh-.01,'set $f(l_2) = l_3$',NVTextL{:},'color',CPap(3,:));
text(0,.275,'$\vdots$',NVTextH{:},'color',o*gr);
text(.02,.275,'join nodes marked by dashed lines',NVTextL{:},'color',o*gr);

% Label 3 unit k
text(XCa3(1,14),XCa3(2,1)+.1,'unit $k$',NVTexthv{:});
text(XCa3(1,15),XCa3(2,15)+.06,'3',NVTexth{:},'color',o*gr,'fontsize',8);
text(XCa3(1,16)+.025,XCa3(2,16)+.015,'4',NVTexth{:},'color',o*gr,'fontsize',8);
% Label 3 unit k+1
text(XCa3(1,14)+sh3(1)-.05,XCa3(2,1)+.1,'unit $k+1$',NVTextlv{:});
text(XCa3(1,15)+sh3(1),XCa3(2,15)-0.035,'1',NVTexth{:},'color',o*gr,'fontsize',8);
text(XCa3(1,16)+0.035+sh3(1),XCa3(2,16),'2',NVTexth{:},'color',o*gr,'fontsize',8);
% Label 2 combine nodes
text(.45,.07,'set $f(l_k) = l_{k+1}$',NVTextL{:},'color',CPap(8,:));

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Size and Save Figure
fName = 'fig_combinea';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');