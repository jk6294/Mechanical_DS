% Supplementary Figure: edge lengths
%% Prepare Space
clear; clc;
params_fig;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 10.1];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  0.00  7.25  4.50]];
% Adjust Position
subp = subp + [fMarg(1) fMarg(3) -sum(fMarg(1:2)) -sum(fMarg(3:4))];
sRat = subp(:,3) ./ subp(:,4);
% Normalize Position
subpN = subp ./ [fSize(1) fSize(2) fSize(1) fSize(2)];
% Label Position in cm from top
labX = -fMarg(1);
labY = fMarg(4)-.18;
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],...
    'Units','centimeters');


%% Quadrifolium
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;

% Load
load superstable.mat
C3a = [197 066 085]/255;

% Correct rotation
Ns = 12;
XCa = XCa - XCa(:,1,:);
XCa(1,:,:) = -XCa(1,:,:);
for i = 1:size(XCa,3)
    Rz = rotz(atan2d(diff(XCa(2,[-2 0]+Ns,i)),diff(XCa(1,[-2 0]+Ns,i))));
    XCa(:,:,i) = Rz(1:2,1:2)'*XCa(:,:,i);
end

% Plot 1: zoomed in subset
plRange = 1:(Ns-2);
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
sc = 0.25; sh = [0.15;2.1];
for i = 1:length(plRange)
    plsInd = (1:3) + plRange(i) - 1;
    pluInd = (1:2) + Ns + 2*(plRange(i)-1);
    Xspl = XCa(:,plsInd,1)*sc+sh;
    Xupl = XCa(:,pluInd,1)*sc+sh;
    
    % Draw units in start configuration
    visualize_network(Xspl,Xupl,conn,'ucolor',C3a,'ms',12);
                  
    % Node labels
    for j = 1:3
        text(Xspl(1,j),Xspl(2,j)-.004,num2str(plsInd(j)),NVTexth{:},'fontsize',FS2);
    end
    for j = 1:2
        text(Xupl(1,j),Xupl(2,j)-.004,num2str(pluInd(j)),NVTexth{:},'fontsize',FS2);
    end
end

% Zoom in lines
plot([0.03 -0.03], [1.4 0.95], 'color', 'k', 'clipping', 0);
plot([1.6 4.55], [1.4 0.95], 'color', 'k', 'clipping', 0);

% Plot 2: zoomed in subset decomposed
sh = [0.15;0.65]; sh2 = [0.76;0.0];
shLx = [-0.10  0.07  0.15 -0.10  0.07  0.05;...
        -0.10  0.07  0.15 -0.10  0.07  0.05;...
        -0.10  0.07  0.15 -0.10  0.07  0.05;...
        -0.10  0.07  0.15 -0.10  0.07  0.05;...
        -0.10  0.07  0.15 -0.10  0.07  0.05];
shLy = [ 0.00 -0.07  0.03  0.00  0.07  0.05;...
         0.00  0.07 -0.03  0.00 -0.07 -0.05;...
         0.00 -0.07  0.03  0.00  0.07  0.05;...
         0.00  0.07 -0.03  0.00 -0.07 -0.05;...
         0.00 -0.07  0.03  0.00  0.07  0.05];
       
for i = 1:5
    plsInd = (1:3) + plRange(i) - 1;
    pluInd = (1:2) + Ns + 2*(plRange(i)-1);
    Xspl = XCa(:,plsInd,1)*sc+sh + sh2*(i-1);
    Xupl = XCa(:,pluInd,1)*sc+sh + sh2*(i-1);
    
    % Draw units in start configuration
    visualize_network(Xspl,Xupl,conn,'ucolor',C3a,'ms',12);
                  
    % Node labels
    for j = 1:3
        text(Xspl(1,j),Xspl(2,j)-.004,num2str(plsInd(j)),NVTexth{:},'fontsize',FS2);
    end
    for j = 1:2
        text(Xupl(1,j),Xupl(2,j)-.004,num2str(pluInd(j)),NVTexth{:},'fontsize',FS2);
    end
    
    % Edge labels
    D = sqrt(sum((Xspl(:,conn(:,1)) - Xupl(:,conn(:,2)-3)).^2))/sc;
    mp = (Xspl(:,conn(:,1)) + Xupl(:,conn(:,2)-3))/2;
    mp = mp + [shLx(i,:); shLy(i,:)];
    for j = 1:6
        text(mp(1,j),mp(2,j),sprintf('%0.2f',D(j)),NVTexth{:},'fontsize',FS);
    end
end

% Text
text(.3,2.95,'zoom in', NVTextR{:});
text(.6,1.25,'decompose', NVTextR{:});

axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;



%% Save
fName = 'suppfig_decompose_superstable';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');