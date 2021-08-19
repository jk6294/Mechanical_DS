% Supplementary Figure: edge lengths
%% Prepare Space
clear; clc;
params_fig;
set(groot,'defaulttextinterpreter','latex');


%% Figure Dimensions
fig = figure(1); clf;
delete(findall(gcf,'type','annotation'));

% Figure Size in cm  [w,h]
fSize = [19 7.5];
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
load 4bar.mat
% Correct rotation
Ns = 114;
XCa = XCa - XCa(:,1,:);

% Plot 1: zoomed in subset
plRange = (1:12);
conn = [1 3; 2 3; 1 4; 2 4];
sc = 0.25; sh = [0.15;1.5];
for i = 1:length(plRange)
    plsInd = (1:4) + 2*(plRange(i) - 1);
    Xspl = XCa(:,plsInd,1)*sc+sh;
    
    % Draw units in start configuration
    visualize_network(Xspl,[],conn,'ms',12);
                  
    % Node labels
    for j = 1:3
        text(Xspl(1,j),Xspl(2,j)-.004,num2str(plsInd(j)),NVTexth{:},'fontsize',FS2);
    end
end

% Zoom in lines
plot([.05 0], [1.15 0.45], 'color', 'k', 'clipping', 0);
plot([1.4 4.55], [1.1 0.65], 'color', 'k', 'clipping', 0);

% Plot 3: zoomed in subset decomposed
sh = [0.15;0.4]; sh2 = [0.76;0.0];
shLx = [-0.10 -0.10 -0.05 -0.08;...
        -0.10 -0.10 -0.05 -0.08;...
        -0.10 -0.10 -0.05 -0.08;...
        -0.10 -0.10 -0.05 -0.08;...
        -0.10 -0.10 -0.05 -0.08;];
shLy = [-0.05 -0.05  0.05  0.05;...
         0.05  0.05 -0.05 -0.05;...
        -0.05 -0.05  0.05  0.05;...
         0.05  0.05 -0.05 -0.05;...
        -0.05 -0.05  0.05  0.05;];
       
for i = 1:5
    plsInd = (1:4) + 2*(plRange(i) - 1);
    Xspl = XCa(:,plsInd,1)*sc+sh + sh2*(i-1);
    
    % Draw units in start configuration
    visualize_network(Xspl,[],conn,'ms',12);
                  
    % Node labels
    for j = 1:4
        text(Xspl(1,j),Xspl(2,j)-.004,num2str(plsInd(j)),NVTexth{:},'fontsize',FS2);
    end
    
    % Edge labels
    D = sqrt(sum((Xspl(:,conn(:,1)) - Xspl(:,conn(:,2))).^2))/sc;
    mp = (Xspl(:,conn(:,1)) + Xspl(:,conn(:,2)))/2;
    mp = mp + [shLx(i,:); shLy(i,:)];
    for j = 1:4
        text(mp(1,j),mp(2,j),sprintf('%0.2f',D(j)),NVTexth{:},'fontsize',FS);
    end
end

% Text
text(.55,0.9,'decompose', NVTextR{:});

axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;



%% Save
fName = 'suppfig_decompose_4bar';
set(gcf, 'Renderer', 'painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');