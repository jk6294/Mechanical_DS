% Figure 1: Motivation and Conformational Motions
%% Prepare Space
clear; clc;
fig = figure(9); clf;
params_fig;
suppfig_disorder_quad_code;


%% Figure Dimensions
% Figure Size in cm  [w,h]
fSize = [17.8 10.0];
% Margins in cm, [l,r,d,u]
fMarg = [.2 .2 .2 .2];
% Subplot position in cm [x,y,w,h]
subp = [[ 0.00  7.00 17.80  3.00];...
        [ 0.00  0.00 17.80  7.00]];
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


%% Compare node placements
pInd = 1;
subplot('position',subpN(pInd,:)); cla; hold on;

sc = .2;              % Scale drawing
sh = [7;0.15];       % Shift drawing

Xucca = Xucca - Xscca(:,Ns);
Xscca = Xscca - Xscca(:,Ns);

visualize_network(XCcc(:,:,1)*sc + sh,[],[1 1]);
visualize_network(Xscca*sc+sh,Xucca*sc+sh,[1 1]);

% Axes
axis([0 sRat(pInd) 0 1]);
set(gca,'visible',0,'xtick',[],'ytick',[]);
drawnow;


%% Distribution of bond lengths
% Center trace
xPc = (xP(:,2:end) + xP(:,1:end-1))/2;
xPc = xPc - mean(xPc,2);

% Iterate
XCccc = zeros(2,Ns-1,size(XCcc,3));
for i = 1:size(XCccc,3)
    xprox = (XCcc(:,2:Ns,i) + XCcc(:,1:Ns-1,i))/2;
    xprox = xprox - mean(xprox,2);
    [u,s,v] = svd(xPc * xprox');
    Rp = u*v'; 
    thp = atan2d(Rp(2,1),Rp(2,2));
    Rz = rotz(thp); Rz = Rz(1:2,1:2);
    XCccc(:,:,i) = Rz * xprox;
end

dAll = squeeze(sum(sqrt(sum((XCccc - xPc).^2))));
[~,sInd] = min(dAll);
sInd = sInd + 0;
clf; hold on;
plot(xPc(1,:),xPc(2,:));
plot(XCccc(1,:,sInd),XCccc(2,:,sInd));
hold off;





BLena = sqrt(sum((Xscca(:,conncc(:,1)) - Xucca(:,conncc(:,2)-Ns)).^2));
BLenb = sqrt(sum((Xscc(:,conncc(:,1)) - Xucc(:,conncc(:,2)-Ns)).^2));

DAll = squeeze((diff(XCcc(:,1:Ns,:),1,2)));


%% Plot
pInd = 2;
subplot('position',subpN(pInd,:)); cla; hold on;

% Parameters
sc = .026;              % Scale drawing
sh = [2.57;0.15];       % Shift drawing

% Draw networks
pI =  [001 050 085 130 180 590];
pIL = [200 200 350 350 350 size(conncc,1)];
for i = 1:length(pI)
    la = .1 + .1*i;
    visualize_network(XCcc(:,unique(conncc(1:pIL(i),1)),pI(i))*sc+sh,...
                      XCcc(:,unique(conncc(1:pIL(i),2)),pI(i))*sc+sh,...
                      conncc(1:pIL(i),:) + [0 max(conncc(1:pIL(i)))-size(Xscc,2)],...
                      'lalpha',la,'msize',3,'ucolor',CSSc(1:length(unique(conncc(1:pIL(i),2))),:).^(.05*i));
end

% Start conformation
% visualize_network(XCcc(:,1:size(Xscc,2),1)*sc+sh,...
%                   XCcc(:,(1:size(Xucc,2))+size(Xscc,2),1)*sc+sh,...
%                   conncc,'msize',2,'ucolor',CSSc);

% End conformation
sh = [2.5; 0.5];
Rz = rotz(47.5); Rz = Rz(1:2,1:2);
Rz = eye(2);
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
fName = 'suppfig_disorder_quad';
set(gcf, 'Renderer', 'painters'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 fSize];
fig.PaperSize = fSize;
saveas(fig, ['Figures/' fName], 'pdf');
set(gcf, 'Renderer', 'opengl');