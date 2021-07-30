%% Quadrifolium
load quadrifolium.mat;

% Plot
figure(1); clf;
visualize_network(Xscc,Xucc,conncc);
for i = 1:size(Xscc,2)
    text(Xscc(1,i),Xscc(2,i),num2str(i));
end
for i = 1:size(Xucc,2)
    text(Xucc(1,i),Xucc(2,i),num2str(i+size(Xscc,2)));
end

Rz = rotz(90);
Rz = Rz(1:2,1:2);

for i = 1:size(XCcc,3)
    XCcc(:,:,i) = Rz * XCcc(:,:,i);
end
XCcc(2,:,:) = -XCcc(2,:,:);
XCcc = XCcc - XCcc(:,1,:);


%% Track
N = size(XCcc,2);
aV = zeros(size(conncc,1),N,size(XCcc,3));
for i = 1:size(XCcc,3)
    disp(i);
% for i = 1:10
    XM = repmat(XCcc(:,conncc(:,1),i) - XCcc(:,conncc(:,2),i),[1 1 N]);
    XMP = [XM(2,:,:); -XM(1,:,:)];
    XD = reshape(XCcc(:,1:N,i),[2,1,N]) - XCcc(:,conncc(:,2),i);
    aV(:,:,i) = squeeze(real(acosd(sum(XD.*XM,1)./(sqrt(sum(XD.^2,1)).*sqrt(sum(XM.^2,1)))))) .* sign(squeeze(sum(XD.*XMP,1)));
end


%% Track2
% find zero crossings
zc1 = abs(sign(aV(:,:,2:end)) - sign(aV(:,:,1:end-1))) > 0;
% clf;
% plot(squeeze(aV(:,2,:))');



%% Track node movement through edges
figure(1); clf;
visualize_network(XCcc(:,:,1),[],conncc);



%% Prepare Space
clear; clc;



%%
% Parameters
sq = sqrt(3)/2;
L = 1;                                          % Initial length
LF = 1.7;                                       % Final length

conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
% Module: increase c
Xs1 = [-sq 0  sq;...
      -.5 1 -.5]*L;
Xs2 = [-sq   0.00   sq;...
       -0.5  1.00  -0.5]*L;
Xf2 = Xs2 * LF - [0;LF-L];
ds2 = sqrt(sum((Xf2(:,1)-Xf2(:,2)).^2));        % Final distance


% Parameters
a = 16.2;                                   % Quadrifolium radius
f = @(tho)[a*sin(2*tho) .* cos(tho);...     % Quadrifolium equation
           a*sin(2*tho) .* sin(tho)];

% Trace quadrifolium
nv = 1000;
thv = linspace(0,2.1*pi,nv) + pi/4;         % Sampled angles
xv = f(thv);                                % Sampled points

% Initialize 
xP = zeros(2,113);                          % Triangle tesselation points
xP(:,1) = f(thv(1))*1.084543;               % Initial triangle point
thop = thv(1);

% Tesselate quadrifolium with isoscelese triangles
options = optimset('TolFun',1e-15,'TolX',1e-15);
for i = 1:size(xP,2)
    [thop, fV] = fzero(@(t) (sqrt(sum((f(t)-xP(:,i)).^2))-ds2/2),...
                             thop+[.01,.2],options);
    xP(:,i+1) = xP(:,i) + 2*(f(thop)-xP(:,i));
end

% Bin triangle distances
Dc = sqrt(sum((xP(:,1:end-2)-xP(:,3:end)).^2));     % Triangle distances: c
[C,~,IC] = unique(round(Dc,3),'stable');            % Bin to 3rd decimal
Xuc = zeros(2,2,length(C));                         % Unique node positions

% Construct units for each unique bin
Xu0a = [-3.5*sq 3.5*sq; 2.2 2.2];                   % Initial conditions
Xu0b = [-sq+.1 -sq+.1; -2 2.5];
Xu0c = [-sq*2 sq*2;  -3.3 -3.3];
for i = 1:length(C)
    % Different initial conditions based on distance c
    if(C(i)<2.2)
        Xu0 = Xu0a;
    elseif(C(i)>4)
        Xu0 = Xu0c;
    else
        Xu0 = Xu0b;
    end
    Xfp = [[-C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)] [0;0]...
           [C(i)/2; -sqrt(ds2^2 - (C(i)/2)^2)]] + [0;Xf2(2,2)];
    Xup = construct_network(Xs1,Xfp,Xu0,conn,0,1);
    Xuc(:,:,i) = Xup(1:2,:);
end

% % Construct chain for quadrifolium network at the start position
[Xscc,Xucc,conncc,CSSc] = network_chain_x(Xs1(:,1:2),Xuc,IC');


%% Lengths
lm = sqrt(sum((Xscc(:,conncc(:,1)) - Xucc(:,conncc(:,2)-size(Xscc,2))).^2));

d_marg = .25;
% d_g = .08;
d_g = .12;
% sW = .17;
sW = .25;

colSh = .25;
rowSh = .25;
coW = 1/8.5;
coH = 1/11;

fig = figure(1); clf;
fSize = [8.5 11];

set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','inches');
set(gcf,'renderer','opengl','Position',[fig.Position(1:2) fSize],'Units','inches');
subplot('position',[0 0 1 1]);
axis([0 1 0 1]);
set(gca,'visible',0);

% Vertical lines
annotation('line',([1 1]*colSh)*coW, [d_marg 11-d_marg]*coH, 'linestyle', '-');
fInd = 1;
hold on;
for i = 1:size(conncc,1)
    annotation('line',([0 sW]+colSh)*coW, [1 1]*rowSh*coH,'linestyle','-');
    annotation('line',([0 sW]+colSh)*coW, ([1 1]*rowSh + 2*d_g + lm(i))*coH,'linestyle','-');
    
    plot(mean(([0 sW]+colSh)*coW),([1 1]*rowSh + d_g)*coH,'k.','markersize',7);
    plot(mean(([0 sW]+colSh)*coW),([1 1]*rowSh + d_g + lm(i))*coH,'k.','markersize',7);
    
%     rShv = .009;
    rShv = .014;
    text(colSh*coW+ rShv, (rowSh + d_g  + lm(i)/2)*coH, num2str(i),...
         'fontsize',8,'verticalalignment','middle','horizontalalignment','center','rotation',90);
    text(colSh*coW+ rShv, (rowSh + d_g/2  + lm(i))*coH, num2str(conncc(i,1)),...
         'fontsize',8,'verticalalignment','middle','horizontalalignment','right','rotation',90);
    text(colSh*coW+ rShv, (rowSh + 3*d_g/2)*coH, num2str(conncc(i,2)),...
         'fontsize',8,'verticalalignment','middle','horizontalalignment','left','rotation',90);
    
    rowSh = rowSh + 2*d_g + lm(i);
    if(rowSh + 2*d_g + lm(i+1) + d_marg > 11)
        rowSh = d_marg;
        colSh = colSh + sW;
        annotation('line',([1 1]*colSh)*coW, [d_marg 11-d_marg]*coH, 'linestyle', '-');
        annotation('line',([1 1]*(colSh+sW))*coW, [d_marg 11-d_marg]*coH, 'linestyle', '-');
        if(colSh > 8)
%             % Save
%             fName = ['draw_quad2_' num2str(fInd)];
%             set(gcf, 'Renderer', 'painters');
%             fig.PaperPositionMode = 'manual';
%             fig.PaperUnits = 'inches';
%             fig.PaperPosition = [0 0 fSize];
%             fig.PaperSize = fSize;
%             saveas(fig, ['Figures/' fName], 'pdf');
%             fInd = fInd + 1;
%             
%             % New figure
%             cla;
%             delete(findall(gcf,'type','annotation'));
%             colSh = d_marg;
%             annotation('line',([1 1]*colSh)*coW, [d_marg 11-d_marg]*coH, 'linestyle', '-');
%             set(gca,'visible',0);
        end
    end
    drawnow;
end
% Xu = Xuc(:,:,IC(1));
% l = sqrt(sum((Xs(:,conn(:,1)) - Xu(:,conn(:,2)-3)).^2));


%% Save
% fName = 'figureq5';
% set(gcf, 'Renderer', 'painters');
% fig.PaperPositionMode = 'manual';
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 fSize];
% fig.PaperSize = fSize;
% saveas(fig, ['Figures/' fName], 'pdf');