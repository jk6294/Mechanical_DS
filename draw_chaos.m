%% Prepare Space
clear; clc;


%% Visualize
figure(2); clf; 
load chaos_network.mat
% XCc5(2,:,:) = -XCc5(2,:,:);
nInt = 1;
nI = 10850;
visualize_network(XCc5(:,:,(nI-100)*nInt),[],connc5,'lcolor',[.8 .8 .8]);
visualize_network(XCc5(:,:,nI*nInt),[],connc5);
axis([0,15,-2,2.5]);


%% Construct Network
% Construct
L = 1;
th = 45;
Rz = rotz(th); Rz = Rz(1:2,1:2);
sh2 = -Rz*[.012;0];
nVal = linspace(-1.6,-0.6,4);
d1dot = -Rz*[1;0];
Xs = [-L*cosd(th)  0  L*cosd(th);...
      -L*sind(th)  0 -L*sind(th)];

% Added node positions
s = .5/sind(th);
Xupd = @(vV) [0 vV; -2*s -L*sind(th)-vV];

% Period double
nValpd = .214659;
nValpdc = .217;

[Xscc,Xucc,conncc] = network_chain_x(Xs(:,1:2),Xupd(nValpdc),ones(1,30));


%% Chaos
% Plot
figure(2); clf;
visualize_network(Xscc,Xucc,conncc);
for i = 1:size(Xscc,2)
    text(Xscc(1,i),Xscc(2,i),num2str(i));
end
for i = 1:size(Xucc,2)
    text(Xucc(1,i),Xucc(2,i),num2str(i+size(Xscc,2)));
end


%% Lengths
lm = 2.5*sqrt(sum((Xscc(:,conncc(:,1)) - Xucc(:,conncc(:,2)-size(Xscc,2))).^2));

d_marg = .25;
d_g = .08;
sW = .17;
% d_g = .12;
% sW = .25;

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
    
    rShv = .009;
%     rShv = .014;
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
            % Save
            fName = ['draw_chaos_' num2str(fInd)];
            set(gcf, 'Renderer', 'painters');
            fig.PaperPositionMode = 'manual';
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 fSize];
            fig.PaperSize = fSize;
            saveas(fig, ['Figures/' fName], 'pdf');
            fInd = fInd + 1;
            
            % New figure
            cla;
            delete(findall(gcf,'type','annotation'));
            colSh = d_marg;
            annotation('line',([1 1]*colSh)*coW, [d_marg 11-d_marg]*coH, 'linestyle', '-');
            set(gca,'visible',0);
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