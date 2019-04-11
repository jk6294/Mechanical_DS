function [] = visualize_network(Xs, Xu, conn, sC)
% This function is to visualize a constructed network.
%
% Inputs
% Xs:       d x n  matrix of specified node positions
% Xu:       d x m  matrix of unspecified node positions
% conn:     k x 2  matrix of edges from node i to node j
% sC:       1 x 1  optional scalar for scaling network
%
% Outputs
% Figure. Assumes the desired figure placement is already selected

if(nargin == 3)
    sC = 1;
end

% Plot Parameters
C_SN = [255 100 100]/255;       % Color of Specified Node   
C_UN = [100 100 255]/255;       % Color of Unspecified Node
d = size(Xs,1);                 % Number of Dimensions
ms = 4*sC;                      % Marker Size
lw = 2*sC;                      % Line Width
ea = .5;                        % Edge Transparency

hold on
X = [Xs Xu];
if(d==2)
    % Edges
    line([X(1,conn(:,1)); X(1,conn(:,2))],...
         [X(2,conn(:,1)); X(2,conn(:,2))],...
         'linewidth', lw, 'color', [0 0 0 ea]);
    % Specified Nodes
    plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
    % Unspecified Nodes
    if(size(Xu,2) > 0)
        plot(Xu(1,:), Xu(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN(1,:));
    end
    set(gca,'visible',0);
    set(gcf,'color','w');
elseif(d==3)
    % Spherical point
    [xSp, ySp, zSp] = sphere(20);
    xSp = xSp/10; 
    ySp = ySp/10; 
    zSp = zSp/10; 
    % Edges
    line([X(1,conn(:,1)); X(1,conn(:,2))],...
         [X(2,conn(:,1)); X(2,conn(:,2))],...
         [X(3,conn(:,1)); X(3,conn(:,2))],...
         'linewidth', lw, 'color', [0 0 0 ea]);
    % Unspecified Nodes
    for i = 1:size(Xu,2)
        s = surf(xSp+Xu(1,i), ySp+Xu(2,i), zSp+Xu(3,i));
        s.FaceColor = C_UN(1,:);
        s.EdgeColor = 'none';
    end
    % Specified Nodes
    for i = 1:size(Xs,2)
        s = surf(xSp+Xs(1,i), ySp+Xs(2,i), zSp+Xs(3,i));
        s.FaceColor = C_SN;
        s.EdgeColor = 'none';
    end
    hold off;
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
            'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
end
hold off;
end