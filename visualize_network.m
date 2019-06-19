function [] = visualize_network(Xs, Xu, conn, sC, C_SN, C_UN, Tr)
% This function is to visualize a constructed network.
%
% Inputs
% Xs:       d x n  matrix of specified node positions
% Xu:       d x m  matrix of unspecified node positions
% conn:     k x 2  matrix of edges from node i to node j
% sC:       1 x 1  optional scalar for scaling network
% C_SN:     n x 4  optional vector of colors
% C_UN:     m x 4  optional vector of colors
% Tr:       1 x 1  optional scalar for transparency: 0 to 1
%
% Outputs
% Figure. Assumes the desired figure placement is already selected

if(nargin == 3)
    sC = 1;
    C_SN = repmat([255 100 100]/255,size(Xs,2),1); 
    C_UN = repmat([100 100 255]/255,size(Xu,2),1);
    Tr = 1;
elseif(nargin == 4)
    C_SN = repmat([255 100 100]/255,size(Xs,2),1); 
    C_UN = repmat([100 100 255]/255,size(Xu,2),1);
    Tr = 1;
elseif(nargin == 5)
    if(size(C_SN,1)==0)
        C_SN = [255 100 100]/255; 
    end
    C_SN = repmat(C_SN,size(Xs,2)/size(C_SN,1),1);
    C_UN = repmat([100 100 255]/255,size(Xu,2),1);
    Tr = 1;
elseif(nargin >= 6)
    if(size(C_SN,1)==0)
        C_SN = [255 100 100]/255; 
    end
    if(size(C_UN,1)==0)
        C_UN = [100 100 255]/255;
    end
    C_SN = repmat(C_SN,size(Xs,2)/size(C_SN,1),1);
    C_UN = repmat(C_UN,size(Xu,2)/size(C_UN,1),1);
    if(nargin == 6)
        Tr = 1;
    end
end
if(size(C_SN,2) == 3)
    C_SN = (C_SN-1)*Tr + 1;
end
if(size(C_UN,2) == 3)
    C_UN = (C_UN-1)*Tr + 1;
end


% Plot Parameters
% C_UN = [100 100 255]/255;         % Color of Unspecified Node
d = size(Xs,1);                     % Number of Dimensions
ms = 3*sC;                          % Marker Size
lw = 1.2*sC;                        % Line Width
bw = 0.5*sC;                        % Length of boarder
eC = [0 0 0 .5*Tr];

hold on
X = [Xs Xu];
if(d==2)
    % Edges
    line([X(1,conn(:,1)); X(1,conn(:,2))],...
         [X(2,conn(:,1)); X(2,conn(:,2))],...
         'linewidth', lw, 'color', eC);
    % Specified Nodes
    for i = 1:size(Xs,2)
        plot(Xs(1,i), Xs(2,i), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN(i,:))
        plot(Xs(1,i), Xs(2,i), 'ko', 'linewidth', bw, 'markersize', ms*2, 'color', [1 1 1]*(1-Tr))
    end
    % Unspecified Nodes
    if(size(Xu,2) > 0)
        for i = 1:size(Xu,2)
            plot(Xu(1,i), Xu(2,i), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN(i,:));
            plot(Xu(1,i), Xu(2,i), 'ko', 'linewidth', bw, 'markersize', 2*ms, 'color', [1 1 1]*(1-Tr));
        end
    end
    set(gca,'visible',0);
    set(gcf,'color','w');
elseif(d==3)
    % Spherical point
    [xSp, ySp, zSp] = sphere(20);
    xSp = xSp/10*sC; 
    ySp = ySp/10*sC; 
    zSp = zSp/10*sC; 
    % Edges
    line([X(1,conn(:,1)); X(1,conn(:,2))],...
         [X(2,conn(:,1)); X(2,conn(:,2))],...
         [X(3,conn(:,1)); X(3,conn(:,2))],...
         'linewidth', lw, 'color', eC);
    % Unspecified Nodes
    for i = 1:size(Xu,2)
        s = surf(xSp+Xu(1,i), ySp+Xu(2,i), zSp+Xu(3,i));
        s.FaceColor = C_UN(i,:);
        s.EdgeColor = 'none';
    end
    % Specified Nodes
    for i = 1:size(Xs,2)
        s = surf(xSp+Xs(1,i), ySp+Xs(2,i), zSp+Xs(3,i));
        s.FaceColor = C_SN(i,:);
        s.EdgeColor = 'none';
    end
    hold off;
    set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
            'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
end
hold off;
end