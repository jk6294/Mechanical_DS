function L = line_coordinates(X,lSh,nW,lw_d,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for computing and plotting coordinates for lines measuring 
% distances between nodes.
% 
% Inputs
% X:    2 x 2 Matrix of node pair coordinates, x top row, y bottom row
% lSh:  1 x 1 Scalar for shift distance from nodes
% nW:   1 x 1 Scalar for width of end nibs
% lw_d: 1 x 1 Scalar for linewidth
%
% Outputs
% L:    2 x 2 x 3 Matrix of line coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parse Inputs
p = inputParser;
addParameter(p,'color',[0 0 0]);
parse(p,varargin{:});

C = p.Results.color;


% Midpoint
pt = mean(X,2);
% Rotation
Rz = rotz(atan2d(diff(X(2,:)),diff(X(1,:))));
Rz = Rz(1:2,1:2);
% Length
l = sqrt(sum(diff(X,[],2).^2));
% Matrix of Rotated Line Coordinates
L = Rz*[[-l l -l -l l l]/2; ([0 0 -nW nW -nW nW]+lSh)];

% Plot
hold on;
line(L(1,1:2)+pt(1),L(2,1:2)+pt(2),'color',C,'LineWidth',lw_d);
line(L(1,3:4)+pt(1),L(2,3:4)+pt(2),'color',C,'LineWidth',lw_d);
line(L(1,5:6)+pt(1),L(2,5:6)+pt(2),'color',C,'LineWidth',lw_d);
