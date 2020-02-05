function [] = arrow(X,Y,ratio,varargin)
% Function for drawing arrows
%
% Inputs
% 
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p,'HeadLength',3);
addParameter(p,'HeadWidth',3);
addParameter(p,'color',[0 0 0]);
addParameter(p,'LineWidth',.5);
parse(p,varargin{:});

HL = p.Results.HeadLength;
HW = p.Results.HeadWidth;
CL = p.Results.color;
LW = p.Results.LineWidth;

% Main line
line(X,Y,'color',CL,'linewidth',LW);
% Perpendicular line
sl = diff(Y)/diff(X);
r = rotz(atan2d(-diff(X),diff(Y))); r = r(1:2,1:2);
x = r*[-HW HW; 0 0]/72 + [X(1);Y(1)];
line(x(1,:),x(2,:),'color',CL,'linewidth',LW);
ah = annotation('arrow','HeadLength',HL,'HeadWidth',HW,...
                'color',CL,'linewidth',LW);
set(ah,'parent',gca,'position',[X(2),Y(2),diff(X)/100,diff(Y)/100*ratio]);
end