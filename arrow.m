function [] = arrow(X,Y,ratio,varargin)
% Function for drawing arrows
%
% Inputs
% X
% Y
% ratio
%
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p,'HeadLength',3);
addParameter(p,'HeadWidth',3);
addParameter(p,'HeadStyle','vback2');
addParameter(p,'color',[0 0 0]);
addParameter(p,'LineWidth',.5);
addParameter(p,'LineStyle','-');
addParameter(p,'lp',1);
parse(p,varargin{:});

HL = p.Results.HeadLength;
HW = p.Results.HeadWidth;
HS = p.Results.HeadStyle;
CL = p.Results.color;
LW = p.Results.LineWidth;
LS = p.Results.LineStyle;
lp = p.Results.lp;

% Plot arrow as annotation
line([0 diff(X)*lp]+X(1),[0 diff(Y)*lp]+Y(1),'color',CL,'linewidth',LW);
ah = annotation('arrow','HeadLength',HL,'HeadWidth',HW,'HeadStyle',HS,...
                'color',CL,'LineStyle','none');
set(ah,'parent',gca,'position',[X(2),Y(2),diff(X)/100,diff(Y)/100*ratio]);
end