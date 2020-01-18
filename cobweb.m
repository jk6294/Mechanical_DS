function [] = cobweb(d,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots a cobweb plot from a set of coordinates
%
% Inputs
% d:        n x 1 vector of cobweb plot coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parse Inputs
p = inputParser;
addParameter(p,'lcolor',[0 0 0]);
addParameter(p,'ncolor',[0 0 0]);
addParameter(p,'msize',3);
addParameter(p,'mdecay',.5);
addParameter(p,'nalpha',1);
addParameter(p,'lalpha',1);
addParameter(p,'linewidth',.8);
addParameter(p,'headlength',3);
addParameter(p,'headwidth',4);
addParameter(p,'headscale',.8);
addParameter(p,'arrowind',[1 1]);
parse(p,varargin{:});

C_L = p.Results.lcolor;
C_N = p.Results.ncolor;
ms = p.Results.msize;
md = p.Results.mdecay;
nTr = p.Results.nalpha;
lTr = p.Results.lalpha;
lw = p.Results.linewidth;
hl = p.Results.headlength;
hw = p.Results.headwidth;
hs = p.Results.headscale;
aI = p.Results.arrowind;

% Replicate colors if size of color list is less than number of nodes
C_N = repmat(C_N,ceil((size(d,1)-1)/size(C_N,1)),1);

% Artificially implement transparency
C_N = (C_N-1)*nTr + 1;
C_L = (C_L-1)*lTr + 1;


%% Draw Arrows
hold on;
dL = [d d]'; dL = dL(:);
line(dL(2:end-2),[dL(3:end-1)],'color',C_L,'linewidth',lw);

for i = (2*aI(1)):(2*aI(2)-1)
    ah = annotation('arrow','HeadLength',hl,'HeadWidth',hw,'color',C_L);
    set(ah,'parent',gca);
    set(ah,'position',[dL(i) dL(i+1) diff(dL(i:i+1))*hs diff(dL(i+1:i+2))*hs]);
end

for i = 1:length(d)-1
    plot(d(i),d(i+1),'o','markersize',ms-(i-1)*md,...
         'linewidth',ms-(i-1)*md,'color',C_N(i,:));
end


end