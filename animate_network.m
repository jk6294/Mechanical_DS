function [] = animate_network(XC,conn,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to animate networks
%
% Inputs:
% XC:       d x dn x T  matrix of simulated position data
% conn:     2 x k       matrix of node connections
%
% Name/Value Pairs:
% scale     1 x 1       size: network node and edge     [1]
% nalpha    1 x 1       transparency: node              [1]
% lalpha    1 x 1       transparency: edge              [1]
% scolor    1 x 3       color: specified node           [255 100 100]/255
% ucolor    1 x 3       color: unspecified node         [100 100 255]/255
% bcolor    1 x 3       color: node boarder             [0 0 0]
% lcolor    1(k) x 3    color: edge                     [0 0 0]
% lwidth    1 x 1       width: edge                     [1.5]
% bwidth    1 x 1       width: node boarder             [0.6]
% figwidth  1 x 1       width: figure (cm)              [10]
% msize     1 x 1       size: nodes                     [4]
% nu        1 x 1       number: unspecified nodes       [0]
% nframe    1 x 1       number: frames in animation     [100]
% save      1 x 1       1 to save, 0 to not save        [0]
% fname     string      file name without ending        ['animation']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parse Inputs
p = inputParser;
addParameter(p,'scale',1);
addParameter(p,'nalpha',1);
addParameter(p,'lalpha',1);
addParameter(p,'scolor',[255 100 100]/255);
addParameter(p,'ucolor',[100 100 255]/255);
addParameter(p,'bcolor',[0 0 0]);
addParameter(p,'lcolor',repmat([0 0 0],size(conn,1),1));
addParameter(p,'lwidth',1.5);
addParameter(p,'bwidth',.6);
addParameter(p,'figwidth',10);
addParameter(p,'msize',4);
addParameter(p,'nu',0);
addParameter(p,'nframe',100);
addParameter(p,'save',0);
addParameter(p,'fname','animation');
addParameter(p,'distance',[]);
parse(p,varargin{:});

sC = p.Results.scale;
nTr = p.Results.nalpha;
lTr = p.Results.lalpha;
C_SN = p.Results.scolor;
C_UN = p.Results.ucolor;
C_L = p.Results.lcolor;
C_B = p.Results.bcolor;
lw = p.Results.lwidth*sC;
ms = p.Results.msize*sC;
bw = p.Results.bwidth*sC;
nu = p.Results.nu;
nF = p.Results.nframe;
fw = p.Results.figwidth;
s = p.Results.save;
fname = p.Results.fname;
d = p.Results.distance;


fig = gcf; clf;
nSV = 1;
fig.Units = 'centimeters';
fName = [fname '.gif'];
dT = 0.03;
pInds = floor(linspace(1,size(XC,3),nF));
ns = size(XC,2)-nu;
a = [min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))];

if(size(d,2) > 0)
    fig.Position(3:4) = [diff(a(1:2))+diff(a(3:4)) diff(a(3:4))] /...
                        (diff(a(1:2))+diff(a(3:4))) * fw;
    sPos1 = [0 0 fig.Position(4)/ fig.Position(3) 1];
    sPos2 = [fig.Position(4)/ fig.Position(3) 0,...
             -diff(fig.Position(3:4))/fig.Position(3) 1];
    dc = squeeze(sqrt(sum(diff(XC,1,2).^2)));
    dc = dc([1:ns-1],:);
else
    fig.Position(3:4) = [diff(a(1:2)) diff(a(3:4))] / diff(a(1:2)) * fw;
    sPos2 = [0 0 1 1];
end

for i = pInds
    if(size(d,2) > 0)
        subplot('Position', sPos1);
        cla;
        plot(d(1,:), d(2,:));
        mind = min(d(:)); 
        maxd = max(d(:));
        diffd = maxd - mind;
        axis([mind-diffd/10 maxd+diffd/10 mind-diffd/10 maxd+diffd/10]);
        hold on;
        plot([min(d(:)) max(d(:))], [min(d(:)) max(d(:))], '--', 'color', [1 1 1]*.7);
        cobweb(dc(:,i),'mdecay',0);
        hold off;
    end
    subplot('Position', sPos2, 'Units', 'normalized');
    cla;
    visualize_network(XC(:,1:ns,i),XC(:,[1:nu]+ns,i),conn,...
                      'scale',sC,'nalpha',nTr,'lalpha',lTr,...
                      'scolor',C_SN,'ucolor',C_UN,'bcolor',C_B,...
                      'lcolor',C_L,'lwidth',lw,'bwidth',bw,'msize',ms);
    axis(a + [[-1 1]*diff(a(1:2)) [-1 1]*diff(a(3:4))]/10);
    set(gca,'visible',0);
    drawnow;

    % Capture the plot as an image
    if(s == 1)
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if nSV == 1
          imwrite(imind,cm,fName,'gif', 'Loopcount',inf,'DelayTime',dT);
          nSV = 0;
        else
          imwrite(imind,cm,fName,'gif','WriteMode','append','DelayTime',dT);
        end
    end
end

end