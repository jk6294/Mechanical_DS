%% Global parameters for figures
% Latex interpreter
set(groot,'defaulttextinterpreter','latex');

% Clear annotations
delete(findall(gcf,'type','annotation'));

% Plot parameters
FS = 10;                % Fontsize
gr = 0.8;               % Ratio for gray color
al = 0.2;               % Value of alpha transparency

% Text parameters
NVTitle  = {'Units','centimeters','fontsize',FS};
NVTitleR  = {'Units','centimeters','fontsize',FS,...
             'HorizontalAlignment','right'};
NVTitleH  = {'Units','centimeters','fontsize',FS,...
             'HorizontalAlignment','center'};
NVTextH  = {'Units','Normalized','fontsize',FS,...
            'HorizontalAlignment','center'};
NVTextR  = {'Units','Normalized','fontsize',FS,...
            'HorizontalAlignment','right'};
NVTextL  = {'Units','Normalized','fontsize',FS};
NVTexth  = {'fontsize',FS,'HorizontalAlignment','center'};
NVTextr  = {'fontsize',FS,'HorizontalAlignment','right'};
NVTextl  = {'fontsize',FS};
NVTexthv  = {'fontsize',FS,'HorizontalAlignment','center',...
             'VerticalAlignment','middle'};
NVTextrv  = {'fontsize',FS,'HorizontalAlignment','right',...
             'VerticalAlignment','middle'};
NVTextlv  = {'fontsize',FS,'HorizontalAlignment','left',...
             'VerticalAlignment','middle'};

% Colors
nT = 1000;              % Color interpolation count
o = [1 1 1];            % Vector of ones
% Inter-node distances
C1a = [047 086 151]/255;
C1b = [140 181 063]/255;
C1c = [231 178 072]/255;
CP1 = interp1([0 .5 1],[C1a;C1b;C1c],linspace(0,1,nT));
% Gradient: Distance c. Interpolate lightness of one color
C2a = [133 37 83]/255;
C2b = C2a.^.1;
CP2 = interp1([0 1],[C2a;C2b],linspace(0,1,nT));
C3a = [197 066 085]/255;

% Rotations
R45 = rotz(45); R45 = R45(1:2,1:2);
R90 = rotz(90); R90 = R90(1:2,1:2);