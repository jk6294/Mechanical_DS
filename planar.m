%% Prepare Space
clear; clc;


%% Points
xs = sym('xs', [4 1]); ys = sym('ys', [4 1]);
xf = sym('xf', [4 1]); yf = sym('yf', [4 1]);
dx = sym('dx', [4 1]); dy = sym('dy', [4 1]);

assume(xs,'real'); assume(ys,'real'); 
assume(xf,'real'); assume(yf,'real');
assume(dx,'real'); assume(dy,'real');

Xs = [-2 -2 2 2]' + rand(4,1)/10; Ys = [1 -1 1 -1]' + rand(4,1)/10;
Xf = [-2 -2 2 2]' + rand(4,1)/10; Yf = [1 -1 1 -1]'*1.2 + rand(4,1)/10;
dX = [-.2 -.2 .3 .4]'; dY = [.2 -.2 .2 -.2]';`

connS = [1 1; 1 2; 2 1; 2 4; 3 2; 3 3; 4 3; 4 4];
connU = [1 2; 2 3; 3 4; 4 1];

fSs = sum([Xs(connS(:,1))-xs(connS(:,2)) Ys(connS(:,1))-ys(connS(:,2))] .* [dX(connS(:,1))-dx(connS(:,2)) dY(connS(:,1))-dy(connS(:,2))],2);
fSf = sum([Xs(connS(:,1))-xs(connS(:,2)) Ys(connS(:,1))-ys(connS(:,2))].^2 - [Xf(connS(:,1))-xf(connS(:,2)) Yf(connS(:,1))-yf(connS(:,2))].^2,2);
fUs = sum([xs(connU(:,1))-xs(connU(:,2)) ys(connU(:,1))-ys(connU(:,2))] .* [dx(connU(:,1))-dx(connU(:,2)) dy(connU(:,1))-dy(connU(:,2))],2);
fUf = sum([xs(connU(:,1))-xs(connU(:,2)) ys(connU(:,1))-ys(connU(:,2))].^2 - [xf(connU(:,1))-xf(connU(:,2)) yf(connU(:,1))-yf(connU(:,2))].^2,2);

f = [fSs; fSf; fUs; fUf];
A = vpasolve(f==0)

% Plot
XS = [A.xs1 A.xs2 A.xs3 A.xs4]';
YS = [A.ys1 A.ys2 A.ys3 A.ys4]';
XF = [A.xf1 A.xf2 A.xf3 A.xf4]';
YF = [A.yf1 A.yf2 A.yf3 A.yf4]';
DX = [A.dx1 A.dx2 A.dx3 A.dx4]';
DY = [A.dy1 A.dy2 A.dy3 A.dy4]';

figure(1); clf;
plot(Xs,Ys,'o');
hold on;
plot(XS,YS,'o');
line([Xs(connS(:,1)) XS(connS(:,2))]', [Ys(connS(:,1)) YS(connS(:,2))]');
line([XS(connU(:,1)) XS(connU(:,2))]', [YS(connU(:,1)) YS(connU(:,2))]');
hold off;


%% Planar
% Planar 1D
% Xs = [ 0  0  0  1.2  1.2  4  8.2  8.2;...
%        2  0 -2  .02 -.02  0  2 -2];
% conn = [1 2; 2 3; 1 4; 2 4; 2 5; 3 5; 4 6; 4 7; 5 6; 5 8; 6 7; 6 8];

% Planar 2D
% Xs = [-4.0  4.0  0.0  -2*sqrt(3)  2*sqrt(3) -2*sqrt(3)  2*sqrt(3) -4*sqrt(3) 4*sqrt(3)  0.0;...
%        0 0 6 6 6 10 10 8 8 12];
% conn = [1 3; 2 3; 1 4; 2 5; 3 4; 3 5; 3 6; 3 7; 4 8; 4 6; 5 7; 5 9; 8 6; 9 7; 6 10; 7 10];
%
% Other
Xs = [-sqrt(3)/2 0 sqrt(3)/2;...
      -1/2       1 -1/2];
s = 5/8 * sqrt(3)/2;
Xu = [-s  s  ;...
      1-sqrt(3)*s 1-sqrt(3)*s];
Xs = [Xs Xu];
conn = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
LVal = sqrt(sum((Xs(:,conn(:,1))-Xs(:,conn(:,2))).^2));
      


figure(1); clf;
subplot(2,2,1);                     % Single Module
visualize_network(Xs,[],conn);
[Xsa, conna] = tesselate_network(Xs,conn,[8.2;0],[10;1]);
drawnow;
subplot(2,2,2);
[XC, fV] = sim_motion(Xs, [], conn, 0.001, 1200, -Xs, 1);
subplot(2,2,3);
d1 = squeeze(sqrt(sum((XC(:,1,:)-XC(:,3,:)).^2)));
d2 = squeeze(sqrt(sum((XC(:,1,:)-XC(:,2,:)).^2)));
plot(d1, d2);
hold on;
plot([0 2], [0 2]);
hold off;

nO = 100;
ad = abs(d1(nO:end)-d2(nO:end));
nV = find(ad == min(ad))+nO-1;

Xs = XC(:,:,nV);

% subplot(2,2,3);                     % Combined Module
% visualize_network(Xsa,[],conna);
% drawnow;
% subplot(2,2,4);                     % Simulation
% LVal = sqrt(sum((Xsa(:,conna(:,1))-Xsa(:,conna(:,2))).^2))';
% [XC, fV] = sim_motion(Xsa,[],conna,.05,2000,Xsa,0);


%% Triangle?
S = zeros(max(max(conna)));
for i = 1:size(conna,1)
    S(conna(i,1),conna(i,2)) = 1;
end
S = S + S';


%% Animate
fig = figure(2); clf;
fName = 'planar_ss.gif';
nSV = 1;
dT = 0.03;
nV = [1:3];
nS = 2;
for i = 1:1:size(XC,3)
    cla;
    visualize_network(XC(:,:,i), [], conn);
    axis([min(min(min(XC(1,:)))) max(max(max(XC(1,:)))) min(min(min(XC(2,:))))  max(max(max(XC(2,:))))]);
    drawnow;

%     % Capture the plot as an image
%     frame = getframe(fig);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if nSV == 1
%       imwrite(imind,cm,fName,'gif', 'Loopcount',inf,'DelayTime',dT);
%       nSV = 0;
%     else
%       imwrite(imind,cm,fName,'gif','WriteMode','append','DelayTime',dT);
%     end
end

