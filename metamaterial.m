%% Prepare Space
clear; clc;


%% Design Network
s3 = sqrt(3)/2;
Xs = [-s3   0.0  s3;...
      -0.5  1.0 -0.5];
Xf1 = Xs*1.5;
Xf2 = [[-s3   0.0  s3]*0.8;...
       -0.5  1.0 -0.5];
% -0.6266 -1.1580
%  1.0280 -0.2957
for th1 = 210
    for th2 = 0:10:350
        if(th1~=th2)
% Xu = [-0.7870  0.5870;...
%        1.8173 -0.5173];
R1 = rotz(-th1); R1 = R1(1:2,1:2);
R2 = rotz(-th2); R2 = R2(1:2,1:2);
v = [0;1];
Xu = [R1*v*1.8 R2*v*1.5];

% Construct
figure(1); clf;
% visualize_conic_finite(Xs, Xf1, [-2 2; -2 2], [400;400], 0, 1, 0);
% visualize_conic_finite(Xs, Xf2, [-2 2; -2 2], [400;400], 0, 1, 0);
conn = [1 4; 1 5; 2 4; 2 5; 3 4; 3 5];
visualize_network(Xs,Xu,conn);

% Simulate
% figure(2); clf;
[XC, fC] = sim_motion(Xs, Xu, conn, 0.03, 500, [Xs Xu], 0);

% Map
figure(2); clf;
d1 = sqrt(squeeze(sum((XC(:,1,:)-XC(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC(:,2,:)-XC(:,3,:)).^2)));
plot(d2,d1);
hold on;
plot([min([d1;d2]) max([d1;d2])], [min([d1;d2]) max([d1;d2])]);
hold off;
axis([1.2 5 1.2 5]);

figure(3); clf;
plot(d2,d1-d2);
hold on;
plot([min(d2) max(d2)], [0 0]);
hold off;
axis([min(d2) max(d2) -.1 .1]);
pause(.1);
        end
    end
end
        


%% Tesselate
Xsp = [Xs [2*s3; 1]];
Xup = [Xu [Xu(1,:)+s3; -Xu(2,:)+.5]];
conn = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
figure(4); clf;
visualize_network(Xsp,Xup,conn);

[Xsa, conna] = tesselate_network([Xsp Xup], conn, [2*s3; 0], [10; 1]);
figure(5); clf;
visualize_network(Xsa, [], conna);


%% Simulate
figure(6); clf;
[XC, fC] = sim_motion(Xsa, [], conna, 0.03, 1000, [Xsa], 0);
% Rigidity
r = zeros(size(XC,3),1);
for i = 1:size(XC,3)
    RP = rigidity(XC(:,:,i),conna);
    r(i) = det(RP*RP');
end
plot(log10(r));
d1 = sqrt(squeeze(sum((XC(:,1,:)-XC(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC(:,2,:)-XC(:,4,:)).^2)));
d3 = sqrt(squeeze(sum((XC(:,4,:)-XC(:,7,:)).^2)));
% plot(d3, d2);
% hold on;
% plot(d2, d1);
% hold off;



%% Animation
fig = figure(7); clf;

% XCP = permute(XCd, [2 1 3]);
XCP = XC;
fName = 'animation.gif';
nSV = 1;
dT = 0.03;
nV = [1:3];
nS = 2;
for i = 1:1:size(XCP,3)
    cla;
    visualize_network(XCP(:,:,i), [], conn);
    axis([min(min(min(XCP(1,:)))) max(max(max(XCP(1,:)))) min(min(min(XCP(2,:))))  max(max(max(XCP(2,:))))]);
    drawnow;

    % Capture the plot as an image
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


