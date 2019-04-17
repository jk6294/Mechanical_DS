% Figure 3: Designing Folding Sequence
%% Prepare Space
clear; clc;

% Subplot Indices
nRow = [6 6 6];     NRow = sum(nRow);   nR = [0 cumsum(nRow)];
nCol = [5 5 5 5];   NCol = sum(nCol);   nC = [0 cumsum(nCol)];
figM = reshape(1:(NRow*NCol), [NCol, NRow])';
cellM = cell(length(nRow)*length(nCol),1);
for i = 1:length(nRow)
    for j = 1:length(nCol)
        cP = figM((nR(i)+1):(nR(i+1)-1), (nC(j)+1):(nC(j+1)-1));
        cellM(i+(j-1)*length(nRow)) = {cP(:)};
    end
end

fig = figure(3); clf;
s = sqrt(3);


% a: Single Network: Amplify
subplot(NRow,NCol,cellM{1});
Xs10 = [-s/2  0    s/2;...
        -0.5  1.0 -0.5];
Xs1T = [-s/2  0    s/2;...
        -1.0  1.5 -1.0];
dXs1 = [ 0.1  0.0 -0.1;...
        -0.2  0.0 -0.4];
visualize_conic(Xs10,dXs1,[-2 2;-2 2],[400;400],0,1,1);
visualize_conic_finite(Xs10,Xs1T,[-2 2;-2 2],[400;400],5,1,1);


%% b: Simulate Single
subplot(NRow,NCol,cellM{2});
Xu1 = [ 0.5865 -1.4190;...
        1.5090 -1.0000];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[XC,fC] = sim_motion(Xs10,Xu1,conn1,.01,1000,[Xs10 Xu1],0);
d1 = sqrt(squeeze(sum((diff(XC(:,1:2,:),[],2)).^2)));
d2 = sqrt(squeeze(sum((diff(XC(:,2:3,:),[],2)).^2)));
plot(d1,d2);
hold on;
plot([0 3],[0 3],'--');
hold off;


   
%% d: Combined Motion
subplot(NRow,NCol,cellM{3});
Xs1c = [Xs10(:,1:2)  Xs10(:,1:2)+[s;0]];
Xu1c = [Xu1 [Xu1(1,:)+s/2;-Xu1(2,:)+.5]];
conn1c = [1 5; 1 6; 2 5; 2 6; 2 7; 2 8; 3 5; 3 6; 3 7; 3 8; 4 7; 4 8];
[Xs1a,Xu1a,conn1a] = tesselate_network_old(Xs1c,Xu1c,conn1c,[s;0],[10;1]);
visualize_network(Xs1a,Xu1a,conn1a);


%% e: Propagate Combined Motion
[XC,fC] = sim_motion(Xs1a,Xu1a,conn1a,.01,2000,-[Xs1a Xu1a],0);
D1 = sqrt(squeeze(sum(diff(XC,1,2).^2)));
D1 = D1(1:size(Xs1a,2)-1,:);





%% Animation
fig = figure(4); clf;
fName = 'fixed_point_sim.gif';
% dT = 0.03;
dT = 0.1;
for i = 1:20:size(XC,3)
    cla;
    dP = [D1(:,i)';D1(:,i)']; dP = dP(:);
    dPa = dP(1:end-1); dPb = [1;dP(3:end)];
    plot(d1,d2);
    hold on;
    plot([0 3],[0 3]);
    line(dP(1:end-1),[1;dP(3:end)]);
    visualize_network(XC(:,1:size(Xs1a,2),i)/5,...
                      XC(:,[1:size(Xu1a,2)]+size(Xs1a,2),i)/5,conn1a,.5);
%     construct_motion(XC(:,1:size(Xsa,2),i), XC(:,1:size(Xsa,2),i+1)-XC(:,1:size(Xsa,2),i), XC(:,[1:size(Xua,2)]+size(Xsa,2),i), conn, 20, 20);
%     axis([min(min(XC(1,:,:))) max(max(XC(1,:,:))) min(min(XC(2,:,:))) max(max(XC(2,:,:)))]);
    set(gca,'visible',0);
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


