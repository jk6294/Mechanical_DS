% Figure 4: Large-Scale Conformational Changes
%% Prepare Space
clear; clc;

% Subplot Indices
NCol = 14;
colBal = {[1 1 3], [1 1 3], [1 1 3]};
rowBal = [9 9 9]; 
NRow = sum(rowBal);
ru = cumsum(rowBal)-1;
rd = [0 ru(1:end-1)+1]+1;
cellL = [0 cumsum(cellfun('length',colBal))];
cellM = cell(1,cellL(end));
for i = 1:length(colBal)
    s = colBal{i}; 
    su = cumsum(s/sum(s)*(NCol+1))-1;
    sd = [1 su(1:end-1)+2];
    for j = 1:length(s)
        cellM(j+cellL(i)) = {[[sd(j) su(j)]+rd(i)*NCol,...
                              [sd(j) su(j)]+ru(i)*NCol]};
    end
end

% Figure Axis Bounds
axM = [-1.4 1.4 -0.8 2];
labX = -.1;
labY = 0.90;
labColY = 1.3;

fig = figure(4); clf;

% Sizing
s = sqrt(3);
pSc = 0.6;

% Node Colors for different curvatures
cTr1 = [126 240 240]/255;
cTr2 = [115 164 211]/255;
cTr3 = [015 082 186]/255;

% Module Parameters
l1 = [0;-s]; l2 = [0;-1.5*s];
a1 = 17; a2 = 38; a3 = 60;
R1 = rotz(-a1/2); R2 = rotz(-a2/2); R3 = rotz(-a3/2);
R1 = R1(1:2,1:2); R2 = R2(1:2,1:2); R3 = R3(1:2,1:2);


%% a: Small Motion
subplot(NRow,NCol,cellM{1}); cla;
Xs10 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs1T = [R1*l2 [0;0] R1\l2] + [0;1];
Xu1 = [-1.8  1.8;...
       3/3 3/3];
conn1 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu1,fV] = construct_network(Xs10,Xs1T,Xu1,conn1,0,1);
Xu1 = Xu1(1:2,:);
% visualize_conic_finite(Xs10,Xs1T,[-2 2;-2 2],[100;100],0,1,1);
visualize_network(Xs10,Xu1,conn1,1,[255 100 100]/255, cTr1);
[XC1,fC] = sim_motion(Xs10,Xu1,conn1,.01,200,-[Xs10 Xu1],0);
axis(axM);
subplot(NRow,NCol,cellM{2}); cla;
d1 = sqrt(squeeze(sum((XC1(:,1,:)-XC1(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC1(:,2,:)-XC1(:,3,:)).^2)));
d3 = sqrt(squeeze(sum((XC1(:,1,:)-XC1(:,3,:)).^2)));
distv = sqrt(sum(diff([Xs1T Xs1T(:,1)],1,2).^2))
distL = sum(abs(distv - [d1 d2 d3]),2);
pInd = find(distL==min(distL),1);
visualize_network(XC1(:,1:3,pInd),XC1(:,4:5,pInd),conn1,...
                  1,[255 100 100]/255, cTr1);
axis(axM);
% plot(d1); hold on; plot(d2); plot(d3); hold off;




%% b
subplot(NRow,NCol,cellM{4}); cla;
Xs20 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs2T = [R2*l2 [0;0] R2\l2] + [0;1];
Xu2 = [-0.8  0.8;...
        1.5  1.5];
conn2 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu2,fV] = construct_network(Xs20,Xs2T,Xu2,conn2,0,1);
Xu2 = Xu2(1:2,:);
% visualize_conic_finite(Xs20,Xs2T,[-2 2;-2 2],[100;100],4,1,1);
visualize_network(Xs20,Xu2,conn2,1,[255 100 100]/255, cTr2);
[XC2,fC] = sim_motion(Xs20,Xu2,conn2,.01,200,[Xs20 Xu2],0);
axis(axM);

subplot(NRow,NCol,cellM{5}); cla;
d1 = sqrt(squeeze(sum((XC2(:,1,:)-XC2(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC2(:,2,:)-XC2(:,3,:)).^2)));
d3 = sqrt(squeeze(sum((XC2(:,1,:)-XC2(:,3,:)).^2)));
distv = sqrt(sum(diff([Xs2T Xs2T(:,1)],1,2).^2))
distL = sum(abs(distv - [d1 d2 d3]),2);
pInd = find(distL==min(distL),1);
visualize_network(XC2(:,1:3,pInd),XC2(:,4:5,pInd),conn1,...
                  1,[255 100 100]/255, cTr2);
axis(axM);
% plot(d1); hold on; plot(d2); plot(d3); hold off;


%% c
subplot(NRow,NCol,cellM{7}); cla;
Xs30 = [-s/2  0    s/2;...
       -0.5  1.0 -0.5];
Xs3T = [R3*l2 [0;0] R3\l2] + [0;1];
Xu3 = [-0.8  0.8;...
        1.0  1.0];
conn3 = [1 4; 2 4; 3 4; 1 5; 2 5; 3 5];
[Xu3,fV] = construct_network(Xs30,Xs3T,Xu3,conn3,0,1);
Xu3 = Xu3(1:2,:);
% visualize_conic_finite(Xs30,Xs3T,[-2 2;-2 2],[100;100],4,1,1);
visualize_network(Xs10,Xu1,conn1,1,[255 100 100]/255, cTr3);
[XC3,fC] = sim_motion(Xs30,Xu3,conn3,.01,200,[Xs30 Xu3],0);

subplot(NRow,NCol,cellM{8}); cla;
d1 = sqrt(squeeze(sum((XC3(:,1,:)-XC3(:,2,:)).^2)));
d2 = sqrt(squeeze(sum((XC3(:,2,:)-XC3(:,3,:)).^2)));
d3 = sqrt(squeeze(sum((XC3(:,1,:)-XC3(:,3,:)).^2)));
distv = sqrt(sum(diff([Xs3T Xs3T(:,1)],1,2).^2))
distL = sum(abs(distv - [d1 d2 d3]),2);
pInd = find(distL==min(distL),1);
visualize_network(XC3(:,1:3,pInd),XC3(:,4:5,pInd),conn1,...
                  1,[255 100 100]/255, cTr3);
% plot(d1); hold on; plot(d2); plot(d3); hold off;
axis(axM);

