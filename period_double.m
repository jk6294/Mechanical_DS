function [cV,lc,LC,dlc] = period_double(Xs,Xui,Xuf,nS,x0,nW)
% Function to compute a period doubling map route to chaos.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameterize added node positions with single variable
syms('c','d1','real');
Xud = Xuf - Xui;
Xu = sym(Xui) + Xud*c;


%% Create Maps
d2 = maps2(Xs,Xu,c,.5);
d2F = matlabFunction(d2, 'Vars', {c,d1});
dd2 = diff(d2,d1);
dd2F = matlabFunction(dd2, 'Vars', {c,d1});
nR = 6;
cV = fliplr((10^nR-logspace(0,nR,nS))/(10^nR-1));
% nW = 50000;

% figure(1); clf;
% d1v = linspace(sqrt(2)/2,2,1000);
d1v = linspace(.5,2,1000);

lc = zeros(1,nS);
LC = cell(1,nS);
dlc = zeros(1,nS);
% Collect Maps
hold on;

fprintf([repmat('.',1,nS) '\n\n']);
tic
parfor i = 1:nS
%     plot(d1v,d2F(cV(i),d1v));
%     drawnow;
    fprintf('\b=\n')
    % Initial runs
    xiv = x0;
    for j = 1:nW(i)
        xiv = d2F(cV(i),xiv);
    end
    
    % Limit cycle runs
    cLim = 2^8;
    xh = zeros(1,cLim);
    nI = 1;
    while(sum(abs(xiv-xh)<1e-8)==0)
        xh(nI) = xiv;
        xiv = d2F(cV(i),xiv);
        nI = nI + 1;
        if(nI > cLim+1)
            disp('greater');
            break;
        end
    end
    lc(i) = nI - find(abs(xiv-xh)<1e-8,1);
    LC{i} = sort(xh((nI-lc(i)):(nI-1)))';
    dlc(i) = prod(dd2F(cV(i),LC{i}));
end
toc

end