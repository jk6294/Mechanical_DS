function lamb = lyap_exp(d2,d0,n,nP,c)
% Function for computing numerical Lyapunov exponent.
%
% Inputs
% d2:           Symbolic expression for map
% d0:   k x 1   Initial point of map
% n:    1 x 2   Number of throwaway and actual iterations
% nP:           Number of iterations to plot
% c:    k x 3   Matrix of colors. Omit if no plot desired
%
% Outputs
% lamb:         Lyapunov Exponent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Functions
f = matlabFunction(d2);
fp = diff(d2);
fpf = matlabFunction(fp);

% Throwaway Iterations
fprintf('Throwaway Iterations\n');
fprintf([repmat('.',1,100) '\n\n']);
tic
nV = [];
nTh = 0;
for i = 1:n(1)
    % Indicate Progress
    if(i/n(1) >= nTh)
        fprintf(repmat('\b',1,length(nV)));
        fprintf('\b=\n');
        nV = num2str(floor(toc*(n(1)-i)/i));
        fprintf(nV);
        nTh = nTh + .01;
    end
    d0 = f(d0);
end
fprintf(repmat('\b',1,length(nV)));

% Compute Trajectory
k = length(d0);
L = zeros(k,n(2));
L(:,1) = d0;
dL = zeros(k,n(2));
dL(:,1) = fpf(d0);
fprintf('Actual Iterations\n');
fprintf([repmat('.',1,100) '\n\n']);
nTh = 0;
nV = [];
tic
for i = 2:n(2)
    if(i/n(2) >= nTh)
        fprintf(repmat('\b',1,length(nV)));
        fprintf('\b=\n');
        nV = num2str(floor(toc*(n(2)-i)/i));
        fprintf(nV);
        nTh = nTh + .01;
    end
    L(:,i) = f(L(:,i-1));
    dL(:,i) = fpf(L(:,i));
end
fprintf(repmat('\b',1,length(nV)));

% Compute Exponent
lamb = sum(log(abs(dL)),2)/n(2);

% Plot
if(size(c,1) > 0)
    dV = linspace(min(d0),max(d0),1000);
    hold on;
    for i = 1:k
        D = [L(i,1:nP);L(i,1:nP)]; D = D(:); D = D(1:end-1);
        line(D(1:end-1),[0;D(3:end)],'color',[c(i,:) 1],'linewidth',.3);
    end
%     plot(dV,f(dV),'k-','linewidth',.7);
    hold off;
end

end

