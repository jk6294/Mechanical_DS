function lamb = lyap_exp(d2,d0,n,c)
% Function for computing numerical Lyapunov exponent.
%
% Inputs
% d2:           Symbolic expression for map
% d0:           Initial point of map
% n:            Number of iterations; should be large.
%
% Outputs
% lamb:         Lyapunov Exponent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary Functions
f = matlabFunction(d2);
fp = diff(d2);
fpf = matlabFunction(fp);

% Compute Trajectory
L = zeros(1,n);
L(1) = d0;
dL = zeros(1,n);
dL(1) = fpf(d0);
for i = 2:n
    L(i) = f(L(i-1));
    dL(i) = fpf(L(i));
end

% Compute Exponent
lamb = sum(abs(log(dL)))/n;

figure(1);
dV = 0.742:.0005:1.872;
D = [L;L]; D = D(1:end-1)';
hold on;
line(D(1:end-1),[.5;D(3:end)],'color',c);
plot(dV,f(dV),'k-','linewidth',2);
hold off;

figure(2); clf;
histogram(L);

