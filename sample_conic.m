function [C,k] = sample_conic(Q,r,n,o)
% Function for regularly sampling a conic from its canonical form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% Q             3 x 3 matrix of conic in matrix form
% r             1 x 4 vector for sample window [xmin,xmax,ymin,ymax]
% n             1 x 1 scalar of number of points to sample
% o             1 x 1 scalar optional for rotating initial sample point
%
% Output
% C             2 x k matrix of sampled coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 3
    o = 0;
end

% Find center
xc = Q(1:2,1:2)\Q(1:2,3);

% Diagonalize
S = det(Q);
[T,l] = eig(Q(1:2,1:2));
l = diag(l);
a2 = -S/(l(1)^2*l(2));
b2 = -S/(l(1)*l(2)^2);

% Draw conic in ab space
thet = linspace(0,360,n+1); thet = thet(1:end-1) + o;
sQ = det(Q(1:2,1:2));
% Ellipse
if sQ > eps
    % Degenerate: Single point
    if abs(S) < eps
        alphthet = 0;
        betthet = 0;
    % Not degenerate
    else 
        % Real ellipse
        if (Q(1,1) + Q(2,2))*S < 0
            a = sqrt(a2);
            b = sqrt(b2);
            alphthet = a*cosd(thet);
            betthet = b*sind(thet);
        % Imaginary
        else
        end
    end
% Hyperbola
elseif sQ < -eps
    if abs(S) < eps
        slo = sqrt(-l(2)/l(1));
        v = T*[slo;1];
        bet = ([r(1:2); r(3:4)] + xc) ./ v;
        betthetF = linspace(min(bet(:)), max(bet(:)),n);
        betthetC = linspace(min(bet(:)), max(bet(:)),n);
        alphthetF =  slo*betthetF;
        alphthetC = -slo*betthetC;
        betthet = [betthetF betthetC];
        alphthet = [alphthetF alphthetC];
    elseif(a2 < b2)
        a = sqrt(-a2);
        b = sqrt(b2);
        alphthet = a*tand(thet);
        betthet = b*secd(thet);
    else
        a = sqrt(a2);
        b = sqrt(-b2);
        alphthet = a*secd(thet);
        betthet = b*tand(thet);
    end
% Parabola
elseif abs(sQ) < eps
end

% Undo linear transformation
ab = [alphthet;betthet];
C = T*ab - xc;
cInds = find(C(1,:) >= r(1) & C(1,:) <= r(2) & C(2,:) >= r(3) & C(2,:) <= r(4));
C = C(:,cInds);
k = sum(cInds < n);
end