function d2 = maps2(Xs,Xu,c,c0)
% Function for generating symbolic map of the distance d_2 as a function
% of distance d_1
%
% Inputs
% Xs        2 x 3:  Matrix of specified node positions. 
% Xu        2 x 2:  Matrix of unspecified node positions
%
% Outputs
% d2        4 x 1:  Symbolic vector of potential functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 1e-5;

% Assign edge lengths
conn = [1 1; 2 1; 3 1; 1 2; 2 2; 3 2];
L = sqrt(sum((Xs(:,conn(:,1))-Xu(:,conn(:,2))).^2));
d1 = sqrt(sum((Xs(:,1)-Xs(:,2)).^2));

% Orient About Pivot Node
Xu = Xu - Xs(:,2);
Xs = Xs - Xs(:,2);
thet = atan2d(Xs(2,1),Xs(1,1));
R = [cosd(-thet) -sind(-thet);...
     sind(-thet)  cosd(-thet)];
% R = rotz(-thet); R = R(1:2,1:2);
Xs = R*Xs; Xu = R*Xu;

% Calculate Variable Node Positions
x4 = (d1^2 + L(2)^2 - L(1)^2)/(2*d1);
x5 = (d1^2 + L(5)^2 - L(4)^2)/(2*d1);
y4v = 0; y5v = 0; y4 = []; y5 = [];
y4p = sqrt(L(2)^2 - x4^2);
y5p = sqrt(L(5)^2 - x5^2);
% Choose positive or negative y coordinates based on initial position
if(abs(y4p - Xu(2,1)) < tol)
    y4 = y4p;
    y4v = 1;
elseif(abs(y4p + Xu(2,1)) < tol)
    y4 = -y4p;
    y4v = -1;
end
if(abs(double(subs(y5p - Xu(2,2),c,c0))) < tol)
    y5 = y5p;
    y5v = 1;
elseif(abs(double(subs(y5p + Xu(2,2),c,c0))) < tol)
    y5 = -y5p;
    y5v = -1;
end

% Calculate Final Specified Node Position
syms({'x3p', 'y3p', 'x4p', 'x5p', 'y4p', 'y5p'});
E1 = L(3)^2 == (x3p-x4p)^2 + (y3p-y4p)^2;
E2 = L(6)^2 == (x3p-x5p)^2 + (y3p-y5p)^2;
s = solve(E1,E2,[x3p y3p]);
x3c = subs(subs(subs(subs(s.x3p,x4p,x4),y4p,y4),x5p,x5),y5p,y5);
y3c = subs(subs(subs(subs(s.y3p,x4p,x4),y4p,y4),x5p,x5),y5p,y5);
sv = sum([x3c y3c]'-Xs(:,3));
sInd = find(abs(double(subs(sv,c,.25))) < tol);

syms d1;
x4p = (d1^2 + L(2)^2 - L(1)^2)/(2*d1);
x5p = (d1^2 + L(5)^2 - L(4)^2)/(2*d1);
y4p = y4v*sqrt(L(2)^2 - x4p^2);
y5p = y5v*sqrt(L(5)^2 - x5p^2);
x3 = subs(s.x3p(sInd));
y3 = subs(s.y3p(sInd));

d2 = simplify(sqrt(x3^2 + y3^2));

end