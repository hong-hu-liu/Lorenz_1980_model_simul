function PAR_Lorenz9D = get_par_Lorenz9D(alpha, F1)
% set up parameters for the 9D Lorenz model
% alpha is the free parameter "a3" in Lorenz'1980 paper.

a1 = 1; a2 = 1; a3 = alpha;
F2 = 0; F3 = 0; 
b1 = 0.5*(a1-a2-a3);
b2 = 0.5*(a2-a3-a1);
b3 = 0.5*(a3-a1-a2);
c = sqrt(b1*b2 + b2*b3 + b3*b1); 
kappa0 = 1/48;
nu0 = 1/48;
g0 = 8;
h1 = -1; h2 = 0; h3 = 0;


% The following "table" shows how the 9 variables are ordered: 
    %  variable: x1     x2     x3    y1    y2    y3    z1    z2   z3
    % index:       1       2      3      4       5     6     7      8      9

%---------------------------------------------
% Constructing the linear part:     
A = zeros(9,9);

A(1,[1 4 7]) = [-nu0*a1, 1, -1];  % for x1-eqn 
A(2,[2 5 8]) = [-nu0*a2, 1, -1];  % for x2-eqn
A(3,[3 6 9]) = [-nu0*a3, 1, -1];  % for x3-eqn

A(4,[1 4]) = [-1, -nu0*a1];  % for y1-eqn
A(5,[2 5]) = [-1, -nu0*a2];  % for y2-eqn
A(6,[3 6]) = [-1, -nu0*a3];  % for y3-eqn

A(7,[1 2 3 5 6 7]) = [g0*a1, b3*h3, b2*h2, -c*h3,  c*h2, -kappa0*a1];  % for z1-eqn
A(8,[1 2 3 4 6 8]) = [b3*h3, g0*a2, b1*h1,  c*h3, -c*h1, -kappa0*a2];  % for z2-eqn
A(9,[1 2 3 4 5 9]) = [b2*h2, b1*h1, g0*a3, -c*h2,  c*h1, -kappa0*a3];  % for z3-eqn

%---------------------------------------------
% Constant forcing part
Forcing = zeros(9,1);
Forcing([7 8 9]) = [F1;F2;F3];

%---------------------------------------------
% Coefficients for the quadratic part
Gx1 = [ b1      % x2x3
        -c*(a1 - a3)/a1 % x2y3
         c*(a1 - a2)/a1 % y2x3
        -2*c^2/a1 % y2y3
         ];

Gx2 = [ b2  % x3x1
        -c*(a2 - a1)/a2 % x3y1
         c*(a2 - a3)/a2 % y3x1
        -2*c^2/a2 % y3y1
         ];

Gx3 = [ b3  % x1x2
        -c*(a3 - a2)/a3 % x1y2
         c*(a3 - a1)/a3 % y1x2
        -2*c^2/a3 % y1y2
         ];
     
Gy1 = [ -a3*b3/a1  % x2y3
        -a2*b2/a1  % y2x3
        c*(a3 - a2)/a1 % y2y3
         ];

Gy2 = [ -a1*b1/a2  % x3y1
        -a3*b3/a2  % y3x1
        c*(a1 - a3)/a2 % y3y1
         ];

Gy3 = [ -a2*b2/a3  % x1y2
        -a1*b1/a3  % y1x2
        c*(a2 - a1)/a3 % y1y2
         ];
     
Gz1 = [ -b3  % x2z3
        -b2  % z2x3
         c   % y2z3
        -c   % y3z2
         ];

Gz2 = [ -b1  % x3z1
        -b3  % z3x1
         c   % y3z1
        -c   % y1z3
         ];

Gz3 = [ -b2  % x1z2
        -b1  % z1x2
         c   % y1z2
        -c   % y2z1
         ];

   
G = cell(9,1);
G{1} = Gx1;
G{2} = Gx2;
G{3} = Gx3;
G{4} = Gy1;
G{5} = Gy2;
G{6} = Gy3;
G{7} = Gz1;
G{8} = Gz2;
G{9} = Gz3;

%---------------------------------------------
% Build the indices for the terms in the quadratic monomials, 
% where the index corresponding to x1, x2, x3, y1, y2, y3, z1, z2, z3 are
% given in line 18-19 above: 

Idx1 = cell(9,1);
Idx2 = cell(9,1);

%------------------------------------------------
% For x1-eqn:
% Note that the indices need to be recorded in the order that matches the order of the corresponding 
% coefficients recorded in Gx1 above. Recall that
%
% Gx1 = [ b1                 % x2x3
%           -c*(a1 - a3)/a1 % x2y3
%            c*(a1 - a2)/a1 % y2x3
%           -2*c^2/a1         % y2y3
%           ];
% 
% So the first term corresponds to x2x3. We store the index of the first variable (which is x2) in Idx1{1}, 
% and store the index of the second variable (which is x3) in Idx2{1}. 
% Thus we set 
% Idx1{1}(1) = 2, % because x2 is assigned with the index 2 as listed in lines 18--19 above.
% Idx2{1}(1) = 3; % because x3 is assigned with the index 3 as listed in lines 18--19 above.

% Likewise, the second term in Gx1 corresponds to x2y3, thus we set 
% Idx1{1}(2) = 2; 
% Idx2{1}(2) = 6; % because y3 is assigned with the index 6 as listed in lines 18--19 above.

% Following this procedure, the completed index listing for the x1-equation
% are given by Idx1{1} and Idx2{1} below. 

% The index listing of the quadratic terms in the other 8 equations is done
% in the same way.

Idx1{1} = [2 2 3 5];    
Idx2{1} = [3 6 5 6];
%------------------------------------------------

%------------------------------------------------
% For x2-eqn
Idx1{2} = [1 3 1 4]; 
Idx2{2} = [3 4 6 6]; 
%------------------------------------------------

%------------------------------------------------
% for x3-eqn
Idx1{3} = [1 1 2 4];
Idx2{3} = [2 5 4 5];
%------------------------------------------------

%------------------------------------------------
% for y1-eqn
Idx1{4} = [2 3 5];    
Idx2{4} = [6 5 6];
%------------------------------------------------

%------------------------------------------------
% for y2-eqn
Idx1{5} = [3 1 4]; 
Idx2{5} = [4 6 6]; 
%------------------------------------------------

%------------------------------------------------
% for y3-eqn
Idx1{6} = [1 2 4];
Idx2{6} = [5 4 5];
%------------------------------------------------

%------------------------------------------------
% for z1-eqn
Idx1{7} = [2 3 5 6];    
Idx2{7} = [9 8 9 8];
%------------------------------------------------

%------------------------------------------------
% for z2-eqn
Idx1{8} = [3 1 6 4]; 
Idx2{8} = [7 9 7 9]; 
%------------------------------------------------

%------------------------------------------------
% for z3-eqn
Idx1{9} = [1 2 4 5];
Idx2{9} = [8 7 8 7];
%------------------------------------------------
  
PAR.A = A;
PAR.Forcing = Forcing;
PAR.G = G;
PAR.Idx1 = Idx1;
PAR.Idx2 = Idx2;
PAR.alpha = alpha;
PAR.F1 = F1;

PAR_Lorenz9D = PAR;


