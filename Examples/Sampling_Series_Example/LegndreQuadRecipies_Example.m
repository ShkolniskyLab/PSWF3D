function [ QuadRule ] = LegndreQuadRecipies_Example( n,x1,x2 )
% This function computes the Gauss-Legendre nodes and weights over [-1,1].
% The construction is based on the corresponding code in "Numerical
% Recipies"
%   
%  Input - n: Degree of ledendre polynomial to be used
%          x1,x2 : End points of the desired segments (x1<x2)
%   
%  Output - QuadRule : The desired quadrature rule. The nodes appear in the  
%                      first row, while the weights appear in the second.

QuadRule = zeros(2,n);

xm = (x1+x2)/2;
xl = (x2-x1)/2;
m = floor((n+1)/2);
EPS = 10^(-14);

for i = 0:(m-1)
    
    z = cos(pi*(i+0.75)/(n+0.5));
    z1 = 10;
    itermax = 300;
    count = 0;
    
    while and((abs(z-z1)>EPS),count<=itermax)
        
        count = count+1;
        p1 = EvalLeg_Example(n,z);
        p2 = EvalLeg_Example(n-1,z);
        pp = n*((z*p1-p2)/(z^2-1));
        z1 = z;
        z = z1 - p1/pp;
        
    end
    
    QuadRule(1,i+1) = xm -xl*z;
    QuadRule(1,n-i) = xm+xl*z;
    QuadRule(2,i+1) = 2*xl/(( 1-z^2)*pp^2);
    QuadRule(2,n-i) =  QuadRule(2,i+1);
    
end



end

