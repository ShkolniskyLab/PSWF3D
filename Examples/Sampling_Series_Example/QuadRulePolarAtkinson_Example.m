function [ QuadRulePolar ] = QuadRulePolarAtkinson_Example( n )
% This function computes the Gauss-Legendre nodes and weights in theta.
% ,where theta is in [0,pi]. The final rule corresponds to the discussion in
% Atkinson's book.
%   
%  Input - n: Degree of ledendre polynomial to be used
%   
%  Output - QuadRule : The desired quadrature rule. The nodes appear in the  
%                      first row, while the weights appear in the second.

QuadRulePolar = zeros(2,n);

QuadRulePolar = LegndreQuadRecipies_Example(n,0,pi);

end

