function [ InnerP ] = GaussInnerP_Example( N1,n1,N2,n2,phiarray,thetaarray )
% This function approximates the inner product of spherical harmonics with
% indices N1,n1 and N2,n2 by using Gauss-Legendre quadrature. Thr number of
% points is M in theta \in [0,pi] and 2M in phi \in [0,2*pi].
%
%   Input - N1,n1,N2,n2 : Indices of the desired spherical harmonics
%           phiarray,thetaarray : Quarature nodes and weights
%
%   Output - InnerP : approximation of inner function.

InnerP = 0;

for j = 1:length(phiarray)
    
   for k = 1:length(thetaarray)
       
      InnerP = InnerP +phiarray(2,j)*thetaarray(2,k)* sin(thetaarray(1,k))*EvalSpherHarmPoint_Example(N1,n1,phiarray(1,j),thetaarray(1,k))*conj(EvalSpherHarmPoint_Example(N2,n2,phiarray(1,j),thetaarray(1,k)));
       
   end    
    
end


end

