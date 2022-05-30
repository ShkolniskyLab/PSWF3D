function [ InnerPtable ] = GaussInnerPTable_Example( Basis,M )
% This function computes the inner product table for all spherical
% harmonics in a given list. The computation is base of Gauss-Legendre
% quadrature with parameter M.
%
%   Input - Basis : List of indices of spherical harmonics of the form N|n
%           M : Gauss-Legendre quadrature parameter
%
%   Output - InnerPtable : Table of inner products

InnerPtable = zeros(length(Basis(:,1)));
phiarray = TrapzAtkinson_Example(2*M);
thetaarray = QuadRulePolarAtkinson_Example(M);

for j = 1:length(Basis)
   
   N1 = Basis(j,1);
   n1 = Basis(j,2);
   
   for k = 1:j
       
       N2 = Basis(k,1);
       n2 = Basis(k,2);
       
       InnerPtable(j,k) = GaussInnerP_Example( N1,n1,N2,n2,phiarray,thetaarray ) ;
       
   end    
    
end


end

