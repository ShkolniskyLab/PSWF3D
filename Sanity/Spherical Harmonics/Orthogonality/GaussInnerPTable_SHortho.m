function [ InnerPtable ] = GaussInnerPTable_SHortho( M,N )
% This function computes the inner product table for all spherical
% harmonics with degree S<=N. The computation is based on Gauss-Legendre
% quadrature with parameter M.
%
%   Input - M : Gauss-Legendre quadrature parameter
%           N : Maximum order of Spherical Harmonics
%
%   Output - InnerPtable : Table of inner products

Basis  = GenSHBasis_SHortho(N);

InnerPtable = zeros(length(Basis(:,1)));
phiarray = TrapzAtkinson_SHortho(2*M);
thetaarray = QuadRulePolarAtkinson_SHortho(M);

for j = 1:length(Basis)
   
   N1 = Basis(j,1);
   n1 = Basis(j,2);
   
   for k = 1:j
       
       N2 = Basis(k,1);
       n2 = Basis(k,2);
       
       InnerPtable(j,k) = GaussInnerP_SHortho( N1,n1,N2,n2,phiarray,thetaarray ) ;
       
   end    
    
end


end

