function [ Basis ] = GenSHBasis_Example( S )
% This function generates a table of indices for generating spherical
% harmonics. The table is of the form N|m . The indices are every possible
% combination of N,m where N<=S
%
%   Input - S : Maximal order of N
%
%   Output - Basis : Table of indices

Basis = [0,0];

for k = 1:S
    
   for l = (-k):k
       
        Basis = vertcat(Basis,[k,l]);
       
   end
    
end


end

