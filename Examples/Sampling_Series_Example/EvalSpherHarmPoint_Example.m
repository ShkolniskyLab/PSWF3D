function [ Val ] = EvalSpherHarmPoint_Example( N,m,phi,theta)
% This function evaluates the spherical harmonic Y_{N,m}(phi,theta) where
% phi is the azimuth (between 0 and 2*pi) and theta is the elevation
% (between 0 and pi). The construction is base on the book 'Numerical Recipies'

% Input - N : Degree of spherical harmonic.
%         m : Index of spherical harmonic.
%         phi : Azimuth.
%         theta : Elevation.
%
% Output - Val : Value of spherical harmonic at the desired point.

assert(abs(m)<=N,'Invalid indices')

if m >= 0
   
    Val = EvalPTilde_Example(N,m,cos(theta)).*exp(1i.*m.*phi);

else
    
    m1 = abs(m);
    Val = ((-1)^m1)*EvalPTilde_Example(N,m1,cos(theta)).*exp((-1)*1i.*m1.*phi);
    
end

end

