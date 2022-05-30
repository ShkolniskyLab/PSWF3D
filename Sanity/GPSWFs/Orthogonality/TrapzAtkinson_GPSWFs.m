function [ QuadRule ] = TrapzAtkinson_GPSWFs( m )
% Thid function generates composite trapezoid rule for itegrating periodic
% functions with period of 2pi. The notation corresponds to that which
% appears on page 167 in Atkinson's book.
%
%   Input - m: Number of segments
%
%   Output - QuadRule : The quadrature rule. The nodes appear in the first
%                       row and the corresponding weights in row2.

h = (2*pi)/m ;
QuadRule = zeros(2,m);

QuadRule(1,:) = [0:h:(m*h-h)];
QuadRule(2,:) = h;


end

