function [ Val ] = EvalLeg_SH( n,x )
% This function evaluates the legendre polynomial P_n of degree n at the
% point given by x
%
%   Input - n : Degree of the desired Legendre polynomial
%           x : Evaluation point
%
%   Output - Val : The value P_n(x)

val0 = 1;
val1 = x;

for i = 1:(n-1)
    
    valtemp = ((2*i+1)*x*val1 - i*val0)/(i+1);
    val0=val1;
    val1=valtemp;
end

Val = val1;

end

