function [ Coeff ] = SamplingCoeff_N( L,fVals,Points,N, prolate_dat, iserr , prolate_dat_tmp,truncationparam,flag )
% This function computes the coeeficients of the sampling series
% corresponding to the indices N,m,n (see the paper), Where N is
% fixed and m,n are all possible values.
%
%   Input - L : Sampling frequency.
%           fVals : Sampling values of the function f on the sampling points.
%           Points : A list of point that were sampled (x|y|z) .          
%           N : Prolate degree .
%           prolate_dat, iserr , prolate_dat_tmp : Data structure from
%           Lederman's code.
%           truncationparam : Truncation index of radial GPSWFs.
%           flag : 0 - if the code is used for construction of coefficients
%                  1 - if the code is used for approximating a value.
%
%   Output - Coeff : Array of Coefficients.

fVals1 = ((prolate_dat.c/(2*pi*L))^3)*fVals;
% Normalize the function values by the sampling frequency.

[ProlatesEval,alphavect] = EvalProlGrid_fixedN( Points, N,prolate_dat, iserr , prolate_dat_tmp,truncationparam,flag,L );
% Evaluate the GPSWFs corresponding to degree N at the sampling points.

Coeff = alphavect.*(conj(ProlatesEval.*repmat(alphavect,1,length(ProlatesEval(2,:))))*fVals1);
% Compute the coefficients.

end

