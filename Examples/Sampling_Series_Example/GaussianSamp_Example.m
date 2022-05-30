function [ fVals ] = GaussianSamp_Example( SamplePoints,sig,mu1,mu2,mu3 )
% This function generates the samples of a gaussian of the form
% f(x1,x2,x3) = 1/(2pi)^(3/2) * 1/sqrt(sig_1*sig_2*sig_3) * exp(-0.5*\sum (x_i-mu_i)^2/sig_i)
% on the points given in SamplePoints
%
%   Input - SamplePoints : Evaluation points for the Gaussian
%           sig1,sig2,sig3 : Variance parameters
%           mu1,mu2,mu3 : Expectation parameters
%
%   Output - fVals : Evaluations of the Gaussian on the points

fVals = zeros(length(SamplePoints),1);

for j = 1:length(SamplePoints)
    
   fact = (((SamplePoints(j,1)-mu1)^2)/sig)+(((SamplePoints(j,2)-mu2)^2)/sig)+(((SamplePoints(j,3)-mu3)^2)/sig);
   
   fVals(j,1) = (1/(2*pi)^(3/2))*(1/sqrt(sig^3))*exp(-0.5*fact); 
    
end


end

