function [  ] = TestSamplingSeries_Example( )

% This is an example of a complete run of the sampling series approximation
% The approximated function is a Gaussian. The appromation grid is a set of
% Gaussian quadrature nodes. The code estimates the L^2 approximation error
% on the unit ball in R^3. The error is displayed and compared to a
% tolerance.

%%%% Part 1 - Assign input variables %%%%

c = 24*pi; %Bandlimit
L = 24; %Sampling frequency
CartGrid = GenSamplePoints3D_V2_Example(L); %Generate Cartesian grid of samples
T = 10; %Truncation Parameter
sig = 0.0136; % Variance of Gaussian
matdim = 300; % Matrix truncation parameter (Lederman's code)
minEigenvalRatio = 10^(-16); % Eigenvalue truncation parameter (Lederman's code)
[quadpoints,quadweights]=GenGridPoints_Example(30,30); % Generate quadrature nodes and weights
Gaussquad = GaussianSamp_Example( quadpoints,sig,10^(-5),10^(-5),0.001 ); % Sample the Gaussian on the quadrature grid
GaussSamp = GaussianSamp_Example( CartGrid,sig,10^(-5),10^(-5),0.001 ); % Sample the Gaussian on the Cartesian grid.
tol = 10^(-13);

%%%% Part 2 - Construct coefficients of sampling series %%%%

[ SampleSeriesCoeff,truncvec,Proldat ] = ConstructSeriesCoeff_Example( T,L,GaussSamp,CartGrid,c, minEigenvalRatio, matdim );

%%%% Part 3 - Evaluate the sampling series on the quadrature grid %%%%

SeriesEval = EvalSeries_Example( quadpoints,SampleSeriesCoeff,truncvec,Proldat,L );

%%%% Part 4 - Compute and display the L^2 error on the unit ball %%%%

norm = 0;

for j = 1:length(quadweights)

    norm = norm + quadweights(j,1)*abs(SeriesEval(1,j)-Gaussquad(j,1))^2;
    
end

if sqrt(norm)>tol
    
    disp(strcat('The approximation error :',num2str(sqrt(norm))));
    disp('Test failed - Error is too large')
    
else
   
    disp(strcat('The approximation error :',num2str(sqrt(norm))));
    disp('Test passed - Error is ok')
    
end

end

