function [err,bound,delt,eps,ro] = CompVal_AllFigs( c,L,T,sig,mu1,mu2,mu3,Smax )
% This function computes the approximation error, error
% bound, sigma ,epsilon, ro(samples outside the unit ball) for a given
% choice of parameters (see Eq. 25)
%
% Input - c: Bandlimit
%         L : Sampling frequency
%         T : Truncation Parameter
%         sig : sigma
%         mu1,mu2,mu3 : Coordinates of the Gaussian mean
%         Smax : Multiple of L (how far to sample ro outside the unit Ball
%
% Output - err : Approximation error
%          bound : Approximation bound
%          delt : Delta
%          eps : Epsilon
%          ro : Ro (sample of Gaussian outside the ball)

matdim = 400;
minEigenvalRatio = 10^(-16);
tol = 10^(-13);

%%%% Part 1 - Delta %%%%
delt = sqrt(((2*pi)^(-3))* (((pi/sig)^(3/2))*erfc(c*sqrt(sig)) + (2*pi/sig)*c*exp(-sig*c^2)));


%%%% Part 2 - Epsilon %%%%

normu = sqrt(mu1^2 + mu2^2 + mu3^2);
eps = sqrt(((pi*sig)^(3/2))*erfc((1-normu)/sqrt(sig)) + (2*pi*sig)*(1-normu)*exp(-(1/sig)*(1-normu)^2));

%%%% Part 3 - Ro %%%%

outerpoints = GenSamplePoints3D_AllFigs( Smax,L );
fVals = GaussianSamp_AllFigs( outerpoints,sig,mu1,mu2,mu3 );
ro = sqrt(sum(abs(fVals).^2));


%%%% Part 4 - Bound %%%%
eta = (4*pi/3)*(c*L)^(3/2);
bound = (eps+delt)*T + (eta/L^3)*ro + 4*delt;


%%%% Part 5 - Approximation Error %%%%

CartGrid = GenSamplePoints3D_V2(L);
[quadpoints,quadweights]=GenGridPoints_AllFigs(36,36);
Gaussquad = GaussianSamp_AllFigs( quadpoints,sig,mu1,mu2,mu3 );
GaussSamp = GaussianSamp_AllFigs( CartGrid,sig,mu1,mu2,mu3 );


[ SampleSeriesCoeff,truncvec,Proldat ] = ConstructSeriesCoeff( T,L,GaussSamp,CartGrid,c, minEigenvalRatio, matdim );

SeriesEval = EvalSeries( quadpoints,SampleSeriesCoeff,truncvec,Proldat,L );

norm = 0;

for j = 1:length(quadweights)

    norm = norm + quadweights(j,1)*abs(SeriesEval(1,j)-Gaussquad(j,1))^2;
    
end

err = sqrt(norm);

end

