function [ SampleSeriesCoeff,truncvec,Proldat ] = ConstructSeriesCoeff_Example( T,L,fVals,SamplingPoints,c, minEigenvalRatio, matdim )
% This function computes the coefficients of the sampling series and the
% truncation vector which corresponds to the parameter T
%
%Input - T : Truncation parameter
%        L : Sampling parameter
%        fVals : Values of function on the Cartesian grid
%        Sampling Points : Cartesian grid
%        c : Bandlimit
%        minEigenvalRatio,matdim : Parameters from Lederman's code
%
%Output - SampleSeriesCoeff: Array of coefficients of the sampling series
%                            ordered according to N
%         truncvec : Number of radial GPSWFs that correspond to each value
%                    of N between 0 and Nmax
%         Proldat : Data array for radial GPSWFs evaluation (generated by
%                   Lederman's code)

%%%%%%%%%%%%%%%% Part1 - Transform points into polar coordinates%%%%%%%%%%%%%%%%%%%%%%
[azimuth,elevation,radius] = cart2sph(SamplingPoints(:,1),SamplingPoints(:,2),SamplingPoints(:,3));
SamplingPointsPolar = horzcat(azimuth,elevation,radius);
SamplingPointsPolar(:,2) = (-1)*SamplingPointsPolar(:,2)+pi/2;

%%%%%%%%%%%%%%%% Part2 - Construct Proldat and truncvec %%%%%%%%%%%%%%%%%%%%%%
Proldat = cell(3,1);
alphfact = (c/(2*pi))^3;
[prolate_dat, iserr , prolate_dat_tmp] = prolate_crea_Example(c, 3, 0, minEigenvalRatio, matdim , 1); 
modalph = alphfact*abs(prolate_dat.alp).^2;
modalph = abs(sqrt(modalph./(1-modalph)));
initpos = find(modalph>T,1,'last');
truncvec = [initpos];

Proldat{1,1} = prolate_dat;
Proldat{2,1} = iserr;
Proldat{3,1} = prolate_dat_tmp;

counter = 1;

while (counter <= 550)
    
    [prolate_dat, iserr , prolate_dat_tmp] = prolate_crea_Example(c, 3, counter, minEigenvalRatio, matdim , 1); 
    
    modalph = alphfact*abs(prolate_dat.alp).^2;
    modalph = abs(sqrt(modalph./(1-modalph)));
    initpos = find(modalph>T,1,'last');
    
    if isempty(initpos) == 1
        
        break
        
    else
      
        Proldat{1,counter+1} = prolate_dat;
        Proldat{2,counter+1} = iserr;
        Proldat{3,counter+1} = prolate_dat_tmp;
        truncvec = horzcat(truncvec,[initpos]);
        counter= counter+1;
    end
      
end

Nmax = length(truncvec) - 1;
assert(Nmax>0,'Truncation parameter is too large ');

%%%%%%%%%%%%%%%% Part3 - Construct Coefficients %%%%%%%%%%%%%%%%%%%%
SampleSeriesCoeff = SamplingCoeff_AllN_Example( L,fVals,SamplingPointsPolar,Nmax,Proldat,truncvec,0);
end

