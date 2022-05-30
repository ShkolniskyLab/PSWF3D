function [ fVals ] = GenEvalsIndic_Example(  SamplePoints,c )
% This function evaluates the inverse trandform of the indicator function of the cube
% given by [-c/3,c/3]^3 on the sampling points.
%
% Input - SamplePoints : Vector of points to be sampled;
%         c : Bandlimit
%
% Output - fVals : Evaluations of the function on the points

fVals = zeros(length(SamplePoints(:,1)),1);

for j=1:length(SamplePoints)
    
    if SamplePoints(j,1) == 0
        
       m1 = c/(3*pi);
                
    else
        
       m1 = sin((c*SamplePoints(j,1))/3)/SamplePoints(j,1);
        
    end
    
    if SamplePoints(j,2) == 0
      
        m2 = c/(3*pi) ;
        
    else
        
        m2 = sin((c*SamplePoints(j,2))/3)/SamplePoints(j,2);
        
    end
    
    if SamplePoints(j,3) == 0
      
        m3 = c/(3*pi) ;
        
    else
        
        m3 = sin((c*SamplePoints(j,3))/3)/SamplePoints(j,3);
        
    end
    
    fVals(j,:) = m1*m2*m3 ;
    
end


end

