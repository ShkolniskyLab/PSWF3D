function [ fVals ] = GenEvalsSombrero_Example( SamplePoints,c )
% This function generates evaluations of the sombrero function given by 
% h_c(x) = (c^(3/2))/(sqrt(2)*pi^(3/2))*((J_{3/2}(c||x||))/(||x||^(3/2))).
% The computation appears in our paper
%   
% Input - SamplePoints : Vector of points to be sampled;
%         c : Bandlimit
%
% Output - fVals : Evaluations of the sombrero function on the points

fVals = zeros(length(SamplePoints(:,1)),1);

for j = 1:length(SamplePoints(:,1))
    
    norm2 = norm(SamplePoints(j,:));
    
    if norm2 == 0
        
        fVals(j,1) = (c^3)/(6*(pi^2));
        
    else
        
        coeff = (1/sqrt(2))*(c/(pi*norm2))^(3/2) ;
        
        bess = sqrt(2/(pi*c*norm2))*((sin(c*norm2)/(c*norm2))-cos(c*norm2));
        
        fVals(j,1) = coeff*bess;
    
    end
    
end

end

