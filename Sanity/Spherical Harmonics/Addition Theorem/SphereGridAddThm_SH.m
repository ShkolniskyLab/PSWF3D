function [ MaxAddThm ] = SphereGridAddThm_SH( p,S )
% This function verifies the addition theorem for a given set of Spherical Harmonics
% for a grid on the unit sphere. The function calls the EvalSpherHarmPoint
% function to evaluate the Spherical Harmonics.
%
%   Input  - p : Number of sub intervals in phi(azimuth),theta(elevation)
%            S : Degree of desired Spherical harmonics
%   
%   Output - MaxAddThm : Maximum error between the RHS and LHS in the 
%                        addition theorem over all grid points

phiarr = [0:((2*pi)/p):(2*pi)];
thetarr = [0:(pi/p):pi];

points = zeros((p+1)^2,2);

for j = 1:length(phiarr)
    
   for m = 1:length(thetarr)
       
       points((p+1)*(j-1)+m,:) = [phiarr(1,j),thetarr(1,m)];
       
   end  
    
end

AddRes = zeros(length(points(:,1)));

for j = 1:length(points(:,1))
    
   for m = 1:j
       
       AddRes(j,m) = AddThmArr_SH(S,points(j,1),points(j,2),points(m,1),points(m,2)); 
       
   end
    
end

MaxAddThm = max(max(AddRes));
   

end

