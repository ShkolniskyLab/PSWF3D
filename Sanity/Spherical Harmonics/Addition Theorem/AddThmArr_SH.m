function [ AddRes ] = AddThmArr_SH( S,phi1,theta1,phi2,theta2 )
% This function verifies the addition theorem for the generated spherical
% harmonics
%
%   Input : S - Degree of spherical harmonics to be verified
%           phi1,phi2,theta1,theta2 - Two points on a sphere

dim = 2*S+1;
surf = 4*pi;
% Generate parameters

indices = zeros(2*S+1,2);
for m = (-S):S
    
   indices(m+S+1,:)=[S,m]; 
    
end
% Generate SphericalHarmonics indices

Rhs=0;

for m = 1:length(indices(:,1))
    
    Rhs = Rhs+EvalSpherHarmPoint(indices(m,1),indices(m,2),phi1,theta1)*conj(EvalSpherHarmPoint(indices(m,1),indices(m,2),phi2,theta2));
    
end

ksi = [cos(phi1)*sin(theta1),sin(phi1)*sin(theta1),cos(theta1)];
eta = [cos(phi2)*sin(theta2),sin(phi2)*sin(theta2),cos(theta2)];

AddRes = abs(Rhs - (dim/surf)*EvalLeg_SH(S,dot(ksi,eta)));

end

