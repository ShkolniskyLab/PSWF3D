function [ Points] = GenSamplePoints3D_AllFigs( Smax,L )

%This function generates a grid of samples outside of the unit ball for the
%Gaussian to be smapled on. 
%
% Input - Smax: multiple of L which determines 'how far' to sample
%         L : Sampling Parameter
%
% Output - Points : Grid of sampling points outside the ball

assert(floor(Smax)==Smax,'Smax is not an integer');
assert(floor(L)==L,'L is not an integer');
assert(Smax>=L,'The resulting points are inside the unit ball');
assert(L>0, 'The parameter L must be positive');

x = (-Smax/L):(1/L):(Smax/L);
y = (-Smax/L):(1/L):(Smax/L);
z = (-Smax/L):(1/L):(Smax/L);
[X,Y,Z] = meshgrid(x,y,z);
M = X.^2+Y.^2+Z.^2;
c = find(M>=1);
[dim1,dim2,dim3] = ind2sub(size(M),c);
Ind = horzcat(horzcat(dim1,dim2),dim3);
Points = zeros(size(Ind));

for n = 1:length(Ind(:,1))
    
    Points(n,:)= [X(Ind(n,1),Ind(n,2),Ind(n,3)),Y(Ind(n,1),Ind(n,2),Ind(n,3)),Z(Ind(n,1),Ind(n,2),Ind(n,3))];
    
end

end

