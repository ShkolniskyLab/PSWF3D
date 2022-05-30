function [  ] = TestAddTheoremMain_SH(  )
% This function verifies the addition theorem for a set of spherical
% harmonics on a grid of spherical points. Both sides of the addition
% theorem are computed on the grid and the maximal error between the RHS
% and LHS is compared to error tolerance.

p = 20;
S = randi([0,12],1,1);
errortol = 10^(-12);

MaxAddThm = SphereGridAddThm_SH( p,S );

if MaxAddThm>errortol
   
    disp(strcat('Chosen degree of Spherical Harmonics :',num2str(S)));
    disp(strcat('The error is :',num2str(MaxAddThm)));
    disp('The error is too large - Test failed');
    
else
    
    disp(strcat('Chosen degree of Spherical Harmonics :',num2str(S)));
    disp(strcat('The error is :',num2str(MaxAddThm)));
    disp('The error is ok - Test passed');
    
end


end

