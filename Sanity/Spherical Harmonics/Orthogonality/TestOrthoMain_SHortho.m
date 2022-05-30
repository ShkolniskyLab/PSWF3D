function [  ] = TestOrthoMain_SHortho(  )

% This function verifies orthonormality of generated Spherical Harmonics.
% An inner product table is computed by Gaussian quadrature for Spherical
% Harmonics of degree less\equal to N and the error
% is compared to an error tolerance.

M = 28;
N = 6;
errortol = 10^(-13);

InnerPTable = GaussInnerPTable_SHortho( M,N );

k = length(InnerPTable(1,:));
ondiagerr = max(max(abs(eye(k)-diag(diag(InnerPTable)))));
offdiagerr = max(max(abs(InnerPTable-diag(diag(InnerPTable)))));

if (ondiagerr>errortol)||(offdiagerr>errortol)
    
    disp(strcat('On diag error :',num2str(ondiagerr)));
    disp(strcat('Off diag error :',num2str(offdiagerr)));
    disp('Test failed - Error is too large')
    
else
    
    disp(strcat('On diag error :',num2str(ondiagerr)));
    disp(strcat('Off diag error :',num2str(offdiagerr)));
    disp('Test passed - Error is ok')
    
end


end

