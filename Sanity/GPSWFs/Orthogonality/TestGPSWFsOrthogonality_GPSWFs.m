function [  ] = TestGPSWFsOrthogonality_GPSWFs(  )
% This test verifies that the GPSWFs of degrees N1 and N2 are orthonormal on the
% unit sphere. The code computes an inner product table and gives one of
% the following
% 1. if N1=N2 - on-diagonal and off-diagonal maximal error and whether it is bigger
%    than the tolerance
% 2 . if N1 !=N2 - the maximal error and whether it is bigger than the
%     tolerance
%
% The values N1 and N2 are chosen in random.

c = pi;
matdim = 600;
phiquad = TrapzAtkinson_GPSWFs(60);
rquad = LegndreQuadRecipies_GPSWFs(30,0,1);
thetaquad = QuadRulePolarAtkinson_GPSWFs( 30 );
tol = 10^(-12);
N1=randi([0,5]);
N2=randi([0,5]);

InnerProdTable = GenInnerProdTable_GPSWFs( matdim,phiquad,rquad,thetaquad,N1,N2,c );

if N1==N2
    
    ondiagerror = max(max(abs(eye(length(InnerProdTable(:,1)))-diag(diag(InnerProdTable)))));
    offdiagerror =  max(max(abs(InnerProdTable-diag(diag(InnerProdTable)))));
    
    if (ondiagerror>tol)||(ondiagerror>tol)
        
        disp(strcat('N1 :',num2str(N1)));
        disp(strcat('N2 :',num2str(N2)));
        disp(strcat('On diagonal error :',num2str(ondiagerror)));
        disp(strcat('Off diagonal error :',num2str(offdiagerror)));
        disp('Test failed - Error is too large');
        
    else

        disp(strcat('N1 :',num2str(N1)));
        disp(strcat('N2 :',num2str(N2)));
        disp(strcat('On diagonal error :',num2str(ondiagerror)));
        disp(strcat('Off diagonal error :',num2str(offdiagerror)));
        disp('Test passed - Error is ok') ;       
        
    end
    
else
    
    totalerr = max(max(abs(InnerProdTable)));
    
    if totalerr>tol
       
       disp(strcat('N1 :',num2str(N1)));
       disp(strcat('N2 :',num2str(N2)));
       disp(strcat('Maximum error is :',num2str(totalerr))); 
       disp('Test failed - Error is too large');
        
    else

       disp(strcat('N1 :',num2str(N1)));
       disp(strcat('N2 :',num2str(N2)));        
       disp(strcat('Maximum error is :',num2str(totalerr))); 
       disp('Test passed - Error is ok');        
        
    end
    
end

end

