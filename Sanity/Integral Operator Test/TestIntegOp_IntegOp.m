function [  ] = TestIntegOp_IntegOp(  )

% This function tests the integral operator on a grid of points. The
% function computes the difference between the RHS and LHS and compares it
% to a given error tolerance.

%%%% Part 1 - Prepare parameters %%%%
tol = 10^(-13);
c= pi;
D = 3;
N = 1;
minEigenvalRatio = 10^(-16);
matdim = 300;
Samples = [0.01:0.01:1];
[prolate_dat, iserr , prolate_dat_tmp] = prolate_crea(c, D, N, minEigenvalRatio, matdim , 1);

%%%% Part 2 - Compute Legendre rule %%%%
LegendreQuad = LegndreQuadRecipies_IntegOp( 200,0,1 );

%%%% Part 3 - Compute LHS in eq (31) in Lederman's code %%%%
evalLHS = bsxfun(@times,bsxfun(@times,prolate_dat.gam,prolate_ev(prolate_dat, [0:prolate_dat.num_prols-1], Samples)),transpose(Samples));

%%%% Part 4 - Compute the RHS is eq (31) in Lederman's code %%%%
evalRHS = zeros(size(evalLHS));
phieval = bsxfun(@times,transpose(LegendreQuad(1,:)),prolate_ev(prolate_dat, [0:prolate_dat.num_prols-1], LegendreQuad(1,:)));
sq = sqrt(c.*transpose(Samples)*LegendreQuad(1,:));
j = besselj(3/2,c.*transpose(Samples)*LegendreQuad(1,:));

for i = 1:prolate_dat.num_prols
    
        evalRHS(:,i) = (sq.*j)*bsxfun(@times,transpose(LegendreQuad(2,:)),phieval(:,i));
    
end

K = max(max(abs(evalLHS-evalRHS)));

if K>=tol
    
    disp(strcat('Maximum error is :',num2str(K)));
    disp('Test failed - error is too large');
    
else
    
    disp(strcat('Maximum error is :',num2str(K)));
    disp('Test passed - error is ok');
    
end


end

