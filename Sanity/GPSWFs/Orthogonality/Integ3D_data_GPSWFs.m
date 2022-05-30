function [ eval1,eval2,weights ] = Integ3D_data_GPSWFs( N1,N2,rquad,phiquad,thetaquad,c,matdim )
% This function calculates the quadrature results for GPSWFs of degree
% N1 and N2
%
%   Input - N1,N2 : Indices of desired GPSWFs
%           rquad : Radial quadrature formula (row 1:nodes, row 2:weights)
%           phiquad : Azimuth quadrature (row 1:nodes, row 2:weights)
%           thetaquad : Elevation quadrature (row 1:nodes, row 2:weights)
%
%   Output - eval1,eval2 : evaluations of the GPSWFs on the quadrature
%                          points
%            weights : quadrature weights

%%%%%%%%% Step 1 - Translation of points into Cartesian Coord %%%%%%%%%
numpoints = length(rquad(1,:))*length(thetaquad(1,:))*length(phiquad(1,:));
points = zeros(numpoints,3);
weights = zeros(numpoints,1);
counter = 1;

for j1 = 1:length(rquad(1,:))
    
   for j2 = 1:length(thetaquad(1,:))
       
      for j3 = 1:length(phiquad(1,:))
          
         xcoor = phiquad(1,j3) ;
         ycoor = thetaquad(1,j2);
         zcoor = rquad(1,j1);
         points(counter,:) = [xcoor,ycoor,zcoor];
         % Transform points into Cartesian form. Needed for evaluation of
         % sphericl harmonics and radial GPSWFs
         
         weights(counter,1) =(rquad(1,j1)^2)*sin(thetaquad(1,j2))*phiquad(2,j3)*rquad(2,j1)*thetaquad(2,j2) ;
         % Generate the quadrature weights. Note that we multiply the
         % product of weights by the jacobian r^2 * sin(theta).
         
         counter = counter+1;
          
      end
       
       
   end
    
    
end

%%%%% Part 2 - Evaluate the GPSWFs on the quadrature nodes %%%%

minEigenvalRatio = 10^(-16);

[prolate_dat1, iserr1 , prolate_dat_tmp1] = prolate_crea(c, 3, N1, minEigenvalRatio, matdim , 1) ;
[prolate_dat2, iserr2 , prolate_dat_tmp2] = prolate_crea(c, 3, N2, minEigenvalRatio, matdim , 1) ;
truncationparam1 = prolate_dat1.num_prols;
truncationparam2 = prolate_dat2.num_prols;

eval1 = EvalProlGrid_fixedN( points, N1,prolate_dat1, iserr1 , prolate_dat_tmp1,truncationparam1,1,80);
eval2 = EvalProlGrid_fixedN( points, N2,prolate_dat2, iserr2 , prolate_dat_tmp2,truncationparam2,1,80);



