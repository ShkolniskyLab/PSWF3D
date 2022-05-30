function [ points,weights ] = GenGridPoints_AllFigs( thetaparam,rparam )
% This function computes the integration nodes and weights for the GPSWFs
% expansion evaluation
%  
%   Input - thetaparam : Number of nodes in theta
%           rparam : Number of nodes in r.
%   
%   Output - points : Evaluation Points for the quadrature
%            EvalWeights : Weights of the gaussian quadrature

phiparam = thetaparam*2;

phiquad = TrapzAtkinson_AllFigs(phiparam);
thetaquad = QuadRulePolarAtkinson_AllFigs(thetaparam);
rquad = LegndreQuadRecipies_AllFigs( rparam,0,1 );

numpoints = length(rquad(1,:))*length(thetaquad(1,:))*length(phiquad(1,:));
points = zeros(numpoints,3);
weights = zeros(numpoints,1);
counter = 1;

for j1 = 1:length(rquad(1,:))
    
   for j2 = 1:length(thetaquad(1,:))
       
      for j3 = 1:length(phiquad(1,:))
          
         [xcoor,ycoor,zcoor] = sph2cart(phiquad(1,j3),(-1)*thetaquad(1,j2)+(pi/2),rquad(1,j1));
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

end

