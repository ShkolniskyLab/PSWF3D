function [ ProlatesEval,alphavect ] = EvalProlGrid_fixedN( PointsPolar, N,prolate_dat, iserr , prolate_dat_tmp,truncationparam,flag,L )
% This fuction evaluates a set of GPSWFs on a list of points
%
%   Input - PointsPolar : A list of point to be sampled in spherical coordinates
%           N : Degree of GPSWFs to be evaluated.
%           prolate_dat, iserr , prolate_dat_tmp : Corresponding data from
%                                                  Lederman's code.
%           truncationparam : Truncation index of radial GPSWFs.
%           flag : 0 - if the code is used for construction of coefficients
%                  1 - if the code is used for approximating a value.
%           L : Sampling parameter
%
%   Output - ProlatesEval : Evaluations of the prolates. The result is a
%                           matrix where each row corresponds to a
%                           different prolate, each column corresponds to a
%                           different point
%            alphavect : Vector of corresponding eigenvalues

%%%%% Part 1 - Evaluate Radial GPSWFs on the radial components %%%%%%%%
assert(or(flag==0,flag==1),'Invalid flag');
radindx = [0:truncationparam-1];
% Generate an array of indices of all radial GPSWFs to be evaluated.

RadProlEv = prolate_ev(prolate_dat, radindx, PointsPolar(:,3));
% Evaluate the radial GPSWFs on the radial components of the points.


%%%%%% Part 2 - Evaluate Spherical Harmonics %%%%%%%%%%%%%%%%%%%%%%%%%%%%
SphHarmEval = zeros(length(PointsPolar(:,1)),2*N+1);
Basis = GenSHBasis(N);
if N > 0
    Basis = Basis(length(Basis)-2*N:length(Basis),:);
end
% Prepare array for evaluations and generate indices of desired spherical
% harmonics

if flag == 1
    
    loc = find(Basis(:,2)==0);

    for j = 1:loc
      
            if j==loc  
           
            SphHarmEval(:,j) = EvalSpherHarmPoint(N,Basis(j,2),PointsPolar(:,1),PointsPolar(:,2));
            
            else
           
            SphHarmEval(:,j) = EvalSpherHarmPoint(N,Basis(j,2),PointsPolar(:,1),PointsPolar(:,2));
            SphHarmEval(:,length(Basis)-j+1) = (-1)^(abs(Basis(j,2))).*conj(SphHarmEval(:,j));
           
            end
          
    end
    
else
    
    loc = find(Basis(:,2)==0); 
    
    for j = 1:loc
        
                
            if j==loc
                
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j) = EvalSpherHarmPoint(N,Basis(j,2),PointsPolar(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,1),PointsPolar(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,2));
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-2,j) = ((1j)^(Basis(j,2))).* SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j);
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-1,j) =((-1)^(Basis(j,2))).* SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j);
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)],j) = ((-1j)^(Basis(j,2))).* SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j) ;
                
                
            else

                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j) = EvalSpherHarmPoint(N,Basis(j,2),PointsPolar(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,1),PointsPolar(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,2));
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,length(Basis)-j+1) = (-1)^(abs(Basis(j,2))).*conj(SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j));
                
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-2,j) = ((1j)^(Basis(j,2))).* SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j);
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-2,length(Basis)-j+1) = (-1)^(abs(Basis(j,2))).*conj(SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-2,j));
                
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-1,j) =((-1)^(Basis(j,2))).* SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j);
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-1,length(Basis)-j+1) = (-1)^(abs(Basis(j,2))).*conj(SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-1,j));
                
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)],j) = ((-1j)^(Basis(j,2))).* SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)]-3,j) ;
                SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)],length(Basis)-j+1) = (-1)^(abs(Basis(j,2))).*conj(SphHarmEval(4.*[1:((length(PointsPolar)-(2*L-1))/4)],j));
                
            end
        
        
    end
    
    
    for j = 1:loc
        
        for k = (length(PointsPolar)-(2*L-1)+1):length(PointsPolar)
        
            if j==loc  
           
            SphHarmEval(k,j) = EvalSpherHarmPoint(N,Basis(j,2),PointsPolar(k,1),PointsPolar(k,2));
            
            else
             
            SphHarmEval(k,j) = EvalSpherHarmPoint(N,Basis(j,2),PointsPolar(k,1),PointsPolar(k,2));
            SphHarmEval(k,length(Basis)-j+1) = (-1)^(abs(Basis(j,2)))*conj(SphHarmEval(k,j));    
                
            end
            
        end
    end
    
end
% The evaluations are stored in SphHarmEval. Every column corresponds to a different spherical
% harmonic and every row corresponds to a different point.


%%%%% Part 3 - Construct GPSWFs Evaluations %%%%%%%%%%%%%%%%%%%%%%
numprols = truncationparam*(2*N+1);
ProlatesEval = zeros(numprols,length(PointsPolar(:,1)));
alphavect = zeros(numprols,1);
% The rows correspond to the GPSWFs and the columns to the evaluation points.
% The rows are ordered according to (1st radial prolate|2nd radial
% prolate|...|last radial prolate)^t


for j = 1:truncationparam
    
    indst = (j-1)*(2*N+1)+1;
    indend = indst+2*N;
    ProlatesEval(indst:indend,:) = bsxfun(@times,transpose(SphHarmEval),transpose(RadProlEv(:,j)));
    alphavect(indst:indend,1) = prolate_dat.alp(j);
    
     
end

end

