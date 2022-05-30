function [ abserr,lenpsimat ] = GenGramianMatrix_AllFigs( T,L,c,minEigenvalRatio,matdim,CartGridPolar)

% This function computes the maximal deviation of eigenvalues of H_c from 1
% for a given parameter of T.
%
% Input - T : Truncation Parameter
%         L : Sampling rate
%         c : Bandlimit
%         minEigenvalRatio,matdim : Parameters from Lederman's code
%         CartGridPolar : Cartesian Grid in Polar representation
%
% Output - abserr : maximum deviation of eigenvalues from 1
%          lenpsimat : size of the matrix \Psi

%%%%%%%%%%%%%%%% Part1 - Construct Proldat and truncvec %%%%%%%%%%%%%%%%%%%%%%
Proldat = cell(3,1);
alphfact = (c/(2*pi))^3;
[prolate_dat, iserr , prolate_dat_tmp] = prolate_crea(c, 3, 0, minEigenvalRatio, matdim , 1); 
modalph = alphfact*abs(prolate_dat.alp).^2;
modalph = abs(sqrt(modalph./(1-modalph)));
initpos = find(modalph>T,1,'last');
truncvec = [initpos];

Proldat{1,1} = prolate_dat;
Proldat{2,1} = iserr;
Proldat{3,1} = prolate_dat_tmp;

counter = 1;

while (counter <= 120)
    
    [prolate_dat, iserr , prolate_dat_tmp] = prolate_crea(c, 3, counter, minEigenvalRatio, matdim , 1); 
    
    modalph = alphfact*abs(prolate_dat.alp).^2;
    modalph = abs(sqrt(modalph./(1-modalph)));
    initpos = find(modalph>T,1,'last');
    
    if isempty(initpos) == 1
        
        break
        
    else
      
        Proldat{1,counter+1} = prolate_dat;
        Proldat{2,counter+1} = iserr;
        Proldat{3,counter+1} = prolate_dat_tmp;
        truncvec = horzcat(truncvec,[initpos]);
        counter= counter+1;
    end
      
end

Nmax = length(truncvec) - 1;
assert(Nmax>0,'Truncation parameter is too large ');

%%%%%%%%%%%%%%%% Part2 - Evaluate the GPSWFs on the cartesian grid %%%%%%%%%%%%%%%%%%%%%%
[ ProlatesEval,alphavect ] = EvalProlGrid_fixedN( CartGridPolar, 0,Proldat{1,1}, Proldat{2,1} , Proldat{3,1},truncvec(1,1),0,L );
Psimat = (1/sqrt(8))*bsxfun(@times,ProlatesEval,alphavect);

for i = 2: min([length(truncvec),7])
    
    N = i-1;
    [ ProlatesEval,alphavect ] = EvalProlGrid_fixedN( CartGridPolar, N,Proldat{1,i}, Proldat{2,i} , Proldat{3,i},truncvec(1,i),0,L );
    Psimat = vertcat(Psimat,(1/sqrt(8))*bsxfun(@times,ProlatesEval,alphavect));
    
end

lenpsimat = length(Psimat(:,1));
eigv = eig(Psimat*Psimat');

abserr = max(abs(eigv-1));

end

