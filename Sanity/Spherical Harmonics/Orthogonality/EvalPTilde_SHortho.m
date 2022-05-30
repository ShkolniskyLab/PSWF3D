function [ Val ] = EvalPTilde_SHortho( N,m,x )
% This function evaluates \tilde{P}_{N,m} at the point x. The code is based
% on the desctiprion in 'Numerical Recipies'.
%
%   Inputs - N : Degree of Legendre polynomials  ( N>=0 )
%            m : Degree of derivative ( N>=m>=0 )
%            x : Evaluation Point
%
%   Output - Val : The value \tilde{P}_{N,m}(x)

if m>N
    
    Val = zeros(size(x));
    
else
    
    Val0 = ((-1).^m).* sqrt((2*m+1)/(4*pi*factorial(2*m))).*prod(1:2:(2*m-1))*sqrt((1-x.^2).^m);
    Val1 = x.*sqrt(2*m+3).*Val0;
    
    if m==N
        
        Val = Val0;
        return
        
    elseif N ==(m+1)
        
        Val = Val1;
        return
        
    else
       
        steps = N-(m+1);
        m0 = m;
        m1 = m+1;
        
        for i = 1:steps
            
            Nrec = m1+1;
            ValTemp = sqrt((4*(Nrec^2)-1)/(Nrec^2-m^2)).*(x.*Val1 - sqrt((((Nrec-1)^2)-m^2)/(4*((Nrec-1)^2)-1)).*Val0);
            Val0 = Val1;
            Val1 = ValTemp ;
            m0 = m0+1;
            m1 = m1+1;
            
        end
        
        Val = Val1;
        
    end
        
    
end

end

