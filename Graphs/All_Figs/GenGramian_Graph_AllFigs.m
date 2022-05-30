function [  ] = GenGramian_Graph_AllFigs( )

% This function Generates the graph of the maximal deviation of the
% eigenvalues of \Psi from 1, as a function of the Truncation parameter
% (see Appendix B) for different values of The truncation Parameter T.

L = 30;
c = pi*L;
minEigenvalRatio = 10^(-16);
matdim = 300;
Tarray = [1,10,100,1000,10000,100000,1000000,10000000];
ErrArray = zeros(1,length(Tarray));
numprolarray = zeros(1,length(Tarray));

%%%%%%%%%%%%%%%% Part1 - Generate Cartesian Grid%%%%%%%%%%%%%%%%%%%%%%
CartGrid = GenSamplePoints3D_V2(L);

%%%%%%%%%%%%%%%% Part2 - Transform points into polar coordinates%%%%%%%%%%%%%%%%%%%%%%
[azimuth,elevation,radius] = cart2sph(CartGrid(:,1),CartGrid(:,2),CartGrid(:,3));
CartGridPolar = horzcat(azimuth,elevation,radius);
CartGridPolar(:,2) = (-1)*CartGridPolar(:,2)+pi/2;

%%%%%%%%%%%%%%%% Part3 - Compute results %%%%%%%%%%%%%%%%%%%%%%

parfor i = 1:length(Tarray)
   
   [ abserr,lenpsimat ] = GenGramianMatrix_AllFigs( Tarray(1,i),L,c,minEigenvalRatio,matdim,CartGridPolar);
   ErrArray(1,i) = abserr;
   numprolarray(1,i) = lenpsimat;
   disp(i);
   
end

m=polyfit(log10(Tarray),log10(ErrArray),1);

%%%%%%%%%%%%%%%% Part4 - Generate Graphs %%%%%%%%%%%%%%%%%%%%%%
figure;
plot(log10(Tarray),log10(ErrArray),'b--o');
hold on
plot(log10(Tarray),m(1,1).*log10(Tarray)+m(1,2),'r-s');
grid on
grid minor
xlim([0 7]);
ylim([-15 1]);
legend('Empirical',strcat('Fit - slope : ',num2str(m(1,1))));
xlabel('$$\log_{10}(T)$$','Interpreter','latex');
ylabel('$$ \log_{10}(\max_k |1-\tau_k|)$$','Interpreter','latex');
title(' L=30 ; c = \pi L ');
hold off;
drawnow

end

