function [  ] = GenGraphs_ErrBound2_AllFigs(  )
% This function Generates the graphs of the approximation error, error
% bound, sigma and epsilon for various choices of sigma and L=34 , c=
% \pi*L, T=1

%%%% Part 1 - Set Parameters %%%%
c = 52*pi;
L = 52;
T = 1;
mu1 = 0.1;
mu2 = 0.1;
mu3 = 0.1;
Smax = 4*L;

sigarray = [0.002:0.001:0.02];
restab = zeros(length(sigarray),5);

%%%% Part 2 - Generate the error,bound,delta,epsilon and ro(sum of samples outside the disk) %%%%

parfor m = 1:length(sigarray)
   
   sig = sigarray(1,m);
   [err,bound,delt,eps,ro] = CompVal_AllFigs( c,L,T,sig,mu1,mu2,mu3,Smax );
   restab(m,:) =  [err,bound,delt,eps,ro]; 
       
end

%%%% Part 3 - Generate Plot %%%%

figure;
plot(sigarray',log10(restab(:,1)),'b--o');
hold on
plot(sigarray,log10(restab(:,2)),'k-.*');
plot(sigarray,log10(restab(:,3)),'r-d');
plot(sigarray,log10(restab(:,4)),'m:+');
grid on
grid minor
xlim([0.002 0.02]);
ylim([-15 2]);
legend('Measured error, T=1','Error bound, T=1','\delta','\epsilon');
xlabel('\sigma');
ylabel('$$ \log_{10}(Error)$$','Interpreter','latex');
title('T=1 ; L=52 ; \mu=(0.1,0.1,0.1)');
drawnow
end

