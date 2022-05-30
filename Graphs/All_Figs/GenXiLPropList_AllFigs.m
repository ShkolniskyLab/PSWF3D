function [  ] = GenXiLPropList_AllFigs(  )

% This function generates the graph of the quotient of the energy of the Gaussian on a
% ball of radius r_1 and and L^6 (see eq A5 in the appendix, under the
% assumption that c =\pi*L) for various values of L.

%%%% Part 1 - Generate data bases %%%%
[quadpoints,quadweights]=GenGridPoints2_AllFigs(36,36,0.95);

Larr = [8,12,16,20,24,28,32,36,40,44,48,52];
Propvec = zeros(1,length(Larr));


%%%% Part 2 - Compute quotients %%%%
for i = 1:length(Larr)
    
    Propvec(1,i) = GenXiLProp_AllFigs(Larr(1,i),quadpoints,quadweights);
    disp(i);
    
end


%%%% Part3 - Generate Plots %%%%
figure;
plot(log10(Larr),log10(Propvec),'b--o');
grid on
xlim([0.7 1.7]);
ylim([-1.2 0]);
xlabel('$$\log_{10}(L)$$','Interpreter','latex');
str = '$$\log_{10}(\frac{1}{c^6}\left|\left| \xi_c \right| \right|_{L^2(r_1R)}^2)$$';
ylabel(str,'Interpreter','latex');
drawnow
end

