function [ IntegTable ] = IntegGPSWFsTable_GPSWFs( eval1,eval2,weights )
% This function generates an inner product table for GPSWFs evlauated in
% eval1 and eval2

numprols1 = length(eval1(:,1));
numprols2 = length(eval2(:,1));

IntegTable = zeros(numprols1,numprols2);

for i = 1:numprols1
        
        IntegTable(i,:) = transpose(eval2*transpose(bsxfun(@times,conj(eval1(i,:)),transpose(weights))));
            
end     


end

