function [ InnerProdTable ] = GenInnerProdTable_GPSWFs( matdim,phiquad,rquad,thetaquad,N1,N2,c )
% This function generates an inner product table of all GPSWFs of degrees
% N1 and N2 which correspond to bandlimit c. The function uses the
% EvalProlGrid_fixedN.m function.
%
% Input - matname : Name of mat file with quadrature&construction
%                   parameters(attached)
%         N1,N2 : Degrees of desired GPSWFs
%         c : bandlimit
%
% Output - InnerProdTable : Table of inner products.


[ eval1,eval2,weights ] = Integ3D_data_GPSWFs( N1,N2,rquad,phiquad,thetaquad,c,matdim );


InnerProdTable = IntegGPSWFsTable_GPSWFs( eval1,eval2,weights );

end

