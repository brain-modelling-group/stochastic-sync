function R = calc_orderParameter(relevant_data)
% Calculating general order parameter function 
%
% Inputs: relevant_data  : node data [mxT]
% Output: R              : order parameter calculation
%
% Original: James Pang, QIMR Berghofer, 2019

%% Main code
R = abs(mean(exp(1i*relevant_data)));
