function coscorr = calc_coscorr(relevant_data)
% Calculating average cosine correlation
%
% Input: relevant_data  : node data [NxT]
% Output: coscorr       : cosine correlation [NxN]
%
% Original: James Pang, QIMR Berghofer, 2019

%%
N = size(relevant_data, 1);
iterations = size(relevant_data, 2);

coscorr = zeros(N, N);

for t = 1:iterations
    coscorr(:,:) = coscorr + cos(relevant_data(:,t)'-relevant_data(:,t));
end
coscorr = coscorr/iterations;

end
