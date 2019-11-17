function [A] = svd_rank_red(U,S,V,max_rank)
%SVD_RANK_RED Summary of this function goes here
%   Detailed explanation goes here
    S(max_rank+1:end,max_rank+1:end) = 0;
    A = U*S*V';
end

