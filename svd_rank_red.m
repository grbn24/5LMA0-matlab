function [A] = svd_rank_red(varargin)
%SVD_RANK_RED Summary of this function goes here
%   Detailed explanation goes here
    switch numel(varargin)
        case 4
            U = varargin{1};
            S = varargin{2};
            V = varargin{3};
            max_rank = varargin{4};
        case 2
            [U,S,V] = svd(varargin{1});
            max_rank = varargin{2};
            
        otherwise
            error("What do you want from me?? ;_;");
    end
    S(max_rank+1:end,max_rank+1:end) = 0;
    A = U*S*V';
end

