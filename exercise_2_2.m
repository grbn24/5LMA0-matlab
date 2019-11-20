[A] = imread('Siep_weiland.jpg');
k = 25;
use_index = false;
if use_index
    [A,cmap]=rgb2ind(A,0);
    A = double(A);
    [U,S,V] = svd(A);
    X = svd_rank_red(U,S,V,k);
    X = round(X);
    A = ind2rgb(A,cmap);
    X = ind2rgb(X,cmap);
else
    X = nan(size(A));
    for ii = 1:size(A,3)
        A_{ii} = double(A(:,:,ii));
        A(:,:,ii) = uint8(A_{ii});
        X_{ii} = svd_rank_red(A_{ii},k);
        X(:,:,ii) = X_{ii};
    end
    X = uint8(X);
end
figure(1);
imshow(A);
figure(2);
imshow(X)
