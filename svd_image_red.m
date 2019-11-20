k = 10;
A = imread('pout.tif');
A = double(A);
[U,S,V] = svd(A);
X = svd_rank_red(U,S,V,k);
figure(1);
imshow(uint8(A));
figure(2);
imshow(uint8(X));