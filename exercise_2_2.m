[A] = imread('Siep_weiland.jpg');
[A,cmap]=rgb2ind(A,1024);
A = double(A);
k = 100;
[U,S,V] = svd(A);
S_x = S;
S_x(k+1:end,k+1:end) = 0;
X = U*S_x*V';
X = round(X);
A = ind2rgb(A,cmap);
X = ind2rgb(X,cmap);
figure(1);
imshow(A);
figure(2);
imshow(X)
