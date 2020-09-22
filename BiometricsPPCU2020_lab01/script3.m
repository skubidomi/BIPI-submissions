close all;
clear;

img = imresize(imread('AlfredoBorba_TuscanLandscape.jpg'), 0.2);
img = rgb2gray(img);
img = im2double(img);

%kernel
K = [1 0 -1; 1 0 -1; 1 0 -1];

%zero padding
s = size(K, 1);
s2 = floor(s/2); %border size
img_padded = padarray(img, [s2, s2], 0, 'both');

%rotate the kernel
K = rot90(K,2);
[rn, cn] = size(img); %row and column size

img_result = zeros(size(img));
for rows=s2+1:s2+rn %
    for cols=s2+1:s2+cn
        img_segment = img_padded(rows-s2:rows+s2, cols-s2:cols+s2);
        img_result(rows-s2, cols-s2) = sum(sum(img_segment.*K));
    end
end

subplot(1,3,1);
imshow(img);
title('original image');

subplot(1,3,2);
imagesc(img_result);
colormap gray;
axis equal;
axis off;
colorbar;
title('with our convolution');

%built in
K_builtin = transpose(fspecial('prewitt'));
img_result_builtin = conv2(img,K_builtin,'same');
subplot(1,3,3);
imagesc(img_result_builtin);
colormap gray;
axis equal;
axis off;
colorbar;
title('with built-in conv2');

