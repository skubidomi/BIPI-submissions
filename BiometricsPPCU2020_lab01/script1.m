close all;
clear;
figure;

%original image
img = imresize(imread('AlfredoBorba_TuscanLandscape.jpg'), 0.1);
fprintf('Size of the image: %d x %d x %d\n', size(img, 1), size(img, 2), size(img, 3));
fprintf('Type of data in the image matrix: uint8: %d, double: %d\n', isa(img, 'uint8'), isa(img, 'double')');
fprintf('Extrema: min=%d, max=%d\n###\n', min(img(:)), max(img(:)));

subplot(2,3,1);
imshow(img);
title('original, color');

%gray image
img_gray = rgb2gray(img);
fprintf('Size of the image: %d x %d x %d\n', size(img_gray, 1), size(img_gray, 2), size(img_gray, 3));
fprintf('Type of data in the image matrix: uint8: %d, double: %d\n', isa(img_gray, 'uint8'), isa(img_gray, 'double')');
fprintf('Extrema: min=%d, max=%d\n###\n', min(img_gray(:)), max(img_gray(:)));

subplot(2,3,2);
imshow(img_gray);
title('gray with rgb2gray');

%double type image
img_double = im2double(img);
fprintf('Size of the image: %d x %d x %d\n', size(img_double, 1), size(img_double, 2), size(img_double, 3));
fprintf('Type of data in the image matrix: uint8: %d, double: %d\n', isa(img_double, 'uint8'), isa(img_double, 'double')');
fprintf('Extrema: min=%d, max=%d\n###\n', min(img_double(:)), max(img_double(:)));

subplot(2,3,3);
imshow(img_double);
title('color with im2double');

%channel intensity
%original image
subplot(2,3,4);
hold on;
plot( 28:48, img(26, 28:48, 1), 'r' );
plot( 28:48, img(26, 28:48, 2), 'g' );
plot( 28:48, img(26, 28:48, 3), 'b' );
xlim([28 48]);
title('pixel channel-intensity values in the 26th row');

%gray image
subplot(2,3,5);
hold on;
plot( 28:48, img_gray(26, 28:48), 'k');
xlim([28 48]);
title('pixel intensity values in the 26th row');

%double type image
subplot(2,3,6);
hold on;
plot( 28:48, img_double(26, 28:48, 1), 'r' );
plot( 28:48, img_double(26, 28:48, 2), 'g' );
plot( 28:48, img_double(26, 28:48, 3), 'b' );
xlim([28 48]);
title('pixel channel-intensity values in the 26th row');
