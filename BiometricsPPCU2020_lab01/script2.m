close all;
clear;
img = rgb2gray(imread('bird.png'));

%original image
H = zeros(1, 256);
for intensity=0:255
    mask = (img==intensity); %2d array
    H(intensity+1) = sum(mask(:)); %vectorize it
end

figure;
subplot(2,2,1);
imshow(img);
title('original grayscale image');

subplot(2,2,3);
title('histogram of the original image');
hold on;
bar(0:255,H);

%stretched image

img_max = double(max(img(:)));%vectorize with (:)
img_min = double(min(img(:)));
img_double = double(img);
img_scaled = (255/(img_max-img_min)) * (img_double-img_min);
img_scaled = uint8(img_scaled);

subplot(2,2,2);
imshow(img_scaled);
title('stretched image');

H_scaled = zeros(1, 256);
for intensity=0:255
    mask = (img_scaled==intensity); %2d array
    H_scaled(intensity+1) = sum(mask(:)); %vectorize it
end

subplot(2,2,4);
title('stretched histogram');
hold on;
bar(0:255,H_scaled);