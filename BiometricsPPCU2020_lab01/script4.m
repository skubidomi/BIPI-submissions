close all;
clear;

I = im2double(rgb2gray(imread('NotreDame.png')));
FT = fft2(I);
shifted_FT = fftshift(FT);
M_FT = log(abs(shifted_FT+1));

M_FT = M_FT./max(M_FT(:));

%4/6
figure;
subplot(2,4,1);
imshow(I);
title('original');
subplot(2,4,5);
imshow(M_FT);
title('magnitude');

%4/7
modified_FT_A = shifted_FT;
modified_FT_A(120:135,50:122) = 0;
modified_FT_A(120:135,136:200) = 0;
M_FT_A = log(abs(modified_FT_A+1));
M_FT_A = M_FT_A./max(M_FT_A(:));
I_A = abs(ifft2(ifftshift(modified_FT_A)));
subplot(2,4,2);
imshow(I_A);
title('result A');
subplot(2,4,6);
imshow(M_FT_A);
title('modified magnitude A');

%4/8
modified_FT_B = shifted_FT;
modified_FT_B(92:97,86:174) = 0;
modified_FT_B(158:164,91:172) = 0;
M_FT_B = log(abs(modified_FT_B+1));
M_FT_B = M_FT_B./max(M_FT_B(:));
I_B = abs(ifft2(ifftshift(modified_FT_B)));
subplot(2,4,3);
imshow(I_B);
title('result B');
subplot(2,4,7);
imshow(M_FT_B);
title('modified magnitude B');

%4/9
modified_FT_C = shifted_FT;
modified_FT_C(124:133,124:133) = 0;
M_FT_C = log(abs(modified_FT_C+1));
M_FT_C = M_FT_C./max(M_FT_C(:));
I_C = abs(ifft2(ifftshift(modified_FT_C)));

subplot(2,4,4);
imshow(I_C);
title('result C');
subplot(2,4,8);
imshow(M_FT_C);
title('modified magnitude C');

