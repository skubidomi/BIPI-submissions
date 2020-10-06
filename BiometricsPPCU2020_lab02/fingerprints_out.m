close all;
clear;
clc;

original_img = imread('102_6.tif');
original_img = im2double(original_img);
original_img = original_img(:, 1:end-2);
img = imgaussfilt(original_img).*0.1; % 0.4 % 0.1
ramp = linspace(0, 0.85, size(img, 2)); % 0.55 % 0.85
img = img + repmat(ramp, size(img, 1), 1);

%%{
figure(1);
subplot(1, 2, 1);
imshow(original_img);
title('original image');
subplot(1, 2, 2);
imshow(img);
title('noisy image, this should be processed')
%%}

% ---
% please do not edit the code above, write your solution under this line,
% please handle original_img as unknown
% ---

figure(2);
hold on;
img_gabor = zeros();
for i = 0:7
    [gb_r, gb_i] = get_gabor_kernel(i/8*pi, 3.8, 1.0, 8.5);
    out_img = apply_gabor_kernel(img, gb_r, gb_i);
    subplot(4, 4, i*2+1); %for the real parts
    imagesc(gb_r);
    title(strcat('real part of Gabor, theta=', num2str(i/8*pi)));
    colormap('gray'), axis equal;
    subplot(4, 4, i*2+2); %for the filtered image
    imagesc(out_img);
    title('result of the specific filtering');
    colormap('gray'), axis equal;
    
    img_gabor = img_gabor+out_img;%sum up
end

figure(3);
subplot(1, 2, 1);
imshow(img);
title('noisy input');
subplot(1, 2, 2);
img_gabor = img_gabor/max(img_gabor(:));
imshow(img_gabor);
title('Gabor filtered, normalized');

%3/5
FT = fftshift(fft2(img_gabor));
%3/6
M_FT = log(abs(FT)+1);
M_FT = M_FT/max(M_FT(:));
%3/7
FT(160:180, 160:180) = 0;
I = abs(ifft2(ifftshift(FT)));
%3/8
I_m = imbinarize(I,0.1);
I_m(:, 1:7) = 0;
I_m(:, end-6:end) = 0;
%3/9
I_m_2 = imdilate(I_m, strel('disk', 11, 8));
%3/10
img_gabor_masked = img_gabor + ~I_m_2;
%3/11
bradley_thresholded = bradley(img_gabor_masked); %default values will be ok.
%3/12
bradley_eroded = bradley_thresholded + ~imerode(I_m_2, strel('disk', 11, 8));

%3/13
figure(4);
subplot(2, 3, 1);
imshow(M_FT);
title('Magnitude of FFT');
subplot(2, 3, 2);
imshow(I_m);
title('baseline mask');
subplot(2, 3, 3);
imshow(I_m_2);
title('dilated mask');
subplot(2, 3, 4);
imshow(img_gabor_masked);
title('masked version of the Gabor-filtered');
subplot(2, 3, 5);
imshow(bradley_thresholded);
title('after Bradley');
subplot(2, 3, 6);
imshow(bradley_eroded);
title('boundaries of Bradley corrected');
%3/14
bw = bradley_eroded;
filled = imfill(~bw, 'hole');
holes = filled & bw;
bigholes = bwareaopen(holes, 20);
smallholes = holes & ~bigholes;
small_holes_removed = bw & ~smallholes;

figure(5);
subplot(2, 3, 1);
imshow(bw);
title('boundaries of Bradley corrected');
subplot(2, 3, 2);
imshow(filled);
title('filled holes');
subplot(2, 3, 3);
imshow(holes);
title('holes');
subplot(2, 3, 4);
imshow(bigholes);
title('big holes');
subplot(2, 3, 5);
imshow(smallholes);
title('small holes');
subplot(2, 3, 6);
imshow(small_holes_removed);
title('small holes removed');

%3/15
skeletonized = ~bwmorph(~small_holes_removed, 'skel', Inf);
%3/16
cleaned = ~skeletonized;
cleaned_struct = regionprops(~skeletonized, 'Area', 'PixelIdxList');

for i = 1:numel(cleaned_struct)
    if cleaned_struct(i).Area < 20
        cleaned(cleaned_struct(i).PixelIdxList) = 0;
    end
end

figure(9);
imshow(cleaned);

%3/17
endpoints = bwmorph(cleaned, 'endpoints');
branchpoints = bwmorph(cleaned, 'branchpoints');
[endy, endx] = ind2sub(size(cleaned) ,find(endpoints));
[branchy, branchx] = ind2sub(size(cleaned) ,find(branchpoints));

%3/18
figure(6);
subplot(1, 3, 1);
imshow(skeletonized);
title('skeletonized image');
subplot(1, 3, 2);
imshow(cleaned);
hold on;
plot(endx, endy, 'ro');
plot(branchx, branchy, 'gs');
title('cleaned from small independent line segment');
legend('end-points', 'branch-points');
subplot(1, 3, 3);
imshow(cleaned);
hold on;

elementcounter = 0;
for i = 1:max(size(endx))
    logicalx = abs(endx-endx(i)) < 4;
    logicaly = abs(endy-endy(i)) < 4; % if the abs value is smaller than the threshold->some logical magic->and the solution of the multiplication is itself, it is alone, it is a real solution
    if sum(logicalx.*logicaly) == 1 %only one which is the center
        if elementcounter == 0
            elementcounter = 1;
            filteredendx = endx(i);
            filteredendy = endy(i);
        else
            filteredendx = [filteredendx endx(i)];
            filteredendy = [filteredendy endy(i)];
        end
    end
end
elementcounter = 0;
for i = 1:max(size(branchx))
    logicalx = abs(branchx-branchx(i)) < 4;
    logicaly = abs(branchy-branchy(i)) < 4;
    if sum(logicalx.*logicaly) == 1 %only one which is the center
        if elementcounter == 0
            elementcounter = 1;
            filteredbranchx = branchx(i);
            filteredbranchy = branchy(i);
        else
            filteredbranchx = [filteredbranchx branchx(i)];
            filteredbranchy = [filteredbranchy branchy(i)];
        end
    end
end
plot(filteredendx, filteredendy, 'ro');
plot(filteredbranchx, filteredbranchy, 'gs');
title('only the real end- and br-points');
legend('real end-points', 'real branch-points');



 