%1,

fprintf('some size: %d x %d x %d\n', size(I, 1), size(I, 2), size(I, 3));
fprintf('type: uint8: %d, double: %d\n', isa(I, 'uint8'), isa(I, 'double'));
fprintf('Extrema: min=%6.3f, max=%d\n', min(I(:)), max(I(:)));

%subplot:
figure;
subplot(2, 3, nmbr); %row-vise
imshow(I);%for the images

subplot(2, 3, 4);%for the plot 
plot( 28:48, I(26, 28:48, 1), 'r' );
hold on;

%2 
I = rgb2gray(imread('asdf.png'));
%histogram is a distribution
H = zeros(1, 256);
for intensity=0:255
    mask = (I==intensity);
    H(intensity+1) = sum(mask(:)); %vectorize it
end

%at the formula we have to convert the number to double
I_max = double(max(I(:)));%vectorize with (:) try it at home
I_min = double(min(I(:)));
I_double = double(I);
I_scaled = (255/(I_max-I_min)) * (I_double-I_min);

%3
I, k

%paddings: zero padding, donut, mirroring
%3x3 kernel esetén elég egy padding sor, 5x5-ösnél 2 additional row and
%column is
%s = size(kernel, 1)
%s2 = floor(s/2)

s = size(k, 1);
s2 = floor(s/2);%alsó egészrész
I_padded = padarray(I, [s2, s2], 0, 'both');
[rn, cn] = size(I);

I_result = zeros(size(I));
for rows=s2+1:s2+rn %
    for cols=s2+1:s2+cn
        I_segment = I_padded(rows-s2:rows+s2, cols-s2:cols+s2);
        I_result(rows-s2, cols-s2) = sum(sum(I_segment.*k_rotated));
    end
end

%4






