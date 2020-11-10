% 1.
close all;
clear;
clc;

% 2.1, 2.2
Gamma = zeros(10304,400);
S = dir('*.jpg');
for i = 1:size(S,1)
    Gamma(:,i) = reshape(imread(S(i).name),[10304,1]);
end

% 3
test_img_column = round(1+399*rand());
test_img = Gamma(:,test_img_column);
Gamma(:,test_img_column) = [];

% 4
Psi = round(mean(Gamma,2));

% 5
A = Gamma - repmat(Psi,1,size(Gamma,2));

% 6.1
L = transpose(A)*A;

% 6.2
[V, D] = eig(L, 'vector');
[D, idx] = sort(D, 'descend');
V = V(:, idx);
Uraw = A*V;

% 6.3
U = Uraw./vecnorm(Uraw);

% 6.4
K = 320;
U = U(:,1:K);

% 7
Y = transpose(U)*A;

% 8
Phi = test_img - Psi;
y = transpose(U)*Phi;

% 9.1
Diff = zeros(1,size(Y,2));
for i = 1:size(Y,2)
    Diff(i) = norm(Y(:,i)-y);
end

% 9.2
[minimum_value, minimum_index] = min(Diff);

% 10
Phi_rec = U*Y(:,minimum_index);
recognised_test_img = Phi_rec + Psi;

% 11
figure(1);
subplot(1,3,1);
imshow(reshape(uint8(test_img),112,92));
subplot(1,3,2);
imshow(reshape(uint8(Psi),112,92));
subplot(1,3,3);
imshow(reshape(uint8(recognised_test_img),112,92));

%12
figure(2);
for i = 1:10
    subplot(2,5,i);
    Uimg = U-
    Uimg = (U(:,i)/max(U(:,i)))*256;
    imshow(reshape(uint8(Uimg),112,92));
end





