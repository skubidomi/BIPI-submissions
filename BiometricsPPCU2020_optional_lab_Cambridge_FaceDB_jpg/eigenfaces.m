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
K = 40;
U = U(:,1:K);






