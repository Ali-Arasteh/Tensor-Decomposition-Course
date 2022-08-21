addpath('tensor_toolbox');
%% Question 1
% HOOI function
I = [5, 3, 2];
R = [3, 2, 2];
G = tensor(rand(R), R);
U = cell(1, length(R));
for i = 1:length(R)
    temp = rand(I(i), R(i));
    for j = 1:R(i)
        for k = 1:j-1
            temp(:, j) = temp(:, j)-temp(:, k)'*temp(:, j)*temp(:, k);
        end
        temp(:, j) = temp(:, j)/norm(temp(:, j));
    end
    U{i} = temp(:, 1:R(i));
end
T = ttm(G, U, 1:length(R));
[recovered_G, recovered_U] = HOOI(T, R);
recovered_T = ttm(recovered_G, recovered_U, 1:length(R));
Error = T-recovered_T;
frobenius_err = norm(Error(:));
%% Question 2
%% Loading dataset
I = [112, 92, 50];
T = zeros(I);
for i = 1:5
    for j = 1:10
        T(:, :, (i-1)*10+j) = imread(['ORL/s', num2str(i), '/', num2str(j), '.pgm']);
    end
end
T = tensor(T, I);
%% HOOI
R = [30, 25, 5];
[~, U] = HOOI(T, R);
figure();
hold on;
for i = 1:5
    plot(1:50, abs(U{3}(:, i)));
    xlim([1, 50]);
end
legend('Class 1', 'Class 2', 'Class 3', 'Class 4', 'Class 5');
title(['HOOI, R = [', num2str(R(1)), ', ', num2str(R(2)), ', ', num2str(R(3)), ']']);
%% Tucker ALS
R = [30, 25, 5];
D = tucker_als(T, R);
figure();
hold on;
for i = 1:5
    plot(1:50, abs(D.U{3}(:, i)));
    xlim([1, 50]);
end
legend('Class 1', 'Class 2', 'Class 3', 'Class 4', 'Class 5');
title(['Tucker ALS, R = [', num2str(R(1)), ', ', num2str(R(2)), ', ', num2str(R(3)), ']']);
%% NTD
R = [30, 25, 5];
opts.maxit = 1000;
opts.tol = 1e-6;
[U, ~, ~] = ntd(T, R, opts);
figure();
hold on;
for i = 1:5
    plot(1:50, abs(U{3}(:, i)));
    xlim([1, 50]);
end
legend('Class 1', 'Class 2', 'Class 3', 'Class 4', 'Class 5');
title(['NTD, R = [', num2str(R(1)), ', ', num2str(R(2)), ', ', num2str(R(3)), ']']);
%% Question 3
%% Part a)
%% Loading dataset
downsample_rate = 4;
I = [10, 9, (192/downsample_rate)*(168/downsample_rate)];
T = zeros(I);
folders = dir('Illumination_Yale');
folders(1:2) = [];
for i = 1:10
    files = dir(['Illumination_Yale/', folders(i).name]);
    files(1:2) = [];
    for j = 1:9
        temp = imresize(imread(['Illumination_Yale/', folders(i).name, '/', files(j).name]), 1/downsample_rate);
        T(i, j, :) = temp(:);
    end
end
T = tensor(T, I);
%% Tucker ALS
R1 = 10;
R2 = [1, 3, 5];
R3 = 90;
recovered_T = cell(1, length(R2));
i = 0;
for r2 = R2
    i = i+1;
    temp = tucker_als(T, [R1, r2, R3]);
    recovered_T{i} = ttm(temp.core, temp.U, 1:3);
end
%% Demonstration
for i = 1:10
    figure();
    for j = 1:4
        if j == 1
            temp = T.data;
        else
            temp = recovered_T{j-1}.data;
        end
        for k = 1:9
            subplot(4, 9, (j-1)*9+k);
            imshow(uint8(reshape(temp(i, k, :), 192/downsample_rate, 168/downsample_rate)));
        end
    end
    sgtitle(['Person ', num2str(i)]);
end
%% Part b)
%% Creating desired matrix
X = reshape(tenmat(T, 2).data, 10*9, (192/downsample_rate)*(168/downsample_rate))';
%% SVD
R = [10, 30, 50];
[U, S, V] = svd(X);
recovered_X = cell(1, length(R2));
i = 0;
for r = R
    i = i+1;
    recovered_X{i} = U(:, 1:r)*S(1:r, 1:r)*V(:, 1:r)';
end
%% Demonstration
for i = 1:10
    figure();
    for j = 1:4
        if j == 1
            temp = X;
        else
            temp = recovered_X{j-1};
        end
        for k = 1:9
            subplot(4, 9, (j-1)*9+k);
            imshow(uint8(reshape(temp(:, (i-1)*9+k), 192/downsample_rate, 168/downsample_rate)));
        end
    end
    sgtitle(['Person ', num2str(i)]);
end
%% Part c)
% in report
%% functions
function [G, U] = HOOI(T, R)
    number_of_itterations = 100;
    U = cell(1, length(R));
    for i = 1:length(R)
        Tn = tenmat(T, i);
        [Un, ~, ~] = svd(Tn.data);
        U{i} = Un(:, 1:R(i));
    end
    for i = 1:number_of_itterations
        for j = 1:length(R)
            temp_U = U;
            for k = 1:length(R)
                temp_U{k} = pinv(temp_U{k});
            end
            temp_U(j) = [];
            temp_R = 1:length(R);
            temp_R(j) = [];
            Z = ttm(T, temp_U, temp_R);
            Zn = tenmat(Z, j);
            [Un, ~, ~] = svd(Zn.data);
            U{j} = Un(:, 1:R(j));
        end
    end
    temp_U = U;
    for k = 1:length(R)
        temp_U{k} = pinv(temp_U{k});
    end
    G = ttm(T, temp_U, 1:length(R));
end