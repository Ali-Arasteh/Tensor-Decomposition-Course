%% Question 1
%% Part a)
% NMF_ALS function
%% Part b)
% NMF_Multiplicative function
%% Part c)
SNR = [-10, 0, 10, 30, 50];
Iter = 1:10;
J = [1, 2, 3, 4];
NMF_ALS_Error = zeros(length(SNR), length(Iter), length(J));
MatLab_NMF_ALS_Error = zeros(length(SNR), length(Iter), length(J));
NMF_Multiplicative_Error = zeros(length(SNR), length(Iter), length(J));
MatLab_NMF_Multiplicative_Error = zeros(length(SNR), length(Iter), length(J));
p = 0;
for snr = SNR
    p = p+1;
    q = 0;
    for iter = Iter
        q = q+1;
        B = rand(6, 3);
        C = rand(3, 4);
        E = rand(6, 4);
        alpha = norm(B*C, 'fro')/(10^(snr/20)*norm(E, 'fro'));
        A = B*C+alpha*E;
        r = 0;
        for j = J
            r = r+1;
            B0 = rand(6, j);
            C0 = rand(j, 4);
            [B_temp, C_temp] = NMF_ALS(A, j, B0, C0);
            NMF_ALS_Error(p, q, r) = norm(A-B_temp*C_temp, 'fro');
            [B_temp, C_temp] = nnmf(A, j, 'algorithm', 'als', 'w0', B0, 'h0', C0);
            MatLab_NMF_ALS_Error(p, q, r) = norm(A-B_temp*C_temp, 'fro');
            [B_temp, C_temp] = NMF_Multiplicative(A, j, B0, C0);
            NMF_Multiplicative_Error(p, q, r) = norm(A-B_temp*C_temp, 'fro');
            [B_temp, C_temp] = nnmf(A, j, 'algorithm', 'mult', 'w0', B0, 'h0', C0);
            MatLab_NMF_Multiplicative_Error(p, q, r) = norm(A-B_temp*C_temp, 'fro');
        end
    end
end
figure();
subplot(2, 2, 1);
hold on;
plot(SNR, mean(NMF_ALS_Error(:, :, 1), 2));
plot(SNR, mean(NMF_ALS_Error(:, :, 2), 2));
plot(SNR, mean(NMF_ALS_Error(:, :, 3), 2));
plot(SNR, mean(NMF_ALS_Error(:, :, 4), 2));
grid on;
legend('j = 1', 'j = 2', 'j = 3', 'j = 4');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('NMF ALS Algorithm');
subplot(2, 2, 2);
hold on;
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 1), 2));
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 2), 2));
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 3), 2));
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 4), 2));
grid on;
legend('j = 1', 'j = 2', 'j = 3', 'j = 4');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('MatLab NMF ALS Algorithm');
subplot(2, 2, 3);
hold on;
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 1), 2));
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 2), 2));
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 3), 2));
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 4), 2));
grid on;
legend('j = 1', 'j = 2', 'j = 3', 'j = 4');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('NMF Multiplicative Algorithm');
subplot(2, 2, 4);
hold on;
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 1), 2));
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 2), 2));
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 3), 2));
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 4), 2));
grid on;
legend('j = 1', 'j = 2', 'j = 3', 'j = 4');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('MatLab NMF Multiplicative Algorithm');
figure();
subplot(2, 2, 1);
hold on;
plot(SNR, mean(NMF_ALS_Error(:, :, 1), 2));
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 1), 2));
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 1), 2));
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 1), 2));
grid on;
legend('ALS', 'MatLab ALS', 'Multiplicative', 'MatLab Multiplicative');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('j = 1');
subplot(2, 2, 2);
hold on;
plot(SNR, mean(NMF_ALS_Error(:, :, 2), 2));
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 2), 2));
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 2), 2));
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 2), 2));
grid on;
legend('ALS', 'MatLab ALS', 'Multiplicative', 'MatLab Multiplicative');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('j = 2');
subplot(2, 2, 3);
hold on;
plot(SNR, mean(NMF_ALS_Error(:, :, 3), 2));
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 3), 2));
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 3), 2));
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 3), 2));
grid on;
legend('ALS', 'MatLab ALS', 'Multiplicative', 'MatLab Multiplicative');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('j = 3');
subplot(2, 2, 4);
hold on;
plot(SNR, mean(NMF_ALS_Error(:, :, 4), 2));
plot(SNR, mean(MatLab_NMF_ALS_Error(:, :, 4), 2));
plot(SNR, mean(NMF_Multiplicative_Error(:, :, 4), 2));
plot(SNR, mean(MatLab_NMF_Multiplicative_Error(:, :, 4), 2));
grid on;
legend('ALS', 'MatLab ALS', 'Multiplicative', 'MatLab Multiplicative');
xlabel('SNR(dB)');
ylabel('mean ||E||_F');
title('j = 4');
%% Question 2)
%%
load('swimmer.mat');
swimmer = A;
A = zeros(length(swimmer), size(swimmer{1}, 1)*size(swimmer{1}, 2));
for i = 1:length(swimmer)
    A(i, :) = reshape(swimmer{i}, 1, []);
end
figure();
imagesc(A);
axis('off');
title('swimmer dataset');
Iter = 1:5;
J = 1:size(swimmer{1}, 1)*size(swimmer{1}, 2);
Swimmer_NMF_ALS_Error = zeros(length(Iter), length(J));
Swimmer_NMF_Multiplicative_Error = zeros(length(Iter), length(J));
p = 0;
for iter = Iter
    p = p+1;
    q = 0;
    for j = J
        q = q+1;
        B0 = rand(size(A, 1), j);
        C0 = rand(j, size(A, 2));
        [B_temp, C_temp] = NMF_ALS(A, j, B0, C0);
        Swimmer_NMF_ALS_Error(p, q) = norm(A-B_temp*C_temp, 'fro');
        [B_temp, C_temp] = NMF_Multiplicative(A, j, B0, C0);
        Swimmer_NMF_Multiplicative_Error(p, q) = norm(A-B_temp*C_temp, 'fro');
    end
end
figure();
plot(J, min(Swimmer_NMF_ALS_Error, [], 1));
xlim([min(J), max(J)]);
xlabel('number of coefficients');
ylabel('||E||_F');
title('ALS Algorithm');
figure();
plot(J, min(Swimmer_NMF_Multiplicative_Error, [], 1));
xlim([min(J), max(J)]);
xlabel('number of coefficients');
ylabel('||E||_F');
title('Multiplicative Algorithm');
%%
% ALS Algorithm
j = 10;
B0 = rand(size(A, 1), j);
C0 = rand(j, size(A, 2));
[B, C] = NMF_ALS(A, j, B0, C0);
Swimmer_NMF_ALS_Error = norm(A-B*C, 'fro');
figure();
for i = 1:10
    subplot(2, 5, i);
    imagesc(reshape(C(i, :), size(swimmer{1}, 1), size(swimmer{1}, 2)));
    axis('off');
end
sgtitle('ALS Algorithm');
% Multiplicative Algorithm
j = 16;
B0 = rand(size(A, 1), j);
C0 = rand(j, size(A, 2));
[B, C] = NMF_Multiplicative(A, j, B0, C0);
Swimmer_NMF_Multiplicative_Error = norm(A-B*C, 'fro');
figure();
for i = 1:16
    subplot(4, 4, i);
    imagesc(reshape(C(i, :), size(swimmer{1}, 1), size(swimmer{1}, 2)));
    axis('off');
end
sgtitle('Multiplicative Algorithm');
%% a
% ALS Algorithm
j = 20;
B0 = rand(size(A, 1), j);
C0 = rand(j, size(A, 2));
[B, C] = NMF_ALS(A, j, B0, C0);
Swimmer_NMF_ALS_Error = norm(A-B*C, 'fro');
figure();
for i = 1:20
    subplot(4, 5, i);
    imagesc(reshape(C(i, :), size(swimmer{1}, 1), size(swimmer{1}, 2)));
    axis('off');
end
sgtitle('ALS Algorithm');
% Multiplicative Algorithm
j = 20;
B0 = rand(size(A, 1), j);
C0 = rand(j, size(A, 2));
[B, C] = NMF_Multiplicative(A, j, B0, C0);
Swimmer_NMF_Multiplicative_Error = norm(A-B*C, 'fro');
figure();
for i = 1:20
    subplot(4, 5, i);
    imagesc(reshape(C(i, :), size(swimmer{1}, 1), size(swimmer{1}, 2)));
    axis('off');
end
sgtitle('Multiplicative Algorithm');
%% functions
function [B, C] = NMF_ALS(A, j, B0, C0)
    if j==size(B0, 2) && j==size(C0, 1)
        B = B0;
        C = C0;
        iter = 100;
        epsilon = 1e-16;
        for i = 1:iter
            B = max(epsilon, A*C'*pinv(C*C'));
            C = max(epsilon, pinv(B'*B)*B'*A);
        end
    else
        printf('Number of components is not right!');
    end
end
function [B, C] = NMF_Multiplicative(A, j, B0, C0)
    if j==size(B0, 2) && j==size(C0, 1)
        B = B0;
        C = C0;
        iter = 100;
        epsilon = 1e-16;
        for i = 1:iter
            B = B.*(A*C')./(B*(C*C')+epsilon);
            C = C.*(B'*A)./((B'*B)*C+epsilon);
        end
    else
        printf('Number of components is not right!');
    end
end