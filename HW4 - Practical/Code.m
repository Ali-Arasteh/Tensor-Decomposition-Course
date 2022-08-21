addpath('./TensorLab');
%% Question 1
% CP_ALS function
I1 = 4;
I2 = 3;
I3 = 2;
R = 2;
T = zeros(I1, I2, I3);
U1_org = randn(I1, R);
U2_org = randn(I2, R);
U3_org = randn(I3, R);
U_org = {U1_org, U2_org, U3_org};
for i = 1:R
    temp = U1_org(:, i)*U2_org(:, i)';
    T = T + reshape(temp(:)*U3_org(:, i)', size(T));
end
U1_0 = randn(I1, R);
U2_0 = randn(I2, R);
U3_0 = randn(I3, R);
U_0 = {U1_0, U2_0, U3_0};
U_ALS = CP_ALS(T, U1_0, U2_0, U3_0);
U1 = U_ALS{1};
U2 = U_ALS{2};
U3 = U_ALS{3};
recovered_T = zeros(I1, I2, I3);
for i = 1:R
    temp = U1(:, i)*U2(:, i)';
    recovered_T = recovered_T + reshape(temp(:)*U3(:, i)', size(recovered_T));
end
Error = T-recovered_T;
forbenius_error = norm(Error(:));
%% Question 2
%% Part a)
Iter = 1:50;
SNR = [0, 20, 40, 60];
ALS_Error = zeros(length(Iter), length(SNR), 4);
I1 = 6;
I2 = 6;
I3 = 6;
R = 3;
for iter = 1:50
    T = zeros(I1, I2, I3);
    U1_org = randn(I1, R);
    U2_org = randn(I2, R);
    U3_org = randn(I3, R);
    U_org = {U1_org, U2_org, U3_org};
    for i = 1:R
        temp = U1_org(:, i)*U2_org(:, i)';
        T = T + reshape(temp(:)*U3_org(:, i)', size(T));
    end
    N = randn(6, 6, 6);
    i = 0;
    for snr = SNR
        i = i+1;
        alpha = norm(T(:))/(norm(N(:))*10^(snr/20));
        noisy_T = T+alpha*N;
        U1_0 = randn(I1, R);
        U2_0 = randn(I2, R);
        U3_0 = randn(I3, R);
        U_0 = {U1_0, U2_0, U3_0};
        U_ALS = CP_ALS(noisy_T, U1_0, U2_0, U3_0);
        ALS_Error(iter, i, 1) = TMSFE(U_org, U_ALS);
        U_MatLab_ALS = cpd_als(noisy_T, U_0);
        ALS_Error(iter, i, 2) = TMSFE(U_org, U_MatLab_ALS);
        U_MatLab_SD = cpd3_sd(noisy_T, U_0);
        ALS_Error(iter, i, 3) = TMSFE(U_org, U_MatLab_SD);
        U_MatLab_minf = cpd_minf(noisy_T, U_0);
        ALS_Error(iter, i, 4) = TMSFE(U_org, U_MatLab_minf);
    end
end
ave = squeeze(mean(ALS_Error, 1));
figure();
hold on;
for i = 1:4
    plot(SNR, ave(:, i));
end
legend('ALS', 'MatLab ALS', 'MatLab SD', 'MatLab MINF');
xlabel('SNR(dB)');
ylabel('TMSFE');
%% Part b)
Iter = 1:50;
SNR = [0, 20, 40, 60];
ALS_with_HOSVD_Error = zeros(length(Iter), length(SNR), 4);
I1 = 6;
I2 = 6;
I3 = 6;
R = 3;
for iter = 1:50
    T = zeros(I1, I2, I3);
    U1_org = randn(I1, R);
    U2_org = randn(I2, R);
    U3_org = randn(I3, R);
    U_org = {U1_org, U2_org, U3_org};
    for i = 1:R
        temp = U1_org(:, i)*U2_org(:, i)';
        T = T + reshape(temp(:)*U3_org(:, i)', size(T));
    end
    N = randn(6, 6, 6);
    i = 0;
    for snr = SNR
        i = i+1;
        alpha = norm(T(:))/(norm(N(:))*10^(snr/20));
        noisy_T = T+alpha*N;
        U_0 = HOSVD(noisy_T, R);
        U_ALS = CP_ALS(noisy_T, U1_0, U2_0, U3_0);
        ALS_with_HOSVD_Error(iter, i, 1) = TMSFE(U_org, U_ALS);
        U_MatLab_ALS = cpd_als(noisy_T, U_0);
        ALS_with_HOSVD_Error(iter, i, 2) = TMSFE(U_org, U_MatLab_ALS);
        U_MatLab_SD = cpd3_sd(noisy_T, U_0);
        ALS_with_HOSVD_Error(iter, i, 3) = TMSFE(U_org, U_MatLab_SD);
        U_MatLab_minf = cpd_minf(noisy_T, U_0);
        ALS_with_HOSVD_Error(iter, i, 4) = TMSFE(U_org, U_MatLab_minf);
    end
end
ave = squeeze(mean(ALS_with_HOSVD_Error, 1));
figure();
hold on;
for i = 1:4
    plot(SNR, ave(:, i));
end
legend('ALS with HOSVD', 'MatLab ALS with HOSVD', 'MatLab SD with HOSVD', 'MatLab MINF with HOSVD');
xlabel('SNR(dB)');
ylabel('TMSFE');
%% Part c)
Iter = 1:50;
SNR = [0, 20, 40, 60];
ALS_Error_correlated = zeros(length(Iter), length(SNR), 4);
I1 = 6;
I2 = 6;
I3 = 6;
R = 3;
for iter = 1:50
    T = zeros(I1, I2, I3);
    U1_org = randn(I1, R);
    U1_org(:, 2) = U1_org(:, 1)+0.5*U1_org(:, 2);
    U2_org = randn(I2, R);
    U2_org(:, 2) = U2_org(:, 1)+0.5*U2_org(:, 2);
    U3_org = randn(I3, R);
    U_org = {U1_org, U2_org, U3_org};
    for i = 1:R
        temp = U1_org(:, i)*U2_org(:, i)';
        T = T + reshape(temp(:)*U3_org(:, i)', size(T));
    end
    N = randn(6, 6, 6);
    i = 0;
    for snr = SNR
        i = i+1;
        alpha = norm(T(:))/(norm(N(:))*10^(snr/20));
        noisy_T = T+alpha*N;
        U1_0 = randn(I1, R);
        U2_0 = randn(I2, R);
        U3_0 = randn(I3, R);
        U_0 = {U1_0, U2_0, U3_0};
        U_ALS = CP_ALS(noisy_T, U1_0, U2_0, U3_0);
        ALS_Error_correlated(iter, i, 1) = TMSFE(U_org, U_ALS);
        U_MatLab_ALS = cpd_als(noisy_T, U_0);
        ALS_Error_correlated(iter, i, 2) = TMSFE(U_org, U_MatLab_ALS);
        U_MatLab_SD = cpd3_sd(noisy_T, U_0);
        ALS_Error_correlated(iter, i, 3) = TMSFE(U_org, U_MatLab_SD);
        U_MatLab_minf = cpd_minf(noisy_T, U_0);
        ALS_Error_correlated(iter, i, 4) = TMSFE(U_org, U_MatLab_minf);
    end
end
ave = squeeze(mean(ALS_Error_correlated, 1));
figure();
hold on;
for i = 1:4
    plot(SNR, ave(:, i));
end
legend('ALS - correlated', 'MatLab ALS - correlated', 'MatLab SD - correlated', 'MatLab MINF - correlated');
xlabel('SNR(dB)');
ylabel('TMSFE');
%% Part d)
Iter = 1:50;
SNR = [0, 20, 40, 60];
ALS_with_HOSVD_Error_correlated = zeros(length(Iter), length(SNR), 4);
I1 = 6;
I2 = 6;
I3 = 6;
R = 3;
for iter = 1:50
    T = zeros(I1, I2, I3);
    U1_org = randn(I1, R);
    U1_org(:, 2) = U1_org(:, 1)+0.5*U1_org(:, 2);
    U2_org = randn(I2, R);
    U2_org(:, 2) = U2_org(:, 1)+0.5*U2_org(:, 2);
    U3_org = randn(I3, R);
    U_org = {U1_org, U2_org, U3_org};
    for i = 1:R
        temp = U1_org(:, i)*U2_org(:, i)';
        T = T + reshape(temp(:)*U3_org(:, i)', size(T));
    end
    N = randn(6, 6, 6);
    i = 0;
    for snr = SNR
        i = i+1;
        alpha = norm(T(:))/(norm(N(:))*10^(snr/20));
        noisy_T = T+alpha*N;
        U_0 = HOSVD(noisy_T, R);
        U_ALS = CP_ALS(noisy_T, U1_0, U2_0, U3_0);
        ALS_with_HOSVD_Error_correlated(iter, i, 1) = TMSFE(U_org, U_ALS);
        U_MatLab_ALS = cpd_als(noisy_T, U_0);
        ALS_with_HOSVD_Error_correlated(iter, i, 2) = TMSFE(U_org, U_MatLab_ALS);
        U_MatLab_SD = cpd3_sd(noisy_T, U_0);
        ALS_with_HOSVD_Error_correlated(iter, i, 3) = TMSFE(U_org, U_MatLab_SD);
        U_MatLab_minf = cpd_minf(noisy_T, U_0);
        ALS_with_HOSVD_Error_correlated(iter, i, 4) = TMSFE(U_org, U_MatLab_minf);
    end
end
ave = squeeze(mean(ALS_with_HOSVD_Error_correlated, 1));
figure();
hold on;
for i = 1:4
    plot(SNR, ave(:, i));
end
legend('ALS with HOSVD - correlated', 'MatLab ALS with HOSVD - correlated', 'MatLab SD with HOSVD - correlated', 'MatLab MINF with HOSVD - correlated');
xlabel('SNR(dB)');
ylabel('TMSFE');
%% Question 3
%% Part a)
load('amino.mat');
I1 = DimX(1);
I2 = DimX(2);
I3 = DimX(3);
T = reshape(X, DimX);
R = [2, 3, 4, 5];
i = 0;
U = cell(6, length(R));
titles = ["ALS", "MatLab ALS", "MatLab SD", "MatLab MINF", "MatLab NLS", "MatLab SGSD"];
for r = R
    i = i+1;
    U1_0 = randn(I1, r);
    U2_0 = randn(I2, r);
    U3_0 = randn(I3, r);
    U_0 = {U1_0, U2_0, U3_0};
    U{1, i} = CP_ALS(T, U1_0, U2_0, U3_0);
    U{2, i} = cpd_als(T, U_0);
    U{3, i} = cpd3_sd(T, U_0);
    U{4, i} = cpd_minf(T, U_0);
    U{5, i} = cpd_nls(T, U_0);
    U{6, i} = cpd3_sgsd(T, U_0);
    figure();
    for l = 1:6
        subplot(3, 2, l);
        hold on;
        for s = 1:r
           plot(EmAx, U{l, i}{2}(:, s));
        end
        title(titles(l));
    end
    sgtitle(['Emission with ', num2str(r), ' components']);
    figure();
    for l = 1:6
        subplot(3, 2, l);
        hold on;
        for s = 1:r
           plot(ExAx, U{l, i}{3}(:, s));
        end
        title(titles(l));
    end
    sgtitle(['Excitation with ', num2str(r), ' components']);
end
%% Part b)
i = 0;
for r = R
    i = i+1;
    for l = 1:6
        figure();
        [consistency, ~, ~, ~] = corcond(T, U{l, i}, [], 1);
        title(titles(l) + ", consistency = " + num2str(consistency) + ", R = " + num2str(r));
    end
end
%% functions
function U = CP_ALS(T, U1_0, U2_0, U3_0)
    U1 = U1_0;
    U2 = U2_0;
    U3 = U3_0;
    number_of_iterations = 100;
    [I1, I2, I3] = size(T);
    T1 = zeros(I1, I2*I3);
    for i3 = 1:I3
        for i2 = 1:I2
            T1(:, (i3-1)*I2+i2) = T(:, i2, i3);
        end
    end
    T2 = zeros(I2, I1*I3);
    for i3 = 1:I3
        for i1 = 1:I1
            T2(:, (i3-1)*I1+i1) = T(i1, :, i3);
        end
    end
    T3 = zeros(I3, I1*I2);
    for i2 = 1:I2
        for i1 = 1:I1
            T3(:, (i2-1)*I1+i1) = T(i1, i2, :);
        end
    end
    for i = 1:number_of_iterations
        U1 = T1*pinv(kr(U3, U2)');
        U2 = T2*pinv(kr(U3, U1)');
        U3 = T3*pinv(kr(U2, U1)');
    end
    U = {U1, U2, U3};
end
function U = HOSVD(T, R)
    [I1, I2, I3] = size(T);
    T1 = zeros(I1, I2*I3);
    for i3 = 1:I3
        for i2 = 1:I2
            T1(:, (i3-1)*I2+i2) = T(:, i2, i3);
        end
    end
    [U1, ~, ~] = svd(T1);
    U1 = U1(:, 1:R);
    T2 = zeros(I2, I1*I3);
    for i3 = 1:I3
        for i1 = 1:I1
            T2(:, (i3-1)*I1+i1) = T(i1, :, i3);
        end
    end
    [U2, ~, ~] = svd(T1);
    U2 = U2(:, 1:R);
    T3 = zeros(I3, I1*I2);
    for i2 = 1:I2
        for i1 = 1:I1
            T3(:, (i2-1)*I1+i1) = T(i1, i2, :);
        end
    end
    [U3, ~, ~] = svd(T1);
    U3 = U3(:, 1:R);
    U = {U1, U2, U3};
end